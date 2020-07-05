//////////////////////////////////////////////////////////////////////
// Main.cc
//////////////////////////////////////////////////////////////////////

#include <climits>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "Assert.h"
#include "MultiSequence.h"
#include "ProbModel.h"
#include "Matrix.h"
#include "ScoreMatrix.h"
#include "SparseMatrix.h"
#include "GlobalAlign.h"
#include "Tree.h"
#include "Utilities.h"
#include "AlignedFragment.h"
#include "LocalAlign.h"
#include "Block.h"
#include "Consistency.h"
#include "PairAligner.h"
#include "Types.h"

bool verbose = true;
int MINLENGTH = 30; // shortest local alignment
bool enableViterbi = true;
bool enableTransitivity = false;
bool fastaOutput = false;

typedef SparseMatrix *SparseMatrixPtr;

bool GetInteger (char *data, int *val);
void RunEM (ProbModel &hmm, int numFilenames, string *filenames);
void RunLocalAligner (int numFilenames, string *filenames);

//////////////////////////////////////////////////////////////////////
// main program
//////////////////////////////////////////////////////////////////////

int main (int argc, char **argv){
  
  PRECOMPUTE_SCORE_TABLES();

  
  fprintf (stderr, "ProDA version 1.0\n\n");

  // usage
  
  if (argc == 1){
    fprintf (stderr, "Usage: proda [-L length] [-silent] [-posterior] [-tran] [-fasta] filename(s)\n");
    return 0;
  }

  // find all command-line flags

  bool computeEMParams = false;
  int numEMreps = 0;
  string *filenames = new string[argc];
  ASSERT (filenames, "Out of memory.");
  int numFilenames = 0;
  
  int i;
  for (i = 1; i < argc; i++){
    if (argv[i][0] == '-'){
      if(!strcmp (argv[i], "-L")){
		  if( i < argc -1){
			  if(!GetInteger (argv[++i], &MINLENGTH)){
				  fprintf(stderr,"ERROR: Invalid integer following option %s :: %s \n",argv[i-1], argv[i]);
				  exit(1);
			  }
		  }
	  }else if(!strcmp(argv[i],"-posterior")){
		  enableViterbi = false;
	  }
	  else if(!strcmp(argv[i],"-tran")){
		  enableTransitivity = true;
	  }
	  else if(!strcmp(argv[i],"-fasta")){
		  fastaOutput = true;
	  }
	  else if(!strcmp(argv[i],"-silent")){
		  verbose = false;
	  }
	  else {
	fprintf (stderr, "Unknown parameter ignored: %s\n", argv[i]);
      }
    } else {
      filenames[numFilenames++] = StrDup (argv[i]);
    } 
  }

  fprintf (stderr,"Minimal block length = %d\n\n",MINLENGTH);

  // run program

  if (computeEMParams){
    ProbModel hmm;
    for (int i = 0; i < numEMreps; i++)
      RunEM (hmm, numFilenames, filenames);
  } else {
    RunLocalAligner (numFilenames, filenames);
  }

  // free memory

  for (i = 0; i < numFilenames; i++)
	  delete filenames[i];
  delete filenames;

  return 0;
}


bool GetInteger (char *data, int *val){
  char *endPtr;
  long int retVal;

  int errno = 0;
  retVal = strtol (data, &endPtr, 0);
  if (retVal == 0 && (errno != 0 || data == endPtr)) return false;
  if (errno != 0 && (retVal == LONG_MAX || retVal == LONG_MIN)) return false;
  if (retVal < (long) INT_MIN || retVal > (long) INT_MAX) return false;
  *val = (int) retVal;
  return true;
}
//////////////////////////////////////////////////////////////////////
// Run EM algorithm
//////////////////////////////////////////////////////////////////////

void RunEM (ProbModel &hmm, int numFilenames, string *filenames){

  Matrix cts (1, NUM_STATES+1, NUM_STATES+1);
  cts.Fill (0);

  for (int k = 0; k < numFilenames; k++){
    
    // load alignment
    MultiSequence *seqs = new MultiSequence();
    ASSERT (seqs, "Out of memory.");
    fprintf (stderr, "Loading: %s\n", filenames[k]);
    seqs->LoadSequences (filenames[k]);
    int n = seqs->GetNumSequences();
    ASSERT (n > 0, "No sequences to align!");
        
    // compute sufficient statistics
    for (int i = 0; i < n; i++){
      for (int j = i+1; j < n; j++){
	Matrix *ep = hmm.ComputeExpectedCounts (seqs->GetSequence (i), seqs->GetSequence (j));
	Matrix &e = *ep;
	
	for (int a = 0; a < NUM_STATES+1; a++)
	  for (int b = 0; b < NUM_STATES+1; b++)
	    cts(0,a,b) += e(0,a,b) / (n*n);
	
	delete ep;
      }
    }
  }

  // compute new params
  hmm.ComputeParams (&cts);
}




//////////////////////////////////////////////////////////////////////
// Run local alignment procedure
//////////////////////////////////////////////////////////////////////

void RunLocalAligner (int numFilenames, string *filenames){
	  

  ProbModel hmm;
  int i,j;

  MultiSequence *seqs = new MultiSequence();
  ASSERT (seqs, "Out of memory.");
  
  // load input sequences
  for (i = 0; i < numFilenames; i++){
    fprintf (stderr, "Loading: %s\n", filenames[i]);
    seqs->LoadSequences (filenames[i]);
  }
  int n = seqs->GetNumSequences();
  ASSERT (n > 0, "No sequences to align!");


  AVECT fragments;


  // compute pairwise similarities and local aligned fragments
  fprintf(stderr,"\nAligning all pairs of sequences. This may take several minutes\n");
  for (i = 0; i < n; i++){
    for (j = i + 1; j < n; j++){

		Sequence seq1 = seqs->GetSequence (i);
		Sequence seq2 = seqs->GetSequence (j);

		PairAligner pAligner(&hmm,&seq1,&seq2);
		
		if(enableViterbi) pAligner.FastPairAlign(fragments);
		else pAligner.PairAlign(fragments);
	}
  }

  char file0[260];
  strcpy(file0,filenames[0]);
  for(i = strlen(file0)-1; i > 0 && file0[i]!='.'; i--);
  if(i > 0) file0[i] = 0;
  FILE *fasta;
  if(fastaOutput){
	  strcat(file0,".fasta");
	  fasta = fopen(file0,"w");
  }
  strcat(file0,".test");
  FILE *output = fopen(file0,"w");
  Block *prohibited = NULL;

  while (fragments.size() > 0) {
        //Form a block
	    fprintf(stderr,"\nForming block\n");
		Block *block = new Block(fragments, seqs, enableTransitivity, prohibited);
		prohibited = NULL;
		int flag = 1;

		fprintf(stderr,"Aligning block\n");
		do{
		  int m = block->size();
		  SparseMatrix **block_posteriors = new SparseMatrixPtr[m*m];
		  for (i = 0; i < m*m; i++) 
			  block_posteriors[i] = NULL;
		  Matrix block_similarity (1, m, m);
		  block_similarity.Fill (1);

		  MultiSequence *block_seqs = new MultiSequence();
		  for (i = 0; i < m; i++){
			  Fragment fm = (*block)[i];
			  Sequence *s = new Sequence(seqs->GetSequence(fm.id));
			  s->Clip(fm.begin, fm.end);
			  s->SetID(i);
			  block_seqs->AddSequence(s);
		  }
		
		  for (i = 0; i < m-1; i++){
			  for (j = i+1; j < m; j++){
				ScoreMatrix *p = hmm.Posterior (block_seqs->GetSequence (i), block_seqs->GetSequence (j),MATCH);
				Matrix *p2 = new Matrix (*p);
				delete p;

				float score;
				int length;
				char *ali = GlobalAlign::ComputeMWTrace(*p2, &score, &length);
				delete ali;
				block_similarity (0,i,j) = block_similarity (0,j,i) = length >= MINLENGTH ? score / length:0;
				block_posteriors[i*m+j] = new SparseMatrix (*p2, 0.01, 0);
				delete p2;
				block_posteriors[j*m+i] = block_posteriors[i*m+j]->ComputeTranspose();
			  }
		  }
		  // compute expected accuracy tree

		  Tree tree (block_similarity, *block_seqs, 0.5);

		  // if not all fragments are related
		  int mm = tree.GetNumSequences();
		  if (mm < m && mm >= 2){
			  block->AdjustAFragmentList(fragments, seqs->GetNumSequences(), &block_similarity, 0.5);
			  IVECT ids;
			  tree.GetIDs(ids);
			  ASSERT(mm == (int)ids.size(), "Wrong tree size\n");
			  //Form new posterior matrix
			  SparseMatrix **newPosteriors = new SparseMatrixPtr[mm*mm];
			  for (i = 0; i < mm*mm; i++)
				  newPosteriors[i] = NULL;
			  for (i = 0; i < mm; i++){
				  for (j = 0; j < mm; j++) if(i != j){
						newPosteriors[i*mm+j] = block_posteriors[ids[i]*m+ids[j]];
				  }
			  }
			  int *used = new int[m];
			  for(i=0; i<m;i++) used[i] = 0;
			  for(i=0;i<mm;i++) used[ids[i]] = 1;
			  for (i = 0; i < m; i++)
				  for (j = 0; j < m; j++)
					  if (i != j && (used[i]==0 || used[j]==0)) delete block_posteriors[i*m+j];
			  delete [] block_posteriors;
			  delete used;
			  Block *newBlock = new Block();
			  for (i=0; i<mm;i++){
					  newBlock->AddFragment((*block)[ids[i]]);
			  }
			  delete block;
			  block = newBlock;
			  block_posteriors = newPosteriors;

			  m = mm;
		  }
		  if (mm < 2) {
			  block->AdjustAFragmentList(fragments, seqs->GetNumSequences(), &block_similarity, 0.5);
			  delete block;
			  delete block_seqs;
			  for (i = 0; i < m*m; i++)
					if (block_posteriors[i]) delete block_posteriors[i];
			  delete [] block_posteriors;
			  break;
		  }
		  // probabilistic consistency
		  
		  for (int c = 0; c < 0; c++){
			  SparseMatrix **newPosteriors = ProbabilisticConsistency (block_posteriors, m);
			  for (i = 0; i < m; i++)
				  for (j = 0; j < m; j++)
					  if (i != j) delete block_posteriors[i*m+j];
			  delete [] block_posteriors;
			  block_posteriors = newPosteriors;
		  }

		  MultiSequence *result = tree.ProgressiveAlignment (m, block_posteriors);
		  result->Sort();

  		  delete block_seqs;
		  for (i = 0; i < m*m; i++)
			if (block_posteriors[i]) delete block_posteriors[i];
		  delete [] block_posteriors;


		  int start, end;
		  result->FindBlock(*block, start, end, int(MINLENGTH*0.9));
		  if(block->size() == m){
			  if(end-start+1 >= MINLENGTH*0.9){
				block->PrintBlock(stdout,seqs);
				block->PrintBlock(output,seqs,1);
				result->WriteCLUSTALW(stdout,start, end);
				if(fastaOutput)
					seqs->WriteFASTA(fasta,block,result,start,end);
			  }
			int deleted = block->AdjustAFragmentList(fragments, seqs->GetNumSequences());
			if(deleted) delete block;
			else prohibited = block;
			flag = 0;
		  }
		  delete result;
		  
		}while(flag);
	  }

	  fclose(output);
	  if(fastaOutput)
		  fclose(fasta);

	delete seqs;
}

