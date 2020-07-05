//////////////////////////////////////////////////////////////////////
// GlobalAlign.cc
//////////////////////////////////////////////////////////////////////

#include <string.h>
#include "Assert.h"
#include "GlobalAlign.h"
#include "Utilities.h"

//////////////////////////////////////////////////////////////////////
// Compute maximum weight trace
//////////////////////////////////////////////////////////////////////

char *GlobalAlign::ComputeMWTrace (const Matrix &m, float *score, int *length){
  const int rows = m.GetNumRows();
  const int cols = m.GetNumCols();
  int i,j;
  int beginx, beginy, endx, endy;
  beginx = beginy = endx = endy = -1;
  // memory allocation
  Matrix *bestp = new Matrix (1, rows, cols);
  ASSERT (bestp, "Out of memory.");
  Matrix &best = *bestp;
  TracebackType *traceback = new TracebackType[rows * cols];
  ASSERT (traceback, "Out of memory.");
  
  // compute best path
  best.Fill (LOG_ZERO_FLOAT);
  best(0,0,0) = 0;
  for (i = 0; i < rows; i++){
    for (j = 0; j < cols; j++){
      traceback[i * cols + j] = NONE;
      if (i > 0 && best(0,i-1,j) > best(0,i,j)){
	best(0,i,j) = best(0,i-1,j);
	traceback[i * cols + j] = UP;
      }
      if (j > 0 && best(0,i,j-1) > best(0,i,j)){
	best(0,i,j) = best(0,i,j-1);
	traceback[i * cols + j] = LEFT;
      }
      if (i > 0 && j > 0 && best(0,i-1,j-1) + m(0,i,j) > best(0,i,j)){
	best(0,i,j) = best(0,i-1,j-1) + m(0,i,j);
	traceback[i * cols + j] = UP_LEFT;
	if(beginx == -1) {beginx = i;beginy = j;}
	endx = i; endy = j;
      }
    }
  }

  if (score){
	  *score = best(0,rows-1,cols-1);
	  *length = min(endy - beginy, endx-beginx) + 1;
  }
  
  delete bestp;
  
  // follow tracebacks
  char *buffer = new char[rows * cols];
  ASSERT (buffer, "Out of memory.");
  int r = rows-1, c = cols-1, len = 0;
  while (traceback[r * cols + c] != NONE){
    switch (traceback[r * cols + c]){
    case UP: r--; buffer[len++] = 'X'; break;
    case LEFT: c--; buffer[len++] = 'Y'; break;
    case UP_LEFT: r--; c--; buffer[len++] = 'B'; break;
    default: ASSERT (false, "Unexpected value found in traceback matrix!");
    }
  }
  
  delete[] traceback;
  
  // reverse alignment path
  char *ret = new char[len+1];
  ASSERT (ret, "Out of memory.");
  for (i = 0; i < len; i++)
    ret[i] = buffer[len - 1 - i];
  ret[len] = '\0';
  
  delete[] buffer;
  
  return ret;
}


//////////////////////////////////////////////////////////////////////
// Insert gaps into aligned sequence
//////////////////////////////////////////////////////////////////////

Sequence *GlobalAlign::InsertGaps (const Sequence &seq, const char *alignmentPath, char ch){
  int len = strlen (alignmentPath);
  const char *data = seq.GetData();
  char *newData = new char[len+2];
  ASSERT (newData, "Out of memory.");

  char *newName = new char[strlen(seq.GetName())+1];
  ASSERT (newName, "Out of memory.");

  newData[0] = '@';
  newData[len+1] = '\0';

  int j = 1;
  for (int i = 0; i < len; i++){
    if (alignmentPath[i] == ch || alignmentPath[i] == 'B')
      newData[i+1] = data[j++];
    else
      newData[i+1] = '-';
  }

  memcpy (newName, seq.GetName(), strlen(seq.GetName())+1);
  
  return new Sequence (newData, newName, len, seq.GetID());
}

//////////////////////////////////////////////////////////////////////
// Build alignment from alignment path
//////////////////////////////////////////////////////////////////////

MultiSequence *GlobalAlign::BuildAlignment (const MultiSequence &group1, 
					    const MultiSequence &group2, 
					    const char *alignmentPath){
  MultiSequence *ret = new MultiSequence();
  ASSERT (ret, "Out of memory.");
  
  for (int i = 0; i < group1.GetNumSequences(); i++)
    ret->AddSequence (InsertGaps (group1.GetSequence(i), alignmentPath, 'X'));
  {for (int i = 0; i < group2.GetNumSequences(); i++)
	  ret->AddSequence (InsertGaps (group2.GetSequence(i), alignmentPath, 'Y'));}
  
  return ret;
}

//////////////////////////////////////////////////////////////////////
// Align two groups of sequences
//////////////////////////////////////////////////////////////////////

MultiSequence *GlobalAlign::AlignGroups (int n,
					 SparseMatrix **posteriors, 
					 const MultiSequence &group1, 
					 const MultiSequence &group2){
  
  int groupLen1 = group1.GetLength();
  int groupLen2 = group2.GetLength();
  Matrix *mp = new Matrix (1, groupLen1+1, groupLen2+1);
  Matrix &m = *mp;
  m.Fill (0);

  for (int s = 0; s < group1.GetNumSequences(); s++){
    for (int t = 0; t < group2.GetNumSequences(); t++){
      int id1 = group1.GetSequence(s).GetID();
      int id2 = group2.GetSequence(t).GetID();
      const SparseMatrix &sm = *posteriors[id1*n+id2];

      int *mapping1 = group1.GetSequence(s).ComputeMapping();
      int *mapping2 = group2.GetSequence(t).ComputeMapping();
      
      for (int i = 0; i < sm.GetNumRows(); i++){
	const SparseMatrix::SparseMatrixEntry *p = sm.GetRowPtr (0, i);
	for (int j = 0; j < sm.GetRowSize(0,i); j++){
	  int gr = mapping1[i];
	  int gc = mapping2[p->column];
	  m(0,gr,gc) += p->value;
	  ++p;
	}
      }
      
      delete[] mapping1;
      delete[] mapping2;
    }
  }
  const char *path = ComputeMWTrace (m);
  delete mp;

  MultiSequence *res = BuildAlignment (group1, group2, path);
  delete[] (char *)path;
  

  return res;
}
				
//////////////////////////////////////////////////////////////////////
// Compute alignment score
//////////////////////////////////////////////////////////////////////

float GlobalAlign::ComputeAlignmentScore (const MultiSequence &seqs, 
					  int n, 
					  SparseMatrix **posteriors){

  int i;
  float score = 0;
  int length = seqs.GetSequence(0).GetLength()+1;
  float *scores = new float[length];
  for(i =0; i < length; i++)
	  scores[i] = 0;
  for (int s = 0; s < n-1; s++){
    for (int t = s+1; t < n; t++){
      int id1 = seqs.GetSequence(s).GetID();
      int id2 = seqs.GetSequence(t).GetID();
      const SparseMatrix &sm = *posteriors[id1*n+id2];
      int *mapping1 = seqs.GetSequence(s).ComputeMapping();
      int *mapping2 = seqs.GetSequence(t).ComputeMapping();
      int len1 = sm.GetNumRows()-1;
      int len2 = sm.GetNumCols()-1;
      
      int pos1 = 1;
      int pos2 = 1;
      while (pos1 <= len1 && pos2 <= len2){
	if (mapping1[pos1] < mapping2[pos2]) pos1++;
	else if (mapping1[pos1] > mapping2[pos2]) pos2++;
	else {
	  score += sm(0,pos1,pos2);
	  scores[mapping1[pos1]] += sm(0,pos1,pos2);
	  pos1++;
	  pos2++;
	}
    }
      
      delete [] mapping1;
      delete [] mapping2;      
    }
  }

  for(i = 1; i < length; i++){
	  fprintf(stderr, "%4.2f ", scores[i]); 
  }
  fprintf(stderr,"\n");
  delete scores;
  
  return score;
}
