////////////////////////////////////////////////////////////////////
// PairAligner.cc
// 
// Implementation of PairAligner class
////////////////////////////////////////////////////////////////////

#include <cstring>
#include <cstdlib>

#include "PairAligner.h"
#include "Utilities.h"
#include "LocalAlign.h"

extern bool verbose;
extern int MINLENGTH;
extern bool enableViterbi;

PairAligner::PairAligner(ProbModel *v_hmm, Sequence *s1, Sequence *s2)
{
	hmm = v_hmm;
	seq1 = s1;
	seq2 = s2;
	xLen = seq1->GetLength();
	yLen = seq2->GetLength();
	map = new int[(xLen+1)* (yLen + 1)];
	ASSERT (map,"Out of memory");
	memset(map,0,sizeof(int)*(xLen + 1) * (yLen + 1));
}

///////////////////////////////////////////////////////////////////////////////
// Update Map to disallow match state
///////////////////////////////////////////////////////////////////////////////
void PairAligner::UpdateMap(AlignedFragment *frag, int self)
{
	int i,k;
	for (k = frag->begin[0],i=0; k <= frag->end[0]; k++,i++){
		if(frag->seq[0][i]!=-1){
			int y = frag->seq[0][i];
			int x;
			for(x = y-MINLENGTH+1; x <y+MINLENGTH;x++){
				if(x >0 && x <= yLen)
					map[k*(yLen+1)+x] = 1;
				if(self && x > 0 && x <= xLen)
					map[x*(yLen+1)+k] = 1;
			}
			for (int z = k-MINLENGTH+1; z < k +MINLENGTH;z++){
				if(z > 0 && z < xLen)
					map[z*(yLen+1)+y] = 1;
				if(self && z >0 && z <= yLen)
					map[y*(yLen+1)+z] = 1;
			}
		}
	}

}

void PairAligner::ConsistencyCheck(AVECT &pair_frags)
{
	int i,k;
	int second;
	int bound[2][6];
	AVECT::iterator it,jt;
	int flag = 0;//no change
	while(!flag){
		flag = 1;
	for (it = pair_frags.begin() ; it != pair_frags.end() && flag; it++){
		AlignedFragment afi = *it;
		bound[0][0] = bound[0][1] = bound[0][2] =afi.begin[0];
		bound[0][3] = bound[0][4] = bound[0][5] =afi.end[0];
		bound[1][0] = bound[1][1] = bound[1][2] = afi.begin[1];
		bound[1][3] = bound[1][4] = bound[1][5] = afi.end[1];
		
		for ( jt = it+1; jt != pair_frags.end() && flag; jt++){
			AlignedFragment afj = *jt;
			for (k = 0; k < 2; k++){
				if(afi.id[0] == afj.id[k] && afi.id[1] == afj.id[1-k] &&
				abs(afi.begin[0] - afj.begin[k]) <=2 &&
				abs(afi.begin[1] - afj.begin[1-k]) <=2 &&
				abs(afi.end[0] - afj.end[k] <=2) &&
				abs(afi.end[1] - afj.end[1-k]) <=2)
				{
					pair_frags.erase(jt);
					flag = 0;
					continue;
				}
				int o_flag = 0;
				if(afi.id[0] == afj.id[k] && afi.id[1] == afj.id[1-k] &&
					Overlap(afi.begin[0],afi.end[0],afj.begin[k],afj.end[k]) > 0 &&
					Overlap(afi.begin[1],afi.end[1],afj.begin[1-k],afj.end[1-k]) > 0){
					int start[2],finish[2];
					PAIRI *pa;
					int startx = max(afi.begin[0],afj.begin[k]);
					PAIRI *xi = afi.GetAlignPos(afi.id[0],startx,second);
					PAIRI *xj = afj.GetAlignPos(afj.id[k],startx,second);
					if(xi->second < xj->second){
						start[0] = max(afi.begin[0],afj.begin[k]);
						pa = afj.GetAlignPos(afj.id[k],start[0],second);
						start[1] = pa->second;
						if(start[1] >= afi.end[1]) {
							pair_frags.erase(jt);
							flag = 0;
							continue;
						}
						finish[1] = min(afi.end[1],afj.end[1-k]);
						pa = afj.GetAlignPos(afj.id[1-k],finish[1],second);
						finish[0] = pa->second;
					}
					else{
						start[1] = max(afi.begin[1],afj.begin[1-k]);
						pa = afj.GetAlignPos(afj.id[1-k],start[1],second);
						start[0] = pa->second;
						if(start[0] >= afi.end[0]) {
							pair_frags.erase(jt);
							flag = 0;
							continue;
						}
						finish[0] = min(afi.end[0],afj.end[k]);
						pa = afj.GetAlignPos(afj.id[k],finish[0],second);
						finish[1] = pa->second;
					}
						PAIRI *imap = afi.GetAlignPos(afi.id[0],start[0],second);
						if(imap != NULL){
							bound[0][1] = start[0];
							bound[1][1] = imap->second;
							delete imap;
						}
						imap = afi.GetAlignPos(afi.id[0],finish[0],second);
						if(imap != NULL){
							bound[0][2] = finish[0];
							bound[1][2] = imap->second;
							delete imap;
						}
						imap = afi.GetAlignPos(afi.id[1],start[1],second);
						if(imap != NULL){
							bound[1][3] = start[1];
							bound[0][3] = imap->second;
							delete imap;
						}
						imap = afi.GetAlignPos(afi.id[1],finish[1],second);
						if(imap != NULL){
							bound[1][4] = finish[1];
							bound[0][4] = imap->second;
							delete imap;
						}
						if(bound[0][3] < bound[0][1]){
							for(int c = 0; c < 2; c++){
								swap(bound[c][1],bound[c][3]);
								swap(bound[c][2],bound[c][4]);
							}
						}
						if(bound[0][3] < bound[0][2]){
							swap(bound[0][2],bound[0][3]);
							swap(bound[1][2],bound[1][3]);
							o_flag = 1;
						}
					it = pair_frags.erase(it);
					AlignedFragment *saf;
					int a,b;
					for (i = 0; i < 5; i++){
						if( i == 2 && o_flag) continue;
						if(i==2 || i ==4)  a = 1;
						else a = 0;
						if(i == 0 || i == 2) b = 1;
						else b = 0;
						if(bound[0][i+1] - bound[0][i] >= MINLENGTH && bound[1][i+1] - bound[1][i] >= MINLENGTH){
							Sequence newSeq1 = *seq1;
							Sequence newSeq2 = *seq2;
							newSeq1.SubStr(bound[0][i]+a,bound[0][i+1]-b);
							newSeq2.SubStr(bound[1][i]+a,bound[1][i+1]-b);
							int *l_map = new int[(xLen+1)* (yLen + 1)];
							ASSERT (l_map,"Out of memory");
							memset(l_map,0,sizeof(int)*(xLen + 1) * (yLen + 1));
							saf = hmm->Viterbi(newSeq1,newSeq2,l_map);
							if(saf != NULL && saf->GetLength() >= MINLENGTH){
								saf->ShiftRight(bound[0][i]+a,bound[1][i]+a);
									pair_frags.push_back(*saf);
							}
							if(saf) delete saf;
							delete l_map;
						}
					}
					flag = 0;
					break;
				}				
			}
		}
	}
	}
}

void PairAligner::FastPairAlign(AVECT &fragments)
{
	AVECT pair_frags;
	int i;
	Sequence &wseq1 = *seq1;
	Sequence &wseq2 = *seq2;
	int self = 0;

	if(wseq1.GetID() == wseq2.GetID()){//Align against itself, disallow diagonal
		self = 1;
		int *s1 = new int[xLen];
		int *s2 = new int[xLen];
		for(i=0;i<xLen;i++)
			s1[i] = s2[i] = i;
		
		AlignedFragment sfrag(wseq1.GetID(),wseq1.GetID(),1,1,xLen,xLen,s1,s2);
		delete s1;delete s2;
		UpdateMap(&sfrag);
		
	}

	AlignedFragment *frag;
	SCORE_PAIR *score = hmm->ViterbiInitialize(wseq1,wseq2,map);
	ScoreMatrix * mp = score->first;
	Matrix *pTrace = score->second;
	ScoreMatrix& m = *mp;
	Matrix& trace = *pTrace;
	frag = hmm->OneAligment(wseq1, wseq2, trace, m);

	while(1){
		
		int len;
		if ((len = frag->GetLength()) < MINLENGTH) { //too short
			delete frag;
			break;
		}

		if(self){
			AVECT one_pair;
			frag->ProcessRepeat(one_pair,MINLENGTH);
			for(i=0;i< (int)one_pair.size();i++){
				AlignedFragment af = one_pair[i];
				pair_frags.push_back(af);
				UpdateMap(&af,1);
			}
			if(one_pair.size()==0)
				UpdateMap(frag,1);
			one_pair.clear();
		}
		else{
			if (verbose)
				frag->Print(stderr);
			pair_frags.push_back(*frag);
			UpdateMap(frag);
		}

		ConsistencyCheck(pair_frags);

		hmm->ViterbiUpdate(mp,pTrace,wseq1,wseq2,map,frag,MINLENGTH);
		delete frag;
		frag = hmm->OneAligment(wseq1, wseq2,trace,m);
	}
	delete mp;
	delete pTrace;
	delete score;
	for(i=0;i< (int)pair_frags.size();i++)
		fragments.push_back(pair_frags[i]);
	
}

////////////////////////////////////////////////////////////////////////////////
// Find all pairwise alignments between two sequences
////////////////////////////////////////////////////////////////////////////////

void PairAligner::PairAlign( AVECT &fragments)
{
	AVECT pair_frags;
	int i;
	int self = 0;

	Sequence &wseq1 = *seq1;
	Sequence &wseq2 = *seq2;

	if(wseq1.GetID() == wseq2.GetID()){//Align against itself, disallow diagonal
		self = 1;
		int *s1 = new int[xLen];
		int *s2 = new int[xLen];
		for(i=0;i<xLen;i++)
			s1[i] = s2[i] = i;
		
		AlignedFragment sfrag(wseq1.GetID(),wseq1.GetID(),1,1,xLen,xLen,s1,s2);
		delete s1;delete s2;
		UpdateMap(&sfrag);
		
	}

	while(1){
		AlignedFragment *frag;
		if(!enableViterbi){
			ScoreMatrix *p = hmm->Posterior (wseq1, wseq2 , map, NUM_STATES);
			Matrix *p2 = new Matrix(*p);
			delete p;

			frag = LocalAlign::ComputeLocalAlignment(wseq1, wseq2, *p2);
			delete p2;
		
		}
		else{
			frag = hmm->Viterbi(wseq1,wseq2,map);
		}
		int len;
		if ((len = frag->GetLength()) < MINLENGTH) { //too short
			delete frag;
			break;
		}

		if(self){
			AVECT one_pair;
			frag->ProcessRepeat(one_pair,MINLENGTH);
			for(i=0;i<(int)one_pair.size();i++){
				AlignedFragment af = one_pair[i];
				pair_frags.push_back(af);
				UpdateMap(&af,1);
			}
			if(one_pair.size()==0)
				UpdateMap(frag,1);
			one_pair.clear();
		}
		else{
			if (verbose)
				frag->Print(stderr);
			pair_frags.push_back(*frag);
			UpdateMap(frag);
		}
		delete frag;
		ConsistencyCheck(pair_frags);
	}
	for(i=0;i<(int)pair_frags.size();i++)
		fragments.push_back(pair_frags[i]);
}

