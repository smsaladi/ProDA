////////////////////////////////////////////////////////////////////////////
// PairAligner.h
//
// Find all pairwise alignments between two sequences
////////////////////////////////////////////////////////////////////////////

#ifndef PAIRALIGNER_H
#define PAIRALIGNER_H

#include "Sequence.h"
#include "ProbModel.h"
#include "AlignedFragment.h"

class PairAligner {
private:
	Sequence *seq1, *seq2;
	int *map;
	ProbModel *hmm;
	int xLen, yLen;

private:
	void UpdateMap(AlignedFragment *frag, int self = 0);
	void ConsistencyCheck(AVECT &pair_frags);
	
public:
	PairAligner(ProbModel *v_hmm, Sequence *s1, Sequence *s2);
	~PairAligner(){if(map) delete map;}
	void FastPairAlign(AVECT &fragments);
	void PairAlign( AVECT &fragments);
};

#endif
