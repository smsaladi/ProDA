//////////////////////////////////////////////////////////////////////
// AlignedFragment.h 
//
// Class for a pairwise local alignment
//////////////////////////////////////////////////////////////////////

#ifndef ALIGNFRAGMENT_H
#define ALIGNFRAGMENT_H

#include <vector>

class Fragment;
class AlignedFragment;

typedef std::vector<AlignedFragment> AVECT;

class AlignedFragment  
{

public:
	void Print(FILE *file);
	int ProcessRepeat(AVECT &fragments, int minlength);
	
	//Getters
	int GetEnd(int i);
	int GetBegin(int i);
	int GetID(int i);
	std::pair<int, int> *GetAlignPos(int seq, int  pos, int &second);
	int GetLength();
	Fragment * GetFragment(int i);
	Fragment * GetAlignFragment(Fragment &fr);
	
	//
	void Adjust(Fragment& fr1, Fragment& fr2, AlignedFragment &rfr1, AlignedFragment &rfr2);
	void Prune();
	AlignedFragment * SubFragment(int begin0,int end0, int begin1, int end1);
	void ShiftRight(int offset1, int offset2);

	
	//Constructors
	AlignedFragment(int d1, int d2, int beg1, int beg2, int e1, int e2, int *s1, int *s2);
	AlignedFragment& operator =(const AlignedFragment af);
	AlignedFragment(const AlignedFragment& af);
	AlignedFragment();
	virtual ~AlignedFragment();

	int id[2];
	int begin[2];
	int end[2]; 
	float similarity;
	int *seq[2];

};

class Fragment{
public:
	Fragment();
	Fragment(int start, int e, int mul, int sid) 
		{ begin = start; end = e; multiply = mul; id = sid; length = end - begin + 1;}
	Fragment (const Fragment& fr)
		{ begin = fr.begin; end = fr.end; length = fr.length; multiply = fr.multiply; id = fr.id;}

	Fragment& operator =(const Fragment f)
		{begin = f.begin; length = f.length; multiply = f.multiply;end = f.end;id = f.id;return *this;}
	int Overlap(Fragment& fr);
	
	int begin;
	int length;
	int multiply; // number of sequences aligned to the fragment
	int id;
	int end;
};

#endif 
