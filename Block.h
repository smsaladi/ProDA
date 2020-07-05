//////////////////////////////////////////////////////////////////////
// Block.h
//
// Interface for block class
//////////////////////////////////////////////////////////////////////

#ifndef BLOCK_H
#define BLOCK_H
#include <vector>
#include "Types.h"

class MultiSequence;
class Block  
{
	std::vector<Fragment> frags;
public:
	int part;
	Fragment seed;

	//Remove fragments already present in blocks
	int AdjustAFragmentList(AVECT &fragments, Matrix *similarity=NULL, float threshold=0);
	int AdjustAFragmentList(AVECT &fragments, int numSeq, Matrix *similarity=NULL, float threshold=0);

	//Add a fragment to block
	void AddFragment(Fragment &fr);

	//Getters
	int GetLength();
	int size();

	//Printing utility
	void PrintBlock(FILE *f, MultiSequence *seqs, int compare = 0);

	//Block operator
	Fragment& operator [](int i);
	Block& operator=(const Block &bl);
	
	//Constructors
	Block();
	Block(std::vector<AlignedFragment>& afrags, MultiSequence *seqs, 
					bool enableTransitivity, Block *prohibited = NULL);
	virtual ~Block();
	
};

#endif 
