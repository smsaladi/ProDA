// LinkTable.h: interface for the LinkTable class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_LINKTABLE_H__DE420983_21EB_4E38_A778_6F31F2460741__INCLUDED_)
#define AFX_LINKTABLE_H__DE420983_21EB_4E38_A778_6F31F2460741__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "Link.h"
#include <vector>
class MultiSequence;
class AlignedFragment;
class Block;
typedef std::vector<AlignedFragment> AVECT;
class LinkTable  
{
public:
	void PrintPlain(FILE *file,int m);
	void Print(FILE *file, int m);
	Block * OneBlock();
	int ExtendLeft(int& orgPos, int orgLen);
	void Print(FILE *file);
	int Entry(int sid, int off);
	LinkTable(MultiSequence *seqs, AVECT& fragments);
	Link **entries;
	int *seqBound;
	int seqNum;
	int entNum;
	LinkTable();
	LinkTable(int num, int cap);
	virtual ~LinkTable();

};

#endif // !defined(AFX_LINKTABLE_H__DE420983_21EB_4E38_A778_6F31F2460741__INCLUDED_)
