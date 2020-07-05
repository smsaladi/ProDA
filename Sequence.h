//////////////////////////////////////////////////////////////////////
// Sequence.h
//
// Class for manipulating single sequences.
//////////////////////////////////////////////////////////////////////

#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <map>
#include "AlignedFragment.h"


//////////////////////////////////////////////////////////////////////
// Sequence object 
//////////////////////////////////////////////////////////////////////

class Sequence {
  char *data;
  char *name;
  int length;
  int id;
  int *align; //number of sequences aligned to a given position
  int *position; //original position used for tracking after erasing fragments

 public:
	 void ClearAlignPosition();
	 int GetAlign(int i) const;
	 void SubStr(int begin, int end);
	 int SetID(int newid);
	 void Clip(int begin, int end);
	 void PrintAlign(FILE *f);
	 void EraseCluster(int start, int end, int n);
	 int OriginPosition(int current);
	 void EraseFragment(int begin, int end);
	 void AddAlignPosition(int begin, int end);

  // constructors
  Sequence (char *data, char *name, int length, int id);
  Sequence (const Sequence &rhs);

  // assignment operator
  const Sequence& operator= (const Sequence &rhs);

  // destructor
  ~Sequence ();
  
  // getters
  const char *GetData () const;
  const char *GetName () const;
  const int GetLength () const;
  const int GetID () const;

  // setters
  void SetData (char *data);

  // compute mapping from letter to positions in sequence
  int *ComputeMapping () const;
};

#endif
