//////////////////////////////////////////////////////////////////////
// MultiSequence.h
//   
// This file contains the routines needed for the creation, 
// maintenance, and use of a MultiSequence object which contains all 
// of the data associated with a set of sequences.  
//////////////////////////////////////////////////////////////////////

#ifndef MULTISEQUENCE_H
#define MULTISEQUENCE_H

#include "Sequence.h"
#include "Block.h"
#include <stdio.h>

//////////////////////////////////////////////////////////////////////
// MultiSequence object
//////////////////////////////////////////////////////////////////////
class MultiSequence {
  Sequence **sequences;
  int numSequences;
  
  // I/O helper routines
  const int AutoDetectFileFormat (const char *filename) const;
  void LoadMFA (const char *filename, bool compressGaps);
  void LoadPILEUP (const char *filename, bool compressGaps);
  void LoadData (const char *filename, bool compressGaps);
  const char ComputeAnnotation (const char *data, const int size) const;

 public:

  //Block operations
  void FindBlock(Block &block, int &start, int &end, int minlength = 20);
  void AddAlignPosition(Fragment * frag);
  void ClearAlignPosition();
	 
  // constructors
  MultiSequence();
  MultiSequence (const MultiSequence &rhs);

  // assignment operator
  const MultiSequence& operator= (const MultiSequence &rhs);
  
  // destructor
  ~MultiSequence();

  // getters
  const int GetNumSequences() const;
  const int GetLength() const;
  const Sequence &GetSequence (int index) const;
  Sequence * GetSequencePtr(int index);

  // add sequences
  void AddSequence (Sequence *seq);

  // sort sequences by id
  void Sort();

  // input
  void LoadSequences (const char *filename);
  void LoadAlignment (const char *filename);

  // output
  void WriteMFA (FILE *file) const;
  void WriteCLUSTALW (FILE *file) const;

  //Block output
  void WriteFASTA(FILE *file,Block *block, MultiSequence *result, int start, int end);
  void WriteCLUSTALW(FILE *file, int start, int end);

};


#endif
