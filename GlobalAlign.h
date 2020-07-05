//////////////////////////////////////////////////////////////////////
// GlobalAlign.h
//
// Global alignment procedures using maximum weight trace.
//////////////////////////////////////////////////////////////////////

#ifndef GLOBALALIGN_H
#define GLOBALALIGN_H

#include "Score.h"
#include "Matrix.h"
#include "MultiSequence.h"
#include "Sequence.h"
#include "ProbModel.h"

//////////////////////////////////////////////////////////////////////
// Global alignment class
//////////////////////////////////////////////////////////////////////

class GlobalAlign {
  
  enum TracebackType { NONE, UP, LEFT, UP_LEFT };

  // insert gaps into aligned sequence
  static Sequence *InsertGaps (const Sequence &seq, const char *alignmentPath, char ch);

 public:

  // maximum weight trace
  static char *ComputeMWTrace (const Matrix &m, float *score = NULL, int *length = NULL);
  
  // convert alignment path into MultiSequence
  static MultiSequence *BuildAlignment (const MultiSequence &group1, 
					const MultiSequence &group2, 
					const char *alignmentPath);

  // align two groups of sequences
  static MultiSequence *AlignGroups (int n,
				     SparseMatrix **posteriors, 
				     const MultiSequence &group1, 
				     const MultiSequence &group2);
  static float ComputeAlignmentScore (const MultiSequence &seqs, 
				      int n, 
				      SparseMatrix **posteriors);
    
};

#endif
