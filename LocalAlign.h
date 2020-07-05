//////////////////////////////////////////////////////////////////////
// LocalAlign.h
//
// Local alignment procedures using maximum weight trace.
//////////////////////////////////////////////////////////////////////


#ifndef LOCALALIGN_H
#define LOCALALIGN_H

#include "Score.h"
#include "Matrix.h"
#include "Sequence.h"
#include "ProbModel.h"
#include "AlignedFragment.h"

//////////////////////////////////////////////////////////////////////
// Local alignment class
//////////////////////////////////////////////////////////////////////

class LocalAlign {
  
  enum TracebackType { NONE, UP, LEFT, UP_LEFT };

 public:
	 static AlignedFragment * ComputeLocalAlignment (const Sequence &seq1, const Sequence &seq2, 
													const Matrix &m, float *score = NULL);
    
};

#endif
