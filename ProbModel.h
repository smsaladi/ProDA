//////////////////////////////////////////////////////////////////////
// ProbModel.h
//
// Probabilistic model routines
//////////////////////////////////////////////////////////////////////

#ifndef PROBMODEL_H
#define PROBMODEL_H

#include "Assert.h"
#include "Types.h"
#include "AlignedFragment.h"
#include "Score.h"
#include "Sequence.h"

enum STATES {
  BEF_X,
  BEF_Y,
  MATCH,
  INS_X,
  INS_Y,
  AFT_X,
  AFT_Y,
  NUM_STATES
};

//////////////////////////////////////////////////////////////////////
// Probabilistic model object
//////////////////////////////////////////////////////////////////////

class ProbModel {

  double A;
  double D;
  double E;
  double T;

  int NUM_TRANS_X;
  int NUM_TRANS_Y;
  int NUM_TRANS_BOTH;

  int TRANSITIONS_EMIT_X[NUM_STATES * NUM_STATES][2];
  int TRANSITIONS_EMIT_Y[NUM_STATES * NUM_STATES][2];
  int TRANSITIONS_EMIT_BOTH[NUM_STATES * NUM_STATES][2];  
  
  SCORE LOG_START[NUM_STATES];
  SCORE LOG_FINAL[NUM_STATES];
  SCORE LOG_TRANS[NUM_STATES][NUM_STATES];
  SCORE LOG_EMIT_2[256][256];
  SCORE LOG_EMIT_1[256];

  // computing forward/backward recurrences
  ScoreMatrix *Forward (const Sequence &sx, const Sequence &sy) const;
  ScoreMatrix *Backward (const Sequence &sx, const Sequence &sy) const;

  // computing partition coefficient (total probability)
  SCORE ComputeTotalProb (const Sequence &sx, const Sequence &sy, const ScoreMatrix &forward, const ScoreMatrix &backward) const;

 public:
	 AlignedFragment * OneAligment(const Sequence &sx, const Sequence &sy, Matrix& trace, ScoreMatrix& m);
	 
	 ScoreMatrix * Backward (const Sequence &sx, const Sequence &sy, int *map) const;
	 ScoreMatrix * Forward (const Sequence &sx, const Sequence &sy, int *map) const;
	 

  // constructor
  ProbModel ();

  // posterior probability computation
  ScoreMatrix *Posterior (const Sequence &sx, const Sequence &sy, STATES state = NUM_STATES) const;
  ScoreMatrix * Posterior (const Sequence &sx, const Sequence &sy, int *map, STATES state = NUM_STATES) const;

  //Viterbi decoding 
  AlignedFragment * Viterbi(const Sequence &sx, const Sequence &sy, int *map);
  void ViterbiUpdate(ScoreMatrix *mp, Matrix *pTrace, const Sequence &sx, const Sequence &sy, 
		 int *map, AlignedFragment *frag, int minlength);
  SCORE_PAIR * ViterbiInitialize(const Sequence &sx, const Sequence &sy, int *map);

  // compute expected sufficient statistics
  Matrix *ComputeExpectedCounts (const Sequence &sx, const Sequence &sy) const;
  
  // compute new parameters
  void ComputeParams (const Matrix *cts);


};

#endif
