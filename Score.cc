//////////////////////////////////////////////////////////////////////
// Score.cc
//////////////////////////////////////////////////////////////////////

#include "Score.h"
#include "Assert.h"

float *EXP_SCORE_TO_FLOAT_TABLE;
SCORE *LOG_EXP_PLUS_1_TABLE;

//////////////////////////////////////////////////////////////////////
// Precompute math tables
//////////////////////////////////////////////////////////////////////

void PRECOMPUTE_SCORE_TABLES(){
  LOG_EXP_PLUS_1_TABLE = new SCORE[TABLE_SIZE];
  ASSERT (LOG_EXP_PLUS_1_TABLE, "Out of memory.");
  
  EXP_SCORE_TO_FLOAT_TABLE = new float[TABLE_SIZE];
  ASSERT (EXP_SCORE_TO_FLOAT_TABLE, "Out of memory.");
  
  for (int i = 0; i < TABLE_SIZE; i++){
    LOG_EXP_PLUS_1_TABLE[i] = (SCORE)(log(exp((double)(i / SCALE)) + 1) * SCALE);
    EXP_SCORE_TO_FLOAT_TABLE[i] = (float) exp((double)(-i / SCALE));
  }

}
