//////////////////////////////////////////////////////////////////////
// ProbModel.cc
//////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <ctype.h>
#include "ProbModel.h"
#include "Utilities.h"

int EMISSIONS[NUM_STATES][2] = { 
  {1, 0},
  {0, 1},
  {1, 1},
  {1, 0},
  {0, 1},
  {1, 0},
  {0, 1}
};

const int START_STATES[NUM_STATES] = { 1, 1, 1, 0, 0, 0, 0 };
const int FINAL_STATES[NUM_STATES] = { 0, 0, 1, 0, 0, 1, 1 };

const int TRANSITIONS[NUM_STATES][NUM_STATES] = {
  { 1, 1, 1, 0, 0, 0, 0 },
  { 0, 1, 1, 0, 0, 0, 0 },
  { 0, 0, 1, 1, 1, 1, 1 },
  { 0, 0, 1, 1, 0, 0, 0 },
  { 0, 0, 1, 0, 1, 0, 0 },
  { 0, 0, 0, 0, 0, 1, 1 },
  { 0, 0, 0, 0, 0, 0, 1 }
};

char *ALPHABET = "ARNDCQEGHILKMFPSTWYV";


float INIT_A = 0.9860202074;
float INIT_D = 0.0207951729;
float INIT_E = 0.6397492290;
float INIT_T = 0.0078469915;

float PROB_SINGLE[20] = {
  0.07831005, 0.05246024, 0.04433257, 0.05130349, 0.02189704,
  0.03585766, 0.05615771, 0.07783433, 0.02601093, 0.06511648,
  0.09716489, 0.05877077, 0.02438117, 0.04463228, 0.03940142,
  0.05849916, 0.05115306, 0.01203523, 0.03124726, 0.07343426
};

float PROB_PAIR[20][20] = {
  {0.02373072, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  {0.00244502, 0.01775118, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  {0.00210228, 0.00207782, 0.01281864, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  {0.00223549, 0.00161657, 0.00353540, 0.01911178, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  {0.00145515, 0.00044701, 0.00042479, 0.00036798, 0.01013470, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  {0.00219102, 0.00253532, 0.00158223, 0.00176784, 0.00032102, 0.00756604, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  {0.00332218, 0.00268865, 0.00224738, 0.00496800, 0.00037956, 0.00345128, 0.01676565, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  {0.00597898, 0.00194865, 0.00288882, 0.00235249, 0.00071206, 0.00142432, 0.00214860, 0.04062876, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  {0.00114353, 0.00132105, 0.00141205, 0.00097077, 0.00026421, 0.00113901, 0.00131767, 0.00103704, 0.00867996, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  {0.00318853, 0.00138145, 0.00104273, 0.00105355, 0.00094040, 0.00100883, 0.00124207, 0.00142520, 0.00059716, 0.01778263, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  {0.00449576, 0.00246811, 0.00160275, 0.00161966, 0.00138494, 0.00180553, 0.00222063, 0.00212853, 0.00111754, 0.01071834, 0.03583921, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  {0.00331693, 0.00595650, 0.00257310, 0.00252518, 0.00046951, 0.00312308, 0.00428420, 0.00259311, 0.00121376, 0.00157852, 0.00259626, 0.01612228, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  {0.00148878, 0.00076734, 0.00063401, 0.00047808, 0.00037421, 0.00075546, 0.00076105, 0.00066504, 0.00042237, 0.00224097, 0.00461939, 0.00096120, 0.00409522, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  {0.00165004, 0.00090768, 0.00084658, 0.00069041, 0.00052274, 0.00059248, 0.00078814, 0.00115204, 0.00072545, 0.00279948, 0.00533369, 0.00087222, 0.00116111, 0.01661038, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  {0.00230618, 0.00106268, 0.00100282, 0.00125381, 0.00034766, 0.00090111, 0.00151550, 0.00155601, 0.00049078, 0.00103767, 0.00157310, 0.00154836, 0.00046718, 0.00060701, 0.01846071, 0.0, 0.0, 0.0, 0.0, 0.0},
  {0.00631752, 0.00224540, 0.00301397, 0.00285226, 0.00094867, 0.00191155, 0.00293898, 0.00381962, 0.00116422, 0.00173565, 0.00250962, 0.00312633, 0.00087787, 0.00119036, 0.00180037, 0.01346609, 0.0, 0.0, 0.0, 0.0},
  {0.00389995, 0.00186053, 0.00220144, 0.00180488, 0.00073798, 0.00154526, 0.00216760, 0.00214841, 0.00077747, 0.00248968, 0.00302273, 0.00250862, 0.00093371, 0.00107595, 0.00147982, 0.00487295, 0.01299436, 0.0, 0.0, 0.0},
  {0.00039119, 0.00029139, 0.00021006, 0.00016015, 0.00010666, 0.00020592, 0.00023815, 0.00038786, 0.00019097, 0.00039549, 0.00076736, 0.00028448, 0.00016253, 0.00085751, 0.00015674, 0.00026525, 0.00024961, 0.00563625, 0.0, 0.0},
  {0.00131840, 0.00099430, 0.00074960, 0.00066005, 0.00036626, 0.00070192, 0.00092548, 0.00089301, 0.00131038, 0.00127857, 0.00219713, 0.00100817, 0.00054105, 0.00368739, 0.00047608, 0.00102648, 0.00094759, 0.00069226, 0.00999315, 0.0},
  {0.00533241, 0.00169359, 0.00136609, 0.00127915, 0.00119152, 0.00132844, 0.00178697, 0.00194579, 0.00071553, 0.01117956, 0.00914460, 0.00210897, 0.00197461, 0.00256159, 0.00135781, 0.00241601, 0.00343452, 0.00038538, 0.00148001, 0.02075171}
};

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////

ProbModel::ProbModel () : A(INIT_A), D(INIT_D), E(INIT_E), T(INIT_T){

  // build transition model

  NUM_TRANS_X = NUM_TRANS_Y = NUM_TRANS_BOTH = 0;
  for (int i = 0; i < NUM_STATES; i++){
    for (int j = 0; j < NUM_STATES; j++){
      if (TRANSITIONS[i][j]){
	if (EMISSIONS[j][0] == 1 && EMISSIONS[j][1] == 0){
	  TRANSITIONS_EMIT_X[NUM_TRANS_X][0] = i;
	  TRANSITIONS_EMIT_X[NUM_TRANS_X][1] = j;
	  NUM_TRANS_X++;
	} else if (EMISSIONS[j][0] == 0 && EMISSIONS[j][1] == 1){
	  TRANSITIONS_EMIT_Y[NUM_TRANS_Y][0] = i;
	  TRANSITIONS_EMIT_Y[NUM_TRANS_Y][1] = j;
	  NUM_TRANS_Y++;
	} else if (EMISSIONS[j][0] == 1 && EMISSIONS[j][1] == 1){

	  TRANSITIONS_EMIT_BOTH[NUM_TRANS_BOTH][0] = i;
	  TRANSITIONS_EMIT_BOTH[NUM_TRANS_BOTH][1] = j;
	  NUM_TRANS_BOTH++;
	}
      }
    }
  }

  // use BLOSUM matrix values

  {for (int i = 0; i < 256; i++){
    LOG_EMIT_1[i] = TO_SCORE(LOG_FLOAT(1e-10));
    for (int j = 0; j < 256; j++)
      LOG_EMIT_2[i][j] = TO_SCORE(LOG_FLOAT(1e-10));
  }}

  {for (int i = 0; i < 20; i++){
    LOG_EMIT_1[(unsigned char) ALPHABET[i]] = 
      LOG_EMIT_1[(unsigned char) tolower(ALPHABET[i])] = 
      TO_SCORE(LOG_FLOAT(PROB_SINGLE[i]));
  }}
  
  {for (int i = 0; i < 20; i++)
    for (int j = 0; j < 20; j++)
      LOG_EMIT_2[(unsigned char) ALPHABET[i]][(unsigned char) ALPHABET[j]] =
	LOG_EMIT_2[(unsigned char) tolower(ALPHABET[i])][(unsigned char) ALPHABET[j]] =
	LOG_EMIT_2[(unsigned char) ALPHABET[i]][(unsigned char) tolower(ALPHABET[j])] =
	LOG_EMIT_2[(unsigned char) tolower(ALPHABET[i])][(unsigned char) tolower(ALPHABET[j])] =
	TO_SCORE(LOG_FLOAT(PROB_PAIR[max(i,j)][min(i,j)]));}
  
  // fill in transition probabilities

  ComputeParams(NULL);
}

//////////////////////////////////////////////////////////////////////
// Compute forward probabilities
//////////////////////////////////////////////////////////////////////

ScoreMatrix *ProbModel::Forward (const Sequence &sx, const Sequence &sy) const {
  int xLen = sx.GetLength();
  int yLen = sy.GetLength();
  const char *x = sx.GetData();
  const char *y = sy.GetData();
  
  ScoreMatrix *mp = new ScoreMatrix (NUM_STATES, xLen+1, yLen+1);
  ASSERT (mp, "Out of memory.");
  ScoreMatrix &m = *mp;
  m.Fill (LOG_ZERO_SCORE);
  
  // initialization condition
  for (int s = 0; s < NUM_STATES; s++) if (START_STATES[s]){
    if (EMISSIONS[s][0] == 1 && EMISSIONS[s][1] == 0){
      m(s, EMISSIONS[s][0], EMISSIONS[s][1]) = LOG_START[s] + LOG_EMIT_1[(unsigned char) x[1]];
    } else if (EMISSIONS[s][0] == 0 && EMISSIONS[s][1] == 1){
      m(s, EMISSIONS[s][0], EMISSIONS[s][1]) = LOG_START[s] + LOG_EMIT_1[(unsigned char) y[1]];
    } else if (EMISSIONS[s][0] == 1 && EMISSIONS[s][1] == 1){
      m(s, EMISSIONS[s][0], EMISSIONS[s][1]) = LOG_START[s] + LOG_EMIT_2[(unsigned char) x[1]][(unsigned char) y[1]];
    }
  }

  SCORE *ij = m.GetPtr(0,0,0);
  SCORE *i1j = m.GetPtr(0,-1,0);
  SCORE *ij1 = m.GetPtr(0,0,-1);
  SCORE *i1j1 = m.GetPtr(0,-1,-1);
  
  // recursion
  for (int i = 0; i <= xLen; i++){ 
    const unsigned char xi = (unsigned char) (i == 0 ? '~' : x[i]);
    for (int j = 0; j <= yLen; j++){ 
      const unsigned char yj = (unsigned char) (j == 0 ? '~' : y[j]);

      if (i > 0){
	for (int k = 0; k < NUM_TRANS_X; k++){
	  int s = TRANSITIONS_EMIT_X[k][0];
	  int t = TRANSITIONS_EMIT_X[k][1];
	  ij[t] = LOG_ADD_SCORE (ij[t], i1j[s] + LOG_TRANS[s][t] + LOG_EMIT_1[xi]);
	}
      }

      if (j > 0){
	for (int k = 0; k < NUM_TRANS_Y; k++){
	  int s = TRANSITIONS_EMIT_Y[k][0];
	  int t = TRANSITIONS_EMIT_Y[k][1];
	  ij[t] = LOG_ADD_SCORE (ij[t], ij1[s] + LOG_TRANS[s][t] + LOG_EMIT_1[yj]);
	}
      }
      
      if (i > 0 && j > 0){
	for (int k = 0; k < NUM_TRANS_BOTH; k++){
	  int s = TRANSITIONS_EMIT_BOTH[k][0];
	  int t = TRANSITIONS_EMIT_BOTH[k][1];
	  ij[t] = LOG_ADD_SCORE (ij[t], i1j1[s] + LOG_TRANS[s][t] + LOG_EMIT_2[xi][yj]);
	}
      }
      
      ij += NUM_STATES;
      i1j += NUM_STATES;
      ij1 += NUM_STATES;
      i1j1 += NUM_STATES;
    }
  }

  return mp;
}

//////////////////////////////////////////////////////////////////////
// Compute backward probabilities
//////////////////////////////////////////////////////////////////////

ScoreMatrix *ProbModel::Backward (const Sequence &sx, const Sequence &sy) const {
  int xLen = sx.GetLength();
  int yLen = sy.GetLength();
  const char *x = sx.GetData();
  const char *y = sy.GetData();
  
  ScoreMatrix *mp = new ScoreMatrix (NUM_STATES, xLen+1, yLen+1);
  ASSERT (mp, "Out of memory.");
  ScoreMatrix &m = *mp;
  m.Fill (LOG_ZERO_SCORE);
  
  // initialization condition
  for (int s = 0; s < NUM_STATES; s++) if (FINAL_STATES[s]){
    m(s, xLen, yLen) = LOG_FINAL[s];
  }
  
  SCORE *ij = m.GetPtr(0,xLen,yLen);
  SCORE *i1j = m.GetPtr(0,xLen+1,yLen);
  SCORE *ij1 = m.GetPtr(0,xLen,yLen+1);
  SCORE *i1j1 = m.GetPtr(0,xLen+1,yLen+1);

  // recursion
  for (int i = xLen; i >= 0; i--){ 
    const unsigned char xi1 = (unsigned char) (i == xLen ? '~' : x[i+1]);
    for (int j = yLen; j >= 0; j--){ 
      const unsigned char yj1 = (unsigned char) (j == yLen ? '~' : y[j+1]);

      if (i < xLen){
	for (int k = 0; k < NUM_TRANS_X; k++){
	  int s = TRANSITIONS_EMIT_X[k][0];
	  int t = TRANSITIONS_EMIT_X[k][1];
	  ij[s] = LOG_ADD_SCORE (ij[s], i1j[t] + LOG_TRANS[s][t] + LOG_EMIT_1[xi1]);
	}
      }

      if (j < yLen){
	for (int k = 0; k < NUM_TRANS_Y; k++){
	  int s = TRANSITIONS_EMIT_Y[k][0];
	  int t = TRANSITIONS_EMIT_Y[k][1];
	  ij[s] = LOG_ADD_SCORE (ij[s], ij1[t] + LOG_TRANS[s][t] + LOG_EMIT_1[yj1]);
	}
      }
      
      if (i < xLen && j < yLen){
	for (int k = 0; k < NUM_TRANS_BOTH; k++){
	  int s = TRANSITIONS_EMIT_BOTH[k][0];
	  int t = TRANSITIONS_EMIT_BOTH[k][1];
	  ij[s] = LOG_ADD_SCORE (ij[s], i1j1[t] + LOG_TRANS[s][t] + LOG_EMIT_2[xi1][yj1]);
	}
      }

      ij -= NUM_STATES;
      i1j -= NUM_STATES;
      ij1 -= NUM_STATES;
      i1j1 -= NUM_STATES;
    }
  }
  
  return mp;
}

//////////////////////////////////////////////////////////////////////
// Compute partition coefficient (total probability)
//////////////////////////////////////////////////////////////////////

SCORE ProbModel::ComputeTotalProb (const Sequence &sx, const Sequence &sy, 
				   const ScoreMatrix &forward, const ScoreMatrix &backward) const {
  int xLen = sx.GetLength();
  int yLen = sy.GetLength();

  SCORE fProb = LOG_ZERO_SCORE;
  SCORE bProb = LOG_ZERO_SCORE;

  for (int s = 0; s < NUM_STATES; s++){
    if (START_STATES[s])
      bProb = LOG_ADD_SCORE (bProb,
			     forward(s,EMISSIONS[s][0],EMISSIONS[s][1]) + 
			     backward(s,EMISSIONS[s][0],EMISSIONS[s][1]));
    if (FINAL_STATES[s])
      fProb = LOG_ADD_SCORE (fProb,
			     forward(s,xLen,yLen) + backward(s,xLen,yLen));
  }
  
  return (fProb + bProb)/2;
}

//////////////////////////////////////////////////////////////////////
// Compute posterior probability matrix
//////////////////////////////////////////////////////////////////////

ScoreMatrix *ProbModel::Posterior (const Sequence &sx, const Sequence &sy, STATES state) const {
  int xLen = sx.GetLength();
  int yLen = sy.GetLength();

  // compute forward and backward probs
  ScoreMatrix *fp = Forward (sx, sy);
  ScoreMatrix *bp = Backward (sx, sy);
  
  ScoreMatrix &f = *fp;
  ScoreMatrix &b = *bp;
  
  SCORE totalProb = ComputeTotalProb (sx, sy, f, b);

  // compute posterior matrix

  if (state != NUM_STATES){
    ScoreMatrix *mp = new ScoreMatrix (1,xLen+1,yLen+1);
    ScoreMatrix &m = *mp;

    for (int i = 0; i <= xLen; i++)
      for (int j = 0; j <= yLen; j++)
	m(0,i,j) = f(state,i,j) + b(state,i,j) - totalProb;

    delete fp;
    delete bp;
    
    return mp;
  } 
  
  for (int i = 0; i <= xLen; i++){ 
    for (int j = 0; j <= yLen; j++){ 
      for (int s = 0; s < NUM_STATES; s++){
	f(s,i,j) = f(s,i,j) + b(s,i,j) - totalProb;
      }
    }
  }

  delete bp;  
  return fp;
}

//////////////////////////////////////////////////////////////////////
// Compute expected sufficient statistics
//////////////////////////////////////////////////////////////////////

Matrix *ProbModel::ComputeExpectedCounts (const Sequence &sx, const Sequence &sy) const {
  int xLen = sx.GetLength();
  int yLen = sy.GetLength();
  const char *x = sx.GetData();
  const char *y = sy.GetData();

  ScoreMatrix *fp = Forward (sx, sy);
  ScoreMatrix *bp = Backward (sx, sy);
  ScoreMatrix &f = *fp;
  ScoreMatrix &b = *bp;

  SCORE totalProb = ComputeTotalProb (sx, sy, f, b);
  
  Matrix *expCountsPtr = new Matrix (1, NUM_STATES+1, NUM_STATES+1);
  Matrix &e = *expCountsPtr;
  e.Fill (ZERO_FLOAT);

  // initialization condition
  for (int s = 0; s < NUM_STATES; s++){
    if (START_STATES[s])
      e(0,NUM_STATES,s) += EXP_SCORE_TO_FLOAT(f(s,EMISSIONS[s][0],EMISSIONS[s][1]) + 
					      b(s,EMISSIONS[s][0],EMISSIONS[s][1]) - totalProb);
    if (FINAL_STATES[s])
      e(0,s,NUM_STATES) += EXP_SCORE_TO_FLOAT(f(s,xLen,yLen) + 
					      b(s,xLen,yLen) - totalProb);
  }  
  
  SCORE *fi1j = f.GetPtr(0,-1,0);
  SCORE *fij1 = f.GetPtr(0,0,-1);
  SCORE *fi1j1 = f.GetPtr(0,-1,-1);
  SCORE *bij = b.GetPtr(0,0,0);
  float *ep = e.GetPtr(0,0,0);

  // recursion
  for (int i = 0; i <= xLen; i++){ 
    const unsigned char xi = (unsigned char) (i == 0 ? '~' : x[i]);

    for (int j = 0; j <= yLen; j++){ 
      const unsigned char yj = (unsigned char) (j == 0 ? '~' : y[j]);
      
      if (i > 0){
	const SCORE base = LOG_EMIT_1[xi] - totalProb;
	for (int k = 0; k < NUM_TRANS_X; k++){
	  int s = TRANSITIONS_EMIT_X[k][0];
	  int t = TRANSITIONS_EMIT_X[k][1];
	  ep[s*(NUM_STATES+1)+t] += EXP_SCORE_TO_FLOAT(fi1j[s] + LOG_TRANS[s][t] + bij[t] + base);
	}
      }

      if (j > 0){
	const SCORE base = LOG_EMIT_1[yj] - totalProb;
	for (int k = 0; k < NUM_TRANS_Y; k++){
	  int s = TRANSITIONS_EMIT_Y[k][0];
	  int t = TRANSITIONS_EMIT_Y[k][1];
	  ep[s*(NUM_STATES+1)+t] += EXP_SCORE_TO_FLOAT(fij1[s] + LOG_TRANS[s][t] + bij[t] + base);
	}
      }
      
      if (i > 0 && j > 0){
	const SCORE base = LOG_EMIT_2[xi][yj] - totalProb;
	for (int k = 0; k < NUM_TRANS_BOTH; k++){
	  int s = TRANSITIONS_EMIT_BOTH[k][0];
	  int t = TRANSITIONS_EMIT_BOTH[k][1];
	  ep[s*(NUM_STATES+1)+t] += EXP_SCORE_TO_FLOAT(fi1j1[s] + LOG_TRANS[s][t] + bij[t] + base);
	}
      }

      fi1j += NUM_STATES;
      fij1 += NUM_STATES;
      fi1j1 += NUM_STATES;
      bij += NUM_STATES;
    }
  }

  delete fp;
  delete bp;
  
  return expCountsPtr;
}

void ProbModel::ComputeParams (const Matrix *ctsPtr){
  
  if (ctsPtr){
    const Matrix &cts = *ctsPtr;
    
    ASSERT (cts.GetNumLayers() == 1, "Invalid sufficient statistics matrix.");
    ASSERT (cts.GetNumRows() == NUM_STATES+1, "Invalid sufficient statistics matrix.");
    ASSERT (cts.GetNumCols() == NUM_STATES+1, "Invalid sufficient statistics matrix.");

    float numA =
      cts(0,NUM_STATES,BEF_X) +
      cts(0,NUM_STATES,BEF_Y) +
      cts(0,BEF_X,BEF_X) +
      cts(0,BEF_X,BEF_Y) +
      cts(0,BEF_Y,BEF_Y) +
      cts(0,MATCH,AFT_X) +
      cts(0,MATCH,AFT_Y) +
      cts(0,AFT_X,AFT_X) +
      cts(0,AFT_X,AFT_Y) +
      cts(0,AFT_Y,AFT_Y);
    
    float num1minusA =
      cts(0,NUM_STATES,BEF_Y) + 
      2*cts(0,NUM_STATES,MATCH) + 
      cts(0,BEF_X,BEF_Y) + 
      2*cts(0,BEF_X,MATCH) + 
      cts(0,BEF_Y,MATCH) + 
      cts(0,MATCH,AFT_Y) + 
      2*cts(0,MATCH,NUM_STATES) + 
      cts(0,AFT_X,AFT_Y) + 
      2*cts(0,AFT_X,NUM_STATES) + 
      cts(0,AFT_Y,NUM_STATES);
    
    float numD =
      cts(0,MATCH,INS_X) + 
      cts(0,MATCH,INS_Y);
    
    float num1minus2D = 
      cts(0,MATCH,MATCH);
    
    float numE =
      cts(0,INS_X,INS_X) +
      cts(0,INS_Y,INS_Y);
    
    float num1minusE = 
      cts(0,INS_X,MATCH) +
      cts(0,INS_Y,MATCH);

    float numT =
      cts(0,MATCH,NUM_STATES) + 
      cts(0,MATCH,AFT_X) +
      cts(0,MATCH,AFT_Y);

    float num1minusT =
      cts(0,MATCH,INS_X) + 
      cts(0,MATCH,INS_Y) +
      cts(0,MATCH,MATCH);

    A = numA / (numA + num1minusA);
    D = 0.5 * numD / (numD + num1minus2D);
    E = numE / (numE + num1minusE);
    T = numT / (numT + num1minusT);
  
    fprintf (stderr, "A = %.10f\n", A);
    fprintf (stderr, "D = %.10f\n", D);
    fprintf (stderr, "E = %.10f\n", E);
    fprintf (stderr, "T = %.10f\n", T);
  }
  
  LOG_START[BEF_X] = TO_SCORE(LOG_FLOAT(A));
  LOG_START[BEF_Y] = TO_SCORE(LOG_FLOAT(A*(1-A)));
  LOG_START[MATCH] = TO_SCORE(LOG_FLOAT((1-A)*(1-A)));

  LOG_TRANS[BEF_X][BEF_X] = TO_SCORE(LOG_FLOAT(A));
  LOG_TRANS[BEF_X][BEF_Y] = TO_SCORE(LOG_FLOAT(A*(1-A)));
  LOG_TRANS[BEF_X][MATCH] = TO_SCORE(LOG_FLOAT((1-A)*(1-A)));
  LOG_TRANS[BEF_Y][BEF_Y] = TO_SCORE(LOG_FLOAT(A));
  LOG_TRANS[BEF_Y][MATCH] = TO_SCORE(LOG_FLOAT(1-A));

  LOG_TRANS[MATCH][MATCH] = TO_SCORE(LOG_FLOAT((1-2*D)*(1-T)));
  LOG_TRANS[MATCH][INS_X] = TO_SCORE(LOG_FLOAT(D*(1-T)));
  LOG_TRANS[MATCH][INS_Y] = TO_SCORE(LOG_FLOAT(D*(1-T)));
  LOG_TRANS[INS_X][MATCH] = TO_SCORE(LOG_FLOAT(1-E));
  LOG_TRANS[INS_Y][MATCH] = TO_SCORE(LOG_FLOAT(1-E));
  LOG_TRANS[INS_X][INS_X] = TO_SCORE(LOG_FLOAT(E));
  LOG_TRANS[INS_Y][INS_Y] = TO_SCORE(LOG_FLOAT(E));
  LOG_TRANS[MATCH][AFT_X] = TO_SCORE(LOG_FLOAT(T*A));
  LOG_TRANS[MATCH][AFT_Y] = TO_SCORE(LOG_FLOAT(T*A*(1-A)));
  
  LOG_TRANS[AFT_X][AFT_X] = TO_SCORE(LOG_FLOAT(A));
  LOG_TRANS[AFT_X][AFT_Y] = TO_SCORE(LOG_FLOAT(A*(1-A)));
  LOG_TRANS[AFT_Y][AFT_Y] = TO_SCORE(LOG_FLOAT(A));

  LOG_FINAL[MATCH] = TO_SCORE(LOG_FLOAT(T*(1-A)*(1-A)));
  LOG_FINAL[AFT_X] = TO_SCORE(LOG_FLOAT((1-A)*(1-A)));
  LOG_FINAL[AFT_Y] = TO_SCORE(LOG_FLOAT(1-A));
 
}

/////////////////////////////////////////////////////////////////////////////////////
// Forward probability with prohibited path given by map
/////////////////////////////////////////////////////////////////////////////////////

ScoreMatrix * ProbModel::Forward(const Sequence &sx, const Sequence &sy, int *map) const
{
  int xLen = sx.GetLength();
  int yLen = sy.GetLength();
  const char *x = sx.GetData();
  const char *y = sy.GetData();
  
  ScoreMatrix *mp = new ScoreMatrix (NUM_STATES, xLen+1, yLen+1);
  ASSERT (mp, "Out of memory.");
  ScoreMatrix &m = *mp;
  m.Fill (LOG_ZERO_SCORE);
  
  // initialization condition
  for (int s = 0; s < NUM_STATES; s++) if (START_STATES[s]){
    if (EMISSIONS[s][0] == 1 && EMISSIONS[s][1] == 0){
      m(s, EMISSIONS[s][0], EMISSIONS[s][1]) = LOG_START[s] + LOG_EMIT_1[(unsigned char) x[1]];
    } else if (EMISSIONS[s][0] == 0 && EMISSIONS[s][1] == 1){
      m(s, EMISSIONS[s][0], EMISSIONS[s][1]) = LOG_START[s] + LOG_EMIT_1[(unsigned char) y[1]];
    } else if (EMISSIONS[s][0] == 1 && EMISSIONS[s][1] == 1){
      m(s, EMISSIONS[s][0], EMISSIONS[s][1]) = LOG_START[s] + LOG_EMIT_2[(unsigned char) x[1]][(unsigned char) y[1]];
    }
  }

  SCORE *ij = m.GetPtr(0,0,0);
  SCORE *i1j = m.GetPtr(0,-1,0);
  SCORE *ij1 = m.GetPtr(0,0,-1);
  SCORE *i1j1 = m.GetPtr(0,-1,-1);
  
  // recursion
  for (int i = 0; i <= xLen; i++){ 
    const unsigned char xi = (unsigned char) (i == 0 ? '~' : x[i]);
    for (int j = 0; j <= yLen; j++){ 
      const unsigned char yj = (unsigned char) (j == 0 ? '~' : y[j]);

      if (i > 0){
	for (int k = 0; k < NUM_TRANS_X; k++){
	  int s = TRANSITIONS_EMIT_X[k][0];
	  int t = TRANSITIONS_EMIT_X[k][1];
	  ij[t] = LOG_ADD_SCORE (ij[t], i1j[s] + LOG_TRANS[s][t] + LOG_EMIT_1[xi]);
	}
      }

      if (j > 0){
	for (int k = 0; k < NUM_TRANS_Y; k++){
	  int s = TRANSITIONS_EMIT_Y[k][0];
	  int t = TRANSITIONS_EMIT_Y[k][1];
	  ij[t] = LOG_ADD_SCORE (ij[t], ij1[s] + LOG_TRANS[s][t] + LOG_EMIT_1[yj]);
	}
      }
      
      if (i > 0 && j > 0){
	for (int k = 0; k < NUM_TRANS_BOTH; k++){
	  int s = TRANSITIONS_EMIT_BOTH[k][0];
	  int t = TRANSITIONS_EMIT_BOTH[k][1];
	  ij[t] = LOG_ADD_SCORE (ij[t], i1j1[s] + LOG_TRANS[s][t] + LOG_EMIT_2[xi][yj]);
	}
	if(map[i*(yLen+1)+j] == 1)
		ij[MATCH] = LOG_ZERO_SCORE;
	  }
      
      ij += NUM_STATES;
      i1j += NUM_STATES;
      ij1 += NUM_STATES;
      i1j1 += NUM_STATES;
    }
  }

  return mp;

}

////////////////////////////////////////////////////////////////////////////////////
// Backward probability with prohibited path given by map
/////////////////////////////////////////////////////////////////////////////////////

ScoreMatrix * ProbModel::Backward(const Sequence &sx, const Sequence &sy, int *map) const
{
  int xLen = sx.GetLength();
  int yLen = sy.GetLength();
  const char *x = sx.GetData();
  const char *y = sy.GetData();
  
  ScoreMatrix *mp = new ScoreMatrix (NUM_STATES, xLen+1, yLen+1);
  ASSERT (mp, "Out of memory.");
  ScoreMatrix &m = *mp;
  m.Fill (LOG_ZERO_SCORE);
  
  // initialization condition
  for (int s = 0; s < NUM_STATES; s++) if (FINAL_STATES[s]){
    m(s, xLen, yLen) = LOG_FINAL[s];
  }
  
  SCORE *ij = m.GetPtr(0,xLen,yLen);
  SCORE *i1j = m.GetPtr(0,xLen+1,yLen);
  SCORE *ij1 = m.GetPtr(0,xLen,yLen+1);
  SCORE *i1j1 = m.GetPtr(0,xLen+1,yLen+1);

  // recursion
  for (int i = xLen; i >= 0; i--){ 
    const unsigned char xi1 = (unsigned char) (i == xLen ? '~' : x[i+1]);
    for (int j = yLen; j >= 0; j--){ 
      const unsigned char yj1 = (unsigned char) (j == yLen ? '~' : y[j+1]);

      if (i < xLen){
	for (int k = 0; k < NUM_TRANS_X; k++){
	  int s = TRANSITIONS_EMIT_X[k][0];
	  int t = TRANSITIONS_EMIT_X[k][1];
	  ij[s] = LOG_ADD_SCORE (ij[s], i1j[t] + LOG_TRANS[s][t] + LOG_EMIT_1[xi1]);
	}
      }

      if (j < yLen){
	for (int k = 0; k < NUM_TRANS_Y; k++){
	  int s = TRANSITIONS_EMIT_Y[k][0];
	  int t = TRANSITIONS_EMIT_Y[k][1];
	  ij[s] = LOG_ADD_SCORE (ij[s], ij1[t] + LOG_TRANS[s][t] + LOG_EMIT_1[yj1]);
	}
      }
      
      if (i < xLen && j < yLen){
	for (int k = 0; k < NUM_TRANS_BOTH; k++){
	  int s = TRANSITIONS_EMIT_BOTH[k][0];
	  int t = TRANSITIONS_EMIT_BOTH[k][1];
	  ij[s] = LOG_ADD_SCORE (ij[s], i1j1[t] + LOG_TRANS[s][t] + LOG_EMIT_2[xi1][yj1]);
	}
	
		if(map[i*(yLen+1)+j] == 1)
			ij[MATCH] = LOG_ZERO_SCORE;
      }

      ij -= NUM_STATES;
      i1j -= NUM_STATES;
      ij1 -= NUM_STATES;
      i1j1 -= NUM_STATES;
    }
  }
  
  return mp;

}

ScoreMatrix * ProbModel::Posterior(const Sequence &sx, const Sequence &sy, int *map, STATES state) const
{
  int xLen = sx.GetLength();
  int yLen = sy.GetLength();

  // compute forward and backward probs
  ScoreMatrix *fp = Forward (sx, sy, map);
  ScoreMatrix *bp = Backward (sx, sy, map);
  
  ScoreMatrix &f = *fp;
  ScoreMatrix &b = *bp;
  
  SCORE totalProb = ComputeTotalProb (sx, sy, f, b);

  // compute posterior matrix

  if (state != NUM_STATES){
    ScoreMatrix *mp = new ScoreMatrix (1,xLen+1,yLen+1);
    ScoreMatrix &m = *mp;

    for (int i = 0; i <= xLen; i++)
      for (int j = 0; j <= yLen; j++)
	m(0,i,j) = f(state,i,j) + b(state,i,j) - totalProb;

    delete fp;
    delete bp;
    
    return mp;
  } 
  
  for (int i = 0; i <= xLen; i++){ 
    for (int j = 0; j <= yLen; j++){ 
      for (int s = 0; s < NUM_STATES; s++){
	f(s,i,j) = f(s,i,j) + b(s,i,j) - totalProb;
      }
    }
  }

  delete bp;  
  return fp;

}


void ConsistencyCheck(AVECT &pair_frags);
void UpdateMap(int *map, int xLen, int yLen, AlignedFragment *frag, int self = 0);

////////////////////////////////////////////////////////////////////////////////////
// Perform Vterbi decoding
// Return one local parwise alignment
////////////////////////////////////////////////////////////////////////////////////

AlignedFragment * ProbModel::Viterbi(const Sequence &sx, const Sequence &sy, int *map)
{
  int xLen = sx.GetLength();
  int yLen = sy.GetLength();
  const char *x = sx.GetData();
  const char *y = sy.GetData();
  SCORE tmp;
  
  ScoreMatrix *mp = new ScoreMatrix (NUM_STATES, xLen+1, yLen+1);
  ASSERT (mp, "Out of memory.");
  ScoreMatrix &m = *mp;
  m.Fill (LOG_ZERO_SCORE);
  Matrix *pTrace = new Matrix(NUM_STATES,xLen+1,yLen+1);
  Matrix &trace = *pTrace;
  trace.Fill(-1);
  
  // initialization condition
  for (int s = 0; s < NUM_STATES; s++) if (START_STATES[s]){
    if (EMISSIONS[s][0] == 1 && EMISSIONS[s][1] == 0){
      m(s, EMISSIONS[s][0], EMISSIONS[s][1]) = LOG_START[s] + LOG_EMIT_1[(unsigned char) x[1]];
    } else if (EMISSIONS[s][0] == 0 && EMISSIONS[s][1] == 1){
      m(s, EMISSIONS[s][0], EMISSIONS[s][1]) = LOG_START[s] + LOG_EMIT_1[(unsigned char) y[1]];
    } else if (EMISSIONS[s][0] == 1 && EMISSIONS[s][1] == 1){
      m(s, EMISSIONS[s][0], EMISSIONS[s][1]) = LOG_START[s] + LOG_EMIT_2[(unsigned char) x[1]][(unsigned char) y[1]];
    }
  }

  SCORE *ij = m.GetPtr(0,0,0);
  SCORE *i1j = m.GetPtr(0,-1,0);
  SCORE *ij1 = m.GetPtr(0,0,-1);
  SCORE *i1j1 = m.GetPtr(0,-1,-1);
  
  // recursion
  int i , j;
  for (i = 0; i <= xLen; i++){ 
    const unsigned char xi = (unsigned char) (i == 0 ? '~' : x[i]);
    for (j = 0; j <= yLen; j++){ 
      const unsigned char yj = (unsigned char) (j == 0 ? '~' : y[j]);

      if (i > 0){
	for (int k = 0; k < NUM_TRANS_X; k++){
	  int s = TRANSITIONS_EMIT_X[k][0];
	  int t = TRANSITIONS_EMIT_X[k][1];
	  tmp = i1j[s] + LOG_TRANS[s][t] + LOG_EMIT_1[xi];
	  if (tmp > ij[t]){
		  ij[t] = tmp;
		  trace(t,i,j) = s;
	  }
	}
      }

      if (j > 0){
	for (int k = 0; k < NUM_TRANS_Y; k++){
	  int s = TRANSITIONS_EMIT_Y[k][0];
	  int t = TRANSITIONS_EMIT_Y[k][1];
	  tmp = ij1[s] + LOG_TRANS[s][t] + LOG_EMIT_1[yj];
	  if (tmp > ij[t]){
		  ij[t] = tmp;
		  trace(t,i,j) = s;
	  }
	}
      }
      
      if (i > 0 && j > 0){
	for (int k = 0; k < NUM_TRANS_BOTH; k++){
	  int s = TRANSITIONS_EMIT_BOTH[k][0];
	  int t = TRANSITIONS_EMIT_BOTH[k][1];
	  tmp = i1j1[s] + LOG_TRANS[s][t] + LOG_EMIT_2[xi][yj];
	  if (tmp > ij[t]){
		  ij[t] = tmp;
		  trace(t,i,j) = s;
	  }
	}
	if(map[i*(yLen+1)+j] == 1){
		ij[MATCH] = LOG_ZERO_SCORE;
		trace(MATCH,i,j) = -1;
	}
	  }
      
      ij += NUM_STATES;
      i1j += NUM_STATES;
      ij1 += NUM_STATES;
      i1j1 += NUM_STATES;
    }
  }
  float bestProb = LOG_ZERO_SCORE;
  int state = -1;
  for(i = 0; i < NUM_STATES; i++){
	  if(FINAL_STATES[i]){
		  float crProb = m(i,xLen,yLen) + LOG_FINAL[i];
		  if(crProb > bestProb){
			  bestProb = crProb;
			  state = i;
		  }
	  }
  }
  delete mp;
  	char *buffer = new char[(xLen+1) * (yLen+1)];
	ASSERT (buffer, "Out of memory.");

	int r = xLen, c = yLen, len = 0;
	int newState;

	while((r != 0 || c != 0) && state != MATCH){
		newState = (int)trace(state,r,c);
		switch(state){
			case AFT_X:r--;break;
			case AFT_Y:c--;break;
			default: ASSERT(false,"Insert state within flanking region");
		}
		state = newState;
	}
	int bestx = r;
	int besty = c;
	while((r != 0 || c != 0) && (state != BEF_X && state !=BEF_Y)){
		newState = (int)trace(state,r,c);
		switch(state){
			case MATCH: c--;r--;buffer[len++] = 'B';break;
			case INS_X: r--;buffer[len++] = 'X';break;
			case INS_Y: c--;buffer[len++] = 'Y';break;
			default: ASSERT(false,"Flanking state inside");
		}
		state = newState;
	}
	  char *ret = new char[len+1];
  ASSERT (ret, "Out of memory.");
  int *s1 = new int[bestx -r +1];
  ASSERT (s1, "Out of memory");
  int *s2 = new int[besty - c + 1];
  ASSERT (s2, "Out of memory");
  int p1, p2;

  for (i = 0, p1 = p2 =0; i < len; i++){
	  ret[i] = buffer[len - 1 - i];
	  if(ret[i] == 'X') s1[p1++] = -1;
	  else if(ret[i] == 'Y') s2[p2++] = -1;
	  else {
		  s1[p1] = p2 + c + 1; 
		  s2[p2] = p1 + r + 1;
		  p1++; p2++;
	  }
  }
  ret[len] = '\0';
  
  delete pTrace;
  delete[] buffer;
  
  delete[] ret; 

  AlignedFragment *res = new AlignedFragment (sx.GetID(), sy.GetID(), r+1, c+1, bestx, besty, s1, s2);
  delete s1;
  delete s2;
  return res;
}

SCORE_PAIR * ProbModel::ViterbiInitialize(const Sequence &sx, const Sequence &sy, int *map)
{
  int xLen = sx.GetLength();
  int yLen = sy.GetLength();
  const char *x = sx.GetData();
  const char *y = sy.GetData();
  SCORE tmp;
  
  ScoreMatrix *mp = new ScoreMatrix (NUM_STATES, xLen+1, yLen+1);
  ASSERT (mp, "Out of memory.");
  ScoreMatrix &m = *mp;
  m.Fill (LOG_ZERO_SCORE);
  Matrix *pTrace = new Matrix(NUM_STATES,xLen+1,yLen+1);
  Matrix &trace = *pTrace;
  trace.Fill(-1);
  
  // initialization condition
  for (int s = 0; s < NUM_STATES; s++) if (START_STATES[s]){
    if (EMISSIONS[s][0] == 1 && EMISSIONS[s][1] == 0){
      m(s, EMISSIONS[s][0], EMISSIONS[s][1]) = LOG_START[s] + LOG_EMIT_1[(unsigned char) x[1]];
    } else if (EMISSIONS[s][0] == 0 && EMISSIONS[s][1] == 1){
      m(s, EMISSIONS[s][0], EMISSIONS[s][1]) = LOG_START[s] + LOG_EMIT_1[(unsigned char) y[1]];
    } else if (EMISSIONS[s][0] == 1 && EMISSIONS[s][1] == 1){
      m(s, EMISSIONS[s][0], EMISSIONS[s][1]) = LOG_START[s] + LOG_EMIT_2[(unsigned char) x[1]][(unsigned char) y[1]];
    }
  }

  SCORE *ij = m.GetPtr(0,0,0);
  SCORE *i1j = m.GetPtr(0,-1,0);
  SCORE *ij1 = m.GetPtr(0,0,-1);
  SCORE *i1j1 = m.GetPtr(0,-1,-1);
  
  // recursion
  int i, j;
  for (i = 0; i <= xLen; i++){ 
    const unsigned char xi = (unsigned char) (i == 0 ? '~' : x[i]);
    for (j = 0; j <= yLen; j++){ 
      const unsigned char yj = (unsigned char) (j == 0 ? '~' : y[j]);

      if (i > 0){
	for (int k = 0; k < NUM_TRANS_X; k++){
	  int s = TRANSITIONS_EMIT_X[k][0];
	  int t = TRANSITIONS_EMIT_X[k][1];
	  tmp = i1j[s] + LOG_TRANS[s][t] + LOG_EMIT_1[xi];
	  if (tmp > ij[t]){
		  ij[t] = tmp;
		  trace(t,i,j) = s;
	  }
	}
      }

      if (j > 0){
	for (int k = 0; k < NUM_TRANS_Y; k++){
	  int s = TRANSITIONS_EMIT_Y[k][0];
	  int t = TRANSITIONS_EMIT_Y[k][1];
	  tmp = ij1[s] + LOG_TRANS[s][t] + LOG_EMIT_1[yj];
	  if (tmp > ij[t]){
		  ij[t] = tmp;
		  trace(t,i,j) = s;
	  }
	}
      }
      
      if (i > 0 && j > 0){
	for (int k = 0; k < NUM_TRANS_BOTH; k++){
	  int s = TRANSITIONS_EMIT_BOTH[k][0];
	  int t = TRANSITIONS_EMIT_BOTH[k][1];
	  tmp = i1j1[s] + LOG_TRANS[s][t] + LOG_EMIT_2[xi][yj];
	  if (tmp > ij[t]){
		  ij[t] = tmp;
		  trace(t,i,j) = s;
	  }
//	  ij[t] = LOG_ADD_SCORE (ij[t], i1j1[s] + LOG_TRANS[s][t] + LOG_EMIT_2[xi][yj]);
	}
	if(map[i*(yLen+1)+j] == 1){
		ij[MATCH] = LOG_ZERO_SCORE;
		trace(MATCH,i,j) = -1;
	}
	  }
      
      ij += NUM_STATES;
      i1j += NUM_STATES;
      ij1 += NUM_STATES;
      i1j1 += NUM_STATES;
    }
  }
  SCORE_PAIR * result = new SCORE_PAIR(mp, pTrace);
	return result;
}

void ProbModel::ViterbiUpdate(ScoreMatrix *mp, Matrix *pTrace, const Sequence &sx, const Sequence &sy, int *map, 
							  AlignedFragment *frag, int minlength)
{
	int beginx, beginy; 
	int xLen = sx.GetLength();
	int yLen = sy.GetLength();
	const char *x = sx.GetData();
	const char *y = sy.GetData();
	SCORE tmp;
	SCORE *buf = new SCORE[NUM_STATES];
  
	ScoreMatrix &m = *mp;
	Matrix &trace = *pTrace;

	int i,j,k,a,b;
	SCORE *ij, *i1j, *ij1, *i1j1 ;

	int *align = frag->seq[0];
	//For each gap-free alignment update Viterbi tables
	for(a = frag->begin[0], b = 0; a <= frag->end[0]; a++, b++){ 
		if (b > 0 && (align[b] == -1 || align[b] - align[b-1] == 1)) continue;
		beginx = a;
		beginy = align[b];

		int startx = max(beginx - minlength + 1, 1);
		int starty = beginy;
  
		// recursion
		for (i = startx; i < beginx; i++){ 
			ij = m.GetPtr(0,i,starty);
			i1j = m.GetPtr(0,i - 1,starty);
			ij1 = m.GetPtr(0,i,starty - 1);
			i1j1 = m.GetPtr(0,i - 1,starty - 1);
			const unsigned char xi = (unsigned char)x[i];
			for (j = starty; j <= yLen; j++){ 
				const unsigned char yj = (unsigned char)y[j];
	
				for(k = 0; k < NUM_STATES; k++)
					buf[k] = LOG_ZERO_SCORE;

				for (k = 0; k < NUM_TRANS_X; k++){
					int s = TRANSITIONS_EMIT_X[k][0];
					int t = TRANSITIONS_EMIT_X[k][1];
					tmp = i1j[s] + LOG_TRANS[s][t] + LOG_EMIT_1[xi];
					if (tmp > buf[t]){
						buf[t] = tmp;
						trace(t,i,j) = s;
					}
				}

				for (k = 0; k < NUM_TRANS_Y; k++){
					int s = TRANSITIONS_EMIT_Y[k][0];
					int t = TRANSITIONS_EMIT_Y[k][1];
					tmp = ij1[s] + LOG_TRANS[s][t] + LOG_EMIT_1[yj];
					if (tmp > buf[t]){
						buf[t] = tmp;
						trace(t,i,j) = s;
					}
				}
      
				for (k = 0; k < NUM_TRANS_BOTH; k++){
					int s = TRANSITIONS_EMIT_BOTH[k][0];
					int t = TRANSITIONS_EMIT_BOTH[k][1];
					tmp = i1j1[s] + LOG_TRANS[s][t] + LOG_EMIT_2[xi][yj];
					if (tmp > buf[t]){
						buf[t] = tmp;
						trace(t,i,j) = s;
					}
				}

				if(map[i*(yLen+1)+j] == 1){
					buf[MATCH] = LOG_ZERO_SCORE;
					trace(MATCH,i,j) = -1;
				}

				if(memcmp(buf,ij,sizeof(SCORE)*NUM_STATES) != 0)
					memcpy(ij,buf,sizeof(SCORE)*NUM_STATES);
				else
					break;		  
      
				ij += NUM_STATES;
				i1j += NUM_STATES;
				ij1 += NUM_STATES;
				i1j1 += NUM_STATES;
			}
		}

		startx = beginx;
		starty = max(beginy - minlength + 1, 1);
		int lastUpdate, currentUpdate = starty-1;
		//For each cell from the left top
		for(bool update = true; update && startx <= xLen && starty <= yLen; startx++,starty++){
		//Update the respective row
			lastUpdate = currentUpdate;
			currentUpdate = starty;
			update = false;
			{
			ij = m.GetPtr(0,startx,starty);
			i1j = m.GetPtr(0,startx - 1,starty);
			ij1 = m.GetPtr(0,startx,starty - 1);
			i1j1 = m.GetPtr(0,startx - 1,starty - 1);
			
			const unsigned char xi = (unsigned char)x[startx];
			for (j = starty; j <= yLen; j++){ 
				const unsigned char yj = (unsigned char)y[j];
	
				for(k = 0; k < NUM_STATES; k++)
					buf[k] = LOG_ZERO_SCORE;

				for (k = 0; k < NUM_TRANS_X; k++){
					int s = TRANSITIONS_EMIT_X[k][0];
					int t = TRANSITIONS_EMIT_X[k][1];
					tmp = i1j[s] + LOG_TRANS[s][t] + LOG_EMIT_1[xi];
					if (tmp > buf[t]){
						buf[t] = tmp;
						trace(t,startx,j) = s;
					}
				}

				for (k = 0; k < NUM_TRANS_Y; k++){
					int s = TRANSITIONS_EMIT_Y[k][0];
					int t = TRANSITIONS_EMIT_Y[k][1];
					tmp = ij1[s] + LOG_TRANS[s][t] + LOG_EMIT_1[yj];
					if (tmp > buf[t]){
						buf[t] = tmp;
						trace(t,startx,j) = s;
					}
				}
      
				for (k = 0; k < NUM_TRANS_BOTH; k++){
					int s = TRANSITIONS_EMIT_BOTH[k][0];
					int t = TRANSITIONS_EMIT_BOTH[k][1];
					tmp = i1j1[s] + LOG_TRANS[s][t] + LOG_EMIT_2[xi][yj];
					if (tmp > buf[t]){
						buf[t] = tmp;
						trace(t,startx,j) = s;
					}
				}

				if(map[startx*(yLen+1)+j] == 1){
					buf[MATCH] = LOG_ZERO_SCORE;
					trace(MATCH,startx,j) = -1;
				}

				if(map[startx*(yLen+1)+j] == 1 || memcmp(buf,ij,sizeof(SCORE)*NUM_STATES) != 0){
					memcpy(ij,buf,sizeof(SCORE)*NUM_STATES);
					update = true;
					currentUpdate = j;
				}
				else if ( j >= lastUpdate+1)
						break;
      
				ij += NUM_STATES;
				i1j += NUM_STATES;
				ij1 += NUM_STATES;
				i1j1 += NUM_STATES;
			}
			}

		//Now update the respective column
			const unsigned char yj = (unsigned char)y[starty];
			for (i = startx+1; i <= xLen; i++){
				ij = m.GetPtr(0,i,starty);
				i1j = m.GetPtr(0,i - 1,starty);
				ij1 = m.GetPtr(0,i,starty - 1);
				i1j1 = m.GetPtr(0,i - 1,starty - 1);
				const unsigned char xi = (unsigned char)x[i];

				for(k = 0; k < NUM_STATES; k++)
					buf[k] = LOG_ZERO_SCORE;

				for (k = 0; k < NUM_TRANS_X; k++){
					int s = TRANSITIONS_EMIT_X[k][0];
					int t = TRANSITIONS_EMIT_X[k][1];
					tmp = i1j[s] + LOG_TRANS[s][t] + LOG_EMIT_1[xi];
					if (tmp > buf[t]){
						buf[t] = tmp;
						trace(t,i,starty) = s;
					}
				}

				for (k = 0; k < NUM_TRANS_Y; k++){
					int s = TRANSITIONS_EMIT_Y[k][0];
					int t = TRANSITIONS_EMIT_Y[k][1];
					tmp = ij1[s] + LOG_TRANS[s][t] + LOG_EMIT_1[yj];
					if (tmp > buf[t]){
						buf[t] = tmp;
						trace(t,i,starty) = s;
					}
				}
      
				for (k = 0; k < NUM_TRANS_BOTH; k++){
					int s = TRANSITIONS_EMIT_BOTH[k][0];
					int t = TRANSITIONS_EMIT_BOTH[k][1];
					tmp = i1j1[s] + LOG_TRANS[s][t] + LOG_EMIT_2[xi][yj];
					if (tmp > buf[t]){
						buf[t] = tmp;
						trace(t,i,starty) = s;
					}
				}

				if(map[i*(yLen+1)+starty] == 1){
					buf[MATCH] = LOG_ZERO_SCORE;
					trace(MATCH,i,starty) = -1;
				}

				if(map[i*(yLen+1)+starty] == 1 || memcmp(buf,ij,sizeof(SCORE)*NUM_STATES) != 0){
					memcpy(ij,buf,sizeof(SCORE)*NUM_STATES);
					update = true;
				}
				else
					break;
			}
		}
		//If the gap-free fragment is to short
		if(startx - beginx < minlength){
			starty = beginy;
		for (i = startx; i <= xLen && i < beginx + minlength; i++){ 
			ij = m.GetPtr(0,i,starty);
			i1j = m.GetPtr(0,i - 1,starty);
			ij1 = m.GetPtr(0,i,starty - 1);
			i1j1 = m.GetPtr(0,i - 1,starty - 1);
			const unsigned char xi = (unsigned char)x[i];
			for (j = starty; j <= yLen; j++){ 
				const unsigned char yj = (unsigned char)y[j];
	
				for(k = 0; k < NUM_STATES; k++)
					buf[k] = LOG_ZERO_SCORE;

				for (k = 0; k < NUM_TRANS_X; k++){
					int s = TRANSITIONS_EMIT_X[k][0];
					int t = TRANSITIONS_EMIT_X[k][1];
					tmp = i1j[s] + LOG_TRANS[s][t] + LOG_EMIT_1[xi];
					if (tmp > buf[t]){
						buf[t] = tmp;
						trace(t,i,j) = s;
					}
				}

				for (k = 0; k < NUM_TRANS_Y; k++){
					int s = TRANSITIONS_EMIT_Y[k][0];
					int t = TRANSITIONS_EMIT_Y[k][1];
					tmp = ij1[s] + LOG_TRANS[s][t] + LOG_EMIT_1[yj];
					if (tmp > buf[t]){
						buf[t] = tmp;
						trace(t,i,j) = s;
					}
				}
      
				for (k = 0; k < NUM_TRANS_BOTH; k++){
					int s = TRANSITIONS_EMIT_BOTH[k][0];
					int t = TRANSITIONS_EMIT_BOTH[k][1];
					tmp = i1j1[s] + LOG_TRANS[s][t] + LOG_EMIT_2[xi][yj];
					if (tmp > buf[t]){
						buf[t] = tmp;
						trace(t,i,j) = s;
					}
				}

				if(map[i*(yLen+1)+j] == 1){
					buf[MATCH] = LOG_ZERO_SCORE;
					trace(MATCH,i,j) = -1;
				}

				if(memcmp(buf,ij,sizeof(SCORE)*NUM_STATES) != 0)
					memcpy(ij,buf,sizeof(SCORE)*NUM_STATES);
				else
					break;		  
      
				ij += NUM_STATES;
				i1j += NUM_STATES;
				ij1 += NUM_STATES;
				i1j1 += NUM_STATES;
			}
		}

		}
	}
		
  delete buf;
}


AlignedFragment * ProbModel::OneAligment(const Sequence &sx, const Sequence &sy, Matrix &trace, ScoreMatrix &m)
{
	int xLen = sx.GetLength();
    int yLen = sy.GetLength();
	float bestProb = LOG_ZERO_SCORE;
    int state = -1;
	int i;
    for(i = 0; i < NUM_STATES; i++){
	  if(FINAL_STATES[i]){
		  float crProb = m(i,xLen,yLen) + LOG_FINAL[i];
		  if(crProb > bestProb){
			  bestProb = crProb;
			  state = i;
		  }
	  }
	}
  	char *buffer = new char[(xLen+1) * (yLen+1)];
	ASSERT (buffer, "Out of memory.");

	int r = xLen, c = yLen, len = 0;
	int newState;
	while((r != 0 || c != 0) && state != MATCH){
		newState = (int)trace(state,r,c);
		switch(state){
			case AFT_X:r--;break;
			case AFT_Y:c--;break;
			default: ASSERT(false,"Insert state within flanking region");
		}
		state = newState;
	}
	int bestx = r;
	int besty = c;
	while((r != 0 || c != 0) && (state != BEF_X && state !=BEF_Y)){
		newState = (int)trace(state,r,c);
		switch(state){
			case MATCH: c--;r--;buffer[len++] = 'B';break;
			case INS_X: r--;buffer[len++] = 'X';break;
			case INS_Y: c--;buffer[len++] = 'Y';break;
			default: ASSERT(false,"Flanking state inside");
		}
		state = newState;
	}
	  char *ret = new char[len+1];
  ASSERT (ret, "Out of memory.");
  int *s1 = new int[bestx -r +1];
  ASSERT (s1, "Out of memory");
  int *s2 = new int[besty - c + 1];
  ASSERT (s2, "Out of memory");
  int p1, p2;

  for (i = 0, p1 = p2 =0; i < len; i++){
	  ret[i] = buffer[len - 1 - i];
	  if(ret[i] == 'X') s1[p1++] = -1;
	  else if(ret[i] == 'Y') s2[p2++] = -1;
	  else {
		  s1[p1] = p2 + c + 1; 
		  s2[p2] = p1 + r + 1;
		  p1++; p2++;
	  }
  }
  ret[len] = '\0';
  
  delete[] buffer;
  
  delete[] ret;

  AlignedFragment *res = new AlignedFragment (sx.GetID(), sy.GetID(), r+1, c+1, bestx, besty, s1, s2);
  delete s1;
  delete s2;
  return res;
}
