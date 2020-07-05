// LocalAlign.cpp: implementation of the LocalAlign class.
//
//////////////////////////////////////////////////////////////////////

#include "LocalAlign.h"


AlignedFragment * LocalAlign::ComputeLocalAlignment(const Sequence &seq1, const Sequence &seq2, 
													const Matrix &m, float *score)
{
	int i, j;
	const int rows = m.GetNumRows();
	const int cols = m.GetNumCols();
  
  // memory allocation
	Matrix *bestp = new Matrix (1, rows, cols);
	ASSERT (bestp, "Out of memory.");
	Matrix &best = *bestp;
	TracebackType *traceback = new TracebackType[rows * cols];
	ASSERT (traceback, "Out of memory.");

	float *nullX = new float[rows];
	float *nullY = new float[cols];
	ASSERT(nullX, "Out of memory.");
	ASSERT(nullY, "Out of memory.");

	for (i = 0; i < rows; i++) 
		nullX[i] = m.SumOfRow(BEF_X, i) + m.SumOfRow(AFT_X, i);
	for (i = 0; i < cols; i++) 
		nullY[i] = m.SumOfColumn(BEF_Y, i) + m.SumOfColumn(AFT_Y, i);

	float bestScore = 0;
	float tmp;
	int bestx, besty;
	bestx = besty = 0;

	best.Fill (LOG_ZERO_FLOAT);
	best(0,0,0) = 0;

	
	for (i = 0; i < rows; i ++){
		for (j = 0; j < cols; j++){
			if (i > 0 && best(0,i-1,j) > best(0,i,j)){
				best(0,i,j) = best(0,i-1,j);
				traceback[i * cols + j] = UP;
			}
			if (j > 0 && best(0,i,j-1) > best(0,i,j)){
				best(0,i,j) = best(0,i,j-1);
				traceback[i * cols + j] = LEFT;
			}
			
			if (i > 0 && j > 0 && (tmp = best(0,i-1,j-1)+m(MATCH,i,j)-(nullX[i]+nullY[j]-nullX[i]*nullY[j])/2) > best(0,i,j)){
				best(0,i,j) = tmp;
				traceback[i * cols + j] = UP_LEFT;
			}
			if (best(0,i,j) <= 0){ 
				best(0,i,j) = 0;
				traceback[i * cols + j] = NONE;
			}
			if (best(0,i,j) > bestScore) {
				bestScore = best(0,i,j);
				bestx = i;
				besty = j;
			}
		}
	}
	 delete bestp;
	 delete nullX;
	 delete nullY;
     
	 if(score) {
		 *score = bestScore;
		 delete traceback;
		 return NULL;
	 }

  // follow tracebacks
	char *buffer = new char[rows * cols];
	ASSERT (buffer, "Out of memory.");

	int r = bestx, c = besty, len = 0;
	while (traceback[r * cols + c] != NONE){
		switch (traceback[r * cols + c]){
		case UP: r--; buffer[len++] = 'X'; break;
		case LEFT: c--; buffer[len++] = 'Y'; break;
		case UP_LEFT: r--; c--; buffer[len++] = 'B'; break;
		default: ASSERT (false, "Unexpected value found in traceback matrix!");
    }
	}
  
  
  delete[] traceback;

  // reverse alignment path
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
		  s1[p1] = p2+c+1; 
		  s2[p2] = p1+r+1;
		  p1++; p2++;
	  }
  }
  ret[len] = '\0';
  
  delete[] buffer;
  
  delete[] ret; 

  AlignedFragment *res = new AlignedFragment (seq1.GetID(), seq2.GetID(), r+1, c+1, bestx, besty, s1, s2);
  delete s1;
  delete s2;
  return res;
}
