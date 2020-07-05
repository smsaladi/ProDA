//////////////////////////////////////////////////////////////////////
// Consistency.cc
//////////////////////////////////////////////////////////////////////

#include "Consistency.h"

void AccumulateConsistencyInfo (SparseMatrix **posteriors, Matrix &p, int x, int y, int z, int n){
  int lenX = posteriors[x*n+y]->GetNumRows()-1;
  SparseMatrix *XZ = posteriors[x*n+z];
  SparseMatrix *ZY = posteriors[z*n+y];
  
  for (int i = 1; i <= lenX; i++){

    const SparseMatrix::SparseMatrixEntry *XZptr = XZ->GetRowPtr(0,i);
    for (int a = 0; a < XZ->GetRowSize(0,i); a++, ++XZptr){
      int k = XZptr->column;
      float val1 = XZptr->value;

      const SparseMatrix::SparseMatrixEntry *ZYptr = ZY->GetRowPtr(0,k);
      for (int b = 0; b < ZY->GetRowSize(0,k); b++, ++ZYptr){
	int j = ZYptr->column;
	float val2 = ZYptr->value;
	
	p(0,i,j) += val1 * val2;
      }
    }
  }
}

typedef SparseMatrix *SparseMatrixPtr;

SparseMatrix **ProbabilisticConsistency (SparseMatrix **posteriors, int n){
  SparseMatrix **newPosteriors = new SparseMatrixPtr[n*n];
  ASSERT (newPosteriors, "Out of memory.");

  int i,j,r,c;
  for (i = 0; i < n*n; i++) newPosteriors[i] = NULL;
 
  for (i = 0; i < n-1; i++){
    for (j = i+1; j < n; j++){
      Matrix *p = new Matrix (*posteriors[i*n+j]);
      int rows = p->GetNumRows();
      int cols = p->GetNumCols();

      for (r = 0; r < rows; r++)
	for (c = 0; c < cols; c++)
	  (*p)(0,r,c) *= 2;

      for (int k = 0; k < n; k++) if (k != i && k != j)
	AccumulateConsistencyInfo (posteriors, *p, i, j, k, n);
      
      for (r = 0; r < rows; r++)
	for (c = 0; c < cols; c++)
	  (*p)(0,r,c) /= n;
      
      newPosteriors[i*n+j] = new SparseMatrix (*p, 0.01, 0);
      delete p;
      newPosteriors[j*n+i] = newPosteriors[i*n+j]->ComputeTranspose();
    }
  }
  
  return newPosteriors;
}
