//////////////////////////////////////////////////////////////////////
// SparseMatrix.h
//
// Sparse matrix storage class for storing a set of two-dimensional
// arrays.
//////////////////////////////////////////////////////////////////////

#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <stdio.h>
#include "Matrix.h"

//////////////////////////////////////////////////////////////////////
// Sparse matrix object
//////////////////////////////////////////////////////////////////////

class Matrix;

class SparseMatrix {
  friend class Matrix;

 public:
  struct SparseMatrixEntry {
    int column;
    float value;
  };
  
 private:  
  float threshold, missing;
  int layers;
  int rows;
  int cols;
  int numEntries;

  SparseMatrixEntry *data;
  int *rowSize;
  SparseMatrixEntry **rowPtrs;

  // default constructor (used for ComputeTranspose())
  SparseMatrix(){}
  
  // printing utility function
  void PrintVal (FILE *file, const float &value) const;

 public:

  // constructor and destructor
  SparseMatrix (const Matrix &matrix, const float &threshold, const float &missing);
  ~SparseMatrix ();

  // compute transpose
  SparseMatrix *ComputeTranspose() const;

  // row accessors
  const SparseMatrixEntry *GetRowPtr (int layer, int row) const;
  const int GetRowSize (int layer, int row) const;

  // matrix dimension accessors
  const int GetNumLayers() const;
  const int GetNumRows() const;
  const int GetNumCols() const;
  

  // printing functions
  void PrintLayer (FILE *file, int layer) const;
  void Print (FILE *file) const;

  // accessor
  const float &operator() (int layer, int row, int col) const;

};

#endif
