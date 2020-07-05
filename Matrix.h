//////////////////////////////////////////////////////////////////////
// Matrix.h
//
// Matrix storage class for storing a set of two-dimensional arrays.
//////////////////////////////////////////////////////////////////////

#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include "Score.h"
#include "SparseMatrix.h"
#include "ScoreMatrix.h"

//////////////////////////////////////////////////////////////////////
// Matrix object
//////////////////////////////////////////////////////////////////////

class SparseMatrix;

class Matrix {
  friend class SparseMatrix;

  int layers;
  int rows;
  int cols;

  float *data;
  
  // printing utility function
  void PrintVal (FILE *file, const float &value) const;
  void PrintRange(FILE *file, int layer, int beginy,int endy,int beginx,int endx);

 public:
	 
  // constructors and destructor
  Matrix (int layers, int rows, int cols);
  Matrix (const Matrix &m);
  Matrix (const ScoreMatrix &m);
  Matrix (const SparseMatrix &sm);
  ~Matrix ();

  // fill all entries with value
  void Fill (const float &value);

  // printing functions
  void PrintLayer (FILE *file, int layer) const;
  void Print (FILE *file) const;

  // compute sum of all entries
  float ComputeSum() const;

  //Computing sum of a row and a column
  float SumOfColumn(int layer, int column) const;
  float SumOfRow(int layer, int row) const;


  //////////////////////////////////////////////////////////////////////
  // Access matrix element
  //////////////////////////////////////////////////////////////////////
  
  float &operator() (int layer, int row, int col){
    ASSERT (0 <= layer && layer < layers, "Requested layer out-of-bounds.");
    ASSERT (0 <= row && row < rows, "Requested row out-of-bounds.");
    ASSERT (0 <= col && col < cols, "Requested column out-of-bounds.");
    return data[(row * cols + col) * layers + layer];
  }

  //////////////////////////////////////////////////////////////////////
  // Access matrix element (const version)
  //////////////////////////////////////////////////////////////////////
  
  const float &operator() (int layer, int row, int col) const {
    ASSERT (0 <= layer && layer < layers, "Requested layer out-of-bounds.");
    ASSERT (0 <= row && row < rows, "Requested row out-of-bounds.");
    ASSERT (0 <= col && col < cols, "Requested column out-of-bounds.");
    return data[(row * cols + col) * layers + layer];
  }
    
  //////////////////////////////////////////////////////////////////////
  // Access matrix element
  //////////////////////////////////////////////////////////////////////
  
  float *GetPtr (int layer, int row, int col){
    return data + (row * cols + col) * layers + layer;
  }

  //////////////////////////////////////////////////////////////////////
  // Access matrix element (const version)
  //////////////////////////////////////////////////////////////////////
  
  const float *GetPtr (int layer, int row, int col) const {
    return data + (row * cols + col) * layers + layer;
  }

  //////////////////////////////////////////////////////////////////////
  // Return number of matrix layers
  //////////////////////////////////////////////////////////////////////
  
  const int GetNumLayers() const {
    return layers;
  }
  
  //////////////////////////////////////////////////////////////////////
  // Return number of matrix rows
  //////////////////////////////////////////////////////////////////////
  
  const int GetNumRows() const {
    return rows;
  }
  
  //////////////////////////////////////////////////////////////////////
  // Return number of matrix columns
  //////////////////////////////////////////////////////////////////////
  
  const int GetNumCols() const {
    return cols;
  }
};

#endif
