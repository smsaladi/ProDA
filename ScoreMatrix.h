//////////////////////////////////////////////////////////////////////
// ScoreMatrix.h
//
// Matrix storage class for storing a set of two-dimensional arrays.
//////////////////////////////////////////////////////////////////////

#ifndef SCOREMATRIX_H
#define SCOREMATRIX_H

#include <stdio.h>
#include "Score.h"
#include "Sequence.h"

//////////////////////////////////////////////////////////////////////
// ScoreMatrix object
//////////////////////////////////////////////////////////////////////

class ScoreMatrix {
  friend class Matrix;

  int layers;
  int rows;
  int cols;

  SCORE *data;
  
  // printing utility function
  void PrintVal (FILE *file, const SCORE &value) const;

 public:
	 
  // constructors and destructor
  ScoreMatrix (int layers, int rows, int cols);
  ScoreMatrix (const ScoreMatrix &m);
  ~ScoreMatrix ();

  // fill all entries with value
  void Fill (const SCORE &value);

  // printing functions
  void PrintLayer (FILE *file, int layer) const;
  void Print (FILE *file) const;
  void PrintSumRange(FILE *file, int beginy, int endy, int beginx, int endx);


  //////////////////////////////////////////////////////////////////////
  // Access matrix element
  //////////////////////////////////////////////////////////////////////
  
  SCORE &operator() (int layer, int row, int col){
/*    ASSERT (0 <= layer && layer < layers, "Requested layer out-of-bounds.");
    ASSERT (0 <= row && row < rows, "Requested row out-of-bounds.");
    ASSERT (0 <= col && col < cols, "Requested column out-of-bounds.");*/
    return data[(row * cols + col) * layers + layer];
  }

  //////////////////////////////////////////////////////////////////////
  // Access matrix element (const version)
  //////////////////////////////////////////////////////////////////////
  
  const SCORE &operator() (int layer, int row, int col) const {
    /*ASSERT (0 <= layer && layer < layers, "Requested layer out-of-bounds.");
    ASSERT (0 <= row && row < rows, "Requested row out-of-bounds.");
    ASSERT (0 <= col && col < cols, "Requested column out-of-bounds.");*/
    return data[(row * cols + col) * layers + layer];
  }
    
  //////////////////////////////////////////////////////////////////////
  // Access matrix element
  //////////////////////////////////////////////////////////////////////
  
  SCORE *GetPtr (int layer, int row, int col){
    return data + (row * cols + col) * layers + layer;
  }

  //////////////////////////////////////////////////////////////////////
  // Access matrix element (const version)
  //////////////////////////////////////////////////////////////////////
  
  const SCORE *GetPtr (int layer, int row, int col) const {
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
