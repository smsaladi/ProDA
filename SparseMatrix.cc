//////////////////////////////////////////////////////////////////////
// SparseMatrix.cc
//////////////////////////////////////////////////////////////////////

#include "Assert.h"
#include "SparseMatrix.h"

typedef SparseMatrix::SparseMatrixEntry *SparseMatrixEntryPtr;

//////////////////////////////////////////////////////////////////////
// Constructor (using Matrix)
//////////////////////////////////////////////////////////////////////

SparseMatrix::SparseMatrix (const Matrix &matrix, const float &threshold, const float &missing) :
  threshold(threshold), missing(missing), 
  layers(matrix.layers), rows(matrix.rows), cols(matrix.cols), numEntries (0) {
  
  ASSERT (layers >= 0, "Number of layers in matrix must be positive.");
  ASSERT (rows >= 0, "Number of rows in matrix must be positive.");
  ASSERT (cols >= 0, "Number of columns in matrix must be positive.");

  // count number of entries needed

  for (int i = 0; i < layers * rows * cols; i++)
    numEntries += (matrix.data[i] >= threshold);
  
  // allocate memory

  data = new SparseMatrixEntry[numEntries];
  ASSERT (data, "Out of memory.");

  rowSize = new int[layers * rows];
  ASSERT (rowSize, "Out of memory.");

  rowPtrs = new SparseMatrixEntryPtr[layers * rows];
  ASSERT (rowPtrs, "Out of memory.");

  // build sparse matrices, layer-by-layer
  
  SparseMatrixEntry *dataPtr = data;
  for (int k = 0; k < layers; k++){
    for (int i = 0; i < rows; i++){
      int numColsUsed = 0;
      for (int j = 0; j < cols; j++){
	if (matrix(k,i,j) >= threshold){
	  dataPtr->column = j;
	  dataPtr->value = matrix(k,i,j);
	  ++dataPtr;
	  numColsUsed++;
	}
      }
      rowSize[k * rows + i] = numColsUsed;
    }
  }

  // compute pointers to beginning of each row

  rowPtrs[0] = data;
  {for (int i = 1; i < layers * rows; i++)
	  rowPtrs[i] = rowPtrs[i-1] + rowSize[i-1];}
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////

SparseMatrix::~SparseMatrix (){
  delete[] data;
  delete[] rowSize;
  delete[] rowPtrs;
}

//////////////////////////////////////////////////////////////////////
// Compute transpose of sparse matrix
//////////////////////////////////////////////////////////////////////

SparseMatrix *SparseMatrix::ComputeTranspose() const {
  SparseMatrix *sm = new SparseMatrix();
  ASSERT (sm, "Out of memory.");
  
  // fill in basic information

  sm->threshold = threshold;
  sm->missing = missing;
  sm->layers = layers;
  sm->rows = cols;
  sm->cols = rows;
  sm->numEntries = numEntries;
  
  // allocate memory

  sm->data = new SparseMatrixEntry[sm->numEntries];
  ASSERT (sm->data, "Out of memory.");

  sm->rowSize = new int[sm->layers * sm->rows];
  ASSERT (sm->rowSize, "Out of memory.");

  sm->rowPtrs = new SparseMatrixEntryPtr[sm->layers * sm->rows];
  ASSERT (sm->rowPtrs, "Out of memory.");

  // compute row sizes
  
  SparseMatrixEntry *dataPtr = data;
  
  for (int k = 0; k < layers; k++){
    for (int j = 0; j < sm->rows; j++)
      sm->rowSize[k * sm->rows + j] = 0;
    
    for (int i = 0; i < rows; i++)
      for (int j = 0; j < rowSize[k * rows + i]; j++)
	sm->rowSize[k * sm->rows + (dataPtr++)->column]++;
  }

  // compute pointers to beginning of each row

  sm->rowPtrs[0] = sm->data;
  for (int i = 1; i < sm->layers * sm->rows; i++)
    sm->rowPtrs[i] = sm->rowPtrs[i-1] + sm->rowSize[i-1];

  // initialize pointers for writing data to new sparse matrix

  SparseMatrixEntry **writePtrs = new SparseMatrixEntryPtr[sm->layers * sm->rows];
  ASSERT (writePtrs, "Out of memory.");
  {for (int i = 0; i < sm->layers * sm->rows; i++)
	  writePtrs[i] = sm->rowPtrs[i];}

  // now find transpose of data
  
  dataPtr = data;
  
  {for (int k = 0; k < layers; k++){
    for (int i = 0; i < rows; i++){
      for (int j = 0; j < rowSize[k * rows + i]; j++){
	writePtrs[k * sm->rows + dataPtr->column]->column = i;
	writePtrs[k * sm->rows + dataPtr->column]->value = dataPtr->value;
	++writePtrs[k * sm->rows + dataPtr->column];
	++dataPtr;
      }
    }
  }}
  
  delete[] writePtrs;
  
  return sm;
}

//////////////////////////////////////////////////////////////////////
// Printing utility function
//////////////////////////////////////////////////////////////////////

void SparseMatrix::PrintVal (FILE *file, const float &value) const {
  if (value == LOG_ZERO_FLOAT)
    fprintf (file, " -inf");
  else
    fprintf (file, "%5.2f", value);
}

//////////////////////////////////////////////////////////////////////
// Print a single matrix layer
//////////////////////////////////////////////////////////////////////

void SparseMatrix::PrintLayer (FILE *file, int layer) const {
  for (int i = 0; i < rows; i++){
    fprintf (file, "%s[", (i == 0 ? "[" : " "));
    for (int j = 0; j < rowSize[layer * rows + i]; j++){
      if (j > 0) fprintf (file, ", ");
      fprintf (file, "(%2d,", rowPtrs[layer * rows + i][j].column);
      PrintVal (file, rowPtrs[layer * rows + i][j].value);
      fprintf (file, ")");
    }
    fprintf (file, "]%s\n", (i == rows-1 ? "]" : ""));
  }
  fprintf (file, "\n");
}

//////////////////////////////////////////////////////////////////////
// Print all matrix layers
//////////////////////////////////////////////////////////////////////

void SparseMatrix::Print (FILE *file) const {
  for (int i = 0; i < layers; i++)
    PrintLayer (file, i);
}

//////////////////////////////////////////////////////////////////////
// Return pointer to row
//////////////////////////////////////////////////////////////////////

const SparseMatrix::SparseMatrixEntry *SparseMatrix::GetRowPtr (int layer, int row) const {
  ASSERT (0 <= layer && layer < layers, "Requested layer out-of-bounds.");
  ASSERT (0 <= row && row < rows, "Requested row out-of-bounds.");
  return rowPtrs[layer * rows + row];
}

//////////////////////////////////////////////////////////////////////
// Return size of row
//////////////////////////////////////////////////////////////////////

const int SparseMatrix::GetRowSize (int layer, int row) const {
  ASSERT (0 <= layer && layer < layers, "Requested layer out-of-bounds.");
  ASSERT (0 <= row && row < rows, "Requested row out-of-bounds.");
  return rowSize[layer * rows + row];
}

//////////////////////////////////////////////////////////////////////
// Return number of matrix layers
//////////////////////////////////////////////////////////////////////

const int SparseMatrix::GetNumLayers() const {
  return layers;
}

//////////////////////////////////////////////////////////////////////
// Return number of matrix rows
//////////////////////////////////////////////////////////////////////

const int SparseMatrix::GetNumRows() const {
  return rows;
}

//////////////////////////////////////////////////////////////////////
// Return number of matrix columns
//////////////////////////////////////////////////////////////////////

const int SparseMatrix::GetNumCols() const {
  return cols;
}

//////////////////////////////////////////////////////////////////////
// Access matrix element (const version)
//////////////////////////////////////////////////////////////////////

const float &SparseMatrix::operator() (int layer, int row, int col) const {
  ASSERT (0 <= layer && layer < layers, "Requested layer out-of-bounds.");
  ASSERT (0 <= row && row < rows, "Requested row out-of-bounds.");
  ASSERT (0 <= col && col < cols, "Requested column out-of-bounds.");
  for (int i = 0; i < rowSize[layer * rows + row]; i++)
    if (rowPtrs[layer * rows + row][i].column == col)
      return rowPtrs[layer * rows + row][i].value;
  return missing;
}

