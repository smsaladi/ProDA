//////////////////////////////////////////////////////////////////////
// Matrix.cc
//////////////////////////////////////////////////////////////////////

#include <string.h>
#include "Assert.h"
#include "Matrix.h"

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////

Matrix::Matrix (int layers, int rows, int cols) : 
  layers(layers), rows(rows), cols(cols) {
  
  ASSERT (layers >= 0, "Number of layers in matrix must be positive.");
  ASSERT (rows >= 0, "Number of rows in matrix must be positive.");
  ASSERT (cols >= 0, "Number of columns in matrix must be positive.");
    
  data = new float[layers * rows * cols];
  ASSERT (data, "Out of memory.");
}

//////////////////////////////////////////////////////////////////////
// Copy constructor
//////////////////////////////////////////////////////////////////////

Matrix::Matrix (const Matrix &m) :
  layers (m.layers), rows(m.rows), cols(m.cols) {
  
  ASSERT (layers >= 0, "Number of layers in matrix must be positive.");
  ASSERT (rows >= 0, "Number of rows in matrix must be positive.");
  ASSERT (cols >= 0, "Number of columns in matrix must be positive.");
    
  data = new float[layers * rows * cols];
  ASSERT (data, "Out of memory.");

  memcpy (data, m.data, sizeof(float) * (layers * rows * cols));
}

//////////////////////////////////////////////////////////////////////
// Copy constructor (using ScoreMatrix)
//////////////////////////////////////////////////////////////////////

Matrix::Matrix (const ScoreMatrix &m) :
  layers (m.layers), rows(m.rows), cols(m.cols) {
  
  ASSERT (layers >= 0, "Number of layers in matrix must be positive.");
  ASSERT (rows >= 0, "Number of rows in matrix must be positive.");
  ASSERT (cols >= 0, "Number of columns in matrix must be positive.");
    
  data = new float[layers * rows * cols];
  ASSERT (data, "Out of memory.");

  for (int i = 0; i < layers * rows * cols; i++)
    data[i] = EXP_FLOAT(TO_FLOAT(m.data[i]));
}

//////////////////////////////////////////////////////////////////////
// Constructor (using SparseMatrix)
//////////////////////////////////////////////////////////////////////

Matrix::Matrix (const SparseMatrix &sm) : 
  layers (sm.layers), rows (sm.rows), cols (sm.cols){

  ASSERT (layers >= 0, "Number of layers in matrix must be positive.");
  ASSERT (rows >= 0, "Number of rows in matrix must be positive.");
  ASSERT (cols >= 0, "Number of columns in matrix must be positive.");
    
  data = new float[layers * rows * cols];
  ASSERT (data, "Out of memory.");

  for (int i = 0; i < layers * rows * cols; i++)
    data[i] = sm.missing;

  {for (int i = 0; i < layers; i++){
    float *dataPtr = data + i;
    for (int j = 0; j < rows; j++){
      for (int k = 0; k < sm.rowSize[i * rows + j]; k++){
	SparseMatrix::SparseMatrixEntry *smCell = &sm.rowPtrs[i * rows + j][k];
	dataPtr[layers * smCell->column] = smCell->value;
      }
      dataPtr += cols * layers;
    }
  }}
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////

Matrix::~Matrix (){
  delete[] data;
}

//////////////////////////////////////////////////////////////////////
// Fill all matrix values
//////////////////////////////////////////////////////////////////////

void Matrix::Fill (const float &value){
  for (int i = 0; i < layers * rows * cols; i++)
    data[i] = value;
}

//////////////////////////////////////////////////////////////////////
// Printing utility function
//////////////////////////////////////////////////////////////////////

void Matrix::PrintVal (FILE *file, const float &value) const {
  if (value == LOG_ZERO_FLOAT)
    fprintf (file, "\t-inf");
  else
    fprintf (file, "\t%d", (int)value);
}

//////////////////////////////////////////////////////////////////////
// Print a single matrix layer
//////////////////////////////////////////////////////////////////////

void Matrix::PrintLayer (FILE *file, int layer) const {
	int i, j;
	for (j = 0; j < cols; j++)
		fprintf (file,"eee\t");
	fprintf(file, "\n");
  for (i = 0; i < rows; i++){
    //fprintf (file, "%s[", (i == 0 ? "[" : " "));
	  fprintf (file, "first");
    for (j = 0; j < cols; j++){
      //if (j > 0) fprintf (file, ", ");
		if (j > 0) fprintf (file, "\t");
      PrintVal (file, operator()(layer,i,j));
    }
    //fprintf (file, "]%s\n", (i == rows-1 ? "]" : ""));
	fprintf(file, "\n");
  }
}

//////////////////////////////////////////////////////////////////////
// Print all matrix layers
//////////////////////////////////////////////////////////////////////

void Matrix::Print (FILE *file) const {
  for (int i = 0; i < layers; i++)
    PrintLayer (file, i);
}

//////////////////////////////////////////////////////////////////////
// Compute sum of all entries in matrix
//////////////////////////////////////////////////////////////////////

float Matrix::ComputeSum() const {
  float total = 0;
  for (int i = 0; i < layers * rows * cols; i++)
    total += data[i];
  return total;
}

//////////////////////////////////////////////////////////////////////
// Compute sum of a row in matrix
//////////////////////////////////////////////////////////////////////

float Matrix::SumOfRow(int layer, int row) const
{
	float max = 0;
	int pos = 0;
	float *ptr;

	ptr = (float *)GetPtr(layer, row, 0);
	for(int i = 0; i < cols; i++){
		if(max < operator()(layer,row,i)){
			max = operator()(layer,row,i);
			pos = i;
		}
	
	}
	return max;
}

//////////////////////////////////////////////////////////////////////
// Compute sum of a column in matrix
//////////////////////////////////////////////////////////////////////

float Matrix::SumOfColumn(int layer, int column) const
{
	float sum = 0;
	int pos = 0;
	for (int i = 0; i < rows; i++){
		
		if(sum < operator()(layer, i, column)){
			sum = operator()(layer, i, column);
			pos = i;
		}
	}

	return sum;
}


void Matrix::PrintRange(FILE *file, int layer, int beginy, int endy, int beginx, int endx)
{
	ASSERT(beginx > 0 && endx < cols && beginy > 0 && endy < rows, "Out of range in PrintRange");
	fprintf(file, " \t");
	for (int k = beginx; k <= endx; k++)
		fprintf(file, "\t%d", k);
	fprintf (file, "\n");
	for (int i = beginy; i <= endy; i++){
		fprintf(file, "\t%d ", i);
    for (int j = beginx; j <= endx; j++){
      //if (j > beginx) fprintf (file, ", ");
      PrintVal (file, operator()(layer,i,j));
    }
    fprintf (file, "\n");
  }
  fprintf (file, "\n");

}
