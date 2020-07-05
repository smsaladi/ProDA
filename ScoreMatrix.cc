//////////////////////////////////////////////////////////////////////
// ScoreMatrix.cc
//////////////////////////////////////////////////////////////////////

#include <string.h>
#include "Assert.h"
#include "ScoreMatrix.h"
#include "ProbModel.h"

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////

ScoreMatrix::ScoreMatrix (int layers, int rows, int cols) : 
  layers(layers), rows(rows), cols(cols) {
  
  ASSERT (layers >= 0, "Number of layers in matrix must be positive.");
  ASSERT (rows >= 0, "Number of rows in matrix must be positive.");
  ASSERT (cols >= 0, "Number of columns in matrix must be positive.");
    
  data = new SCORE[layers * rows * cols];
  ASSERT (data, "Out of memory.");
}

//////////////////////////////////////////////////////////////////////
// Copy constructor
//////////////////////////////////////////////////////////////////////

ScoreMatrix::ScoreMatrix (const ScoreMatrix &m) :
  layers (m.layers), rows(m.rows), cols(m.cols) {
  
  ASSERT (layers >= 0, "Number of layers in matrix must be positive.");
  ASSERT (rows >= 0, "Number of rows in matrix must be positive.");
  ASSERT (cols >= 0, "Number of columns in matrix must be positive.");
    
  data = new SCORE[layers * rows * cols];
  ASSERT (data, "Out of memory.");

  memcpy (data, m.data, sizeof(SCORE) * (layers * rows * cols));
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////

ScoreMatrix::~ScoreMatrix (){
  delete[] data;
}

//////////////////////////////////////////////////////////////////////
// Fill all matrix values
//////////////////////////////////////////////////////////////////////

void ScoreMatrix::Fill (const SCORE &value){
  SCORE *p = data;
  SCORE *pEnd = data + layers * rows * cols;
  while (p != pEnd) *p++ = value;
}

//////////////////////////////////////////////////////////////////////
// Printing utility function
//////////////////////////////////////////////////////////////////////

void ScoreMatrix::PrintVal (FILE *file, const SCORE &value) const {
  if (value == LOG_ZERO_SCORE)
    fprintf (file, "      -inf");
  else
	  fprintf (file, "%10d", value);
 //   fprintf (file, "%7.3f", exp(TO_FLOAT(value)));
}

//////////////////////////////////////////////////////////////////////
// Print a single matrix layer
//////////////////////////////////////////////////////////////////////

void ScoreMatrix::PrintLayer (FILE *file, int layer) const {
	int i, j;
	for (j = 0; j < cols; j++)
		fprintf (file,"eee\t");
	fprintf(file, "\n");
  for (i = 0; i < rows; i++){
	  fprintf (file, "first");
    for (j = 0; j < cols; j++){
		if (j > 0) fprintf (file, "\t");
      PrintVal (file, operator()(layer,i,j));
    }
	fprintf(file, "\n");
  }
}

//////////////////////////////////////////////////////////////////////
// Print all matrix layers
//////////////////////////////////////////////////////////////////////

void ScoreMatrix::Print (FILE *file) const {
  for (int i = 0; i < layers; i++)
    PrintLayer (file, i);
}


void ScoreMatrix::PrintSumRange(FILE *file, int beginy, int endy, int beginx, int endx)
{
	ASSERT(beginx > 0 && endx < cols && beginy > 0 && endy < rows, "Out of range in PrintRange");
	fprintf(file, "    ");
	for (int k = beginx; k <= endx; k++)
		fprintf(file, "%9d", k);
	fprintf (file, "\n");
	for (int i = beginy; i <= endy; i++){
		fprintf(file, "%3d ", i);
    fprintf (file, "%s[", (i == beginy ? "[" : " "));
    for (int j = beginx; j <= endx; j++){
      if (j > beginx) fprintf (file, ", ");
	  fprintf (file, "%7.3f", exp(TO_FLOAT(operator()(MATCH,i,j))) - 
		  exp(TO_FLOAT(operator()(BEF_X,i,j))) - 
		  exp(TO_FLOAT(operator()(BEF_Y,i,j))) -
		  exp(TO_FLOAT(operator()(AFT_X,i,j))) -
		  exp(TO_FLOAT(operator()(AFT_Y,i,j)))/* -
		  exp(TO_FLOAT(operator()(INS_X,i,j))) -
		  exp(TO_FLOAT(operator()(INS_Y,i,j)))*/);
      
    }
    fprintf (file, "]%s\n", (i == endy ? "]" : ""));
  }
  fprintf (file, "\n");

}
