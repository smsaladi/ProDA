//////////////////////////////////////////////////////////////////////
// Consistency.h
//
// Routines for probabilistic consistency.
//////////////////////////////////////////////////////////////////////

#ifndef CONSISTENCY_H
#define CONSISTENCY_H

#include "SparseMatrix.h"
#include "Matrix.h"

void AccumulateConsistencyInfo (const SparseMatrix **posteriors, Matrix &p, int x, int y, int z, int n);
SparseMatrix **ProbabilisticConsistency (SparseMatrix **posteriors, int n);

#endif
