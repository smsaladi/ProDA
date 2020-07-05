//////////////////////////////////////////////////////////////////////
// Types.h
//
// Data types for ProDA
//////////////////////////////////////////////////////////////////////

#ifndef TYPES_H
#define TYPES_H

#include <vector>
#include <list>
#include "Matrix.h"
#include "ScoreMatrix.h"

typedef char *string;
typedef std::pair<Fragment, Fragment> FPAIR;
typedef std::vector<FPAIR> PVECT ;
typedef std::pair<int , int> PAIRI;
typedef std::vector<PAIRI> VECT;
typedef std::list<Sequence> SEQLIST;

typedef std::vector<int> IVECT;

typedef std::pair<ScoreMatrix *, Matrix *> SCORE_PAIR;


#endif
