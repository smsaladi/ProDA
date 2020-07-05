//////////////////////////////////////////////////////////////////////
// Utilities.h
//
// Miscellaneous utility routines for ProDA.
//////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdio.h>

// reading character buffer of arbitrary size from a file

int GetData (FILE *file, char *&buffer, 
	     const char *terminatingChars,
	     const char *skipChars);

// duplicate string

char *StrDup (const char *s);

char *SubString (const char *s, int i, int j);

// math utility functions

inline int    min (int    a, int    b){ if (a < b) return a; return b; }
inline float  min (float  a, float  b){ if (a < b) return a; return b; }
inline double min (double a, double b){ if (a < b) return a; return b; }
inline int    max (int    a, int    b){ if (a > b) return a; return b; }
inline float  max (float  a, float  b){ if (a > b) return a; return b; }
inline double max (double a, double b){ if (a > b) return a; return b; }
inline void swap(int &a, int &b){int tmp = a;a = b; b= tmp;}

//
int Overlap(int b1, int e1, int b2, int e2);

#endif
