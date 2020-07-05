//////////////////////////////////////////////////////////////////////
// Assert.h
//
// Extension of C assert() that allows for error messages.
//////////////////////////////////////////////////////////////////////

#ifndef ASSERT_H
#define ASSERT_H

int _ASSERT_FAILED (char *filename, int line_number, const char *error_msg);

#ifdef NDEBUG
#define ASSERT(test,error_msg)
#else
#define ASSERT(test,error_msg) (test ? 0 : _ASSERT_FAILED(__FILE__, __LINE__, error_msg))
#endif
  
#endif
