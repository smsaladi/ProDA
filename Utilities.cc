//////////////////////////////////////////////////////////////////////
// Utilities.cc
//////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <string.h>
#include "Assert.h"
#include "Utilities.h"

//////////////////////////////////////////////////////////////////////
// Read data from file.
//   
// This routine will read characters from a file object and 
// store them in the resulting buffer.  
//
// -- characters from the "terminatingChars" string are used
//    to signify when reading should stop; these characters
//    are not included in the buffer
//     
// -- the "skipChars" string denotes any other characters 
//    that should be skipped
//
// -- a NULL character is appended to the end of the read
//    string
//	
// The result returned is the length of the read buffer,
// excluding the NULL character.  If a character appears in
// both the terminating characters and the skipped characters
// strings, the former takes precedence.
//////////////////////////////////////////////////////////////////////

int GetData (FILE *file, char *&buffer, 
	     const char *terminatingChars,
	     const char *skipChars){
  
  bool isTerm[256];
  bool isSkip[256];
  
  int length = 0, capacity = 1;
  char *temp = new char[capacity];
  char ch;
  
  ASSERT (temp, "Out of memory.");
  
  // precompute character detection flags
  
  for (int i = 0; i < 256; i++) isTerm[i] = isSkip[i] = false;
  {for (int i = strlen(terminatingChars) - 1; i >= 0; i--) 
	  isTerm[(unsigned char) terminatingChars[i]] = true;}
  {for (int i = strlen(skipChars) - 1; i >= 0; i--)
	  isSkip[(unsigned char) skipChars[i]] = true;}

  // read buffer
  
  while ((ch = fgetc (file)) != EOF){
    
    if (isTerm[(unsigned char) ch]) break;
    if (isSkip[(unsigned char) ch]) continue;

    if (length == capacity){
      buffer = new char[capacity *= 2];
      ASSERT (buffer, "Out of memory.");
      
      memcpy (buffer, temp, sizeof(char) * length);
      delete[] temp;
      temp = buffer;	
    }
    
    temp[length++] = ch;
  }
  
  // trim buffer to correct length
  
  buffer = new char[length+1];
  ASSERT (buffer, "Out of memory.");
  
  memcpy (buffer, temp, sizeof(char) * length);
  buffer[length] = '\0';
  delete[] temp;

  return length;
}

//////////////////////////////////////////////////////////////////////
// Duplicate string
//////////////////////////////////////////////////////////////////////

char *StrDup (const char *s){
  int len = strlen(s);
  char *ret = new char[len+1];
  ASSERT (ret, "Out of memory.");
  memcpy (ret, s, len+1);
  return ret;
}

//////////////////////////////////////////////////////////////////////
// Substring
//////////////////////////////////////////////////////////////////////

char *SubString (const char *s, int i, int j){
  ASSERT (i >= 0 && i <= (int) strlen(s), "Invalid index.");
  ASSERT (j >= i && j <= (int) strlen(s), "Invalid index.");
  char *ret = new char[j - i + 1];
  ASSERT (ret, "Out of memory.");
  memcpy (ret, s + i, j - i);
  ret[j - i] = '\0';
  return ret;
}

///////////////////////////////////////////////////////////////////////////////
// Returns overlap length of two intervals
///////////////////////////////////////////////////////////////////////////////
int Overlap(int b1, int e1, int b2, int e2)
{
	int b = b1 > b2 ? b1 : b2;
	int e = e1 < e2 ? e1 : e2;
	return e - b + 1;
}


