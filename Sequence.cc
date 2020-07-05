//////////////////////////////////////////////////////////////////////
// Sequence.cc
//////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <string.h>
#include <map>
#include "Sequence.h"
#include "Assert.h"
#include "AlignedFragment.h"

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////

Sequence::Sequence (char *data, char *name, int length, int id) :
  data (data), name (name), length (length), id (id) {
  align = new int[length+1];
  ASSERT (align, "Out of memory.");
  memset(align, 0, (length + 1)*sizeof(int));
  position = new int[length+1];
  ASSERT (align, "Out of memory.");
  for (int i = 0; i < length+1; i++)
		position[i] = i;
  }

//////////////////////////////////////////////////////////////////////
// Copy constructor
//////////////////////////////////////////////////////////////////////

Sequence::Sequence (const Sequence &rhs) : 
  data (NULL), name (NULL), length (rhs.length), id (rhs.id) {

  if (length > 0){
    data = new char[length+2];
    ASSERT (data, "Out of memory.");
    ASSERT (length + 1 == (int) strlen(rhs.data), "Sequence of incorrect length.");
    memcpy (data, rhs.data, sizeof(char) * (length+2));
	align = new int[length+1];
	ASSERT (align, "Out of memory.");
	memcpy (align, rhs.align, sizeof(int) * (length+1));
	position = new int[length+1];
	ASSERT (position, "Out of memory.");
	memcpy (position, rhs.position, sizeof(int) * (length+1));
  }
  
  if (rhs.name){
    name = new char[strlen(rhs.name)+1];
    ASSERT (name, "Out of memory.");
    memcpy (name, rhs.name, sizeof(char) * (strlen(rhs.name)+1));
  }
}
  
//////////////////////////////////////////////////////////////////////
// Assignment operator
//////////////////////////////////////////////////////////////////////

const Sequence& Sequence::operator= (const Sequence &rhs){
  if (this != &rhs){
    delete[] data; data = NULL;
    delete[] name; name = NULL;
	delete[] align;
	delete[] position;
    
    length = rhs.length;
    id = rhs.id;
    
    if (length > 0){
      data = new char[length+2];
      ASSERT (data, "Out of memory.");
      memcpy (data, rhs.data, sizeof(char) * (length+2));
	  align = new int[length+1];
	  ASSERT (align, "Out of memory.");
	  memcpy (align, rhs.align, sizeof(int) * (length+1));
	  position = new int[length+1];
	  ASSERT (position, "Out of memory.");
	  memcpy (position, rhs.position, sizeof(int) * (length+1));
    } 
    
    if (rhs.name){
      name = new char[strlen(rhs.name)+1];
      ASSERT (name, "Out of memory.");
      memcpy (name, rhs.name, sizeof(char) * (strlen(rhs.name)+1));
    }      
  }
  
  return *this;
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////

Sequence::~Sequence (){
  delete[] data;
  delete[] name;
  delete[] align;
  delete[] position;
}

//////////////////////////////////////////////////////////////////////
// Retrieve sequence data
//////////////////////////////////////////////////////////////////////

const char *Sequence::GetData () const {
  return data;
}

//////////////////////////////////////////////////////////////////////
// Retrieve sequence name 
//////////////////////////////////////////////////////////////////////

const char *Sequence::GetName () const {
  return name;
}

//////////////////////////////////////////////////////////////////////
// Retrieve sequence length 
//////////////////////////////////////////////////////////////////////

const int Sequence::GetLength () const {
  return length;
}

//////////////////////////////////////////////////////////////////////
// Retrieve sequence id 
//////////////////////////////////////////////////////////////////////

const int Sequence::GetID () const {
  return id;
}

//////////////////////////////////////////////////////////////////////
// Store sequence data 
//////////////////////////////////////////////////////////////////////

void Sequence::SetData (char *data) {
  ASSERT (data[0] == '@', "Invalid sequence data format.");
  this->data = data;
  this->length = strlen(data) - 1;
}

//////////////////////////////////////////////////////////////////////
// Compute mapping from letter to positions in sequence
//////////////////////////////////////////////////////////////////////

int *Sequence::ComputeMapping() const {
  int numLetters = 0;

  for (int i = 1; i <= length; i++)
    numLetters += (data[i] != '-');

  int *ret = new int[numLetters+1];
  ASSERT (ret, "Out of memory.");
  
  int j = 0;
  {for (int i = 0; i <= length; i++)
    if (i == 0 || data[i] != '-')
		ret[j++] = i;}

  return ret;
}

//////////////////////////////////////////////////////////////////////
// Increases the number of sequences aligned to begin-end fragment
//////////////////////////////////////////////////////////////////////

void Sequence::AddAlignPosition(int begin, int end)
{
	for (int i = begin; i <= end; i++)
		align[i]++;
}

//////////////////////////////////////////////////////////////////////
// Erase a frament for next local alignment
//////////////////////////////////////////////////////////////////////

void Sequence::EraseFragment(int begin, int end)
{
	int len = end - begin + 1;
	if ( len >= length) return;
	for (int i = begin; i+len<length; i++){
		data[i] = data[i+len];
		position[i] = position[i+len];
	}
	length -= len;
}

//////////////////////////////////////////////////////////////////////
// Return position in the origianl sequence
//////////////////////////////////////////////////////////////////////

int Sequence::OriginPosition(int current)
{
	return position[current];
}


void Sequence::EraseCluster(int start, int end, int n)
{
	for (int i = start; i <= end; i++)
		align[i] -= n;
}

void Sequence::PrintAlign(FILE *f)
{
	for (int i = 1; i <= length; i++)
		fprintf(f,"%d ", align[i]);
	fprintf(f,"\n");
}

////////////////////////////////////////////////////////////////////////////////////
// Clips the sequence
////////////////////////////////////////////////////////////////////////////////////

void Sequence::Clip(int begin, int end)
{
	length = end-begin+1;
	char *tmp = new char[length+2];
	tmp[0] = '@';
	memcpy (tmp+1, data+begin, sizeof(char)*length);
	tmp[length+1] = 0;
	delete data;
	data = tmp;
}

int Sequence::SetID(int newid)
{
	int oldid = id;
	id = newid;
	return oldid;
}

///////////////////////////////////////////////////////////////////////////////
// Shrink the sequence to substring from begin to end inclusively
///////////////////////////////////////////////////////////////////////////////

void Sequence::SubStr(int begin, int end)
{
	ASSERT ( begin > 0 && end <= length, "Wrong substring index");
	length = end - begin + 1;
	char *d = new char[length+2];
	memcpy (d+1, data+begin, sizeof(char)*(length+1));
	delete data;
	data = d;
	data[0] = '@';data[length+1] = 0;
	int *p = new int [length + 1];
	memcpy (p+1, position + begin,sizeof(int)*(length));
	delete position;
	position = p;
	int *a = new int [length + 1];
	memcpy (a, align + begin,sizeof(int)*(length+1));
	delete align;
	align = a;
	
}

int Sequence::GetAlign(int i) const
{
	if (i <0 || i >length) return -1;
	return align[i];
}

void Sequence::ClearAlignPosition()
{
	for (int i = 0 ; i <= length; i++)
		align[i] = 0;
}

