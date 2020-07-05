//////////////////////////////////////////////////////////////////////
// MultiSequence.cc
//////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "Assert.h"
#include "Utilities.h"
#include "MultiSequence.h"
#include "Types.h"

const int NUM_COLUMNS = 60;
const char *WHITE_SPACE_PLUS_GAPS = " \n\r\t\v\f-.";
const char *WHITE_SPACE = " \n\r\t\v\f";
const char *END_OF_LINE = "\n\r\f";

typedef Sequence *SequencePtr;
typedef char *charPtr;


//////////////////////////////////////////////////////////////////////
// Default constructor
//////////////////////////////////////////////////////////////////////

MultiSequence::MultiSequence() : sequences (NULL), numSequences (0) {}
  
//////////////////////////////////////////////////////////////////////
// Copy constructor
//////////////////////////////////////////////////////////////////////

MultiSequence::MultiSequence (const MultiSequence &rhs) : 
  sequences (NULL), numSequences (rhs.numSequences){
  
  if (rhs.sequences){
    sequences = new SequencePtr[numSequences];      
    ASSERT (sequences, "Out of memory.");
    
    for (int i = 0; i < numSequences; i++){
      sequences[i] = new Sequence (*rhs.sequences[i]);
      ASSERT (sequences[i], "Out of memory.");
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Assignment operator
//////////////////////////////////////////////////////////////////////

const MultiSequence& MultiSequence::operator= (const MultiSequence &rhs){
  
  if (this != &rhs){
    numSequences = rhs.numSequences;
    sequences = NULL;
    
    if (rhs.sequences){
      sequences = new SequencePtr[numSequences];      
      ASSERT (sequences, "Out of memory.");
      
      for (int i = 0; i < numSequences; i++){
	sequences[i] = new Sequence (*rhs.sequences[i]);
	ASSERT (sequences[i], "Out of memory.");
      }
    }
  }
  
  return *this;
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////

MultiSequence::~MultiSequence (){
  if (sequences){
    for (int i = 0; i < numSequences; i++)
      delete sequences[i];
    delete[] sequences;
  }
}

//////////////////////////////////////////////////////////////////////
// Return number of sequences
//////////////////////////////////////////////////////////////////////

const int MultiSequence::GetNumSequences() const {
  return numSequences;
}

//////////////////////////////////////////////////////////////////////
// Return length of first sequence
//////////////////////////////////////////////////////////////////////

const int MultiSequence::GetLength() const {
  ASSERT (numSequences > 0, "MultiSequence must have at least one sequence to retrieve length.");
  return sequences[0]->GetLength();
}

//////////////////////////////////////////////////////////////////////
// Rerieve sequence
//////////////////////////////////////////////////////////////////////
  
const Sequence &MultiSequence::GetSequence (int index) const {
  ASSERT (0 <= index && index < numSequences, "Invalid sequence index.");
  return *sequences[index];
}

//////////////////////////////////////////////////////////////////////
// Add new sequence
//////////////////////////////////////////////////////////////////////

void MultiSequence::AddSequence (Sequence *seq){
  Sequence **temp = new SequencePtr[numSequences+1];
  ASSERT (temp, "Out of memory.");
  
  if (sequences) memcpy (temp, sequences, sizeof(SequencePtr) * numSequences);
  temp[numSequences++] = seq;
  
  delete[] sequences;
  sequences = temp;
}

//////////////////////////////////////////////////////////////////////
// Auto-detect file format 
//
// 0 = MFA
// 1 = PILEUP
//////////////////////////////////////////////////////////////////////

const int MultiSequence::AutoDetectFileFormat (const char *filename) const {
  int fileType = -1;
  char *header;
  
  FILE *file = fopen (filename, "r");
  ASSERT (file, "Unable to open input file!");
  
  int length = GetData (file, header, END_OF_LINE, "");
  if (length){
    if (header[0] == '>'){
      fileType = 0;
    } else {
      if (!strncmp (header, "PileUp", 6))
	fileType = 1;
    }
  }
  
  fclose (file);
  
  return fileType;
}

//////////////////////////////////////////////////////////////////////
// Load sequences from MFA file
//////////////////////////////////////////////////////////////////////

void MultiSequence::LoadMFA (const char *filename, bool compressGaps){
  FILE *file = fopen (filename, "r");
  bool firstSequence = true;
  
  ASSERT (file, "Unable to open input file!");
  
  while (true){
    int length;
    
    // read MFA header
    
    char *name;
    length = GetData (file, name, END_OF_LINE, "");
    if (length == 0) break;
    
    if (firstSequence){
      ASSERT (name[0] == '>', "MFA sequence header should begin with '>'.");
      char *temp = new char[length];
      ASSERT (temp, "Out of memory.");
      
      memcpy (temp, name+1, sizeof(char) * length);
      delete[] name;
      name = temp;
      firstSequence = false;
    }
    
    // read MFA character data
    
    char *data;      
    length = GetData (file, data, ">", (compressGaps ? WHITE_SPACE_PLUS_GAPS : WHITE_SPACE));
    for (int i = 0; i < length; i++){
      if (data[i] == '.') data[i] = '-';
      ASSERT (('A' <= data[i] && data[i] <= 'Z') ||
	      ('a' <= data[i] && data[i] <= 'z') ||
	      (data[i] == '*' || data[i] == '-'),
	      "Unknown character encountered in MFA sequence data.");
    }
    
    // insert '@' at the beginning of the sequence
    
    char *temp = new char[length+2];
    ASSERT (temp, "Out of memory.");
    
    memcpy (temp+1, data, sizeof(char) * (length+1));
    temp[0] = '@';
    delete[] data;
    data = temp;	
    
    // add sequence
    
    Sequence *seq = new Sequence (data, name, length, numSequences);
    ASSERT (seq, "Out of memory.");
    
    AddSequence (seq);
  }
  
  ASSERT (!firstSequence, "No sequences read.");
  
  fclose (file);  
}

//////////////////////////////////////////////////////////////////////
// Load sequences from PILEUP file
//////////////////////////////////////////////////////////////////////

void MultiSequence::LoadPILEUP (const char *filename, bool compressGaps){
  FILE *file = fopen (filename, "r");
  ASSERT (file, "Unable to open input file!");
  int numRead = 0;
  
  // process header
  
  while (true){
    char *text;
    int length = GetData (file, text, END_OF_LINE, "");
    if (length == 0 && feof (file)) break;
    if (strstr (text, "//")) break;
    
    // parse sequence description 
    
    char *ptr = strstr (text, "Name:");
    if (ptr){
      int res;
      char *temp = new char[length];
      ASSERT (temp, "Out of memory.");
      
      res = sscanf (ptr + 5, "%s", temp);
      ASSERT (res != EOF && res > 0, "Failed to read sequence name in PILEUP file.");
      char *name = new char[strlen(temp)+1];
      ASSERT (name, "out of memory.");
      memcpy (name, temp, sizeof(char) * (strlen(temp)+1));
      for (int i = 0; i < numSequences; i++)
	ASSERT (strcmp (sequences[i]->GetName(), name), "Duplicate sequence name found.");
      
      
      ptr = strstr (text, "Len:");
      ASSERT (ptr, "Length field expected for sequence in PILEUP file.");
      sscanf (ptr + 4, "%d", &length);
      ASSERT (res != EOF && res > 0, "Failed to read sequence length in PILEUP file.");
      ASSERT (length > 0, "Length of sequence must be positive in PILEUP file.");
      
      delete[] temp;
      
      Sequence *seq = new Sequence (NULL, name, length, numSequences);
      ASSERT (seq, "Out of memory.");
      
      numRead++;
      AddSequence (seq);
    }
  }
  
  // prepare buffers for reading
  
  int *charsRead = new int[numRead];
  ASSERT (charsRead, "Out of memory.");
  char **data = new charPtr[numRead];
  ASSERT (data, "Out of memory.");
  for (int i = 0; i < numRead; i++){
    int index = i + numSequences - numRead;
    charsRead[i] = 0;
    data[i] = new char[sequences[index]->GetLength()+2];
    ASSERT (data[i], "Out of memory.");
    strcpy (data[i], "@");
  }
  
  // read sequences 
  
  while (true){
    char *text;
    int length = GetData (file, text, WHITE_SPACE, "");
    if (length == 0 && feof (file)) break;
    
    // search for sequence name
    
    int foundSequence = -1;
    for (int i = numSequences - numRead; i < numSequences; i++){
      if (!strcmp (text, sequences[i]->GetName())){
	foundSequence = i; 
	break;
      }
    }      
    delete[] text;
    if (foundSequence == -1) continue;
    
    // read sequence data
    
    length = GetData (file, text, END_OF_LINE, (compressGaps ? WHITE_SPACE_PLUS_GAPS : WHITE_SPACE));
    {for (int i = 0; i < length; i++){
      if (text[i] == '.') text[i] = '-';
      ASSERT (('A' <= text[i] && text[i] <= 'Z') ||
	      ('a' <= text[i] && text[i] <= 'z') ||
	      (text[i] == '*' || text[i] == '-'),
	      "Unknown character encountered in PILEUP sequence data.");
    }}
    
    // add to sequence data
    
    charsRead[foundSequence - numSequences + numRead] += length;
    ASSERT (charsRead[foundSequence - numSequences + numRead] <= sequences[foundSequence]->GetLength(),
	    "Sequence longer than reported length in PILEUP file.");
    strcat (data[foundSequence - numSequences + numRead], text);
  }
  
  {for (int i = numSequences - numRead; i < numSequences; i++){
    sequences[i]->SetData (data[i - numSequences + numRead]);
    if (!compressGaps){
      ASSERT (sequences[i]->GetLength() == charsRead[i],
	      "Actual sequence length inconsistent with reported length in PILEUP file");
    }
  }}
  
  ASSERT (numRead > 0, "No sequences read.");
  
  delete[] charsRead;
  delete[] data;
  
  fclose (file);  
}

//////////////////////////////////////////////////////////////////////
// Load data from file.
//////////////////////////////////////////////////////////////////////

void MultiSequence::LoadData (const char *filename, bool compressGaps){
  int fileFormat = AutoDetectFileFormat (filename);
  
  switch (fileFormat){
  case 0: LoadMFA (filename, compressGaps); break;
  case 1: LoadPILEUP (filename, compressGaps); break;
  default: ASSERT (false, "Unrecognized input file type.");
  }
}

//////////////////////////////////////////////////////////////////////
// Load sequences from file
//////////////////////////////////////////////////////////////////////

void MultiSequence::LoadSequences (const char *filename){
  LoadData (filename, true);
}

//////////////////////////////////////////////////////////////////////
// Load alignment from file
//////////////////////////////////////////////////////////////////////
  
void MultiSequence::LoadAlignment (const char *filename){
  LoadData (filename, false);
}

//////////////////////////////////////////////////////////////////////
// Print sequence in MFA format
//////////////////////////////////////////////////////////////////////
  
void MultiSequence::WriteMFA (FILE *file) const {
  
  for (int i = 0; i < numSequences; i++){
    fprintf (file, ">%s\n", sequences[i]->GetName());
    int length = sequences[i]->GetLength();
    const char *data = sequences[i]->GetData();
    
    for (int j = 1; j <= length; j++){
      fprintf (file, "%c", data[j]);
      if (j % NUM_COLUMNS == 0) fprintf (file, "\n");
    }
    if (length % NUM_COLUMNS != 0) fprintf (file, "\n");
  }
}

//////////////////////////////////////////////////////////////////////
// Compute CLUSTALW annotation character for alignment column
//////////////////////////////////////////////////////////////////////

const char MultiSequence::ComputeAnnotation (const char *data, const int size) const {
  static char *groups[47] = {
    
    // Identities
    "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", 
    "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "*",
    
    // Strong groups
    "STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW",
    
    // Weaker groups
    "CSA", "ATV", "SAG", "STNK", "STPA", "SGND", "SNDEQK", "NDEQHK", "NEQHRK", "FVLIM", "HFY"
  };
  
  for (int m = 0; m < 47; m++){
    bool isConserved = true;
    for (int j = 0; isConserved && j < size; j++)
      isConserved = strchr (groups[m], toupper(data[j]));
    if (isConserved){
      if (m < 27) return '*';
      if (m < 36) return ':';
      return '.';
    }
  }
  
  return ' ';
}

//////////////////////////////////////////////////////////////////////
// Print sequence in CLUSTALW format
//////////////////////////////////////////////////////////////////////
  
void MultiSequence::WriteCLUSTALW (FILE *file) const {
  
//  fprintf (file, "PROBCONS version %s multiple sequence alignment\n\n", "2.0");
  
  if (numSequences == 0) return;
  
  // Get sequence length and length of longest sequence name
  
  int length = sequences[0]->GetLength();
  int nameLength = strlen(sequences[0]->GetName());
  for (int i = 1; i < numSequences; i++){
    ASSERT (sequences[i]->GetLength() == length, 
	    "ERROR: Sequences of unequal length in CLUSTALW output.");
    nameLength = max (nameLength, (int) strlen(sequences[i]->GetName()));
  }
  
  // Print out sequences
  
  char *buffer = new char[numSequences];
  ASSERT (buffer, "Out of memory.");
  
  {for (int i = 1; i <= length; i += NUM_COLUMNS){
    for (int j = 0; j < numSequences; j++){
      fprintf (file, "%*s    ", nameLength, sequences[j]->GetName());
      const char *data = sequences[j]->GetData();
      for (int k = i; k <= min (i + NUM_COLUMNS - 1, length); k++)
	fprintf (file, "%c", data[k]);
      fprintf (file, "\n");
    }
    
    // Compute annotation line
    
    fprintf (file, "%*s    ", nameLength, "");
    for (int k = i; k <= min (i + NUM_COLUMNS - 1, length); k++){
      for (int j = 0; j < numSequences; j++)
	buffer[j] = sequences[j]->GetData()[k];
      char ch = ComputeAnnotation (buffer, numSequences);
      fprintf (file, "%c", ch);
    }
    fprintf (file, "\n");
    
    if (i + NUM_COLUMNS <= length)
      fprintf (file, "\n");
  }}
}

//////////////////////////////////////////////////////////////////////
// Sort sequences by ID number
//////////////////////////////////////////////////////////////////////

void MultiSequence::Sort(){
  for (int i = 0; i < numSequences; i++){
    for (int j = i+1; j < numSequences; j++){
      if (sequences[i]->GetID() > sequences[j]->GetID()){
	Sequence *temp = sequences[i];
	sequences[i] = sequences[j];
	sequences[j] = temp;
      }
    }
  }
}


//////////////////////////////////////////////////////////////////////
// Retrieve pointer to a sequence
//////////////////////////////////////////////////////////////////////
  

Sequence * MultiSequence::GetSequencePtr(int index)
{
  ASSERT (0 <= index && index < numSequences, "Invalid sequence index.");
  return sequences[index];
}


void MultiSequence::AddAlignPosition(Fragment * frag)
{
	sequences[frag->id]->AddAlignPosition(frag->begin,frag->end);
}

////////////////////////////////////////////////////////////////////////
// Return block of all aligned sequences
////////////////////////////////////////////////////////////////////////
typedef const char * pchar;
typedef int * pint;
void MultiSequence::FindBlock(Block &block, int &start, int &end, int minlength)
{
	int i,j;
	if (numSequences <=0) {start = end = 0;return;}

	start = end = -1;
	int length = sequences[0]->GetLength();

	const char **data = new pchar[numSequences];
	for ( i = 0; i < numSequences; i++)
		data[i] = sequences[i]->GetData();
	
	for( i = 1; i <= length; i++){
		for ( j = 0; j < numSequences && data[j][i] != '-';j++);
		if(j == numSequences){
			start = i;
			break;
		}
	}
	for( i = length; i > 0; i--){
		for ( j = 0; j < numSequences && data[j][i] != '-';j++);
		if(j == numSequences){
			end = i;
			break;
		}
	}

	int numFrag = numSequences ;//gradually decreases number of fragments

		int **map = new pint[numSequences];
		for (i = 0; i < numSequences; i++)
			map[i] = new int[length+1];
		for (i = 0; i < numSequences; i++){
			int count = 0;
			for (j = 1; j <= length; j++){
				if(data[i][j] != '-') count++;
				map[i][j] = count;
			}
		}
		int ml = end - start +1;
		for(i=0; i < numSequences && ml >= minlength; i++)
			if (ml > map[i][end] - map[i][start]+1) ml = map[i][end] - map[i][start]+1;
		
		if (ml < minlength){
			int *counts = new int[length+1];
			for (i = 1 ; i <= length; i++) counts[i] = 0;
			for (j = 0; j < numSequences; j++){
				for (i = 1 ; i <= length; i++)			
					if(data[j][i] != '-') counts[i]++;

			}
			for (numFrag = numSequences - 1; numFrag >=2 ; numFrag--){
				for (i = 1 ; i <= length - minlength && counts[i] < numFrag; i++);
				start = i;
				if(start > length - minlength) continue;
				for (i = length; i > start && counts[i] < numFrag; i--);
				end = i;
			
				int good = 0;
				for(i=0; i < numSequences; i++)
					if (map[i][end] - map[i][start]+1 >= minlength) good++;
				if (good >= numFrag){
					int match,bindex,eindex,flag = 0;;	
					for (bindex = start; bindex < end-minlength && flag ==0;bindex++){
						for(eindex = end; eindex > bindex+minlength;eindex--){
							for(match=0,j=0;j<numSequences;j++){
								if(data[j][bindex] != '-' && data[j][eindex] != '-' &&
									map[j][eindex] - map[j][bindex]+1 >= minlength) 
									match++;
							}
							if(match == numFrag && eindex - bindex >= minlength){
								start = bindex; end = eindex; flag = 1;
								break;
							}
						}
					}
					if(flag) break;
				}
			}

		}
		
	

	
	if(numFrag >= 2 && end-start+1 >= minlength && data[0][start] != '-' && data[0][end] != '-' &&
			map[0][end] - map[0][start]+1 >= minlength){
		int b,e;
	Block tmp;
	int ok = 1;
	tmp.seed = block.seed;
	for (i = 0,j=0; i < numSequences; i++){
		if(data[i][start] == '-' || data[i][end] == '-' ||
			map[i][end] - map[i][start]+1 < minlength) continue;
		for(b = block[i].begin, j = 0; j < start; j++)
			if (data[i][j] != '-') b++;
		for (e = b,j = start; j < end; j++)
			if (data[i][j] != '-') e++;
		Fragment fr(b-1,e-1,0,block[i].id);//Not the original ID
		tmp.AddFragment(fr);
	}
	if(numFrag < numSequences){ //remove some sequences
		Sequence **temp = new SequencePtr[numFrag];
		ASSERT (temp, "Out of memory.");
		for (i = 0,j=0; i < numSequences; i++){
			if(data[i][start] == '-' || data[i][end] == '-'  ||
				map[i][end] - map[i][start]+1 < minlength)
				delete sequences[i];
			else
				temp[j++] = sequences[i];
		}
		delete[] sequences;
		sequences = temp;
		numSequences = j;
	}
	if(ok)
		block = tmp;
	else
		start = end = 0;
	}
	else{
		start = end = 0;
	}
	
	if(start == 0) block.part = 1;
	
	for (i = 0; i < numSequences; i++)
		delete map[i];
	delete map;
	delete data;
	
}


//////////////////////////////////////////////////////////////////////////
//  Output from start to end
//////////////////////////////////////////////////////////////////////////
void MultiSequence::WriteCLUSTALW(FILE *file, int start, int end)
{
  if (numSequences == 0) return;
  
  // Get sequence length and length of longest sequence name
  
  int length = sequences[0]->GetLength();
  int nameLength = strlen(sequences[0]->GetName());
  for (int i = 1; i < numSequences; i++){
    ASSERT (sequences[i]->GetLength() == length, 
	    "ERROR: Sequences of unequal length in CLUSTALW output.");
    nameLength = max (nameLength, (int) strlen(sequences[i]->GetName()));
  }
  
  // Print out sequences
  
  char *buffer = new char[numSequences];
  ASSERT (buffer, "Out of memory.");
  
  {for (int i = start; i <= end; i += NUM_COLUMNS){
    for (int j = 0; j < numSequences; j++){
      fprintf (file, "%*s    ", nameLength, sequences[j]->GetName());
      const char *data = sequences[j]->GetData();
      for (int k = i; k <= min (i + NUM_COLUMNS - 1, end); k++)
	fprintf (file, "%c", data[k]);
      fprintf (file, "\n");
    }
    
    // Compute annotation line
    
    fprintf (file, "%*s    ", nameLength, "");
    for (int k = i; k <= min (i + NUM_COLUMNS - 1, end); k++){
      for (int j = 0; j < numSequences; j++)
	buffer[j] = sequences[j]->GetData()[k];
      char ch = ComputeAnnotation (buffer, numSequences);
      fprintf (file, "%c", ch);
    }
    fprintf (file, "\n");
    
    if (i + NUM_COLUMNS <= length)
      fprintf (file, "\n");
  }}

}

void MultiSequence::ClearAlignPosition()
{
	for (int i = 0; i < numSequences; i++)
		sequences[i]->ClearAlignPosition();
}



void MultiSequence::WriteFASTA(FILE *file, Block *block, MultiSequence *result, int start, int end)
{
	int i,j;
	if(block->size() > numSequences){
		fprintf(stderr,"\nBlock contains repeats\n\n");
		exit(1);
	}
	int *map = new int[numSequences];
	for (i = 0 ; i < numSequences; i++)
		map[i] = -1;
	for (i =0; i < block->size(); i++)
		map[(*block)[i].id] = i;

	int maxBegin = 0;
	int maxEnd = 0;
	for (i = 0; i < block->size(); i++){
		if(maxBegin < (*block)[i].begin)
			maxBegin = (*block)[i].begin;
	}

	for(i=0;i<numSequences;i++){
		int len;
		if(map[i] != -1)
			len = maxBegin + end-start+sequences[i]->GetLength()-(*block)[map[i]].end;
		else
			len = sequences[i]->GetLength()+end-start+1;
		if(maxEnd < len)
			maxEnd = len;
	}

	for (i = 0; i < numSequences; i++){
		fprintf(file,">%s\n",sequences[i]->GetName());
		int cur = map[i];
		const char *data = sequences[i]->GetData();
		int count = 0;
		if(cur != -1){
			for (j = 1; j < (*block)[cur].begin; j++,count++)
				fprintf(file,"%c",data[j]+'a'-'A');
			for (j = (*block)[cur].begin; j < maxBegin; j++,count++)
				fprintf(file,".");
			const char *alignData = result->GetSequence(cur).GetData();
			for (j = start; j <= end; j++,count++)
				fprintf(file,"%c",alignData[j]);
			for (j = (*block)[cur].end+1; j <= sequences[i]->GetLength();j++,count++)
				fprintf(file,"%c",data[j]+'a'-'A');
		}
		else{
			int beginGap = min (maxBegin, sequences[i]->GetLength()+1);
			for (j = 1; j <beginGap;j++,count++)
				fprintf(file,"%c",data[j]+'a'-'A');
			if (beginGap < sequences[i]->GetLength()+1){
				for (j = start; j <= end; j++,count++)
					fprintf(file,"-");
				for (j = beginGap; j <= sequences[i]->GetLength(); j++,count++)
					fprintf (file,"%c",data[j]+'a'-'A');
			}
		}
		for (; count < maxEnd; count++)
			fprintf(file,".");
		fprintf(file,"\n");
	}
	fprintf(file,"#\n");
	delete map;
}
