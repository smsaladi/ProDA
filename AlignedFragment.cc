// AlignedFragment.cpp: implementation of the AlignedFragment class.
//
//////////////////////////////////////////////////////////////////////
using namespace std;
#include <cstdio>
#include <cstring>

#include "AlignedFragment.h"
#include "Assert.h"
#include "Utilities.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

AlignedFragment::AlignedFragment()
{

}

AlignedFragment::~AlignedFragment()
{
	if(seq[0] != NULL) delete seq[0];
	if(seq[1] != NULL) delete seq[1];
}


AlignedFragment::AlignedFragment(int d1, int d2, int beg1, int beg2, int e1, int e2, int *s1, int *s2)
{
	id[0] = d1; id[1] = d2;
	begin[0] = beg1; begin[1] = beg2;
	end[0] = e1; end[1] = e2;
	int *p[2];
	p[0] = s1;p[1] = s2;
	for (int k = 0; k < 2; k++){
		if (begin[k] >= end[k]) {
			begin[k] = -1;
			end[k] = begin[k] - 1;
			seq[k] = NULL;
			continue;
		}
		seq[k] = new int[end[k] - begin[k] + 1];
		ASSERT(seq[k],"Not enough memory");
		for (int i = 0; i <= end[k] - begin[k]; i++)
			seq[k][i] = p[k][i];
	}
}

AlignedFragment::AlignedFragment(const AlignedFragment& af)
{
	for (int k = 0; k < 2; k++){
		id [k] = af.id[k];
		begin[k] = af.begin[k];
		end[k] = af.end[k];
		if(begin[k]>=end[k]){
			seq[k]=NULL;
			continue;
		}
		seq[k] = new int[end[k] - begin[k] + 1];
		ASSERT(seq[k],"Out of memory");
		for (int i = 0; i <= end[k] - begin[k]; i++)
			seq[k][i] = af.seq[k][i];
	}
	similarity = af.similarity;
}

AlignedFragment& AlignedFragment::operator =(const AlignedFragment af){
	for (int k = 0; k < 2; k++){
		id [k] = af.id[k];
		begin[k] = af.begin[k];
		end[k] = af.end[k];
		if(begin[k]>=end[k]){
			seq[k]=NULL;
			continue;
		}
		seq[k] = new int[end[k] - begin[k] + 1];
		ASSERT(seq[k],"Not enough memory");
		for (int i = 0; i <= end[k] - begin[k]; i++)
			seq[k][i] = af.seq[k][i];
	}
	similarity = af.similarity;
	return *this;
}

////////////////////////////////////////////////////////////////////
// Returns the position aligned to the given  position
////////////////////////////////////////////////////////////////////

std::pair<int, int> *AlignedFragment::GetAlignPos(int sequence, int pos, int &second)
{
	int k;
	for (k = 0; k < 2; k++){
		if(sequence != id[k]) continue;
		if (pos >= begin[k] && pos <= end[k]) {
			while (pos >= begin[k] && seq[k][pos-begin[k]] == -1) pos--;
			second = 1-k;
			if (pos < begin[k]) 
				return NULL;
			return new std::pair<int, int>(id[1-k], seq[k][pos-begin[k]]);
		}
	}
	return NULL;
}


//////////////////////////////////////////////////////////////////
// Truncats to the overlap with another fragment
//////////////////////////////////////////////////////////////////

int Fragment::Overlap(Fragment &fr)
{
	if (id != fr.id) return -1;
	int b = begin > fr.begin ? begin : fr.begin;
	int e = begin < fr.begin ? begin : fr.begin;
	
	if (e - b > 0){
		begin = b;
		end = e;
		length = e - b + 1;
	}
	return e - b + 1;
}

////////////////////////////////////////////////////////////////////////////////
//  Erases a fragment from AlignedFrament, returns two shorter ones
////////////////////////////////////////////////////////////////////////////////

void AlignedFragment::Adjust(Fragment &fr1, Fragment &fr2, AlignedFragment &rfr1, AlignedFragment &rfr2)
{
	ASSERT (fr1.id == id[0] || fr1.id == id[1], "AlignedFragment Adjusting fault");
	ASSERT (fr2.id == id[0] || fr2.id == id[1], "AlignedFragment Adjusting fault");
	if (fr1.id == id[0] && Overlap(fr1.begin,fr1.end,begin[0],end[0]) > 1){
		rfr1 = *SubFragment(begin[0],fr1.begin-1,begin[1],fr2.begin-1);
		rfr2 = *SubFragment(fr1.end+1,end[0],fr2.end+1,end[1]);
	}
	else{
		rfr1 = *SubFragment(begin[0],fr2.begin-1,begin[1],fr1.begin-1);
		rfr2 = *SubFragment(fr2.end+1,end[0],fr1.end+1,end[1]);
	}

}

Fragment::Fragment()
{

}

int AlignedFragment::GetLength()
{
	if(begin[0] <=0 || end[0] <=0 || begin[1] <=0 || end[1] <=0) return 0;
	if(begin[0] == end[0] || begin[1] == end[1]) return 0;
	int i,k,res[2];
	for(k=0; k<2;k++){
		for (res[k] = 0, i = 0; i <= end[k] - begin[k];i++)
			if(seq[k][i] > 0) res[k]++;
	}
	int result = res[0] < res[1] ? res[0] : res[1];
 	if (result <=1) return 0;
 	return result;
}

int AlignedFragment::GetID(int i)
{
	return id[i];
}

int AlignedFragment::GetBegin(int i)
{
	return begin[i];
}

int AlignedFragment::GetEnd(int i)
{
	return end[i];
}


///////////////////////////////////////////////////////////////////////
// Prune unaligned positions at the two ends
///////////////////////////////////////////////////////////////////////
void AlignedFragment::Prune()
{
	if(end[0]<=begin[0] || end[1]<=begin[1]) return;
	int i,b[2],e[2],k;

	for(k=0;k<2;k++){
		for(i=end[k]-begin[k],e[k]=0;i>=0 && seq[k][i] < 0;i--) e[k]++;
		for(i=0,b[k]=0;i<=end[k]-begin[k] && seq[k][i] < 0;i++) b[k]++;
	}
	for(k=0;k<2;k++){
		if(b[k]+e[k]>0){
			int *s = new int[end[k]-begin[k]+1-b[k]-e[k]];
			memcpy(s,seq[k]+b[k],sizeof(int)*(end[k]-begin[k]+1-b[k]-e[k]));
			delete seq[k];
			seq[k]=s;
			begin[k] += b[k];
			end[k] -= e[k];
		}
	}
}

Fragment * AlignedFragment::GetFragment(int i)
{
	Fragment *fr = new Fragment(begin[i],end[i],0,id[i]);
	return fr;
}

int AlignedFragment::ProcessRepeat(AVECT &fragments, int minlength)
{
	if(id[0] != id[1]) return 0;
	int i;
	if(min(end[0],end[1]) >= max(begin[0],begin[1])){//Ovelapping
		int bound[2][4];
		bound[0][0] = begin[0];bound[0][3] = end[0];
		bound[1][0] = begin[1]; bound[1][3] = end[1];
		for(int k =0; k <2; k++){
		if(begin[k] < begin[1-k]){
			for(i=begin[1-k]-begin[k];i>0 && seq[k][i]==-1;i--);
			bound[1-k][1] = seq[k][i]+bound[1-k][0];
			bound[k][1] = i+begin[k];

			for(i=end[k]-begin[1-k];i>0 && seq[1-k][i] == -1; i--);
			bound[1-k][2]=i+begin[1-k];
			bound[k][2] = seq[1-k][i]+bound[k][0];
		}
		}
		if(bound[0][1] >bound[0][2]){
			int tmp;
			tmp = bound[0][1];
			bound[0][1] = bound[0][2];
			bound[0][2] = tmp;
			tmp = bound[1][1];
			bound[1][1] = bound[1][2];
			bound[1][2] = tmp;
		}
		else{
			bound[0][1]--;bound[1][1]--;
			bound[0][2]++;bound[1][2]++;
		}
		for(i=0;i<3;i+=2){
			AlignedFragment * sub = this->SubFragment(bound[0][i],bound[0][i+1],
				bound[1][i],bound[1][i+1]);
			if(sub != NULL && sub->GetLength() >= minlength){
					fragments.push_back(*sub);
					sub->Print(stderr);
			}
			if(sub) delete sub;
		}
		return 1;
	}
	else{
		Print(stderr);
		fragments.push_back(*this);
		return 0;
	}
}

AlignedFragment * AlignedFragment::SubFragment(int begin0, int end0, int begin1, int end1)
{
	if(begin0 < begin[0] || end0 > end[0] || begin0 >= end0 ||
		begin1 < begin[1] || end1 > end[1] || begin1 >= end1)
		return new AlignedFragment(0,0,0,0,0,0,0,0);
	int b0,e0;
	for(b0 = begin0; b0 < end0 && seq[0][b0-begin[0]] <begin1; b0++);
	for(e0 = b0; e0 < end0 && seq[0][e0-begin[0]] < end1; e0++);

	if (b0 >= e0) new AlignedFragment(0,0,0,0,0,0,0,0);
	AlignedFragment *res = new AlignedFragment(id[0],id[1],b0,seq[0][b0-begin[0]],e0,seq[0][e0-begin[0]],
								seq[0]+b0-begin[0],seq[1]+seq[0][b0-begin[0]]-begin[1]);
	res->Prune();
	return res;
}

void AlignedFragment::Print(FILE *file)
{
	fprintf(file,"<%d %d %d><%d %d %d>\n",id[0],begin[0],end[0],id[1],begin[1],end[1]);
}

Fragment * AlignedFragment::GetAlignFragment(Fragment &fr)
{
	Fragment *res;
	int second1, second2;
	std::pair<int, int> *start,*finish;
	start = GetAlignPos(fr.id,fr.begin,second1);
	finish = GetAlignPos(fr.id,fr.end,second2);
	if(start != NULL && finish != NULL && second1 == second2)
		res = new Fragment(start->second,finish->second,0,start->first);
	else 
		res = NULL;
	if(start) delete start;
	if(finish) delete finish;
	return res;
}


void AlignedFragment::ShiftRight(int offset1, int offset2)
{
	int i,j;
	begin[0] += offset1 - 1;
	end[0] += offset1 - 1;
	for(i = 0, j = begin[0]; j <= end[0]; i++,j++)
		if(seq[0][i] >=0 )
			seq[0][i] += offset2 - 1;
	begin[1] += offset2 - 1;
	end[1] += offset2 - 1;
	for( i = 0, j = begin[1]; j <= end[1]; i++,j++)
		if(seq[1][i] >= 0)
			seq[1][i]+= offset1-1;
}
