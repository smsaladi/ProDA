// Block.cpp: implementation of the Block class.
//
//////////////////////////////////////////////////////////////////////
#include <cstdlib>

#include "Block.h"
#include "AlignedFragment.h"
#include "Sequence.h"
#include "MultiSequence.h"
#include "Utilities.h"
#include "Types.h"

extern int MINLENGTH;
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Block::Block()
{
	part = 0;
}

Block::~Block()
{
	frags.clear();
}

int Block::size()
{
	return frags.size();
}

Fragment& Block::operator [](int i){
	return frags.at(i);
}

////////////////////////////////////////////////////////////////////////
// Given aligned fragments find a block of maximal number of fragments
// All the fragments of a block must contain an overlap
////////////////////////////////////////////////////////////////////////

Block::Block(std::vector<AlignedFragment>& afrags, MultiSequence *seqs, 
			 bool enableTransitivity, Block *prohibited)
{
	part = 0;
	int i,a;
	int j, k;

//Make a hash
	int seqNum = seqs->GetNumSequences();
	IVECT *hash = new IVECT[seqNum*seqNum];
	for (i = 0; i < (int)afrags.size(); i++){
		hash[afrags[i].id[0]*seqNum + afrags[i].id[1]].push_back(i);
	}

	std::vector<Fragment> t_frags;//temporal set of fragments
	//First find a seed - fragment aligned to maximal number of other fragment
	int best,best_id, best_pos, multi,end,bfrag,bk,best_end;
	int *used = new int[afrags.size()]; //indicators of used fragmentsem
	int *bad = new int[afrags.size()];//fragments that cannot be a seed
	for(i=0;i<(int)afrags.size();i++)  bad[i] = -1;
	std::vector<int> starts,finishs;
	Fragment *frs;

	do{
		best = 0;
	for( i = 0; i < (int)afrags.size(); i++){ 
		int sec;
		for ( k = 0; k < 2; k++){
			if(k == bad[i]) continue;
			multi = 0;
			int cr_end = 100000;
			for (j = 0; j < seqNum; j++){
				int cId = afrags[i].id[k] * seqNum + j;
				for (a = 0; a < (int)hash[cId].size(); a++){
					int jj = hash[cId][a];
					if (afrags[jj].GetAlignPos(afrags[i].id[k],afrags[i].begin[k],sec) != NULL &&
						afrags[jj].end[1-sec]-afrags[i].begin[k] +1 >= MINLENGTH){
						multi++;
						cr_end = min(cr_end,afrags[jj].end[1-sec]);
					}
				}
			}
			for (j = 0; j < seqNum; j++){
				int cId = j * seqNum + afrags[i].id[k];
				for (a = 0; a < (int)hash[cId].size(); a++){
					int jj = hash[cId][a];
					if (afrags[jj].GetAlignPos(afrags[i].id[k],afrags[i].begin[k],sec) != NULL &&
						afrags[jj].end[1-sec]-afrags[i].begin[k] +1 >= MINLENGTH){
						multi++;
						cr_end = min(cr_end,afrags[jj].end[1-sec]);
					}
				}
			}
			if (multi > best){
				best = multi; best_id = afrags[i].id[k]; best_pos = afrags[i].begin[k];end = afrags[i].end[k];
				bfrag = i;bk=k;
				best_end = cr_end;
			}
		}
	}
	
	multi = best;
	Fragment first(best_pos, best_end, multi, best_id);
	seed = first;
	for (i = 0; i < (int)afrags.size(); i++) used[i] = 0;
	//Grow maximal block from the seed;
	frs = new Fragment[multi+1];
	frs[0] = first;
	k=1;
		for  (i = 0; i < (int)afrags.size(); i++){ 
			if(used[i]) continue;
			std::pair<int, int> *start,*finish;
			int ssecond,fsecond;
			if((start = afrags[i].GetAlignPos(first.id, first.begin, ssecond)) != NULL &&
				(finish = afrags[i].GetAlignPos(first.id, first.end, fsecond)) != NULL &&
				ssecond == fsecond && finish->second - start->second+1 >= first.length*0.6){
				Fragment fr(start->second, finish->second,
					0, start->first);
				starts.push_back(afrags[i].begin[1-ssecond]);
				finishs.push_back(afrags[i].end[1-ssecond]);
				frs[k++] = fr;
				used[i] = 1;
				delete start;delete finish;
			}
		}

	multi = k -1;
	if(multi == 0) {
		bad[bfrag]=bk;
		delete frs;
	}
	}while(multi==0);

	int mean = 0; 
	for (i = 0; i < (int)starts.size(); i++)
		mean += starts[i];
	mean /= starts.size();
	int sd = 0;
	for (i = 0; i < (int)starts.size(); i++)
		sd += (starts[i]-mean) * (starts[i]-mean);
	sd = (int)sqrt(sd/starts.size());
	int newMean = 0, num = 0;
	for(i = 0; i < (int)starts.size(); i++){
		if(abs(starts[i]-mean) <= sd){
			num++;
			newMean += starts[i];
		}
	}
	newMean /= num;
	int extendLeft = seed.begin - newMean;

	mean = 0; 
	for (i = 0; i < (int)finishs.size(); i++)
		mean += finishs[i];
	mean /= finishs.size();
	sd = 0;
	for (i = 0; i < (int)finishs.size(); i++)
		sd += (finishs[i]-mean) * (finishs[i]-mean);
	sd = (int)sqrt(sd/finishs.size());
	newMean = num = 0;
	for(i = 0; i < (int)finishs.size(); i++){
		if(abs(finishs[i]-mean) <= sd){
			num++;
			newMean += finishs[i];
		}
	}
	newMean /= num;
	int extendRight = newMean - seed.end;

	for(i = 0; i < multi+1; i++){
		frs[i].begin -= extendLeft;//extension ;
		if(frs[i].begin < 1)
			frs[i].begin = 1;
		frs[i].end += extendRight;//extension;
		if(frs[i].end > seqs->GetSequence(frs[i].id).GetLength())
			frs[i].end = seqs->GetSequence(frs[i].id).GetLength();
	}

	for(i = 0; i < multi+1; i++)
		t_frags.push_back(frs[i]);
	delete frs;

		//Extend block through transitivity
	int size = t_frags.size();
	if(enableTransitivity){
	for (i = 1; i < size; i++){
		Fragment fr = t_frags[i];
		fr.begin += min(extendLeft+2,10);
		if(fr.begin < 1) fr.begin = 1;
		fr.end -= min(extendRight+2,10);
		if(fr.end > seqs->GetSequence(fr.id).GetLength()) fr.end = seqs->GetSequence(fr.id).GetLength();
		for (j = 0; j < (int)afrags.size(); j++){
			if(used[j]) continue;
			AlignedFragment af = afrags[j];
			for ( k = 0; k < 2; k++){
				if (fr.id == af.id[k] && fr.begin >= af.begin[k] && fr.end <= af.end[k]){
					Fragment *frn = af.GetAlignFragment(fr);
					frn->begin -= min(extendLeft+2,10);
					if(frn->begin < 1) frn->begin = 1;
					frn->end += min(extendRight+2,10);
					if(frn->end > seqs->GetSequence(frn->id).GetLength()) 
						frn->end = seqs->GetSequence(frn->id).GetLength();
					int flag = 1;
					for(int m = 0; m < (int)t_frags.size() && flag; m++){
						Fragment frm = t_frags[m];
						if (frm.id == frn->id && 
							Overlap(frm.begin,frm.end,frn->begin,frn->end) > 10)
							flag = 0;
					}
					for (int a = 1; flag && prohibited && a<(int)prohibited->size(); a++){
						if((*prohibited)[a].id == frn->id  &&
							Overlap((*prohibited)[a].begin,(*prohibited)[a].end,frn->begin,frn->end) > 0)
							flag = 0;
					}
					if(flag){
						t_frags.push_back(*frn);
						used[j] = 1;
					}
					delete frn;
				}
			}
		}
	}
	}

	//Copy vector to an array for speed and changes
	frs = new Fragment[t_frags.size()];
	size = t_frags.size();
	for(i = 0; i < size; i++)
		frs[i] = t_frags[i];
//Get the shortest distance between two fragments in a same sequence
	int closest = 1000000;
	for(i = 0; i < size; i++){
		for (j = i+1; j < size; j++){
			if(frs[i].id == frs[j].id){
				int distance = frs[i].end - frs[j].end > 0? frs[i].begin-frs[j].end : frs[j].begin-frs[i].end;
				if(closest > distance) 
					closest = distance;
			}
		}
	}
	
	closest = closest/2 + closest % 2;

	for( i = 0; i < size; i++){
		for ( j = i+1; j < size; j++){
			if(frs[i].id == frs[j].id && 
				Overlap(frs[i].begin,frs[i].end,frs[j].begin,frs[j].end) > 0){
				if(frs[i].begin > frs[j].begin){
					frs[i].begin = (frs[j].end + frs[i].begin)/2;
					frs[j].end = frs[i].begin - 1;
				}
				else{
					frs[j].begin = (frs[i].end + frs[j].begin)/2;
					frs[i].end = frs[j].begin - 1;
				}
			}
		}
	}
	for(i =0 ; i < size; i++)
		frags.push_back(frs[i]);
	delete frs;
	delete used;
	delete bad;
	delete [] hash;
	if(prohibited) delete prohibited;
}


void Block::PrintBlock(FILE *f, MultiSequence *seqs, int compare)
{
	int unsigned i;
	if(!compare){//Normal output
		fprintf(f,"\n");
		for (i = 0; i< frags.size(); i++){
			fprintf(f, "%s(%d-%d) ", seqs->GetSequence(frags[i].id).GetName(), frags[i].begin, 
				frags[i].end);
		}
		fprintf(f,"\n\n");
	}
	else{//Output for comparison with blast
		fprintf(f,">\n");
		for (i = 0; i< frags.size(); i++){
			fprintf(f, "%s\t%d\t%d\n", seqs->GetSequence(frags[i].id).GetName(), frags[i].begin, 
				frags[i].end);
		}
	}

}

///////////////////////////////////////////////////////////////////////////////////////
// Returns the size of the shortest fragment
///////////////////////////////////////////////////////////////////////////////////////

int Block::GetLength()
{
	if (frags.size() == 0) return -1;
	int shortest = frags[0].length;
	for (int unsigned i = 1; i < frags.size(); i++)
		if (shortest > frags[i].length)
			shortest = frags[i].length;
	return shortest;
}

void Block::AddFragment(Fragment &fr)
{
	frags.push_back(fr);
}


Block& Block::operator =(const Block &bl)
{
	frags.clear();
	frags = bl.frags;
	seed = bl.seed;
	part = bl.part;
	return *this; 
}

/////////////////////////////////////////////////////////////////////////////////
// Remove used AlignedFragments or cut them based on the block
/////////////////////////////////////////////////////////////////////////////////
typedef AlignedFragment *pAF;
int Block::AdjustAFragmentList(AVECT &fragments, Matrix *similarity, float threshold)
{
	int result = 0;
	int i, j, m,k,flag;
	int asize = fragments.size();
	int size = frags.size();
	int *used = new int[asize];
	AVECT tmp = fragments;
	fragments.clear();
	for (i=0; i < asize; i++)
		used[i] = 0;
	for (i = 0; i < size-1; i++){
		Fragment fi = frags[i];
		for (j = i+1; j < size; j++){
			if(similarity && (*similarity)(0,i,j) >= threshold) continue;
			Fragment fj = frags[j];
			flag = 0;
			for (m = 0; m < asize; m++){
				if (used[m]) continue;
				AlignedFragment af = tmp[m];
				for (k=0;k<2;k++){
					if (af.id[k] == fi.id && af.id[1-k] == fj.id &&
						Overlap(af.begin[k],af.end[k],fi.begin,fi.end) > 1 &&
						Overlap(af.begin[1-k],af.end[1-k],fj.begin,fj.end) > 1){
						flag = 1;
						break;
					}
				}
				if (flag){//found overlap
					result++;
					AlignedFragment n1, n2;
					af.Adjust(fi,fj,n1,n2);
					if (n1.GetLength() >= MINLENGTH) fragments.push_back(n1);
					if (n2.GetLength() >= MINLENGTH) fragments.push_back(n2);
					used[m] = 1;
					flag = 0;
				}
			}
		}
	}
	
	for (m = 0; m < asize; m++)
		if(used[m] == 0) fragments.push_back(tmp[m]);
	tmp.clear();
	delete used;
	return result;
}

int Block::AdjustAFragmentList(AVECT &fragments, int seqNum, Matrix *similarity, float threshold)
{
	int result = 0;
	int i, j, m,k,flag;
	int asize = fragments.size();
	int size = frags.size();

//Make a hash
	IVECT *hash = new IVECT[seqNum*seqNum];
	for (i = 0; i < asize; i++){
		hash[fragments[i].id[0]*seqNum + fragments[i].id[1]].push_back(i);
	}


	int *used = new int[asize];
	AVECT tmp = fragments;
	fragments.clear();
	for (i=0; i < asize; i++)
		used[i] = 0;
	for (i = 0; i < size-1; i++){
		Fragment fi = frags[i];
		for (j = i+1; j < size; j++){
			if(similarity && (*similarity)(0,i,j) >= threshold) continue;
			Fragment fj = frags[j];
			flag = 0;
			for(k=0;k<2;k++){
				int cId = k==0? fi.id * seqNum + fj.id : fj.id * seqNum +fi.id;
			for (int a = 0; a < (int) hash[cId].size(); a++){
				m = hash[cId][a];
				if (used[m]) continue;
				AlignedFragment af = tmp[m];
				if (af.id[k] == fi.id && af.id[1-k] == fj.id &&
						Overlap(af.begin[k],af.end[k],fi.begin,fi.end) > 1 &&
						Overlap(af.begin[1-k],af.end[1-k],fj.begin,fj.end) > 1){
						flag = 1;
				}
				if (flag){//found overlap
					result++;
					AlignedFragment n1, n2;
					af.Adjust(fi,fj,n1,n2);
					if (n1.GetLength() >= MINLENGTH) fragments.push_back(n1);
					if (n2.GetLength() >= MINLENGTH) fragments.push_back(n2);
					used[m] = 1;
					flag = 0;
				}
			}
			}
		}
	}
	
	for (m = 0; m < asize; m++)
		if(used[m] == 0) fragments.push_back(tmp[m]);
	tmp.clear();
	delete used;
	delete [] hash;
	return result;
}


