//////////////////////////////////////////////////////////////////////
// Tree.cc
//////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include "Tree.h"
#include "GlobalAlign.h"

//////////////////////////////////////////////////////////////////////
// TreeNode constructor
//////////////////////////////////////////////////////////////////////

Tree::TreeNode::TreeNode (bool isLeaf, int numSequences, MultiSequence *seqs,
			  TreeNode *leftChild, TreeNode *rightChild) : 
  isLeaf(isLeaf), numSequences(numSequences), seqs(seqs),
  leftChild(leftChild), rightChild(rightChild) {}

//////////////////////////////////////////////////////////////////////
// TreeNode destructor
//////////////////////////////////////////////////////////////////////

Tree::TreeNode::~TreeNode (){
  if (seqs) delete seqs;
  if (!isLeaf){
    if (leftChild) delete leftChild;
    if (rightChild) delete rightChild;
  }
}

//////////////////////////////////////////////////////////////////////
// Returns whether this node is a leaf or not
//////////////////////////////////////////////////////////////////////

const bool Tree::TreeNode::GetIsLeaf() const {
  return isLeaf;
}

//////////////////////////////////////////////////////////////////////
// Returns number of sequences in this subtree
//////////////////////////////////////////////////////////////////////

const int Tree::TreeNode::GetNumSequences() const {
  return numSequences;
}

//////////////////////////////////////////////////////////////////////
// Returns the MultiSequence data associated with this node.
//////////////////////////////////////////////////////////////////////

const MultiSequence *Tree::TreeNode::GetSequences() const {
  return seqs;
}

//////////////////////////////////////////////////////////////////////
// Returns the left child TreeNode
//////////////////////////////////////////////////////////////////////

const Tree::TreeNode *Tree::TreeNode::GetLeftChild() const {
  return leftChild;
}

//////////////////////////////////////////////////////////////////////
// Returns the right child TreeNode
//////////////////////////////////////////////////////////////////////

const Tree::TreeNode *Tree::TreeNode::GetRightChild() const {
  return rightChild;
}

//////////////////////////////////////////////////////////////////////
// Print out subtree recursively
//////////////////////////////////////////////////////////////////////

void Tree::TreeNode::Print (FILE *file) const {
  if (isLeaf){
    fprintf (file, "%s", seqs->GetSequence(0).GetName());
  } else {
    fprintf (file, "(");
    leftChild->Print(file);
    if (leftChild->isLeaf || rightChild->isLeaf)
      fprintf (file, " ");
    rightChild->Print(file);
    fprintf (file, ")");
  }
}

//////////////////////////////////////////////////////////////////////
// Perform progressive alignment for subtree
//////////////////////////////////////////////////////////////////////

void Tree::TreeNode::ProgressiveAlignment (int numSequences, SparseMatrix **posteriors){

  if (seqs) return;

  leftChild->ProgressiveAlignment (numSequences, posteriors);
  rightChild->ProgressiveAlignment (numSequences, posteriors);

  seqs = GlobalAlign::AlignGroups (numSequences, posteriors, 
				   *leftChild->seqs, *rightChild->seqs);
  seqs->Sort();
}

void Tree::TreeNode::GetIDs(IVECT &ids)
{
	if(isLeaf){
		ids.push_back(seqs->GetSequence(0).GetID());
		Sequence *ptr = seqs->GetSequencePtr(0);
		ptr->SetID(ids.size()-1);
	}
	else{
		leftChild->GetIDs(ids);
		rightChild->GetIDs(ids);
	}
}

void Tree::TreeNode::UpdateIDs(int *used)
{
		
	if(isLeaf){

		Sequence *ptr = seqs->GetSequencePtr(0);
		ptr->SetID(used[ptr->GetID()]);
	}
	else{
		leftChild->UpdateIDs(used);
		rightChild->UpdateIDs(used);
	}
}

//////////////////////////////////////////////////////////////////////
// Tree constructor
//////////////////////////////////////////////////////////////////////

typedef Tree::TreeNode *TreeNodePtr;

Tree::Tree (Matrix similarity, const MultiSequence &seqs, float threshold){
  // get number of sequences
  
  int n = similarity.GetNumRows();
  ASSERT (n == similarity.GetNumCols(), "Similarity matrix not square.");
  
  int i,j,k;
  // initialize diagonal of distance matrix
  for (i = 0; i < n; i++)
    similarity(0,i,i) = -1;
  
  // build initial set of trees
  
  TreeNode **trees = new TreeNodePtr[n];
  ASSERT (trees, "Out of memory.");
  
  for (i = 0; i < n; i++){
    MultiSequence *tempSeqs = new MultiSequence();
    ASSERT (tempSeqs, "Out of memory.");
    Sequence *seq = new Sequence (seqs.GetSequence (i));
    ASSERT (seq, "Out of memory.");
    tempSeqs->AddSequence (seq);
    trees[i] = new TreeNode (true, 1, tempSeqs, NULL, NULL);
  }
  
  // perform n - 1 merges
  
  for (k = 0; k < n-1; k++){
    
    // find nearest neighbors
    
    int bi = 0, bj = 0;
    for (i = 0; i < n; i++){
      for (j = i+1; j < n; j++) if (i != j){
	if (similarity(0,i,j) > similarity(0,bi,bj)){
	  bi = i;
	  bj = j;
	}
      }
    }
	//stop if similarity drops bellow the threshold
	int numSeq = trees[bi]->GetNumSequences() + trees[bj]->GetNumSequences();
	if((similarity(0,bi,bj) < threshold && numSeq <= 2) ||
	   similarity(0,bi,bj) < threshold - 0.2) break;

    // merge trees
    
    TreeNode *temp = new TreeNode (false, numSeq,
				   NULL, trees[bi], trees[bj]);
    ASSERT (temp, "Out of memory.");
    trees[bi] = temp;
    trees[bj] = NULL;
    
    // update distances
    similarity(0,bi,bj) = similarity(0,bj,bi) = -1;
    for (int m = 0; m < n; m++) if (m != bi && m != bj){
      similarity(0,bi,m) = similarity(0,m,bi) = 
	(similarity(0,bi,m) * temp->GetLeftChild()->GetNumSequences() + 
	 similarity(0,bj,m) * temp->GetRightChild()->GetNumSequences()) /
	temp->GetNumSequences();
      similarity(0,bj,m) = similarity(0,m,bj) = -1;
    }
  }
  
  root = trees[0];
  if( root->GetNumSequences() == 1) {
	  for (i = 1; i < n; i++){
		if (trees[i] != NULL && trees[i]->GetNumSequences() > root->GetNumSequences()){
			root = trees[i];
		}
	  }
  }
  for (i = 0; i < n; i++){
	  if(trees[i] && trees[i] != root)
		  delete trees[i];
  }
  delete[] trees;
}

//////////////////////////////////////////////////////////////////////
// Tree destructor
//////////////////////////////////////////////////////////////////////

Tree::~Tree (){
  delete root;
}

//////////////////////////////////////////////////////////////////////
// Print tree
//////////////////////////////////////////////////////////////////////

void Tree::Print (FILE *file) const {
  ASSERT (root, "Tree not created.");

  root->Print (file);
  fprintf (file, "\n");
}

//////////////////////////////////////////////////////////////////////
// Perform progressive alignment on the tree
//////////////////////////////////////////////////////////////////////

MultiSequence *Tree::ProgressiveAlignment (int numSequences, SparseMatrix **posteriors){
  root->ProgressiveAlignment (numSequences, posteriors);
  return new MultiSequence(*root->GetSequences());
  
}

int Tree::GetNumSequences()
{
	if(root)
		return root->GetNumSequences();
	else return 0;
}


void Tree::GetIDs(IVECT &ids)
{
	root->GetIDs(ids);
}

void Tree::UpdateIDs(int *used)
{
	root->UpdateIDs(used);
}

