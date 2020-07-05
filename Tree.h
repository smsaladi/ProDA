//////////////////////////////////////////////////////////////////////
// Tree.h
//////////////////////////////////////////////////////////////////////

#ifndef TREE_H
#define TREE_H

#include <vector>
#include "Matrix.h"
#include "MultiSequence.h"

//////////////////////////////////////////////////////////////////////
// Expected accuracy tree
//////////////////////////////////////////////////////////////////////

typedef std::vector<int> IVECT;
class Tree {

 public:

  //////////////////////////////////////////////////////////////////////
  // Tree node struct
  //////////////////////////////////////////////////////////////////////
  
  class TreeNode {
    
    bool isLeaf;
    int numSequences;
    MultiSequence *seqs;
    
    TreeNode *leftChild;
    TreeNode *rightChild;

  public:
	  void UpdateIDs(int *used);
	  void GetIDs(IVECT &ids);

    // constructor and destructor
    TreeNode (bool isLeaf, int numSequences, MultiSequence *seqs,
	      TreeNode *leftChild, TreeNode *rightChild);
    ~TreeNode ();
    
    // getters
    const bool GetIsLeaf() const;
    const int GetNumSequences() const;
    const MultiSequence *GetSequences() const;
    const TreeNode *GetLeftChild() const;
    const TreeNode *GetRightChild() const;

    // print subtree starting at this node
    void Print (FILE *file) const;
    
    // progressive alignment
    void ProgressiveAlignment (int numSequences, SparseMatrix **posteriors);
  };
  
 private:

  TreeNode *root;

 public:
	 void UpdateIDs(int *used);
	 void GetIDs(IVECT &ids);
	 int GetNumSequences();
  
  // constructor and destructor
  Tree (Matrix similarity, const MultiSequence &seqs, float threshold = 0);
  ~Tree ();
  
  // print tree
  void Print (FILE *file) const;
  
  // progressive alignment
  MultiSequence *ProgressiveAlignment (int numSequences, SparseMatrix **posteriors);
};

#endif

