/**
*
*
*  \date Jul, 2016
*  \author Peter Wind <peter.wind@uit.no> \n
*  CTCC, University of Tromsø
*
*/

#ifndef SERIALTREE_H_
#define SERIALTREE_H_

#include <Eigen/Core>
#include "parallel.h"
#include "NodeIndex.h"
#include <vector>

template<int D> class MultiResolutionAnalysis;
template<int D> class ProjectedNode;
template<int D> class GenNode;
template<int D> class MWNode;
template<int D> class MWTree;
template<int D> class FunctionTree;
template<int D> class FunctionNode;

template<int D>
class SerialTree  {
public:
    SerialTree(FunctionTree<D> *tree, int max_nodes);
    ~SerialTree();

    FunctionTree<D>* getTree() { return this->tree_p; }

    void allocRoots(FunctionTree<D> &tree);
    void allocChildren(FunctionNode<D> &parent);
    void allocGenChildren(FunctionNode<D> &parent);

    void deallocNodes(int serialIx);
    void deallocGenNodes(int serialIx);

    ProjectedNode<D>* createSnode(const NodeIndex<D> &nIdx);
    void GenS_nodes(MWNode<D>* Node);
    void S_mwTransform(double* coeff_in, double* coeff_out, bool ReadOnlyScalingCoeff, int Children_Stride, bool overwrite=true);
    void S_mwTransformBack(double* coeff_in, double* coeff_out, int Children_Stride);

    void serialTreeAdd(double c, FunctionTree<D>* &TreeB, FunctionTree<D>* &TreeC);
    void serialTreeAdd_Up(double c, FunctionTree<D>* &TreeB, FunctionTree<D>* &TreeC);
    void rewritePointers(int Nchunks);

    friend class MWTree<D>;
    friend class MWNode<D>;
    friend class FunctionNode<D>;
    friend class ProjectedNode<D>;
    friend class GenNode<D>;

protected:
    int maxNodes;               //max number of nodes that can be defined
    int maxGenNodes;            //max number of Gen nodes that can be defined
    int sizeNodeCoeff;          //size of coeff for one node
    int sizeGenNodeCoeff;       //size of coeff for one Gen node

    int nNodes;                 //number of ProjectedNodes already defined
    int nGenNodes;              //number of GenNodes already defined
    int nNodesCoeff;            //number of nodes coeff already defined
    int nGenNodesCoeff;         //number of GenNodes coeff already defined

    double **coeffStack;
    double **genCoeffStack;

    int *nodeStackStatus;
    int *genNodeStackStatus;

    char *cvptr_ProjectedNode;//virtual table pointer for ProjectedNode
    char *cvptr_GenNode;// virtual table pointer for GenNode

    int maxNodesPerChunk;
    std::vector<ProjectedNode<D>*>  nodeChunks;
    std::vector<GenNode<D>*>  genNodeChunks;
    std::vector<double*>  nodeCoeffChunks;
    std::vector<double*>  genNodeCoeffChunks;

    ProjectedNode<D> *sNodes;   //serial ProjectedNodes
    GenNode<D> *sGenNodes;      //serial GenNodes
    double *sNodesCoeff;        //serial ProjectedNodes coefficients
    double *sGenNodesCoeff;     //serial GenNodes coefficients

    FunctionTree<D> *tree_p;
    ProjectedNode<D>* lastNode; //pointer to the last active node
    GenNode<D>* lastGenNode;    //pointer to the last active Gen node

    ProjectedNode<D>* allocNodes(int nAlloc, int* serialIx, double **coefs_p);
    GenNode<D>* allocGenNodes(int nAlloc, int* serialIx, double **coefs_p);

private:
#ifdef HAVE_OPENMP
    omp_lock_t Stree_lock;
#endif
};

#endif /* SERIALTREE_H_*/
