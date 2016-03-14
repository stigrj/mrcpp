#include "GridCleaner.h"
#include "TreeAdaptor.h"
#include "DefaultCalculator.h"
#include "MWTree.h"
#include "Timer.h"

using namespace std;

template<int D>
GridCleaner<D>::GridCleaner(const MultiResolutionAnalysis<D> &mra)
        : TreeBuilder<D>(mra, -1) {
    this->calculator = new DefaultCalculator<D>();
    this->adaptor = new TreeAdaptor<D>();
}

template<int D>
GridCleaner<D>::GridCleaner(const MultiResolutionAnalysis<D> &mra,
                            const TreeAdaptor<D> &a)
        : TreeBuilder<D>(mra, -1) {
    this->calculator = new DefaultCalculator<D>();
    this->adaptor = a.copy();
}

template<int D>
GridCleaner<D>::~GridCleaner() {
    this->clearCalculator();
    this->clearAdaptor();
}

template<int D>
int GridCleaner<D>::clean(MWTree<D> &tree) const {
    Timer clean_t, split_t;
    if (this->calculator == 0) MSG_ERROR("Calculator not initialized");
    println(10, " == Clearing tree");

    split_t.restart();
    MWNodeVector *newVec = new MWNodeVector;
    MWNodeVector *workVec = tree.copyEndNodeTable();
    this->adaptor->splitNodeVector(*newVec, *workVec);
    int nSplit = newVec->size();
    delete workVec;
    delete newVec;
    split_t.stop();

    printout(10, "  -- #  0: Split        ");
    printout(10, setw(6) << nSplit << " nodes\n");

    clean_t.restart();
    MWNodeVector nodeVec;
    tree.makeNodeTable(nodeVec);
    int nClear = nodeVec.size();
    this->calculator->calcNodeVector(nodeVec);//clear all coefficients
    clean_t.stop();

    printout(10, "  -- #  1: Cleared      ");
    printout(10, setw(6) << nClear << " nodes\n");

    println(10, "");
    println(10, "Time split          " << split_t);
    println(10, "Time clean          " << clean_t);
    println(10, endl);

    tree.resetEndNodeTable();
    tree.clearSquareNorm();
    return nSplit;
}

template class GridCleaner<1>;
template class GridCleaner<2>;
template class GridCleaner<3>;
