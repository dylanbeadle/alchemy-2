#ifndef LVGRB_ESTIMATOR_H
#define LVGRB_ESTIMATOR_H

#include <vector>
#include <iostream>
#include <set>
#include <map>
using namespace std;
#include "lvrmln.h"
#include "decomposer.h"
#include "splitter.h"
#include "logdouble.h"
#include "heuristics.h"
#include "normalizer.h"
#include "clusterutil.h"
#include "ptptree.h"
#include "lvrnormpropagation.h"
#include "lvrsingletonnormpropagation.h"
#include "lvrptptreesampling.h"

struct LvgRBEstimator
{
	LvgRBEstimator(LvrMLN& mln);
	~LvgRBEstimator()
	{
		delete ptpTree;
		delete decomposer;
	}
	LogDouble weightedModelCount(vector<WClause*>& CNF,int parentId,int parentBranchNo);
	void addTreeNode(Atom* atom,vector<LogDouble> weights,vector<LogDouble> binCoeffs,
		LogDouble norm,NodeType nodeType,int nodeId, int parentId, int parentBranchNo,
		vector<LogDouble>* branchWeights = NULL,LogDouble nodeValue = LogDouble(0,false));
	LPTPNode* addTreeNode(NodeType nodeType,int nodeId, int parentId, 
		int parentBranchNo,int powerFactor = -1, LogDouble nodeValue=LogDouble(1,false));

	LogDouble CNFWeight(vector<WClause*>& CNF);
	bool decomposeCNF(vector<WClause*>& CNF,int& powerFactor);
	void updateRBEstimates(vector<WClause*>& CNF);
	LPTPTree* ptpTree;
	LDecomposer* decomposer;
	LvrMLN& mln;
	int id;
};


#endif
