#ifndef LIS_APPROX_INFERENCE_H_
#define LIS_APPROX_INFERENCE_H_

#include <vector>
#include <iostream>
#include <set>
#include <map>
using namespace std;
#include "lvrmln.h"
#include "decomposer.h"
#include "propositionalresolution.h"
#include "logdouble.h"
#include "extensions.h"
#include "heuristics.h"
#include "normalizer.h"
#include "clusterutil.h"
#include "proposalstructure.h"
#include "lvrnormpropagation.h"
#include "lvrsingletonnormpropagation.h"
#include "lvparams.h"
#include "lvgrbestimator.h"

struct LISApproxInference
{
	void setSamplingMode(ESamplingMode samplingMode_)
	{
		samplingMode = samplingMode_;
	}
	void setProposalDistribution(LProposalDistribution* distribution_)
	{
		distribution = distribution_;
	}
	LISApproxInference(LvrMLN& mln);
	~LISApproxInference()
	{
		delete decomposer;
		delete heuristics;
		delete lvrNormPropagate;
		if(rbestimator)
			delete rbestimator;
	}
	LogDouble doLvApproxPartitionInformed(vector<WClause*>& CNF);
	LogDouble doLvApproxPartitionInformedV1(vector<WClause*>& CNF);
	LogDouble doLvApproxPartitionBinomial(vector<WClause*>& CNF);
	void estimateApproxMarginals(LvrParams* params, LProposalDistribution* distribution = NULL);
	LogDouble estimatePartitionFunction(LvrParams* params,LProposalDistribution* distribution= NULL);
	Atom* selectAtomToCondition(vector<WClause*> clauses)
	{
		return heuristics->getAtomToSplit(clauses);
	}
	bool decomposeCNF(vector<WClause*>& CNF,int& powerFactor);
	LogDouble doLvApproxPartitionInformedRB(vector<WClause*>& CNF);

private:
	LvrMLN& mln;
	LDecomposer* decomposer;
	LHeuristics* heuristics;	
	LvrNormPropagation* lvrNormPropagate;
	ESamplingMode samplingMode;
	
	LvgRBEstimator* rbestimator;
	//not owned
	LProposalDistribution* distribution;
	bool raoblackwelize;
	void seperateDisjointCNF(vector<WClause*> CNF, vector<vector<WClause*> >& disjointCNFList);
	LogDouble CNFWeight(vector<WClause*>& CNF);
};
#endif
