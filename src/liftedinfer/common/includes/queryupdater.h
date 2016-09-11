#ifndef __LVRQUERY_UPDATER
#define __LVRQUERY_UPDATER
#include "lvrmln.h"
#include "lvratomhashtemplate.h"
#include "lvrpermutations.h"


struct LvrQueryUpdater
{
	static LvrQueryUpdater* Instance();
	static void createInstance(vector<vector<int> > queryIntRep,vector<string> queryStrings, 
		string resultfilename, vector<vector<int> >* evidenceIntRep = NULL);
	void updateQueryValues(Atom* atom,int sampledTrueVal);
	void updateQueryValues(Atom* atom,vector<bool> isolatedTerms,vector<int> sampledTrueVals);
	void updateQueryValuesLVGibbs(Atom* atom,int sampledTrueVal);
	//void updateQueryValuesLVGibbsRB(Atom* atom,LogDouble prob,int sampledvalue);
	void updateQueryValuesLVGibbsRB(Atom* atom,vector<LogDouble> prob);
	void updateDontCare();
	void updateGibbsDontCareRB();
	void writeGibbsRBToFile(int iterations);
	void normalizeRB(int iterations);
	static bool isInstanceCreated();
	void normalize(int iterations);
	void writeToFile(int numIterations);
	int getCurrentSampleValue(vector<int> queryintRep);
	void updateAllImportanceWeights(LogDouble currentIterWeight);
	void writeToFile(LogDouble totalWeight);
	bool isNormIdInQuery(int id)
	{
		if(queryNormIds.find(id)!=queryNormIds.end())
			return true;
		return false;
	}
	void updateExactQueryWeights(set<int> queryHashes,LogDouble value);
	vector<double> getCurrentProbs()
	{
		return currentProbabilities;
	}
	void resetallupdateflags();
	void updateAllImportanceWeightsRB(LogDouble currentIterWeight);
	void updateISRBEstimates(Atom* atom,vector<LogDouble> probs);
private:
	LvrAtomHashTemplate<string>* queryHashTemplate;
	LvrAtomHashTemplate<bool>* lvrAtomHashUpdateFlags;
	LvrQueryUpdater(vector<vector<int> > queryIntRep,vector<string> queryStrings,vector<vector<int> >* evidenceIntRep = NULL);
	static LvrQueryUpdater* m_pInstance;
	vector<double> currentProbabilities;

	//used for gibbs sampling to estimate MLE counts
	LvrAtomHashTemplate<int>* lvrAtomHashTemplate;
	//used for importance sampling MAR, exact inference ptp
	LvrAtomHashTemplate<LogDouble*>* cumulativeWeight;
	LvrAtomHashTemplate<int>* currentSampledValue;
	//Rao-Blackwell estimates
	LvrAtomHashTemplate<LogDouble*>* rbestimates;

	//store the norm ids for fast searching
	map<int,bool> queryNormIds;
	void doUpdateQueryValues(Atom* atom,int sampledTrueVal);
	void doUpdateQueryValues(Atom* atom,vector<bool> isolatedTerms,vector<int> sampledTrueVals);
	

};
#endif
