#include "lisapproxinference.h"
#include "cleanuputils.h"
#include "rulesutil.h"
#include <time.h>
#include "mathutils.h"
#include "samplingalgs.h"
#include "queryupdater.h"
#include <sstream>
#include <fstream>
using namespace std;

LogDouble LISApproxInference::CNFWeight(vector<WClause*>& CNF)
{
	/*for(unsigned int i=0;i<CNF.size();i++)
	{
		if(CNF[i]->satisfied)
			cout<<"Satisfied v";
		CNF[i]->print();
	}
	*/
	vector<Atom*> removedAtoms;
	set<int> completedIds;
	LogDouble totalWeight(1,false);
	vector<WClause*> remainingCNF;
	for(unsigned int i=0;i<CNF.size();i++)
	{
		if(CNF[i]->satisfied)
		{
			for(unsigned int j=0;j<CNF[i]->atoms.size();j++)
			{
				if(completedIds.count(CNF[i]->atoms[j]->symbol->id) == 0)
				{
					removedAtoms.push_back(CNF[i]->atoms[j]);
					completedIds.insert(CNF[i]->atoms[j]->symbol->id);
				}
				
			}
			//get weight of formula
			set<LvrTerm*> terms;
			for(unsigned int kk=0;kk<CNF[i]->atoms.size();kk++)
			{
				for(unsigned int mm=0;mm<CNF[i]->atoms[kk]->terms.size();mm++)
				{
					terms.insert(CNF[i]->atoms[kk]->terms[mm]);
				}
			}
			int sz=1;
			for(set<LvrTerm*>::iterator it=terms.begin();it!=terms.end();it++)
			{
				sz *= (*it)->domain.size();
			}
			//totalWeight += CNF[i]->weight*LogDouble(sz,false);
			//if(!CNF[i]->weight.is_zero)
			{
				LogDouble tmpL = CNF[i]->weight*LogDouble(sz,false);
				double tmp = exp(tmpL.value);
				LogDouble newWt(tmp,true);
				totalWeight *= newWt;
			}
		}
		else
		{
			remainingCNF.push_back(LvrMLN::create_new_clause(CNF[i]));
		}
		
	}
	int totalDontcare = 0;
	for(vector<Atom*>::iterator it=removedAtoms.begin();it!=removedAtoms.end();it++)
	{
		//find the number of groundings
		int sz = 1;
		for(unsigned int i=0;i<(*it)->terms.size();i++)
		{
			sz *= (*it)->terms[i]->domain.size();
		}
		totalDontcare += sz;
	}
	LogDouble d(totalDontcare,false);
	LogDouble out(1,false);
	LogDouble c(2,false);
	LogDouble::LDPower(c,totalDontcare,out);
	LogDouble tVal = totalWeight*out;
	cleanup(CNF);
	CNF = remainingCNF;

	return tVal;
}


bool LISApproxInference::decomposeCNF(vector<WClause*>& CNF,int& powerFactor)
{
	vector<Decomposer*> decomposer_list;
	decomposer->find_decomposer(CNF,decomposer_list);
	
#ifdef __DEBUG_PRINT__
	cout<<"Decomposer={ ";
	for(unsigned int i=0;i<decomposer_list.size();i++)
	{
		decomposer_list[i]->print();
	}
	cout<<"}"<<endl;
#endif
	if(decomposer_list.size()==0)
		return false;
	else
	{
		//replace all decomposer terms with a constant from the domain
		powerFactor = 1;
		for(unsigned int i=0;i<decomposer_list.size();i++)
		{
			powerFactor*=decomposer_list[i]->decomposer_terms[0]->domain.size();
			int domSize = decomposer_list[i]->decomposer_terms[0]->domain.size();
			for(unsigned int j=0;j<decomposer_list[i]->decomposer_terms.size();j++)
			{
				//in all clauses that this term occurs change weight
				/*
				for(unsigned int jj=0;jj<CNF.size();jj++)
				{
					for(unsigned int kk=0;kk<CNF[jj]->atoms.size();kk++)
					{
						bool done = false;
						for(unsigned int mm=0;mm<CNF[jj]->atoms[kk]->terms.size();mm++)
						{
							if(CNF[jj]->atoms[kk]->terms[mm] == decomposer_list[i]->decomposer_terms[j])
							{
								//change weight
								CNF[jj]->weight = CNF[jj]->weight*LogDouble(domSize,false);
								done=true;
							}
						}
						if(done)
							break;
					}
				}
				*/
				//store pre-decomp domain
				decomposer_list[i]->decomposer_terms[j]->origDomain.clear();
				for(unsigned int jj=0;jj<decomposer_list[i]->decomposer_terms[j]->domain.size();jj++)
					decomposer_list[i]->decomposer_terms[j]->origDomain.push_back(decomposer_list[i]->decomposer_terms[j]->domain[jj]);

				//Erase all but one grounding
				decomposer_list[i]->decomposer_terms[j]->domain.erase(decomposer_list[i]->decomposer_terms[j]->domain.begin(),decomposer_list[i]->decomposer_terms[j]->domain.begin()+decomposer_list[i]->decomposer_terms[j]->domain.size()-1);
				//insert the original size into the term
				decomposer_list[i]->decomposer_terms[j]->origDomainSize = domSize;
			}
		}
	}
	cleanup(decomposer_list);
	return true;
}
////IMPLEMENTATION1: NO DISJOINT SEPERATION MORE GROUNDING, faster in some cases
LogDouble LISApproxInference::doLvApproxPartitionInformedV1(vector<WClause*>& CNF)
{
	LogDouble totalVal(1,false);
	if(CNF.size()==0)
		return totalVal;
	int completedCount=0;
	for(unsigned int i=0;i<CNF.size();i++)
	{
		if(CNF[i]->atoms.size()==0 || CNF[i]->satisfied)
		{
			completedCount++;
		}
	}
	if(completedCount == CNF.size())
	{
		LogDouble weight = CNFWeight(CNF);
		cleanup(CNF);
		return weight;
	}
	
	int powerFactor;
	bool isDecomposed = decomposeCNF(CNF,powerFactor);
	if(isDecomposed)
	{
		LogDouble mcnt = doLvApproxPartitionInformedV1(CNF);
		LogDouble val = LogDouble(1,false);
		LogDouble::LDPower(mcnt,powerFactor,val);
		totalVal = totalVal*val;
		return totalVal;
	}
	
	Atom* tmpAtom = heuristics->getAtomToSplit(CNF);	
	if(tmpAtom==NULL)
	{	
		LogDouble wt1 = CNFWeight(CNF);
		cleanup(CNF);
		return totalVal*wt1;
			
	}
	Atom* atom = LvrMLN::create_new_atom(tmpAtom);
	//sample according to mode
	LProposalDistributionElement* lpe = lpe = distribution->findElement(atom);
	//use Isolated terms rule
	vector<bool> isolatedTerms;
	bool isIsolated = LRulesUtil::computeIsolatedTerms(atom,CNF,isolatedTerms);
	double binCoeff = 1;
	LogDouble probOfSample(1,false);
	LogDouble sampleWeight(1,false);
	vector<int> sampledValues;
	if(isIsolated)
	{
		int isolatedSize=1;
		int nonisolatedSize = 1;
		for(unsigned int jj=0;jj<isolatedTerms.size();jj++)
		{
			if(isolatedTerms[jj])
			{
				isolatedSize *= atom->terms[jj]->domain.size();
			}
			else
			{
				nonisolatedSize *=  atom->terms[jj]->domain.size();
			}
		}
		LSamplingAlgs::Instance()->sample(isolatedSize,nonisolatedSize,sampleWeight,sampledValues,samplingMode,lpe);
		if(LvrQueryUpdater::isInstanceCreated())
		{
			LvrQueryUpdater::Instance()->updateQueryValues(atom,isolatedTerms,sampledValues);
		}
	}
	else
	{
		int numGroundings = atom->getNumberOfGroundings();
		LSamplingAlgs::Instance()->sample(numGroundings,1,sampleWeight,sampledValues,samplingMode,lpe);
		if(LvrQueryUpdater::isInstanceCreated())
		{
			LvrQueryUpdater::Instance()->updateQueryValues(atom,sampledValues[0]);
		}
		
	}
	lvrNormPropagate->propagateNormalizedCNF(CNF,atom,isolatedTerms,sampledValues);
	LogDouble mcnt = doLvApproxPartitionInformedV1(CNF);
	totalVal *= sampleWeight*mcnt;
	delete atom;
	cleanup(CNF);
	return totalVal;
}

LogDouble LISApproxInference::doLvApproxPartitionInformed(vector<WClause*>& CNF)
{
	LogDouble totalVal(1,false);
	if(CNF.size()==0)
	{
		return totalVal;
	}	
	int powerFactor;
	bool isDecomposed = decomposeCNF(CNF,powerFactor);
	if(isDecomposed)
	{
		LogDouble mcnt = doLvApproxPartitionInformed(CNF);
		LogDouble val = LogDouble(1,false);
		LogDouble::LDPower(mcnt,powerFactor,val);
		totalVal = totalVal*val;
	}
	else
	{
		int singletonIndex;
		bool singleton=false;
		Atom* constatom=NULL;
		Atom* otheratom=NULL;
		Atom* singatom=NULL;
		Atom* tmpAtom = NULL;
		for(int ii=0;ii<CNF.size();ii++)	
		{
		  	bool foundatom=false;
			for(int jj=0;jj<CNF[ii]->atoms.size();jj++)
			{
				otheratom = CNF[ii]->atoms[jj];
				singleton = CNF[ii]->atoms[jj]->isSingletonAtom(singletonIndex);
				if(singleton)
				{
					singatom = CNF[ii]->atoms[jj];
					foundatom=true;
					break;
				}
				else if(CNF[ii]->atoms[jj]->isConstant())
					constatom = CNF[ii]->atoms[jj];
			}
			if(foundatom)
				break;
		}
		if(singatom)
			tmpAtom = singatom;
		else if(constatom)
			tmpAtom=constatom;
		else
			tmpAtom=otheratom;
		if(tmpAtom==NULL)
		{
			LogDouble wt1 = CNFWeight(CNF);
			cleanup(CNF);
		    totalVal = totalVal*wt1;
			return totalVal;			
		}
		bool nonpoly;
		if(!singleton)
			nonpoly=true;
		else
			nonpoly=false;
		Atom* atom = LvrMLN::create_new_atom(tmpAtom);
		vector<bool> isolatedTerms;
		bool isIsolated = false;
		if(nonpoly)
		{
			if(!atom->isConstant())
			{
				if(samplingMode==EINFORMED)
				{
					isIsolated = LRulesUtil::computeIsolatedTerms(atom,CNF,isolatedTerms);
				}
				if(!isIsolated)
				{
					nonpoly=false;
					/*int maxdomid=0;
					int maxdomsize=0;
					for( int ii=0;ii<atom->terms.size();ii++)
					{
						if(atom->terms[ii]->domain.size()>maxdomsize)
						{
							maxdomid=ii;
							maxdomsize=atom->terms[ii]->domain.size();
						}
					}
					singletonIndex = maxdomid;
					*/
					bool first=true;
					for( int ii=0;ii<atom->terms.size();ii++)
					{
						if(atom->terms[ii]->domain.size()>1)
						{
							if(first)
							{
								singletonIndex = ii;
								first=false;
							}
							double r = LvRandomGenUtil::Instance()->getNormRand();
							if(r < 0.5)
							{
								singletonIndex = ii;
								break;
							}
						}
					}

				}
			}
		}
		LProposalDistributionElement* lpe = NULL;
		if(samplingMode==EINFORMED)
			lpe = distribution->findElement(atom);
		if(!nonpoly)
		{
			LogDouble probOfSample(1,false);
			LogDouble sampleWeight(1,false);
			vector<int> sampledValues;
			int domSize = atom->terms[singletonIndex]->domain.size();
			if(samplingMode==EINFORMED)
			{	
				LSamplingAlgs::Instance()->sample(domSize,1,sampleWeight,sampledValues,samplingMode,lpe);
			}
			else if(samplingMode==EBINOMIAL)
			{	
				LSamplingAlgs::Instance()->sample(domSize,1,sampleWeight,sampledValues,EBINOMIAL);
			}	
			else
				LSamplingAlgs::Instance()->sample(domSize,1,sampleWeight,sampledValues);
			if(LvrQueryUpdater::isInstanceCreated())
			{
				LvrQueryUpdater::Instance()->updateQueryValues(atom,sampledValues[0]);
			}
			LvrSingletonNormPropagation::propagateNormalizedCNF(CNF,atom,singletonIndex,sampledValues[0]);
			LogDouble mcnt = doLvApproxPartitionInformed(CNF);
			totalVal *= sampleWeight*mcnt;
		}
		else if(atom->isConstant())
		{
			LogDouble sampleWeight(1,false);
			vector<int> sampledValues;
			int numGroundings = atom->getNumberOfGroundings();
			bool setTrue=false;
			if(samplingMode==EINFORMED)
			{
				LSamplingAlgs::Instance()->sample(numGroundings,1,sampleWeight,sampledValues,samplingMode,lpe);
				if(sampledValues[0]==1)
					setTrue=true;
			}
			else
			{
				double r = (double)rand()/(double)RAND_MAX;
				if(r > 0.5)
				{
					sampledValues.push_back(1);
					setTrue=true;
				}
				else
					sampledValues.push_back(0);
				sampleWeight = LogDouble(2,false);
			}
			if(LvrQueryUpdater::isInstanceCreated())
			{
				LvrQueryUpdater::Instance()->updateQueryValues(atom,sampledValues[0]);
			}
			for(int ii=0;ii<CNF.size();ii++)
			{
				for(int jj=0;jj<CNF[ii]->atoms.size();jj++)
				{
					if(CNF[ii]->atoms[jj]->symbol->id == atom->symbol->id)
					{
						if(CNF[ii]->sign[jj] && !setTrue)
							CNF[ii]->satisfied= true;
						else if(!CNF[ii]->sign[jj] && setTrue)
							CNF[ii]->satisfied= true;
						CNF[ii]->removeAtom(jj);
						jj--;
					}
				}
				if(CNF[ii]->atoms.size()==0 && !CNF[ii]->satisfied)
				{
					removeItem(CNF,ii);
					ii--;
				}
			}
			LogDouble mcnt = doLvApproxPartitionInformed(CNF);
			totalVal *= sampleWeight*mcnt;
		}
		else if(isIsolated)
		{	
			//use Isolated terms rule
			double binCoeff = 1;
			LogDouble probOfSample(1,false);
			LogDouble sampleWeight(1,false);
			vector<int> sampledValues;
			int isolatedSize=1;
			int nonisolatedSize = 1;
			for(unsigned int jj=0;jj<isolatedTerms.size();jj++)
			{
				if(isolatedTerms[jj])
				{
					isolatedSize *= atom->terms[jj]->domain.size();
				}
				else
				{
					nonisolatedSize *=  atom->terms[jj]->domain.size();
				}
			}
			LSamplingAlgs::Instance()->sample(isolatedSize,nonisolatedSize,sampleWeight,sampledValues,samplingMode,lpe);
			if(LvrQueryUpdater::isInstanceCreated())
			{
				LvrQueryUpdater::Instance()->updateQueryValues(atom,isolatedTerms,sampledValues);
			}
			lvrNormPropagate->propagateNormalizedCNF(CNF,atom,isolatedTerms,sampledValues);
			LogDouble mcnt = doLvApproxPartitionInformed(CNF);
			totalVal *= sampleWeight*mcnt;
		}
		delete atom;
	}
	cleanup(CNF);
	return totalVal;
}



LogDouble LISApproxInference::doLvApproxPartitionInformedRB(vector<WClause*>& CNF)
{
	LogDouble totalVal(1,false);
	if(CNF.size()==0)
	{
		return totalVal;
	}	
	int powerFactor;
	bool isDecomposed = decomposeCNF(CNF,powerFactor);
	if(isDecomposed)
	{
		LogDouble mcnt = doLvApproxPartitionInformedRB(CNF);
		LogDouble val = LogDouble(1,false);
		LogDouble::LDPower(mcnt,powerFactor,val);
		totalVal = totalVal*val;
	}
	else
	{
		int singletonIndex;
		bool singleton=false;
		Atom* constatom=NULL;
		Atom* otheratom=NULL;
		Atom* singatom=NULL;
		Atom* tmpAtom = NULL;
		for(int ii=0;ii<CNF.size();ii++)	
		{
		  	bool foundatom=false;
			for(int jj=0;jj<CNF[ii]->atoms.size();jj++)
			{
				if(LvrQueryUpdater::Instance()->isNormIdInQuery(CNF[ii]->atoms[jj]->symbol->normParentId))
					continue;
				otheratom = CNF[ii]->atoms[jj];
				singleton = CNF[ii]->atoms[jj]->isSingletonAtom(singletonIndex);
				if(singleton)
				{
					singatom = CNF[ii]->atoms[jj];
					foundatom=true;
					break;
				}
				else if(CNF[ii]->atoms[jj]->isConstant())
					constatom = CNF[ii]->atoms[jj];
			}
			if(foundatom)
				break;
		}
		if(singatom)
			tmpAtom = singatom;
		else if(constatom)
			tmpAtom=constatom;
		else
			tmpAtom=otheratom;
		if(tmpAtom==NULL)
		{
			LogDouble wt1 = CNFWeight(CNF);
			//exact inference on the rest of the MLN
			rbestimator->updateRBEstimates(CNF);
			cleanup(CNF);
		    totalVal = totalVal*wt1;
			return totalVal;			
		}
		bool nonpoly;
		if(!singleton)
			nonpoly=true;
		else
			nonpoly=false;
		Atom* atom = LvrMLN::create_new_atom(tmpAtom);
		vector<bool> isolatedTerms;
		bool isIsolated = false;
		if(nonpoly)
		{
			if(!atom->isConstant())
			{
				if(samplingMode==EINFORMED)
				{
					isIsolated = LRulesUtil::computeIsolatedTerms(atom,CNF,isolatedTerms);
				}
				if(!isIsolated)
				{
					nonpoly=false;
					/*int maxdomid=0;
					int maxdomsize=0;
					for( int ii=0;ii<atom->terms.size();ii++)
					{
						if(atom->terms[ii]->domain.size()>maxdomsize)
						{
							maxdomid=ii;
							maxdomsize=atom->terms[ii]->domain.size();
						}
					}
					singletonIndex = maxdomid;
					*/
					bool first=true;
					for( int ii=0;ii<atom->terms.size();ii++)
					{
						if(atom->terms[ii]->domain.size()>1)
						{
							if(first)
							{
								singletonIndex = ii;
								first=false;
							}
							double r = LvRandomGenUtil::Instance()->getNormRand();
							if(r < 0.5)
							{
								singletonIndex = ii;
								break;
							}
						}
					}

				}
			}
		}
		LProposalDistributionElement* lpe = NULL;
		if(samplingMode==EINFORMED)
			lpe = distribution->findElement(atom);
		if(!nonpoly)
		{
			LogDouble probOfSample(1,false);
			LogDouble sampleWeight(1,false);
			vector<int> sampledValues;
			int domSize = atom->terms[singletonIndex]->domain.size();
			if(samplingMode==EINFORMED)
			{	
				LSamplingAlgs::Instance()->sample(domSize,1,sampleWeight,sampledValues,samplingMode,lpe);
			}
			else if(samplingMode==EBINOMIAL)
			{	
				LSamplingAlgs::Instance()->sample(domSize,1,sampleWeight,sampledValues,EBINOMIAL);
			}	
			else
				LSamplingAlgs::Instance()->sample(domSize,1,sampleWeight,sampledValues);
			LvrSingletonNormPropagation::propagateNormalizedCNF(CNF,atom,singletonIndex,sampledValues[0]);
			LogDouble mcnt = doLvApproxPartitionInformedRB(CNF);
			totalVal *= sampleWeight*mcnt;
		}
		else if(atom->isConstant())
		{
			LogDouble sampleWeight(1,false);
			vector<int> sampledValues;
			int numGroundings = atom->getNumberOfGroundings();
			bool setTrue=false;
			if(samplingMode==EINFORMED)
			{
				LSamplingAlgs::Instance()->sample(numGroundings,1,sampleWeight,sampledValues,samplingMode,lpe);
				if(sampledValues[0]==1)
					setTrue=true;
			}
			else
			{
				double r = (double)rand()/(double)RAND_MAX;
				if(r > 0.5)
				{
					sampledValues.push_back(1);
					setTrue=true;
				}
				else
					sampledValues.push_back(0);
				sampleWeight = LogDouble(2,false);
			}
			for(int ii=0;ii<CNF.size();ii++)
			{
				for(int jj=0;jj<CNF[ii]->atoms.size();jj++)
				{
					if(CNF[ii]->atoms[jj]->symbol->id == atom->symbol->id)
					{
						if(CNF[ii]->sign[jj] && !setTrue)
							CNF[ii]->satisfied= true;
						else if(!CNF[ii]->sign[jj] && setTrue)
							CNF[ii]->satisfied= true;
						CNF[ii]->removeAtom(jj);
						jj--;
					}
				}
				if(CNF[ii]->atoms.size()==0 && !CNF[ii]->satisfied)
				{
					removeItem(CNF,ii);
					ii--;
				}
			}
			LogDouble mcnt = doLvApproxPartitionInformedRB(CNF);
			totalVal *= sampleWeight*mcnt;
		}
		else if(isIsolated)
		{	
			//use Isolated terms rule
			double binCoeff = 1;
			LogDouble probOfSample(1,false);
			LogDouble sampleWeight(1,false);
			vector<int> sampledValues;
			int isolatedSize=1;
			int nonisolatedSize = 1;
			for(unsigned int jj=0;jj<isolatedTerms.size();jj++)
			{
				if(isolatedTerms[jj])
				{
					isolatedSize *= atom->terms[jj]->domain.size();
				}
				else
				{
					nonisolatedSize *=  atom->terms[jj]->domain.size();
				}
			}
			LSamplingAlgs::Instance()->sample(isolatedSize,nonisolatedSize,sampleWeight,sampledValues,samplingMode,lpe);
			lvrNormPropagate->propagateNormalizedCNF(CNF,atom,isolatedTerms,sampledValues);
			LogDouble mcnt = doLvApproxPartitionInformedRB(CNF);
			totalVal *= sampleWeight*mcnt;
		}
		delete atom;
	}
	cleanup(CNF);
	return totalVal;
}


/*
LogDouble LISApproxInference::doLvApproxPartitionInformed(vector<WClause*>& CNF1)
{
	LogDouble totalVal(1,false);
	if(CNF1.size()==0)
		return totalVal;
	vector<vector<WClause*> > dCNFs;
	LClusterUtil::seperateDisjointCNF(CNF1,dCNFs);
	cleanup(CNF1);
	for(unsigned int t=0;t<dCNFs.size();t++)
	{
		vector<WClause*> CNF = dCNFs[t];
		int powerFactor;
		bool isDecomposed = decomposeCNF(CNF,powerFactor);
		if(isDecomposed)
		{
			LogDouble mcnt = doLvApproxPartitionInformed(CNF);
			LogDouble val = LogDouble(1,false);
			LogDouble::LDPower(mcnt,powerFactor,val);
			totalVal = totalVal*val;
			continue;
		}		
		Atom* tmpAtom = heuristics->getAtomToSplit(CNF);	
		if(tmpAtom==NULL)
		{	
			LogDouble wt1 = CNFWeight(CNF);
			cleanup(CNF);
		    totalVal = totalVal*wt1;
			continue;
			
		}
		Atom* atom = LvrMLN::create_new_atom(tmpAtom);
		LProposalDistributionElement* lpe = lpe = distribution->findElement(atom);
		//check for nonpoly cases
		//check for self joins
		bool selfJoined = false;
		for(int m=0;m<CNF.size();m++)
		{
			if(CNF.at(m)->isSelfJoinedOnAtom(atom))
			{
				selfJoined = true;
				break;
			}
		}
		int singletonIndex;
		bool singleton = atom->isSingletonAtom(singletonIndex);
		bool nonpoly = false;
		if(!singleton)
			nonpoly = true;
		else if(selfJoined)
		{
			//check for blocking rule
			bool blocked = LRulesUtil::isBlocked(atom,singletonIndex,CNF);
			if(!blocked)
				nonpoly=true;
		}
		if(!nonpoly)
		{
			LogDouble probOfSample(1,false);
			LogDouble sampleWeight(1,false);
			vector<int> sampledValues;
			int domSize = atom->terms[singletonIndex]->domain.size();
			LSamplingAlgs::Instance()->sample(domSize,1,sampleWeight,sampledValues,EINFORMED,lpe);
			if(LvrQueryUpdater::isInstanceCreated())
			{
				LvrQueryUpdater::Instance()->updateQueryValues(atom,sampledValues[0]);
			}
			LvrSingletonNormPropagation::propagateNormalizedCNF(CNF,atom,singletonIndex,sampledValues[0]);
			LogDouble mcnt = doLvApproxPartitionInformed(CNF);
			totalVal *= sampleWeight*mcnt;

		}
		else
		{
			//use Isolated terms rule
			vector<bool> isolatedTerms;
			bool isIsolated = LRulesUtil::computeIsolatedTerms(atom,CNF,isolatedTerms);
			double binCoeff = 1;
			LogDouble probOfSample(1,false);
			LogDouble sampleWeight(1,false);
			vector<int> sampledValues;
			if(isIsolated)
			{
				int isolatedSize=1;
				int nonisolatedSize = 1;
				for(unsigned int jj=0;jj<isolatedTerms.size();jj++)
				{
					if(isolatedTerms[jj])
					{
						isolatedSize *= atom->terms[jj]->domain.size();
					}
					else
					{
						nonisolatedSize *=  atom->terms[jj]->domain.size();
					}
				}
				LSamplingAlgs::Instance()->sample(isolatedSize,nonisolatedSize,sampleWeight,sampledValues,samplingMode,lpe);
				if(LvrQueryUpdater::isInstanceCreated())
				{
					LvrQueryUpdater::Instance()->updateQueryValues(atom,isolatedTerms,sampledValues);
				}
			}
			else
			{
				int numGroundings = atom->getNumberOfGroundings();
				LSamplingAlgs::Instance()->sample(numGroundings,1,sampleWeight,sampledValues,EINFORMED,lpe);
				if(LvrQueryUpdater::isInstanceCreated())
				{
					LvrQueryUpdater::Instance()->updateQueryValues(atom,sampledValues[0]);
				}
		
			}
			lvrNormPropagate->propagateNormalizedCNF(CNF,atom,isolatedTerms,sampledValues);
			LogDouble mcnt = doLvApproxPartitionInformed(CNF);
			totalVal *= sampleWeight*mcnt;
		}

		delete atom;
		cleanup(CNF);
	}
	return totalVal;
}
*/
/*
LogDouble LISApproxInference::doLvApproxPartitionBinomial(vector<WClause*>& CNF1)
{
	LogDouble totalVal(1,false);
	if(CNF1.size()==0)
		return totalVal;
	int completedCount=0;
	for(unsigned int i=0;i<CNF1.size();i++)
	{
		if(CNF1[i]->atoms.size()==0 || CNF1[i]->satisfied)
		{
			completedCount++;
		}
	}
	if(completedCount == CNF1.size())
	{
		LogDouble weight = CNFWeight(CNF1);
		cleanup(CNF1);
		return weight;
	}
	vector<vector<WClause*> > dCNFs;
	LClusterUtil::seperateDisjointCNF(CNF1,dCNFs);
	cleanup(CNF1);
	for(unsigned int t=0;t<dCNFs.size();t++)
	{
		vector<WClause*> CNF = dCNFs[t];
		int powerFactor;
		bool isDecomposed = decomposeCNF(CNF,powerFactor);
		if(isDecomposed)
		{
			LogDouble mcnt = doLvApproxPartitionBinomial(CNF);
			LogDouble val = LogDouble(1,false);
			LogDouble::LDPower(mcnt,powerFactor,val);
			totalVal = totalVal*val;
			continue;
		}		
		Atom* tmpAtom = heuristics->getAtomToSplit(CNF);	
		if(tmpAtom==NULL)
		{	
			LogDouble wt1 = CNFWeight(CNF);
			cleanup(CNF);
		    totalVal*wt1;
			continue;
			
		}
		Atom* atom = LvrMLN::create_new_atom(tmpAtom);
		LogDouble probOfSample(1,false);
		LogDouble sampleWeight(1,false);
		int numGroundings = atom->getNumberOfGroundings();
		//select the largest domain grounding
		int maxDomainTerm = 0;
		int maxDomain = -1;
		for(unsigned int jj=0;jj<atom->terms.size();jj++)
		{
			int up = atom->terms[jj]->domain.size();
			if(up > maxDomain)
			{
				maxDomain = atom->terms[jj]->domain.size();
				maxDomainTerm = jj;
			}
		}
		numGroundings = atom->terms[maxDomainTerm]->domain.size();
		vector<int> sampledValues;
		if(samplingMode==EBINOMIAL)
		{
			LSamplingAlgs::Instance()->sample(numGroundings,1,sampleWeight,sampledValues,EBINOMIAL);
		}
		else
			LSamplingAlgs::Instance()->sample(numGroundings,1,sampleWeight,sampledValues);
		if(LvrQueryUpdater::isInstanceCreated())
		{
			LvrQueryUpdater::Instance()->updateQueryValues(atom,sampledValues[0]);
		}
		LvrSingletonNormPropagation::propagateNormalizedCNF(CNF,atom,maxDomainTerm,sampledValues[0]);
		LogDouble mcnt = doLvApproxPartitionBinomial(CNF);
		totalVal *= sampleWeight*mcnt;
		delete atom;
		cleanup(CNF);
	}
	return totalVal;
}
*/
/*
LogDouble LISApproxInference::estimatePartitionFunction(LvrParams* params,LProposalDistribution* distribution)
{
	cout<<"Lifted Importance Sampling for estimating Z..."<<endl;
	time_t start;
	time(&start);
	cout<<start<<endl;
	cout<<"Time ="<<ctime(&start)<<endl;
	LogDouble ZApprox;
	int iterations = 0;
	if(distribution)
		setProposalDistribution(distribution);
	setSamplingMode(params->samplingMode);
	if(params->learningRate <=0 )
		params->learningRate = LEARNINGRATE;
	if(params->proposalUpdateInterval <= 0)
		params->proposalUpdateInterval = PROPOSALUPDATEINTERVAL;
	time_t autotime;
	time(&autotime);
	int itime=0;
	LogDouble totalZ;
	while(1)
	{
		//make a copy of normalized clauses
		vector<WClause*> clauses;
		LvrMLN::copyAllClauses(mln.clauses,clauses);
		LogDouble currZ;
		currZ = doLvApproxPartitionInformed(clauses);
		totalZ = totalZ + currZ;
		iterations++;
		if(params->samplingMode == EINFORMED)
		{
			if(iterations % params->proposalUpdateInterval == 0)
				distribution->updateDistributions(params->learningRate);
		}
		//take average
		ZApprox = totalZ / LogDouble(iterations,false);
		if(params->autotest)
		{
			time_t curr;
			time(&curr);
			int seconds = difftime(curr,autotime);			
			if(seconds > params->printinterval)
			{
				time(&autotime);
				string outname(params->fileprefix);
				itime += params->iprintinterval;
				stringstream st;
				st<<itime;
				outname.append(st.str());
				outname.append(".dat");
				ofstream out(outname.c_str());
				out<<iterations<<" "<<ZApprox.value<<endl;
				out.close();
				if(itime >= (int)params->endtime)
					return ZApprox;
			}
		}
		else
		{
			time_t curr;
			time(&curr);
			int seconds = difftime(curr,start);
			if(iterations%PRINTRESULTSINTERVAL == 0 ||
					seconds > params->maxSeconds || iterations >= params->maxSteps)
			{
				cout<<"iteration="<<iterations<<", currZ : ";
				currZ.printValue();
				cout<<", ZApprox : ";
				ZApprox.printValue();
				cout<<endl;
				if(seconds > params->maxSeconds || iterations >= params->maxSteps)
				{
					return ZApprox;
				}
			}
		}
	}
}
*/


LogDouble LISApproxInference::estimatePartitionFunction(LvrParams* params,LProposalDistribution* distribution)
{
	cout<<"Lifted Importance Sampling for estimating Z..."<<endl;
	time_t start;
	time(&start);
	cout<<"Time ="<<ctime(&start)<<endl;
	LogDouble ZApprox;
	int iterations = 0;
	if(distribution)
		setProposalDistribution(distribution);
	setSamplingMode(params->samplingMode);
	if(params->learningRate <=0 )
		params->learningRate = LEARNINGRATE;
	if(params->proposalUpdateInterval <= 0)
		params->proposalUpdateInterval = PROPOSALUPDATEINTERVAL;
	time_t autotime;
	int itime=0;
	LogDouble totalZ;
	int samplestothrow=25;
	LogDouble rejectionmax;
	int initialiters=0;
	vector<int> rejarr;
	while(1)
	{
		//make a copy of normalized clauses
		vector<WClause*> clauses;
		LvrMLN::copyAllClauses(mln.clauses,clauses);
		LogDouble currZ;
		currZ = doLvApproxPartitionInformed(clauses);
		initialiters++;
		if(initialiters <= samplestothrow)
		{
			rejarr.push_back((int)currZ.value);
			if(initialiters==samplestothrow)
			{
				sort(rejarr.begin(),rejarr.end());
				int val = rejarr[rejarr.size()/2];
				rejectionmax = LogDouble(val,true); 
			}
			continue;
		}
		//accept/reject sample
		LogDouble p1(1,false);
		if(!rejectionmax.is_zero)
		{
			p1 =  currZ/rejectionmax;
			double r = LvRandomGenUtil::Instance()->getNormRand();
			LogDouble p2(r,false);
			if(p2 > p1)
			{
				//reject sample
				continue;
			}
		}
		if(p1 > LogDouble(1,false))
			p1 = LogDouble(1,false);
		if(iterations==0)
			time(&autotime);
		iterations++;
		currZ = currZ/p1;
		totalZ = totalZ + currZ;
		iterations++;
		if(params->samplingMode == EINFORMED)
		{
			if(iterations % params->proposalUpdateInterval == 0)
				distribution->updateDistributions(params->learningRate);
		}
		//take average
		ZApprox = totalZ/LogDouble(iterations,false);
		if(params->autotest)
		{
			time_t curr;
			time(&curr);
			int seconds = difftime(curr,autotime);			
			if(seconds > params->printinterval)
			{
				time(&autotime);
				string outname(params->fileprefix);
				itime += params->iprintinterval;
				stringstream st;
				st<<itime;
				outname.append(st.str());
				outname.append(".dat");
				ofstream out(outname.c_str());
				out<<iterations<<" "<<ZApprox.value<<endl;
				out.close();
				if(itime >= (int)params->endtime)
					return ZApprox;
			}
		}
		else
		{
			time_t curr;
			time(&curr);
			int seconds = difftime(curr,start);
			if(iterations%PRINTRESULTSINTERVAL == 0 ||
					seconds > params->maxSeconds || iterations >= params->maxSteps)
			{
				cout<<"iteration="<<iterations<<", currZ : ";
				currZ.printValue();
				cout<<", ZApprox : ";
				ZApprox.printValue();
				cout<<endl;
				if(seconds > params->maxSeconds || iterations >= params->maxSteps)
				{
					return ZApprox;
				}
			}
		}
	}
}

void LISApproxInference::estimateApproxMarginals(LvrParams* params,LProposalDistribution* distribution)
{
	if(params->lisRB)
	{
		rbestimator = new LvgRBEstimator(mln);
		cout<<"Running Rao Blackwellized IS"<<endl;
	}
	
	setProposalDistribution(distribution);
	cout<<"Estimating Marginals using Lifted Importance Sampling..."<<endl;
	time_t start;
	time(&start);
	cout<<"Time ="<<ctime(&start)<<endl;
	LogDouble ZApprox;
	int iterations = 0;
	int printInterval = PRINTRESULTSINTERVAL;
	setSamplingMode(params->samplingMode);
	LogDouble totalWeight;
	if(params->learningRate <=0 )
		params->learningRate = LEARNINGRATE;
	if(params->proposalUpdateInterval <= 0)
		params->proposalUpdateInterval = PROPOSALUPDATEINTERVAL;
	int samplestothrow=25;
	LogDouble rejectionmax;
	int initialiters=0;
	vector<int> rejarr;
	while(1)
	{
		//make a copy of normalized clauses
		vector<WClause*> clauses;
		LvrMLN::copyAllClauses(mln.clauses,clauses);
		LogDouble currWt;
		if(params->lisRB)
			currWt = doLvApproxPartitionInformedRB(clauses);
		else
			currWt = doLvApproxPartitionInformed(clauses);
		initialiters++;
		if(initialiters <= samplestothrow)
		{
			LvrQueryUpdater::Instance()->resetallupdateflags();
			rejarr.push_back((int)currWt.value);
			if(initialiters==samplestothrow)
			{
				sort(rejarr.begin(),rejarr.end());
				int val = rejarr[rejarr.size()/2];
				rejectionmax = LogDouble(val,true); 
			}
			continue;
		}
		//accept/reject sample
		LogDouble p1(1,false);
		if(!rejectionmax.is_zero)
		{
			p1 =  currWt/rejectionmax;
			//cout<<currWt.value<<" "<<p1<<" "<<rejectionmax.value<<endl;
			double r = LvRandomGenUtil::Instance()->getNormRand();
			LogDouble p2(r,false);
			if(p2 > p1)
			{
				//reject sample
				//cout<<"reject"<<endl;
				continue;
			}
		}
		if(p1 > LogDouble(1,false))
			p1 = LogDouble(1,false);
		//cout<<"accept"<<endl;
		iterations++;
		currWt = currWt/p1;
		//update the cumulative weight
		totalWeight += currWt;
		//totalWeight += currWt;
		//update approximation
		if(params->lisRB)
		{
			LvrQueryUpdater::Instance()->updateGibbsDontCareRB();
			LvrQueryUpdater::Instance()->updateAllImportanceWeightsRB(currWt);
		}
		else
		{
			LvrQueryUpdater::Instance()->updateDontCare();
			LvrQueryUpdater::Instance()->updateAllImportanceWeights(currWt);
		}
		
		if(params->samplingMode == EINFORMED)
		{
			if(iterations % params->proposalUpdateInterval == 0)
			{
				distribution->updateDistributions(params->learningRate);
			}
		}
		time_t curr;
		time(&curr);
		int seconds = difftime(curr,start);
		if(iterations % PRINTRESULTSINTERVAL == 0 ||
				seconds > params->maxSeconds || iterations >= params->maxSteps)
		{
			LvrQueryUpdater::Instance()->writeToFile(totalWeight);
			cout<<"iteration "<< iterations<<endl;
			cout<<"Z-curr = ";
			currWt.printValue();
			cout<<endl;
			cout<<"cumulative-Z = ";
			totalWeight.printValue();
			cout<<endl;
			if(seconds > params->maxSeconds || iterations >= params->maxSteps)
				break;
		}
	}
}

LISApproxInference::LISApproxInference(LvrMLN& mln_): mln(mln_),rbestimator(NULL),raoblackwelize(false)
{
	decomposer = new LDecomposer(mln);
	heuristics = new LHeuristics(*decomposer,mln);
	lvrNormPropagate = new LvrNormPropagation(mln);
}
