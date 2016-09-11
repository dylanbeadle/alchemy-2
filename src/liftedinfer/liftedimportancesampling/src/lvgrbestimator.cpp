#include "propositionalresolution.h"
#include "randomgenutil.h"
#ifndef _MSC_VER 
#include <unistd.h>
#endif
#include "lvgrbestimator.h"
#include "queryupdater.h"
#include "cleanuputils.h"
LogDouble LvgRBEstimator::CNFWeight(vector<WClause*>& CNF)
{
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


bool LvgRBEstimator::decomposeCNF(vector<WClause*>& CNF,int& powerFactor)
{
	vector<Decomposer*> decomposer_list;
	decomposer->find_decomposer(CNF,decomposer_list);
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
				//make a copy of the pre-decomp domain
				for(unsigned int jj=0;jj<decomposer_list[i]->decomposer_terms[j]->domain.size();jj++)
					decomposer_list[i]->decomposer_terms[j]->origDomain.push_back(decomposer_list[i]->decomposer_terms[j]->domain[jj]);
				//Erase all but one grounding
				decomposer_list[i]->decomposer_terms[j]->domain.erase
					(decomposer_list[i]->decomposer_terms[j]->domain.begin(),
					decomposer_list[i]->decomposer_terms[j]->domain.begin()+decomposer_list[i]->decomposer_terms[j]->domain.size()-1);
				
			}
		}
	}
	cleanup(decomposer_list);
	return true;
}

void LvgRBEstimator::addTreeNode(Atom* atom,vector<LogDouble> weights,vector<LogDouble> binCoeffs,
	LogDouble norm,NodeType nodeType,int nodeId, int parentId, 
	int parentBranchNo,vector<LogDouble>* branchWeights,LogDouble nodeValue)
{
	Atom* newAtom = LvrMLN::create_new_atom(atom);
	LPTPNode* nd = new LPTPNode();
	nd->atom = newAtom;
	nd->zValues = weights;
	nd->norm = norm;
	nd->binCoeff = binCoeffs;
	nd->type = nodeType;
	nd->id = nodeId;
	nd->parentId = parentId;
	nd->parentBranchNo = parentBranchNo;
	if(branchWeights!=NULL)
		nd->branchWeights = (*branchWeights);
	ptpTree->treeNodes.push_back(nd);
	nd->nodeValue = nodeValue;
}

LPTPNode* LvgRBEstimator::addTreeNode(NodeType nodeType,int nodeId, int parentId, int parentBranchNo,int powerFactor,LogDouble nodeValue)
{
	LPTPNode* nd = new LPTPNode();
	nd->type = nodeType;
	nd->id = nodeId;
	nd->parentId = parentId;
	nd->parentBranchNo = parentBranchNo;
	nd->powerFactor = powerFactor;
	nd->nodeValue = nodeValue;
	ptpTree->treeNodes.push_back(nd);
	return nd;
}


LogDouble LvgRBEstimator::weightedModelCount(vector<WClause*>& CNF,int parentId,int parentBranchNo)
{
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
	vector<int> powerFactors;
	vector<bool> isDecomposedCNF;
	vector<vector<WClause*> > decomposedList;
	LClusterUtil::seperateDisjointCNF(CNF,decomposedList);
	//CNF can be cleaned up since split into clusters
	cleanup(CNF);
	LogDouble totalVal(1,false);

	bool isParentChanged=false;
	int changedParent;
	if(decomposedList.size() > 1)
	{
		addTreeNode(EDISJOINT,id,parentId,parentBranchNo);
		id++;
		isParentChanged = true;
		changedParent = id-1;
	}

	for(unsigned int t=0;t<decomposedList.size();t++)
	{
		if(isParentChanged)
		{
		parentId = changedParent;
		parentBranchNo = t;
		}

		int powerFactor;
		vector<map<int,int> > decomposerMappings;
		vector<int> domainPowerFactors;
		bool isDecomposed = decomposeCNF(decomposedList[t],powerFactor);
		if(isDecomposed)
		{
			int idToUse = id++;
			//use power rule
			LogDouble mcnt = weightedModelCount(decomposedList[t],idToUse,0);
			LogDouble val = LogDouble(1,false);
			LogDouble::LDPower(mcnt,powerFactor,val);
			totalVal = totalVal*val;
			LPTPNode* node = addTreeNode(EDECOMPOSER,idToUse,parentId,parentBranchNo,powerFactor);
			node->decomposerMappings = decomposerMappings;
			node->domainPowerFactors = domainPowerFactors;
			continue;
		}
		else
		{
			vector<vector<WClause*> > tempCNFList;
			//choose an atom
			Atom* tmpatom = NULL;
			for(unsigned int ii=0;ii<decomposedList[t].size();ii++)
			{
				if(decomposedList[t][ii]->atoms.size() > 0)
				{
					tmpatom = decomposedList[t][ii]->atoms[0];
					break;
				}
			}
			if(tmpatom==NULL)
			{	
				LogDouble wt1 = CNFWeight(decomposedList[t]);
				cleanup(decomposedList[t]);
				totalVal *= wt1;
				continue;
			}
			Atom* atom = LvrMLN::create_new_atom(tmpatom);
			int singletonIndex;
			bool singleton = atom->isSingletonAtom(singletonIndex);
			/*if(!singleton && !atom->isConstant())
			{
				//choose an index randomly
				for(unsigned int jj=0;jj<atom->terms.size();jj++)
				{
					if(atom->terms[jj]->domain.size() > 1)
						singletonIndex = jj;
					double r = LvRandomGenUtil::Instance()->getNormRand();
					if(r < 0.5)
						break;
				}
			}
			*/
			//check for self joins
			bool selfJoined = false;
			for(unsigned int m=0;m<decomposedList[t].size();m++)
			{
				if(decomposedList[t].at(m)->isSelfJoinedOnAtom(atom))
				{
					selfJoined = true;
					break;
				}
			}
			if(atom->isConstant())
			{
				int idToUse = id++;
				//propositional atom
				if(!selfJoined)
				{
					PropositionalResolution::doPropositionalSplit(decomposedList[t],atom,tempCNFList);
					LogDouble w1 = weightedModelCount(tempCNFList[0],idToUse,0);
					LogDouble w2 = weightedModelCount(tempCNFList[1],idToUse,1);
					totalVal *= (w1+w2);
					vector<LogDouble> wts;
					wts.push_back(w1);
					wts.push_back(w2);
					LogDouble norm = w1+w2;
					vector<LogDouble> bincoefs;
					bincoefs.push_back(LogDouble(1,false));
					bincoefs.push_back(LogDouble(1,false));
					addTreeNode(atom,wts,bincoefs,norm,ESPLIT,idToUse,parentId,parentBranchNo);
					delete atom;
					cleanup(tempCNFList[0]);
					cleanup(tempCNFList[1]);
					cleanup(decomposedList[t]);

				}
				else
				{
					PropositionalResolution::doPropositionalSplitSelfJoins(decomposedList[t],atom,tempCNFList);
					LogDouble w1 = weightedModelCount(tempCNFList[0],idToUse,0);
					LogDouble w2 = weightedModelCount(tempCNFList[1],idToUse,1);
					totalVal *= (w1+w2);
					vector<LogDouble> wts;
					wts.push_back(w1);
					wts.push_back(w2);
					LogDouble norm = w1+w2;
					vector<LogDouble> bincoefs;
					bincoefs.push_back(LogDouble(1,false));
					bincoefs.push_back(LogDouble(1,false));
					addTreeNode(atom,wts,bincoefs,norm,ESPLIT,idToUse,parentId,parentBranchNo);

					cleanup(tempCNFList[0]);
					cleanup(tempCNFList[1]);
					tempCNFList.clear();
					delete atom;
					cleanup(decomposedList[t]);
				}
			}
			else
			{
				int idToUse = id++;
				int totalSize = atom->terms[singletonIndex]->domain.size();
				LogDouble v;
				vector<LogDouble> zvals;
				vector<LogDouble> binceoffs;
				LogDouble norm;
				for(unsigned int i=0;i<=totalSize;i++)
				{
					vector<WClause*> clausesCopy;
					LvrMLN::copyAllClauses(decomposedList[t],clausesCopy);
					LvrSingletonNormPropagation::propagateNormalizedCNF(clausesCopy,
						atom,singletonIndex,i);
					LogDouble binCoeff = LogDouble::Binomial(totalSize,i);
					LogDouble mcnt = weightedModelCount(clausesCopy,idToUse,i);
					cleanup(clausesCopy);
					v = v + binCoeff*mcnt;
					zvals.push_back(mcnt);
					binceoffs.push_back(binCoeff);
					norm = norm + binCoeff*mcnt;
				}

				totalVal *= v;
				addTreeNode(atom,zvals,binceoffs,norm,ESPLIT,idToUse,parentId,parentBranchNo);
				delete atom;
				cleanup(decomposedList[t]);
			}
		}
	}
	return totalVal;
}


LvgRBEstimator::LvgRBEstimator(LvrMLN& mln_): mln(mln_)
{
	decomposer = new LDecomposer(mln);
	ptpTree = new LPTPTree();
}

void LvgRBEstimator::updateRBEstimates(vector<WClause*>& CNF)
{
	ptpTree->cleanTree();
	id=0;
	weightedModelCount(CNF,-1,-1);
	for(unsigned int i=0;i<ptpTree->treeNodes.size();i++)
	{
		if(ptpTree->treeNodes[i]->type==ESPLIT)
		{
			vector<LogDouble> probs;
			for(unsigned int j=0;j<ptpTree->treeNodes[i]->zValues.size();j++)
			{
				probs.push_back(ptpTree->treeNodes[i]->getProbability(j));
			}
			//update the query
			LvrQueryUpdater::Instance()->updateISRBEstimates(ptpTree->treeNodes[i]->atom,probs);
		}
	}
}