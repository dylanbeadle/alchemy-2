#include "liftedalgshandler.h"
#include "wmconvertor.h"
#include <sstream>
using namespace std;

LiftedAlgsHandler::LiftedAlgsHandler(LvrMLN& mln_):mln(mln_)
{
	ptpSearch = new LPTPSearch(mln);
	clusterCreator = new LClusterCreator(&mln);
	proposalConstructor = new LProposalConstructor(mln);
	gibbsProcessHandler = new GibbsProcessHandler(mln);
}

LiftedAlgsHandler::~LiftedAlgsHandler()
{
	delete ptpSearch;
	delete clusterCreator;
	delete proposalConstructor;
	delete gibbsProcessHandler;
}

void LiftedAlgsHandler::LBGApproxMarginals(LvrParams* params)
{
	params->autotest = 0;
	char* autoTestS = getenv("AUTO_TESTING");
	if (autoTestS) {
		params->autotest = atoi(autoTestS);
	}
	char* printervalS = getenv("PRINT_INTERVAL");
	if (printervalS) {
		params->printinterval = atoi(printervalS);
		params->iprintinterval = atoi(printervalS);
	}
	char* endtimeS = getenv("END_TIME");
	if (endtimeS) {
		params->endtime = atoi(endtimeS);
	}
	char* algp = getenv("FILEPREFIX");
	if(algp)
		params->fileprefix = string(algp);

	char* burnS = getenv("BURN_IN");
	if (burnS) {
		params->burnMaxSteps = atoi(burnS);
	}
	params->gibbsRB = true;
	gibbsProcessHandler->startLBGForMar(params);
}

void LiftedAlgsHandler::LISApproxPartition(LvrParams* params)
{
	/*if(!params->isWeightLearning)
	{
		if(params->samplingMode == EINFORMED)
		{
			buildProposal(params);
		}
	}
	*/
	//Auto testing env vars
	params->autotest = 0;
	char* autoTestS = getenv("AUTO_TESTING");
	if (autoTestS) {
		params->autotest = atoi(autoTestS);
	}
	char* printervalS = getenv("PRINT_INTERVAL");
	if (printervalS) {
		params->printinterval = atoi(printervalS);
		params->iprintinterval = atoi(printervalS);
	}
	char* endtimeS = getenv("END_TIME");
	if (endtimeS) {
		params->endtime = atoi(endtimeS);
	}
	char* algp = getenv("FILEPREFIX");
	if(algp)
		params->fileprefix = string(algp);
	char* proposaltypeS = getenv("PROPOSAL_TYPE");
	if (proposaltypeS) {
		string proposaltype(proposaltypeS);
		if(proposaltype.compare("INFORMED")==0)
		{
			params->samplingMode = EINFORMED;
		}
		else if(proposaltype.compare("BINOMIAL")==0)
		{
			params->samplingMode = EBINOMIAL;
		}
		else if(proposaltype.compare("UNIFORM")==0)
		{
			params->samplingMode = EUNIFORM;
		}

	}
	char* learnrateS = getenv("LEARNING_RATE");
	if (learnrateS) {
		string lrate(learnrateS);
		stringstream st(lrate);
		st>>params->learningRate;
	}

	proposalConstructor->startPartitionFunction(params);
}

void LiftedAlgsHandler::LISApproxMarginals(LvrParams* params)
{
	/*
	if(!params->isWeightLearning)
	{
		if(params->samplingMode == EINFORMED)
		{
			buildProposal(params);
		}
	}
	*/
	proposalConstructor->startMARInference(params);
}

void LiftedAlgsHandler::buildProposal(LvrParams* params)
{
	proposalConstructor->startConstruction(params);
	cout<<"done construction"<<endl;
}

void LiftedAlgsHandler::WMPTPZApprox(LvrParams* params)
{
	//convert to WM PTP
	cout<<"Pre-processing..Converting to Weighted Model Counting..."<<endl;
	LWMConvertor* lw = new LWMConvertor(mln);
	lw->convertMLN();
	delete lw;		
	cout<<"WM Converted MLN"<<endl;
	for(int i=0;i<mln.clauses.size();i++)
		mln.clauses[i]->print();
	cout<<endl;
	LFileDump::dumpMLNToFile(&mln);
	cout<<"done writing dump files ("<<MLNDUMPFILENAME<<","<<SYMBOLSDUMPFILENAME<<")"<<endl;
	ptpSearch->startApproxWeightedModelCounting(params);
}


void LiftedAlgsHandler::WMPTPZExact(LvrParams* params)
{
	cout<<"Pre-processing..Converting to Weighted Model Counting..."<<endl;
	//convert to WM PTP
	LWMConvertor* lw = new LWMConvertor(mln);
	lw->convertMLN();
	delete lw;		
	cout<<"WM Converted MLN"<<endl;
	for(int i=0;i<mln.clauses.size();i++)
		mln.clauses[i]->print();
	cout<<endl;
	LFileDump::dumpMLNToFile(&mln);
	cout<<"done writing dump files ("<<MLNDUMPFILENAME<<","<<SYMBOLSDUMPFILENAME<<")"<<endl;
	cout<<"Exact Z = ";
	ptpSearch->startExactWeightedModelCounting(params).printValue();
	cout<<endl;
}


