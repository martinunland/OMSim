#include "OMSimRunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

#include <ctime>
#include <sys/time.h>

#include "OMSimAnalysisManager.hh"
#include <time.h>
#include <sys/time.h>
#include "OMSimCommandArgsTable.hh"
#include "OMSimLogger.hh"


OMSimRunAction::OMSimRunAction(){}
OMSimRunAction::~OMSimRunAction(){}

void OMSimRunAction::BeginOfRunAction(const G4Run*)
{  startingtime = clock() / CLOCKS_PER_SEC;
 log_info("Start of run action");
	
}

void OMSimRunAction::EndOfRunAction(const G4Run*)
{
log_info("End of run action");
std::string lFileName = OMSimCommandArgsTable::getInstance().get<std::string>("output_file");

OMSimAnalysisManager& lAnalysisManager = OMSimAnalysisManager::getInstance();
lAnalysisManager.datafile.open(lFileName.c_str(), std::ios::out|std::ios::app);
lAnalysisManager.WriteAccept();
// 	Close output data file
lAnalysisManager.datafile.close();
lAnalysisManager.Reset();
//double finishtime=clock() / CLOCKS_PER_SEC;
//G4cout << "Computation time: " << finishtime-startingtime << " seconds." << G4endl;
}

