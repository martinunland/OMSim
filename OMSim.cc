#include "OMSimDetectorConstruction.hh"
#include "OMSimPhysicsList.hh"
#include "OMSimPrimaryGeneratorAction.hh"
#include "OMSimRunAction.hh"
#include "OMSimEventAction.hh"
#include "OMSimTrackingAction.hh"
#include "OMSimSteppingAction.hh"
#include "OMSimSteppingVerbose.hh"
#include "OMSimAnalysisManager.hh"
#include "OMSimPMTResponse.hh"
#include "OMSimCommandArgsTable.hh"
#include "OMSimUIinterface.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4ThreeVector.hh"
#include "G4Navigator.hh"

#include "G4UItcsh.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh" //xxx

#include "argtable2.h"
#include <ctime>
#include <sys/time.h>

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>

#include <cmath> // for abs() of doubles
// since Geant4.10: include units manually
#include "G4SystemOfUnits.hh"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <TGraph.h>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

std::vector<int> arange(int start, int stop, int step) {
    std::vector<int> values;
    for (int value = start; value < stop; value += step) {
        values.push_back(value);
    }
    return values;
}

int OMSim()
{
	struct timeval time_for_randy;
	gettimeofday(&time_for_randy, NULL);

	long randseed = time_for_randy.tv_sec + 4294 * time_for_randy.tv_usec;
	CLHEP::HepRandom::setTheEngine(new CLHEP::RanluxEngine(randseed, 3));

	std::stringstream command;

	G4RunManager *runManager = new G4RunManager;

	OMSimDetectorConstruction *detector;
	detector = new OMSimDetectorConstruction();
	runManager->SetUserInitialization(detector);

	G4VUserPhysicsList *physics = new OMSimPhysicsList;
	runManager->SetUserInitialization(physics);

	auto visManager = new G4VisExecutive;
	visManager->Initialize();

	G4VUserPrimaryGeneratorAction *gen_action = new OMSimPrimaryGeneratorAction();
	runManager->SetUserAction(gen_action);

	G4UserRunAction *run_action = new OMSimRunAction();
	runManager->SetUserAction(run_action);

	G4UserEventAction *event_action = new OMSimEventAction();
	runManager->SetUserAction(event_action);

	G4UserTrackingAction *tracking_action = new OMSimTrackingAction();
	runManager->SetUserAction(tracking_action);

	G4UserSteppingAction *stepping_action = new OMSimSteppingAction();
	runManager->SetUserAction(stepping_action);

	runManager->Initialize();

	OMSimPMTResponse &lPhotocathodeResponse = OMSimPMTResponse::getInstance();
	OMSimAnalysisManager &lAnalysisManager = OMSimAnalysisManager::getInstance();
	OMSimUIinterface &lUIinterface = OMSimUIinterface::getInstance();
	lUIinterface.setUI(G4UImanager::GetUIpointer());

	G4Navigator *navigator = new G4Navigator();
	navigator->SetWorldVolume(detector->mWorldPhysical);
	navigator->LocateGlobalPointAndSetup(G4ThreeVector(0., 0., 0.));

	G4TouchableHistory *history = navigator->CreateTouchableHistory();

	OMSimCommandArgsTable &lArgs = OMSimCommandArgsTable::getInstance();
	lUIinterface.applyCommand("/control/execute ", lArgs.get<bool>("visual"));

	double startingtime = clock() / CLOCKS_PER_SEC;

	G4double lDistance = lArgs.get<G4double>("dist");
	G4double cos_phi, sin_phi, cos_theta, sin_theta, rho, posX, posY, posZ;


	std::vector<G4PV2DDataVector> data = detector->mData->loadtxt("theta_phi.txt", false);
	std::vector<G4double> lThetas = data.at(0);
	std::vector<G4double> lPhis = data.at(1);

	std::vector<int> lWavelengths;
	lWavelengths.push_back(400);
	//arange(250, 710, 10);
	for (const auto& lWavelength : lWavelengths) {
	for(std::vector<int>::size_type i = 0; i != lThetas.size(); i++) { 
	G4double lPhi = lThetas.at(i);//lArgs.get<G4double>("phi");
	G4double lTheta = lPhis.at(i);//lArgs.get<G4double>("theta");
	cos_phi = cos(lPhi*deg);
	sin_phi = sin(lPhi*deg);
	cos_theta = cos(lTheta*deg);
	sin_theta = sin(lTheta*deg);

	rho = lDistance*sin_theta;
	posX = rho*cos_phi;
	posY = rho*sin_phi;
	posZ = lDistance*cos_theta;

	lUIinterface.applyCommand("/control/execute OMSim.gps");

	//lUIinterface.applyCommand("/gps/particle opticalphoton");
	lUIinterface.applyCommand("/gps/energy", 1239.84193 / lWavelength,"eV");

	//lUIinterface.applyCommand("/gps/pos/type Plane");
	//lUIinterface.applyCommand("/gps/pos/shape Circle");
	lUIinterface.applyCommand("/gps/pos/centre", posX, posY, posZ, "mm");
	lUIinterface.applyCommand("/gps/pos/radius", lArgs.get<G4double>("diam")*0.5, "mm");
	double x,y,z; // vector entries for plane positioning
	x = -sin_phi;	// d/dphi of positionVector (original divided by sin_theta, because length one not needed)
	y = cos_phi; 
	z = 0;
	lUIinterface.applyCommand("/gps/pos/rot1", x,y,z);
	lUIinterface.applyCommand("/gps/ang/rot1", x,y,z);
	x = -cos_phi * cos_theta;	// -d/dtheta of positionVector (divided by sin_theta, because length one not needed)
	y = -sin_phi * cos_theta;
	z = sin_theta;
	lUIinterface.applyCommand("/gps/pos/rot2", x,y,z);
	lUIinterface.applyCommand("/gps/ang/rot2", x,y,z);
	lUIinterface.runBeamOn();

	std::string lFileName = lArgs.get<std::string>("output_file") + ".txt";
	OMSimAnalysisManager& lAnalysisManager = OMSimAnalysisManager::getInstance();
	lAnalysisManager.datafile.open(lFileName.c_str(), std::ios::out|std::ios::app);
	lAnalysisManager.datafile << std::fixed << std::setprecision(3) << lPhi << "\t" << lTheta << "\t" <<  lArgs.get<G4double>("diam") << "\t" << lWavelength << "\t" << OMSimCommandArgsTable::getInstance().get<G4int>("numevents") <<"\t" ;	

	lAnalysisManager.WriteAccept();
	// 	Close output data file
	lAnalysisManager.datafile.close();
	lAnalysisManager.Reset();

	}
	}
	//detector->mMDOM->setNavigator(navigator);
	//detector->mMDOM->runBeamOnFlasher(1, 9);

	//lUIinterface.runBeamOn();
	// opening user interface prompt and visualization after simulation was run
	if (lArgs.get<bool>("visual"))
	{
		char *argumv[] = {"all", NULL};
		G4UIExecutive *UIEx = new G4UIExecutive(1, argumv);
		lUIinterface.applyCommand("/control/execute ../aux/init_vis.mac");
		UIEx->SessionStart();
		delete UIEx;
	}
	double finishtime = clock() / CLOCKS_PER_SEC;
	G4cout << "Computation time: " << finishtime - startingtime << " seconds." << G4endl;

	delete navigator;
	delete history;
	delete visManager;
	delete runManager;
	return 0;
}

int main(int argc, char *argv[])
{
	try
	{
		// Do not use G4String as type here...
		po::options_description desc("User arguments available");
		desc.add_options()
		("help,h", "produce this help message")("world_radius,w", po::value<G4double>()->default_value(3.0), "radius of world sphere in m")
		("diam,d", po::value<G4double>()->default_value(560.0), "beam diameter in mm")
		("dist,r", po::value<G4double>()->default_value(2000), "emitter distance from origin, in mm")
		("theta,t", po::value<G4double>()->default_value(0.0), "theta (= zenith) in deg")
		("phi,f", po::value<G4double>()->default_value(0.0), "phi (= azimuth) in deg")
		("wavelength,l", po::value<G4double>()->default_value(400.0), "wavelength of incoming light in nm")
		("numevents,n", po::value<G4int>()->default_value(0), "number of photons emitted per angle")
		("gunfile,g", po::value<std::string>()->default_value("OMSim.gps"), "file containing GPS parameters")
		("glass", po::value<G4int>()->default_value(0), "index to select glass type [VITROVEX = 0, Chiba = 1, Kopp = 2, myVitroVex = 3, myChiba = 4, WOMQuartz = 5, fusedSilica = 6]")
		("gel", po::value<G4int>()->default_value(1), "index to select gel type [Wacker = 0, Chiba = 1, IceCube = 2, Wacker_company = 3]")
		("detector_type", po::value<G4int>()->default_value(2), "module type [custom = 0, Single PMT = 1, mDOM = 2, pDDOM = 3, LOM16 = 4]")
		("cad", po::bool_switch(), "activates CAD import for supported modules (LOM16)")("harness", po::bool_switch(), "activates harness [OFF, ON]")
		("environment", po::value<G4int>()->default_value(0), "medium in which the setup is emmersed [AIR = 0, ice = 1, spice = 2]")
		("output_file", po::value<std::string>()->default_value("mdom_testoutput.txt"), "filename for output")
		("reflective_surface", po::value<G4int>()->default_value(0), "index to select reflective surface type [Refl_V95Gel = 0, Refl_V98Gel = 1, Refl_Aluminium = 2, Refl_Total98 = 3]")
		("visual,v", po::bool_switch(), "shows visualization of module after run")
		("pmt_model", po::value<G4int>()->default_value(0), "R15458 (mDOM) = 0,  R7081 (DOM) = 1, 4inch (LOM) = 2, R5912_20_100 (D-Egg)= 3")
		("abs_length", po::value<G4double>()->default_value(1), "Abs length of photocathode")
		("string_pos_angle", po::value<G4double>()->default_value(45), "Polar angle of main data cable (viewed from above)");
		po::variables_map lVariablesMap;
		po::store(po::parse_command_line(argc, argv, desc), lVariablesMap);
		po::notify(lVariablesMap);

		if (lVariablesMap.count("help"))
		{
			std::cout << desc << "\n";
			return 1;
		}

		OMSimCommandArgsTable &lArgs = OMSimCommandArgsTable::getInstance();

		// Now store the parsed parameters in the OMSimCommandArgsTable instance
		for (const auto &option : lVariablesMap)
		{
			lArgs.setParameter(option.first, option.second.value());
		}

		// Now that all parameters are set, "finalize" the OMSimCommandArgsTable instance so that the parameters cannot be modified anymore
		lArgs.finalize();

		OMSim();
	}
	catch (std::exception &e)
	{
		std::cerr << "error: " << e.what() << "\n";
		return 1;
	}
	catch (...)
	{
		std::cerr << "Exception of unknown type!\n";
	}

	return 0;
}
