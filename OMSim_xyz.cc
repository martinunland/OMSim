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

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4ThreeVector.hh"

#include "G4UItcsh.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"	//xxx

#include "argtable2.h"
#include <ctime>
#include <sys/time.h>

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>

#include <cmath>	// for abs() of doubles
//since Geant4.10: include units manually
#include "G4SystemOfUnits.hh"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include </data/root/include/TGraph.h>
#include <TGraph.h>

unsigned int	stats_buffer_max_size = 10;	// how many hits to keep in memory before purging to file in EndOfEventAction
unsigned int	stats_buffer_last_purge_at = 0;	// at what hits count was the hits file last written to
// std::vector<G4int>	stats_PMT_hit;
// std::vector<G4int>	stats_OM_hit;

// double hc_eVnm = 1239.84193; // h*c in eV * nm
double pi = M_PI;
G4double	gworldsize;
G4int		gsimevents;
G4String	ggunfilename;
G4double	gDistance;
G4double	gBeamDiam;
G4double	gposX, gposY, gposZ;
G4double	gtheta;
G4double	gphi;
G4double    gwavelen;
G4String	ghitsfilename;
G4String	gHittype;
G4int		gPMT;
G4int		gEnvironment;
G4bool		gVisual;
G4bool		gCADImport;
G4bool 		gPlaceHarness;
G4bool		gInteractive;
G4bool		gHeader;

G4int		gGlass;
G4int		gGel;
G4double	gRefCone_angle;
G4int		gConeMat;
G4int		gHolderColor;
G4int		gDOM;
G4int       gHarness;
G4int       gRopeNumber;


// G4String	greffilename;

G4bool		gKillAll;
G4long		current_event_id;

struct timeval	gTime_Run_Start;
struct timeval	gTime_Run_End;
long randseed;

OMSimAnalysisManager gAnalysisManager;



OMSimPMTResponse* gPhotocathodeResponse = new OMSimPMTResponse();


// void clearstats() {
// 	stats_PMT_hit.clear();
// 	stats_OM_hit.clear();
// }

std::vector<G4String> explode(G4String s, char d) {
	std::vector<G4String> o;
	int i,j;
	i = s.find_first_of("#");
	if (i == 0) return o;
 	while (s.size() > 0) {
		i = s.find_first_of(d);
		j = s.find_last_of(d);
		o.push_back(s.substr(0, i));
		if (i == j) {
			o.push_back(s.substr(j+1));
			break;
		}
		s.erase(0,i+1);
 	}
	return o;// o beinhaltet s ohne d
}
std::vector<G4String> explode(char* cs, char d) {
	std::vector<G4String> o;
	G4String s = cs;
	return explode(s,d);
}
std::vector<double> readColumnDouble (G4String fn, int col) {
	std::vector<double>	values;
	unsigned int c;
	double	a;
	c = col;
	std::ifstream	infile;
	std::vector<G4String> n;
	char l[256];
	G4String l2;
	infile.open(fn);
	while (infile.good() && !infile.eof()) {
		infile.getline(l,255);
		l2 = l;
		n = explode(l2,'\t');
		if (n.size()>=c) {
			a = atof(n.at(c-1));
			values.push_back(a);
		}
	}
	infile.close();

	return values;//values enthält den c. Wert aus fn (aus jeder Spalte,welche  nach 255 zeichen oder durch \n beendet wird?)
}

int OMSim() {
	struct timeval time_for_randy;
	gettimeofday(&time_for_randy, NULL);

	randseed = time_for_randy.tv_sec+4294*time_for_randy.tv_usec;
	CLHEP::HepRandom::setTheEngine(new CLHEP::RanluxEngine(randseed,3));
	
	std::stringstream command;

	G4RunManager* runManager = new G4RunManager;

	OMSimDetectorConstruction* detector;
	detector = new OMSimDetectorConstruction();
	runManager->SetUserInitialization(detector);


	G4VUserPhysicsList* physics = new OMSimPhysicsList;
	runManager->SetUserInitialization(physics);


 		auto visManager = new G4VisExecutive;
 		visManager->Initialize();
 		

	G4VUserPrimaryGeneratorAction* gen_action = new OMSimPrimaryGeneratorAction();
	runManager->SetUserAction(gen_action);

	G4UserRunAction* run_action = new OMSimRunAction();
	runManager->SetUserAction(run_action);

	G4UserEventAction* event_action = new OMSimEventAction();
	runManager->SetUserAction(event_action);

 	G4UserTrackingAction* tracking_action = new OMSimTrackingAction();
 	runManager->SetUserAction(tracking_action);

	G4UserSteppingAction* stepping_action = new OMSimSteppingAction();
	runManager->SetUserAction(stepping_action);

	runManager->Initialize();

	auto UI = G4UImanager::GetUIpointer();
	double startingtime = clock() / CLOCKS_PER_SEC;
// setting up light source:
	command.str("");
	command << "/control/execute " << ggunfilename;
	UI->ApplyCommand(command.str());
	
	double BeamRad = 0.5*gBeamDiam;
	double cos_phi, sin_phi, cos_theta, sin_theta;
// 	double rho, posX,posY,posZ;
	
	cos_phi = cos(gphi*deg);
	sin_phi = sin(gphi*deg);
	cos_theta = cos(gtheta*deg);
	sin_theta = sin(gtheta*deg);



    TGraph* ZCorr = new TGraph("Zcorrection.txt");
    ZCorr->SetName("Zcorrection");
// 	rho = gDistance*sin_theta;
// 	posX = rho*cos_phi;
// 	posY = rho*sin_phi;
// 	posZ = gDistance*cos_theta;
	for (int xx=-41; xx<42;xx++){

        for (int yy=-41; yy<42;yy++){
            if( (std::sqrt(xx*xx+yy*yy) < 42) ){

    command.str("");
    command << "ScanDataTest/" << xx << "_" << yy << ".txt";
    ghitsfilename = command.str();

	command.str("");
	command << "/control/execute " << ggunfilename;
	UI->ApplyCommand(command.str());

	
	command.str("");
	command << "/gps/energy " << (1239.84193 / gwavelen) << " eV ";
	UI->ApplyCommand(command.str());


    gposX = xx;
    gposY = yy;
    gposZ = ZCorr->Eval(std::sqrt(xx*xx+yy*yy));
	if (gposZ<4.8) { gposZ = 4.8;}
    G4cout << gposX << " " << gposZ << " " << std::sqrt(xx*xx+yy*yy)<<" "<< gposZ << G4endl;
	command.str("");
	command << "/gps/pos/centre "<< gposX <<" "<< gposY <<" "<< gposZ <<" mm";
	UI->ApplyCommand(command.str());


	/*
	command.str("");
	command << "/gps/pos/radius " << BeamRad << " mm";
	UI->ApplyCommand(command.str());

	double x,y,z; // vector entries for plane positioning

	x = -sin_phi;	// d/dphi of positionVector (original divided by sin_theta, because length one not needed)
	y = cos_phi; 
	z = 0;
	command.str("");
	command << "/gps/pos/rot1 " << x <<" "<< y <<" "<< z;
	UI->ApplyCommand(command.str());
	command.str("");
	command << "/gps/ang/rot1 " << x <<" "<< y <<" "<< z;
	UI->ApplyCommand(command.str());

	x = -cos_phi * cos_theta;	// -d/dtheta of positionVector (divided by sin_theta, because length one not needed)
	y = -sin_phi * cos_theta;
	z = sin_theta;
	command.str("");
	command << "/gps/pos/rot2 " << x <<" "<< y <<" "<< z;
	UI->ApplyCommand(command.str());
	command.str("");
	command << "/gps/ang/rot2 " << x <<" "<< y <<" "<< z;
	UI->ApplyCommand(command.str());

	*/


// // preparing file for results.
std::ofstream datafile;
 	datafile.open(ghitsfilename.c_str(), std::ios::out|std::ios::app);
// 	if (gHeader) {
// 		datafile << "# theta\tphi\tposX\tposY\tposZ\twavelength\temitted photons\tbeam radius [mm]\thits per PMT (#0 to #23)\ttotal hits\n";
// 	}
	datafile << std::fixed << std::setprecision(3) << gposX << "\t" << gposY << "\t" << gposZ << "\t" << gwavelen << "\t" ;	
 	datafile.close();

// preparing file for results.
	
	
	
// starting run (analysis manager takes care of results)
	if (gsimevents > 0) {
		command.str("");
		command << "/run/beamOn " << gsimevents;
		UI->ApplyCommand(command.str());
	}
        }

        }
    }
// opening user interface prompt and visualization after simulation was run
	if (gInteractive){
		int argumc = 1;
		char* argumv[] = {"all", NULL};
		G4UIExecutive* UIEx = new G4UIExecutive(argumc, argumv);
		if (gVisual){
			UI->ApplyCommand("/control/execute ../aux/init_vis.mac");
			command.str("");
			command << "/run/beamOn " << gsimevents;
			UI->ApplyCommand(command.str());
		}
		UIEx->SessionStart();
		delete UIEx;
	}
double finishtime=clock() / CLOCKS_PER_SEC;
G4cout << "Computation time: " << finishtime-startingtime << " seconds." << G4endl;

	delete visManager;


	delete runManager;
	return 0;
}


int main(int argc,char *argv[])
{
	struct arg_dbl  *worldsize	= arg_dbl0("wW", "world","<n>","\t\tradius of world sphere in m");
	struct arg_dbl  *diameter	= arg_dbl0("dD", "diam","<n>","\t\tbeam diameter in mm");
    struct arg_dbl  *xpos		= arg_dbl0("xX", "posx, posX","<n>","\t\tx in mm");
	struct arg_dbl  *ypos		= arg_dbl0("yY", "posy, posY","<n>","\t\ty in mm");
	struct arg_dbl  *zpos		= arg_dbl0("zZ", "posz, posZ","<n>","\t\tz in mm");
	struct arg_dbl  *distance	= arg_dbl0("rR", "dist, rad","<n>","\t\temitter distance from origin, in mm");
	struct arg_dbl  *theta		= arg_dbl0("tT", "theta","<n>","\t\ttheta (= zenith) in deg");
 	struct arg_dbl  *phi		= arg_dbl0("fF", "phi","<n>","\t\tphi (= azimuth) in deg");
    struct arg_dbl  *wavelen	= arg_dbl0("lL", "lambda","<n>","\t\twavelength of incoming light in nm");
	struct arg_int  *events		= arg_int0("nN", "numevents,nevents","<n>","\tnumber of photons emitted per angle");
	struct arg_file *gunfile	= arg_file0("gG","gun","<file.txt>","\t\tfile containing GPS parameters");
	struct arg_int  *pmt		= arg_int0("pP", "pmt,PMT","<n>","\t\tPMT type [12199S, etel, 12199e]");

	struct arg_int  *glass		= arg_int0("uU", "glass","<n>","\t\t\tglass type [VITROVEX, Chiba, Kopp, myVitroVex, myChiba, WOMQuartz, fusedSilica]");
	struct arg_int	*gel 		= arg_int0("jJ", "gel", "<n>", "\t\t\tgel type [Wacker, Chiba, IceCube, Wacker_company]");
	struct arg_dbl  *cone_ang   = arg_dbl0("aA", "cone_ang","<n>","\t\t\topening semi-angle of cone; (51 deg)");	
	struct arg_int	*conemat 	= arg_int0("kK", "conemat", "<n>", "\t\t\tcone material [V95, v98, aluminium, total98]");
	struct arg_int	*dom 		= arg_int0("mM", "om, dom", "<n>", "\t\t\tmodule type [Single PMT, MDOM, PDom, Lom16, custom]");
	struct arg_lit	*CADimport 	= arg_lit0("cC", "cad", "\t\t\tactivates CAD import for supported modules (LOM16)");
	struct arg_lit	*PlaceHarness 	= arg_lit0("sS", "harness", "\t\t\tactivates harness [OFF, ON]");
	
	
	struct arg_int  *environment= arg_int0("eE", "environment","<n>","\t\tmedium in which the setup is emmersed [AIR, ice, spice]");
	struct arg_file *outputfile	= arg_file0("oO","output","<file.txt>","\t\tfilename for hits data");
	struct arg_int  *hittype	= arg_int0("hH", "hits","<n>","\t\thit collection [individual, COLLECTIVE]");
	struct arg_lit	*interactive= arg_lit0("iI","interact","\t\topens user interface after run");
	struct arg_lit	*visual		= arg_lit0("vV","visual","\t\tshows visualization of module after run (also calls interactive)");
	struct arg_lit	*nohead		= arg_lit0("qQ","nh, nohead","\t\tno header in outputfile");
	struct arg_lit	*help		= arg_lit0(NULL,"help","\t\tprint this help and exit");
	struct arg_end  *end		= arg_end(20);
	
	void* argtable[] = {worldsize,
						diameter,xpos, ypos, zpos,
						distance,
						theta,
						phi,
						wavelen,
						events,
						gunfile,
						pmt,
						
						glass,
						gel,
						cone_ang,
						conemat,
						dom,
						CADimport,
						PlaceHarness,
						
						environment,
						outputfile,
						hittype,
						interactive,
						visual,
						nohead,
						help, end};
						
	const char* progname = "OMSim";
	int nerrors;
	int exitcode=0;

	// verify the argtable[] entries were allocated sucessfully
	if (arg_nullcheck(argtable) != 0) {
		/* NULL entries were detected, some allocations must have failed */
		printf("%s: insufficient memory\n",progname);
		exitcode=1;
		arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
		return exitcode;
	}

	// set any command line default values prior to parsing
	worldsize->dval[0] = 3.0;	// world diameter in meters
	//diameter->dval[0] = 420.0;	// 400 mm # for 14" sphere, 480 mm # for 17" sphere 
	diameter->dval[0] = 560.0;	// LOM18 has 270*2 = 540 mm pole to pole length + 10mm safety margin
	xpos->dval[0] = 0.0;
	ypos->dval[0] = 0.0;
	zpos->dval[0] = 300.;
	//distance->dval[0] = 0.5 * 356 + 27.5 + 2; // here value for mDOM scan with 2*mm safety margin
	//distance->dval[0] = 270 + 20; // value for LOM18 glass with 20 mm saftety margin
	distance->dval[0] = 2000; // close to world border for harness simulations...no scattering in ice so it doesn't matter how far away it it
	


	theta->dval[0] = 0.0;
	phi->dval[0] = 0.0;
    wavelen->dval[0] = 400.0;	// [nm]
	events->ival[0] = 0;
	gunfile->filename[0] = "OMSim_xyz.gps";
	pmt->ival[0] = 0;			// use new R12199 version as default
	
	glass->ival[0] = 0;	// use VITROVEX as default
	gel->ival[0] = 1;	// use Chiba
	cone_ang->dval[0] = 51.0; // [degrees]	
	conemat->ival[0] = 0;	// use Alemco V95 as default
	dom->ival[0] = 1;	// use mDOM as default
	
	environment->ival[0] = 0;	// use air as default
	outputfile->filename[0] = "mdom_testoutput.txt";
	hittype->ival[0] = 1;		// store information on collective hits as default

	/* Parse the command line as defined by argtable[] */
    nerrors = arg_parse(argc,argv,argtable);

    /* special case: '--help' takes precedence over error reporting */
    if (help->count > 0)
	{
        printf("\nGEANT4 simulation of the mDOM: angular acceptance scan\n");
        printf("\nUsage: %s", progname);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
        printf("\n");
        exitcode=0;
        goto hell;
	}

    /* If the parser returned any errors then display them and exit */
    if (nerrors > 0)
	{
        /* Display the error details contained in the arg_end struct.*/
        arg_print_errors(stdout,end,progname);
        printf("Try '%s --help' for more information.\n",progname);
        exitcode=1;
        goto hell;
	}

    /* special case: uname with no command line options induces brief help */
    if (argc==1)
	{
        printf("Try '%s --help' for more information.\n",progname);
        exitcode=0;
        goto hell;
	}

//	assign command-line arguments to variables:
	gworldsize = worldsize->dval[0];
	gBeamDiam = diameter->dval[0];
    gposX = xpos->dval[0];
	gposY = ypos->dval[0];
	gposZ = zpos->dval[0];
	gDistance = distance->dval[0];
	gtheta 	= theta->dval[0];
	gphi 	= phi->dval[0];
	gwavelen = wavelen->dval[0];
	gsimevents = events->ival[0];
	ggunfilename = gunfile->filename[0];
	gPMT = pmt->ival[0];
	
	gGlass = glass->ival[0];
	gGel = gel->ival[0];
	gRefCone_angle = cone_ang->dval[0];	
	gConeMat = conemat->ival[0];
	gDOM = dom->ival[0];

	
	gEnvironment = environment->ival[0];
	ghitsfilename = outputfile->filename[0];
	if (hittype->ival[0]==0){
		gHittype = "individual";
	}
	if (hittype->ival[0]==1){
		gHittype = "collective";
	}
	if (interactive->count > 0) gInteractive = true; else gInteractive = false;
	if (CADimport->count > 0) gCADImport = true; else gCADImport = false;
	if (PlaceHarness->count > 0) gPlaceHarness = true; else gPlaceHarness = false;
	if (visual->count > 0) {
		gVisual = true;
		gInteractive = true;
	}
	else {
		gVisual = false;
	}
	if (nohead->count > 0) gHeader = false; else gHeader = true;
	
	//	check params for sanity
	OMSim();

hell:
    /* deallocate each non-null entry in argtable[] */
	arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
	//	return exitcode;

}

