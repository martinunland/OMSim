#ifndef OMSimLOM16_h
#define OMSimLOM16_h 1
#include "abcDetectorComponent.hh"
#include "OMSimPMTConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "OMSimInputData.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"

extern G4double gRefCone_angle;


class LOM16 : public abcDetectorComponent
{
public:
    LOM16(OMSimInputData* pData, G4bool pPlaceHarness = false);
    void Construction();
    G4String mDataKey = "om_LOM16";

  

private:
    OMSimPMTConstruction* mPMTManager;

    //functions
    void GetSharedData();
    G4UnionSolid* PressureVessel(const G4double pOutRad, G4String pSuffix);

    //Lom specific functions
    void PlaceCADSupportStructure();

    void AppendEquatorBand();
    
    //for gelpad and PMT creation
    void PlacePMTsAndGelpads(G4VSolid* lGelSolid,G4LogicalVolume* lGelLogical);
    void SetPMTAndGelpadPositions();
    void CreateGelpadLogicalVolumes(G4VSolid* lGelSolid);
    void PlacePMTs(G4LogicalVolume* lInnerVolumeLogical);
    void PlaceGelpads(G4LogicalVolume* lInnerVolumeLogical);

    //selection variables
    G4bool mPlaceHarness = true;
    G4bool mHarnessUnion = true; //it should be true for the first module that you build, and then false

    //vectors for positions and rotations
    std::vector<G4ThreeVector> mPMTPositions;
    std::vector<G4ThreeVector> mGelpadPositions;
    std::vector<G4double> mPMT_theta;
    std::vector<G4double> mPMT_phi;


    //Shared data from jSON file
    G4double mGlassOutRad;
    G4int mNrPolarPMTs;
    G4int mNrEqPMTs;
    G4int mTotalNrPMTs;

    //helper variables
    std::stringstream mConverter;
    std::stringstream mConverter2;

    //for gelpads
    G4double mGelPadDZ;

    //logical of gelpads
    std::vector<G4LogicalVolume*> mGelPad_logical;


};

#endif
//
