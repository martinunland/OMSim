// #ifndef DEGG_Test_h
// #define DEGG_Test_h 1
#include "abcDetectorComponent.hh"
#include "OMSimPMTConstruction.hh"

#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "OMSimInputData.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
class dEGGHarness;

class DEgg : public abcDetectorComponent
{
public:
    
    DEgg(OMSimInputData* pData,G4bool pPlaceHarness=true);
    void Construction();

private:
    OMSimPMTConstruction *mPMTManager;
    dEGGHarness* mHarness;

    G4String mDataKey = "om_DEGG";
    
    void AppendPMTs();
    void AppendInternalComponentsFromCAD();
    void AppendPressureVesselFromCAD();
    
    G4VSolid* CreateEggSolid(G4int pSegments_1,
                            G4double pSphereRmax,
                            G4double pSpheredTheta,
                            G4double pSphereTransformZ,
                            G4double pTorus1R,
                            G4double pCenterOfTorus1R,
                            G4int pSegments_2,
                            G4double pTorus2R,
                            G4double pCenterOfTorus2R,
                            G4double pCenterOfTorus2_z,
                            G4double pTorus2_Zmin,
                            G4double pTorus2_Zmax,
                            G4double pTorus2_Z0,
                            G4double pTorus1TransformZ);

};
