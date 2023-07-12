#ifndef OMSimPMTResponse_h
#define OMSimPMTResponse_h 1

#include "G4SystemOfUnits.hh"
#include "G4Types.hh"
#include <functional>
#include <map>
#include <TGraph.h>
#include <TGraph2D.h>

#include "interpolation.h" //alglib

class OMSimPMTResponse
{
public:

    static OMSimPMTResponse& getInstance() { // Meyers singleton
        static OMSimPMTResponse instance;
        return instance;
    }


    struct PMTPulse{
    G4double PE;
    G4double TransitTime;
    G4double DetectionProbability;
  };

    PMTPulse ProcessPhotocathodeHit(G4double pX, G4double pY, G4double pWavelength);


private:
    std::vector<G4double> mScannedWavelengths{460*nm, 480*nm, 500*nm, 520*nm, 540*nm, 560*nm,  580*nm, 600*nm, 620*nm, 640*nm};

    double mX;
    double mY;
    G4double mWavelength;


    TGraph2D* scan;
    TGraph* mRelativeDetectionEfficiencyInterp;
    // std::map<G4double, TGraph2D*> mGainG2Dmap;
    // std::map<G4double, TGraph2D*> mGainResolutionG2Dmap;
    // std::map<G4double, TGraph2D*> mTransitTimeG2Dmap;
    // std::map<G4double, TGraph2D*> mTransitTimeSpreadG2Dmap;

    void fetchDataForInterpolation(const std::string& pFilePath, alglib::real_1d_array& pX, alglib::real_1d_array& pY, alglib::real_1d_array& pF); 
    void buildSplineFromDataFile(const std::string& pFileName, std::map<double, alglib::spline2dinterpolant>& pMap, const double& pKey);
    std::map<G4double, alglib::spline2dinterpolant> mGainG2Dmap;
    std::map<G4double, alglib::spline2dinterpolant> mGainResolutionG2Dmap;
    std::map<G4double, alglib::spline2dinterpolant> mTransitTimeG2Dmap;
    std::map<G4double, alglib::spline2dinterpolant> mTransitTimeSpreadG2Dmap;

    G4double GetCharge(G4double pWavelengthKey);
    G4double GetCharge(G4double pWavelengthKey1, G4double pWavelengthKey2);

    G4double GetTransitTime(G4double pWavelengthKey);
    G4double GetTransitTime(G4double pWavelengthKey1, G4double pWavelengthKey2);

    PMTPulse GetPulseFromInterpolation(G4double pWavelengthKey1, G4double pWavelengthKey2);
    PMTPulse GetPulseFromKey(G4double pWavelengthKey);

    G4double WavelengthInterpolatedValue(std::map<G4double, alglib::spline2dinterpolant> pMap, G4double pWavelengthKey1, G4double pWavelengthKey2);

    OMSimPMTResponse();
    ~OMSimPMTResponse() = default;
    OMSimPMTResponse(const OMSimPMTResponse&) = delete;
    OMSimPMTResponse& operator=(const OMSimPMTResponse&) = delete;

};

#endif
//
