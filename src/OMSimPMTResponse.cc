/** @file OMSimPMTResponse.cc
 *  @brief PMT response simulator
 *
 *  @author Martin Unland
 *  @date April 2022
 *
 *  @version Geant4 10.7
 */

#include "OMSimPMTResponse.hh"
#include "OMSimLogger.hh"
#include "G4SystemOfUnits.hh"
#include "G4Types.hh"
#include "Randomize.hh"
#include <functional>
#include <TGraph2D.h>
#include <TGraph.h>

#include "interpolation.h" //alglib
#include "OMSimInputData.hh"

void OMSimPMTResponse::fetchDataForInterpolation(const std::string& pFilePath, alglib::real_1d_array& pX, alglib::real_1d_array& pY, alglib::real_1d_array& pF) {
    // Load the text file into a vector of vectors
    std::vector<std::vector<double>> lData = InputDataManager::loadtxt(pFilePath, true, 0, '\t');

    // Ensure the data has at least three columns (X, Y, F)
    if (lData.size() < 3) throw std::runtime_error("Data has less than three columns in file: " + pFilePath);

    // Convert std::vector to alglib::real_1d_array
    pX.setcontent(lData[0].size(), lData[0].data());
    pY.setcontent(lData[1].size(), lData[1].data());
    pF.setcontent(lData[2].size(), lData[2].data());
}

void OMSimPMTResponse::buildSplineFromDataFile(const std::string& pFileName, std::map<double, alglib::spline2dinterpolant>& pMap, const double& pKey)
{
    alglib::real_1d_array lX, lY, lF;
    fetchDataForInterpolation(pFileName, lX, lY, lF);
    alglib::spline2dinterpolant c;
    log_debug("here");
    alglib::spline2dbuildbicubicv(lX, lX.length(), lY, lY.length(), lF, 1, c);
    log_debug("here2");
    pMap[pKey] = c;
}

OMSimPMTResponse::OMSimPMTResponse() {
    log_debug("Opening photocathode scans data...");
    std::string path = "../data/PMT_scans/";

    mRelativeDetectionEfficiencyInterp = new TGraph((path+"weightsVsR_vFit_220nm.txt").c_str());
    mRelativeDetectionEfficiencyInterp->SetName("RelativeDetectionEfficiencyWeight");

    // for(const auto& lKey : mScannedWavelengths) {
    //     std::string lWv = std::to_string((int)(lKey/nm));
    //     mGainG2Dmap[lKey] = new TGraph2D((path+"Gain_PE_"+lWv+".dat").c_str());
    //     mGainG2Dmap[lKey]->SetName(("Gain_PE_"+lWv).c_str());
    //     mGainResolutionG2Dmap[lKey] = new TGraph2D((path+"SPEresolution_"+lWv+".dat").c_str());
    //     mGainResolutionG2Dmap[lKey]->SetName(("SPEresolution_"+lWv).c_str());
    //     mTransitTimeG2Dmap[lKey] = new TGraph2D((path+"TransitTime_"+lWv+".dat").c_str());
    //     mTransitTimeG2Dmap[lKey]->SetName(("TransitTime_"+lWv).c_str());
    //     mTransitTimeSpreadG2Dmap[lKey] = new TGraph2D((path+"TransitTimeSpread_"+lWv+".dat").c_str());
    //     mTransitTimeSpreadG2Dmap[lKey]->SetName(("TransitTimeSpread_"+lWv).c_str());
    // }

    for(const auto& lKey : mScannedWavelengths) {
        std::string lWv = std::to_string((int)(lKey/nm));
        buildSplineFromDataFile(path + "Gain_PE_" + lWv + ".dat", mGainG2Dmap, lKey);
        buildSplineFromDataFile(path + "SPEresolution_" + lWv + ".dat", mGainResolutionG2Dmap, lKey);
        buildSplineFromDataFile(path + "TransitTime_" + lWv + ".dat", mTransitTimeG2Dmap, lKey);
        buildSplineFromDataFile(path + "TransitTimeSpread_" + lWv + ".dat", mTransitTimeSpreadG2Dmap, lKey);
    }

    log_debug("Finished opening photocathode scans data...");
}


/**
 * Get mean charge and SPE resolution from measurements for a given wavelength and sample from gaussian with mean PE and std SPE_res.
 * There are only a limited number of scans, which one can select from the scan map. We use TGraph2D from ROOT to interpolate selected scan.
 * @param pWavelengthKey  key for selecting a scan 
 * @return G4double charge in PE
 */
G4double OMSimPMTResponse::GetCharge(G4double pWavelengthKey) {

    // G4double lMeanPE = mGainG2Dmap[pWavelengthKey]->Interpolate(mX, mY);
    // G4double lSPEResolution = mGainResolutionG2Dmap[pWavelengthKey]->Interpolate(mX, mY);

    G4double lMeanPE = alglib::spline2dcalc(mGainG2Dmap[pWavelengthKey], mX, mY);
    G4double lSPEResolution = alglib::spline2dcalc(mGainResolutionG2Dmap[pWavelengthKey], mX, mY);

    double lToReturn = -1;
    double lCounter = 0;
    while (lToReturn<0){
        lToReturn = G4RandGauss::shoot(lMeanPE, lSPEResolution);
        lCounter ++;
        if (lCounter>10) return 0;
    }

    return lToReturn; 
}


G4double OMSimPMTResponse::GetCharge(G4double pWavelength1, G4double pWavelength2) {

    G4double lMeanPE =  WavelengthInterpolatedValue(mGainG2Dmap, pWavelength1, pWavelength2);
    G4double lSPEResolution = WavelengthInterpolatedValue(mGainResolutionG2Dmap, pWavelength1, pWavelength2);

    double lToReturn = -1;
    double lCounter = 0;
    while (lToReturn<0){
        lToReturn = G4RandGauss::shoot(lMeanPE, lSPEResolution);
        lCounter ++;
        if (lCounter>10) return 0;
    }

    return lToReturn; 
}

/**
 * Get transit time and transit time spread from measurements for a given wavelength and sample from gaussian with mean TT and std TTS.
 * There are only a limited number of scans, which one can select from the scan map. We use TGraph2D from ROOT to interpolate selected scan.
 * @param pWavelengthKey  key for selecting a scan 
 * @return G4double transit time
 */
G4double OMSimPMTResponse::GetTransitTime(G4double pWavelengthKey) {

    // G4double lMeanTransitTime = mTransitTimeG2Dmap[pWavelengthKey]->Interpolate(mX, mY) * ns;
    // G4double lTTS = mTransitTimeSpreadG2Dmap[pWavelengthKey]->Interpolate(mX, mY) * ns;

    G4double lMeanTransitTime = alglib::spline2dcalc(mTransitTimeG2Dmap[pWavelengthKey], mX, mY)* ns;
    G4double lTTS = alglib::spline2dcalc(mTransitTimeSpreadG2Dmap[pWavelengthKey], mX, mY)* ns;

    return G4RandGauss::shoot(lMeanTransitTime, lTTS); 
}


G4double OMSimPMTResponse::GetTransitTime(G4double pWavelength1, G4double pWavelength2) {

    G4double lMeanTransitTime = WavelengthInterpolatedValue(mTransitTimeG2Dmap, pWavelength1, pWavelength2 )*ns; 
    G4double lTTS = WavelengthInterpolatedValue(mTransitTimeSpreadG2Dmap, pWavelength1, pWavelength2 )*ns; 

    return G4RandGauss::shoot(lMeanTransitTime, lTTS); 
}

/**
 * Get interpolated value between two scans of different wavelenths
 * @param pMap  Map with TGraph2D scans to be used
 * @param pWavelength1  first reference wavelength for interpolation 
 * @param pWavelength2  second reference wavelength for interpolation 
 * @return G4double interpolated value
 */
G4double OMSimPMTResponse::WavelengthInterpolatedValue(std::map<G4double, alglib::spline2dinterpolant> pMap ,G4double pWavelength1, G4double pWavelength2){

    // G4double lValue1 = pMap[pWavelength1]->Interpolate(mX, mY);
    // G4double lValue2 = pMap[pWavelength2]->Interpolate(mX, mY);

    G4double lValue1 = alglib::spline2dcalc(pMap[pWavelength1], mX, mY)* ns;
    G4double lValue2 = alglib::spline2dcalc(pMap[pWavelength2], mX, mY)* ns;

    return lValue1 + (mWavelength-pWavelength1) * (lValue2-lValue1) / (pWavelength2 - pWavelength1);
}


/**
 * Get a finished Pulse using a key for the maps directly
 * @param pWavelengthKey  wavelength key to be used 
 * @return PMTPulse Struct with transit time, charge and detection probability
 */
OMSimPMTResponse::PMTPulse OMSimPMTResponse::GetPulseFromKey(G4double pWavelengthKey){
    OMSimPMTResponse::PMTPulse lPulse;
    lPulse.PE = GetCharge(pWavelengthKey);
    lPulse.TransitTime = GetTransitTime(pWavelengthKey);

    return lPulse;
}

/**
 * Get a finished Pulse interpolating scan results between measured wavelengths (simulated photons has a wavelength between these two)
 * @param pWavelengthKey1 first wavelength key to be used 
 * @param pWavelengthKey2 second wavelength key to be used 
 * @return PMTPulse Struct with transit time, charge and detection probability
 */
OMSimPMTResponse::PMTPulse OMSimPMTResponse::GetPulseFromInterpolation(G4double pWavelengthKey1, G4double pWavelengthKey2){
    OMSimPMTResponse::PMTPulse lPulse;
    lPulse.PE = GetCharge(pWavelengthKey1, pWavelengthKey2);
    lPulse.TransitTime = GetTransitTime(pWavelengthKey1, pWavelengthKey2);
    return lPulse;
}



/**
 * Transform a hit on the photocathode to a Pulse, with the transit time, charge and detection probability. 
 * @param pX  x position on photocathode
 * @param pY  x position on photocathode
 * @param pWavelength  Wavelength of photon
 * @return PMTPulse Struct with transit time, charge and detection probability
 */
OMSimPMTResponse::PMTPulse OMSimPMTResponse::ProcessPhotocathodeHit(G4double pX, G4double pY, G4double pWavelength) {
    
    mX = pX/mm;
    mY = pY/mm;
    G4double lR = std::sqrt(mX*mX+mY*mY);
    mWavelength = pWavelength;

    OMSimPMTResponse::PMTPulse lPulse;

    lPulse.DetectionProbability = mRelativeDetectionEfficiencyInterp->Eval(lR);
    
    //If wavelength matches with one of the measured ones
    if ( std::find(mScannedWavelengths.begin(), mScannedWavelengths.end(), pWavelength) != mScannedWavelengths.end() ){
        lPulse = GetPulseFromKey(pWavelength);
    }
    //If wavelength smaller than scanned, return for first measured wavelength 
    else if (pWavelength < mScannedWavelengths.at(0)){
        lPulse = GetPulseFromKey(mScannedWavelengths.at(0));
    } 
    //If wavelength larger than scanned, return for last measured wavelength 
    else if (pWavelength > mScannedWavelengths.back()){
        lPulse = GetPulseFromKey(mScannedWavelengths.back());
    }
    //If wavelength in range of scanned, return interpolated values
    else {
        //get first index of first element larger than seeked wavelength
        auto const lLowerBoundIndex = std::lower_bound(mScannedWavelengths.begin(), mScannedWavelengths.end(), pWavelength) - mScannedWavelengths.begin();

        G4double lW1 = mScannedWavelengths.at(lLowerBoundIndex-1);
        G4double lW2 = mScannedWavelengths.at(lLowerBoundIndex);

        lPulse = GetPulseFromInterpolation(lW1, lW2); 
    }
    return lPulse;
}




