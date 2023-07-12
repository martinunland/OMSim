#include "OMSimSteppingAction.hh"

#include "G4VProcess.hh"
#include "G4RunManager.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
//since Geant4.10: include units manually
#include "G4SystemOfUnits.hh"

#include "OMSimAnalysisManager.hh"

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

OMSimSteppingAction::OMSimSteppingAction()
{ 
}


void OMSimSteppingAction::UserSteppingAction(const G4Step* aStep)
{    G4Track* aTrack = aStep->GetTrack();
    
    //kill particles that are stuck... e.g. doing a loop in the pressure vessel
    if ( aTrack-> GetCurrentStepNumber() > 100000) {
        G4cout << "Particle stuck   " <<  aTrack->GetDefinition()->GetParticleName()  << " " << 1239.84193/(aTrack->GetKineticEnergy()/eV)<< G4endl;
        if ( aTrack->GetTrackStatus() != fStopAndKill ) {
            aTrack->SetTrackStatus(fStopAndKill);
        }
    }    




if ( aTrack->GetDefinition()->GetParticleName() == "opticalphoton" ) {
     if ( aTrack->GetTrackStatus() != fStopAndKill ) {
        if (aStep->GetPostStepPoint()->GetMaterial()->GetName() == "RiAbs_Photocathode" ){

        /*
        size_t idx_rindex1        = 0;
        const G4Material* aMaterial = aTrack->GetMaterial();
        G4MaterialPropertiesTable* aMaterialPropertiesTable = aMaterial->GetMaterialPropertiesTable();

        G4MaterialPropertyVector* RindexMPV = aMaterialPropertiesTable->GetProperty(kRINDEX);
        G4cout << RindexMPV->Value( aTrack->GetDynamicParticle()->GetTotalMomentum(), idx_rindex1)<< G4endl;
        */
        G4double lEkin = aTrack->GetKineticEnergy();

        G4double h = 4.135667696E-15*eV*s;
        G4double c = 2.99792458E17*nm/s;

        std::vector<G4String> n = explode(aTrack->GetStep()->GetPreStepPoint()->GetTouchableHandle()->GetVolume(2)->GetName(),'_');

        //OMSimPMTResponse& lPhotocathodeResponse = OMSimPMTResponse::getInstance();
	    OMSimAnalysisManager& lAnalysisManager = OMSimAnalysisManager::getInstance();

        lAnalysisManager.stats_PMT_hit.push_back(atoi(n.at(1)));	
        G4ThreeVector lGlobalPosition = aTrack->GetPosition();
        G4ThreeVector lLocalPosition = aTrack->GetStep()->GetPostStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(lGlobalPosition);

        G4double x = lLocalPosition.x()/mm;
        G4double y = lLocalPosition.y()/mm;
        G4double lR = std::sqrt(x*x+y*y);

        //lAnalysisManager.lPulses.push_back(lPhotocathodeResponse.ProcessPhotocathodeHit(x, y, h*c/lEkin));

            G4ThreeVector lDeltaPos = aTrack->GetVertexPosition() - lGlobalPosition;
            lAnalysisManager.stats_photon_direction.push_back(aTrack->GetMomentumDirection());
            lAnalysisManager.stats_photon_position.push_back(lLocalPosition);//aTrack->GetPosition());
            lAnalysisManager.stats_event_id.push_back(lAnalysisManager.current_event_id);
            lAnalysisManager.stats_photon_flight_time.push_back(aTrack->GetLocalTime());
            lAnalysisManager.stats_photon_track_length.push_back(aTrack->GetTrackLength()/m);
            lAnalysisManager.stats_hit_time.push_back(aTrack->GetGlobalTime());
            lAnalysisManager.stats_photon_energy.push_back(lEkin/eV);
            lAnalysisManager.stats_event_distance.push_back(lDeltaPos.mag()/m);
            aTrack->SetTrackStatus(fStopAndKill);
    }
    }
}


}
