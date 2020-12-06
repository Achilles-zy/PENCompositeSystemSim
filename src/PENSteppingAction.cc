// PENSteppingAction.cc

#include "PENSteppingAction.hh"

#include "PENDetectorConstruction.hh"
#include "PENEventAction.hh"
#include "PENRunAction.hh"

#include "G4Track.hh"
#include "G4SteppingManager.hh"

#include "G4ThreeVector.hh"
#include "G4Step.hh"
#include "globals.hh"
#include "G4RunManager.hh"
#include "G4HadronicProcessType.hh"
#include "G4EmProcessSubType.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PENSteppingAction::PENSteppingAction(
				     PENEventAction* evt,PENRunAction* run)
:PENEvent(evt),PENRun(run),EnableAcc(false)
{ }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PENSteppingAction::UserSteppingAction(const G4Step* aStep)
{
	EnableAcc = PENRun->GetAccelerate();
	//G4cout << "EnableAcc: " << PENRun->GetAccelerate() << G4endl;
	auto volume = aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
	auto touchable = aStep->GetPostStepPoint()->GetTouchableHandle();
	auto particle_name = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();

	auto edep = aStep->GetTotalEnergyDeposit();

	//**********************For Acceleration**********************//
	if (EnableAcc == true) {
		if (PENEvent->GetAcceletateStatus() == false) {
			SignalSiPMCount = 0;
			G4int rownb = PENEvent->GetRowNb();
			G4int columnnb = PENEvent->GetColumnNb();
			for (int i = 0; i < rownb; i++)
			{
				for (int j = 0; j < columnnb; j++)
				{
					if (PENEvent->GetSiPMPhotonCount(i, j) > 3) {
						SignalSiPMCount++;
					}
				}
			}
			if (SignalSiPMCount >= 2) {
				PENEvent->SetAccelerate(true);
			}
		}

		if (PENEvent->GetAcceletateStatus() == true && particle_name == "opticalphoton") {
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		}
		//if (PENEvent->GetTotalSiPMPhotonCnt() > 2 && particle_name == "opticalphoton") {
		//	aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		//}
	}

	//***********************************************************//

	const PENDetectorConstruction* detectorConstruction
		= static_cast<const PENDetectorConstruction*>
		(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

	if (volume == detectorConstruction->GetBulk() && particle_name != "opticalphoton") {
		//PENEvent->BulkTrue();
		PENEvent->AddBulkEnergy(edep);
	}
	
	if (volume == detectorConstruction->GetEnv() && particle_name == "opticalphoton") {
		PENEvent->DetectableTrue();
	}

	for (int i = 0; i <= 4; i++) {
		if (volume == detectorConstruction->GetSiPM(i) && particle_name == "opticalphoton" && detectorConstruction->GetSiPM(i) != nullptr) {
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			//G4cout << i << G4endl;
			//G4cout << detectorConstruction->GetSiPM(i)->GetName() << G4endl;
			G4int SiPMArrayID = touchable->GetCopyNumber(1);
			G4int SiPMID = touchable->GetCopyNumber(0);
			PENEvent->AddToSiPM(SiPMArrayID, SiPMID);
			PENEvent->CountTotalSiPMPhoton(1);
		}
	}

	G4int processtype = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessSubType();
	G4int creatorprocess = aStep->GetTrack()->GetCreatorProcess()->GetProcessSubType();
	G4int parentID = aStep->GetTrack()->GetParentID();

	if (parentID == 1 ) {
		if (processtype == fScintillation || processtype == fRadioactiveDecay) {
			//aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
		}
	}

	if (processtype == fRadioactiveDecay) {
		//G4cout << "Parent ID =" <<parentID << G4endl;
	}

	if (parentID == 1) {
		//G4cout << "Parent ID =" << parentID << "Particle =" << aStep->GetTrack()->GetParticleDefinition()->GetParticleName() << "Process =" << aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() <<aStep->GetTotalEnergyDeposit()<< G4endl;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

