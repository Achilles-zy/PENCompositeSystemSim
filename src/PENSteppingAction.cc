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
#include "G4SystemOfUnits.hh"

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
	G4int TrackID = aStep->GetTrack()->GetTrackID();
	//**********************For Acceleration**********************//
	if (EnableAcc == true) {
		if (PENEvent->GetAcceletateStatus() == false && TrackID % 200 == 0) {
			SignalSiPMCount = 0;
			ContainerSignalSiPMCount = 0;
			G4int rownb = PENEvent->GetRowNb();
			G4int columnnb = PENEvent->GetColumnNb();
			for (int i = 0; i < rownb; i++)
			{
				for (int j = 0; j < columnnb; j++)
				{
					if (PENEvent->GetSiPMSignalCount(i, j) > 0) {
						SignalSiPMCount++;
					}
				}
			}

			G4int crownb = PENEvent->GetContainerRowNb();
			G4int ccolumnnb = PENEvent->GetContainerColumnNb();
			for (int i = 0; i < crownb; i++)
			{
				for (int j = 0; j < ccolumnnb; j++)
				{
					if (PENEvent->GetContainerSiPMSignalCount(i, j) > 0) {
						ContainerSignalSiPMCount++;
					}
				}
			}
			G4int totalcnt = SignalSiPMCount + ContainerSignalSiPMCount;
			if (totalcnt >= 2) {
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

	//G4cout << aStep->GetPostStepPoint()->GetPosition() << G4endl;
	for (int i = 0; i < 5; i++) {
		if (volume == detectorConstruction->GetSiPM(i) && detectorConstruction->GetSiPM(i) != nullptr && particle_name == "opticalphoton" ) {
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);

			G4int GetCopyNumber0 = touchable->GetCopyNumber(0);
			G4int GetCopyNumber1 = touchable->GetCopyNumber(1);
			G4double Energy = aStep->GetPostStepPoint()->GetKineticEnergy() / (1 * eV);
			G4double WaveLength = 1242 / Energy;//nm
			G4double SiPMEff = GetEfficiency(WaveLength);
			G4double rnd = G4UniformRand();
			if (rnd < SiPMEff) {
				PENEvent->AddToSiPMSignal(GetCopyNumber1, GetCopyNumber0);
			}
			PENEvent->AddToSiPM(GetCopyNumber1, GetCopyNumber0);
			PENEvent->CountTotalSiPMPhoton(1);
		}
	}

	for (int j = 0; j < 5; j++) {
		if (volume == detectorConstruction->GetContainerSiPM(j) && detectorConstruction->GetContainerSiPM(j) != nullptr && particle_name == "opticalphoton") {
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);

			G4int GetCopyNumber0 = touchable->GetCopyNumber(0);
			G4int GetCopyNumber1 = touchable->GetCopyNumber(1);
			G4double Energy = aStep->GetPostStepPoint()->GetKineticEnergy() / (1 * eV);
			G4double WaveLength = 1242 / Energy;//nm
			G4double SiPMEff = GetEfficiency(WaveLength);
			G4double rnd = G4UniformRand();
			if (rnd < SiPMEff) {
				PENEvent->AddToContainerSiPMSignal(GetCopyNumber1, GetCopyNumber0);
			}
			PENEvent->AddToContainerSiPM(GetCopyNumber1, GetCopyNumber0);
			PENEvent->CountTotalSiPMPhoton(1);
		}
	}

	//G4int processtype = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessSubType();
	//G4int creatorprocess = aStep->GetTrack()->GetCreatorProcess()->GetProcessSubType();
	//G4int parentID = aStep->GetTrack()->GetParentID();

	//if (parentID == 1 ) {
	//	if (processtype == fScintillation || processtype == fRadioactiveDecay) {
	//		//aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
	//	}
	//}

	//if (processtype == fRadioactiveDecay) {
	//	//G4cout << "Parent ID =" <<parentID << G4endl;
	//}

	//if (parentID == 1) {
	//	//G4cout << "Parent ID =" << parentID << "Particle =" << aStep->GetTrack()->GetParticleDefinition()->GetParticleName() << "Process =" << aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() <<aStep->GetTotalEnergyDeposit()<< G4endl;
	//}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PENSteppingAction::GetEfficiency(G4double wavelength) {
	G4double eff;
	if (wavelength > 400 && wavelength <= 900) {
		eff = 0.5 - (wavelength - 400 * nm) / 1000;
	}
	else if (wavelength > 100 && wavelength <= 400) {
		eff = (wavelength - 100) * 5 / 3000;
	}
	else {
		eff = 0;
	}
	return eff;
}