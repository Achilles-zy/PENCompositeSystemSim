#include "PENPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"

//#include "TMath.h"
#include "Randomize.hh"

PENPrimaryGeneratorAction::PENPrimaryGeneratorAction(PENDetectorConstruction* det):
	G4VUserPrimaryGeneratorAction(),
	PrimaryE(0),
	PrimaryName(""),
	SrcType("PENShell")
{
    fPENGPS = new G4GeneralParticleSource();
	fPrimaryMessenger = new PENPrimaryGeneratorMessenger(this);
	fDetCons = det;
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName = "e-";
	G4double particleEnergy = 0.1 * MeV;
	fPENGPS->SetParticleDefinition(particleTable->FindParticle(particleName));

	/*
	G4PhysicalVolumeStore* PVStore = G4PhysicalVolumeStore::GetInstance();
	G4int i = 0;
	G4VPhysicalVolume* tempPV = NULL;
	while (i < G4int(PVStore->size())) {
		tempPV = (*PVStore)[i];
		G4cout << i << " " << " " << tempPV->GetName() << " " << G4endl;
		i++;
	}
	*/
}

PENPrimaryGeneratorAction::~PENPrimaryGeneratorAction()
{
    delete fPENGPS;
}

void PENPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	//G4String particleName = "e-";
	//G4double particleEnergy = 0 * MeV;
	G4String mode = fDetCons->GetMode();
	if (mode == "Unit") {
		if (SrcType == "Wire") {
			G4double Radius = fDetCons->GetWireRadius();
			G4double Length = fDetCons->GetWireLength();
			G4ThreeVector WirePos = fDetCons->GetWirePos();
			G4cout << "==========================Primary Info==========================" << G4endl;
			G4cout << "Wire Position: " << WirePos << G4endl;
			G4cout << "Wire Radius: " << Radius << G4endl;
			G4cout << "Wire Length: " << Length << G4endl;
			G4cout << "================================================================" << G4endl;
			fPENGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fPENGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
			fPENGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
			fPENGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
			fPENGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(WirePos);
			fPENGPS->GetCurrentSource()->GetPosDist()->SetRadius(Radius);
			fPENGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(Length / 2);

			fPENGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("Wire");
		}
		else if (SrcType == "PENShell") {
			G4double Radius = fDetCons->GetPENShellRadius();
			G4double Length = fDetCons->GetPENShellLength();
			G4cout << "==========================Primary Info==========================" << G4endl;
			G4cout << "Sample Region Radius: " << Radius << G4endl;
			G4cout << "Sample Region Length: " << Length << G4endl;
			G4cout << "================================================================" << G4endl;
			fPENGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fPENGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
			fPENGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
			fPENGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
			fPENGPS->GetCurrentSource()->GetPosDist()->SetRadius(Radius * 2);
			fPENGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(Length / 2);

			fPENGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("PENShell");
		}

		else {
			G4cout << "Error: Src type not found!" << G4endl;
		}
	}
	else {
		G4double Radius = fDetCons->GetPENShellRadius();
		G4double Length = fDetCons->GetPENShellLength();
		G4cout << "==========================Primary Info==========================" << G4endl;
		G4cout << "Sample Region Radius: " << Radius << G4endl;
		G4cout << "Sample Region Length: " << Length << G4endl;
		G4cout << "================================================================" << G4endl;
		fPENGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
		fPENGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
		fPENGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
		fPENGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
		fPENGPS->GetCurrentSource()->GetPosDist()->SetRadius(Radius * 2);
		fPENGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(Length / 2);

		//fPENGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("PENShell");
	}

    fPENGPS->GeneratePrimaryVertex(anEvent);
	anEvent->GetEventID();

	if (anEvent->GetEventID() == 0) {
		PrimaryE = fPENGPS->GetCurrentSource()->GetParticleEnergy();
		PrimaryName = fPENGPS->GetCurrentSource()->GetParticleDefinition()->GetParticleName();
	}
}