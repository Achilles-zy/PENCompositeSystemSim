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
	InitialE(1 * keV),
	PrimaryName(""),
	SrcType("PENShell")
{
    fPENGPS = new G4GeneralParticleSource();
	fPrimaryMessenger = new PENPrimaryGeneratorMessenger(this);
	fDetCons = det;
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName = "e-";
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

			//G4cout << "==========================Primary Info==========================" << G4endl;
			//G4cout << "Wire Position: " << WirePos << G4endl;
			//G4cout << "Wire Radius: " << Radius << G4endl;
			//G4cout << "Wire Length: " << Length << G4endl;
			//G4cout << "================================================================" << G4endl;

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

			if (G4RunManager::GetRunManager()->GetRunManagerType() == 1) {
				G4cout << "==========================Primary Info==========================" << G4endl;
				G4cout << "Sample Region Radius: " << Radius << G4endl;
				G4cout << "Sample Region Length: " << Length << G4endl;
				G4cout << "================================================================" << G4endl;
			}

			fPENGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fPENGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
			fPENGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
			fPENGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
			fPENGPS->GetCurrentSource()->GetPosDist()->SetRadius(Radius * 2);
			fPENGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(Length / 2);

			fPENGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("PENShell");
		}
		else if (SrcType == "Point") {
			fPENGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Point");
			fPENGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0, 0, 25 * mm));
			fPENGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fPENGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
			fPENGPS->GetCurrentSource()->GetEneDist()->SetEmax(InitialE);
			fPENGPS->GetCurrentSource()->GetEneDist()->SetEmin(InitialE);
		}
		else if (SrcType == "InnerShellSurface") {
			fPENGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Surface");
			fPENGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fPENGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
			fPENGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
			fPENGPS->GetCurrentSource()->GetPosDist()->SetRadius(51 * mm);
			fPENGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(51 * mm / 2);
		}
		else {
			G4cout << "Error: Src type not found! Using Geant4 default settings." << G4endl;
		}
	}
	else if (mode == "SArUnit") {
		if (SrcType == "Crystal") {
			G4double Radius = fDetCons->GetSArBrickRadius();
			G4double Length = fDetCons->GetSArBrickHeight();
			if (G4RunManager::GetRunManager()->GetRunManagerType() == 1) {
				G4cout << "==========================Primary Info==========================" << G4endl;
				G4cout << "Sample Region Radius: " << Radius << G4endl;
				G4cout << "Sample Region Length: " << Length << G4endl;
				G4cout << "================================================================" << G4endl;
			}

			fPENGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fPENGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
			fPENGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
			fPENGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
			fPENGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0, 0, Length / 2));
			fPENGPS->GetCurrentSource()->GetPosDist()->SetRadius(Radius);
			fPENGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(Length / 2);

			fPENGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("SArCrystal");
		}

		//Needs to be fixed
		
		else if (SrcType == "Container") {
			G4double Radius = fDetCons->GetSArBrickRadius();
			G4double Length = fDetCons->GetSArBrickHeight();
			if (G4RunManager::GetRunManager()->GetRunManagerType() == 1) {
				G4cout << "==========================Primary Info==========================" << G4endl;
				G4cout << "Sample Region Radius: " << Radius << G4endl;
				G4cout << "Sample Region Length: " << Length << G4endl;
				G4cout << "================================================================" << G4endl;
			}
			fPENGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fPENGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
			fPENGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
			fPENGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
			fPENGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0, 0, Length / 2));
			fPENGPS->GetCurrentSource()->GetPosDist()->SetRadius(Radius);
			fPENGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(Length / 2);

			fPENGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("SArBrick");
		}
		
		else {
			G4cout << "Error: Src type not found! Using Geant4 default settings." << G4endl;
		}
	}
	else {
		G4cout << "Error: Mode not found! Using Geant4 default settings." << G4endl;
	}

    fPENGPS->GeneratePrimaryVertex(anEvent);
	anEvent->GetEventID();

	if (anEvent->GetEventID() == 0) {
		PrimaryE = fPENGPS->GetCurrentSource()->GetParticleEnergy();
		PrimaryName = fPENGPS->GetCurrentSource()->GetParticleDefinition()->GetParticleName();
	}
}