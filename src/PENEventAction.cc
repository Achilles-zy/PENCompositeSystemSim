//PENEventAction.cc
#include <string.h>
#include "PENEventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "g4root.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
PENEventAction::PENEventAction(PENRunAction* runaction)
	: edepBulk(0.),
    SiPMPhotonCount{ 0 },
	Total(0),
	ID(0),
    EscapedPhotonCount(0),
	ifSiPM(false),
	ifBulk(false),
    ifDetectable(false),
	run(runaction)
    //ResultFile("Distribution_Results_NTuple.root","RECREATE"),
    //Distribution_Results("Distribution_Results","Distribution_Results")
{
    PhotonCut_0 = 1;
    PhotonCut_1 = 3;
    PhotonCut_2 = 5;
    PhotonCut_3 = 10;
    SignalSiPMCount_0 = 0;
    SignalSiPMCount_1 = 0;
    SignalSiPMCount_2 = 0;
    SignalSiPMCount_3 = 0;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
PENEventAction::~PENEventAction()
{
  //Distribution_Results.Write();
  //ResultFile.Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void PENEventAction::BeginOfEventAction(const G4Event* evt)
{
  edepBulk = 0;
  memset(SiPMPhotonCount, 0, sizeof(SiPMPhotonCount));
  Total = 0;
  TotalSiPMPhotonCount = 0;
  SignalSiPMCount_0 = 0;
  SignalSiPMCount_1 = 0;
  SignalSiPMCount_2 = 0;
  SignalSiPMCount_3 = 0;
  MinSignalSiPMCount = 2;
  EscapedPhotonCount = 0;
  ifSiPM = false;
  ifBulk = false;
  ifDetectable = false;
  // G4cout<<ID<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void PENEventAction::EndOfEventAction(const G4Event* evt)
{


  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillH1(0, edepBulk);
  analysisManager->FillH1(1, TotalSiPMPhotonCount);
  analysisManager->FillH1(2, TotalSiPMPhotonCount);

  G4int rownb = sizeof(SiPMPhotonCount) / sizeof(SiPMPhotonCount[0]);
  G4int columnnb = sizeof(SiPMPhotonCount[0]) / sizeof(SiPMPhotonCount[0][0]);
  for (int i = 0; i < rownb; i++)
  {
      for (int j = 0; j < columnnb; j++)
      {
          G4int NtupleColumnID = i * columnnb + j;
          analysisManager->FillNtupleIColumn(1, NtupleColumnID, TotalSiPMPhotonCount);
          if (SiPMPhotonCount[i][j] > PhotonCut_0) {
              SignalSiPMCount_0++;
          }
          if (SiPMPhotonCount[i][j] > PhotonCut_1) {
              SignalSiPMCount_1++;
          }
          if (SiPMPhotonCount[i][j] > PhotonCut_2) {
              SignalSiPMCount_2++;
          }
          if (SiPMPhotonCount[i][j] > PhotonCut_3) {
              SignalSiPMCount_3++;
          }
      }
  }

  analysisManager->FillNtupleIColumn(1, 0, TotalSiPMPhotonCount);
  analysisManager->AddNtupleRow(1);
  //G4cout << SiPMPhotonCount << G4endl;

	if (edepBulk > 0 && SignalSiPMCount_0 >= MinSignalSiPMCount) {
		run->CountVetoEvent();
    }
    if (edepBulk > 0 && SignalSiPMCount_1 >= MinSignalSiPMCount) {
        run->CountVetoEvent_1();
    }
    if (edepBulk > 0 && SignalSiPMCount_2 >= MinSignalSiPMCount) {
        run->CountVetoEvent_2();
    }
    if (edepBulk > 0 && SignalSiPMCount_3 >= MinSignalSiPMCount) {
        run->CountVetoEvent_3();
    }

	if (edepBulk > 0) {
		run->CountBulkEvent();
        analysisManager->FillH1(3, TotalSiPMPhotonCount);
	}

	if (TotalSiPMPhotonCount > 0) {
		run->CountSiPMEvent();
	}

    if (ifDetectable == true) {
        run->CountDetectableEvent();
    }

    if (ifDetectable == true && edepBulk > 0) {
        run->CountVetoPossibleEvent();
    }

    if (edepBulk > 2000 * keV && edepBulk < 2100 * keV) {
        run->CountROIEvent();

        if (SignalSiPMCount_1 >= MinSignalSiPMCount) {
            run->CountROIVetoEvent();
        }

        if (ifDetectable == true) {
            run->CountROIVetoPossibleEvent();
        }
    }

    G4int evtID = evt->GetEventID();

    if (evtID % 5000 == 0) {
        G4cout << evtID << G4endl;
    }

  ID++;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PENEventAction::AddBulkEnergy(G4double de)
{
    edepBulk += de;
}

void PENEventAction::AddToSiPM(G4int i,G4int j){
  Total++;
  SiPMPhotonCount[i][j]++;
}