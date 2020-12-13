#ifndef PENEventAction_h
#define PENEventAction_h 1

#include "G4UserEventAction.hh"
#include "PENRunAction.hh"
#include "globals.hh"
//#include "TROOT.h"
//#include "TFile.h"
//#include "TNtuple.h"
//#include "Rtypes.h"

class TNtuple;
class TFile;
class G4Event;

class PENEventAction : public G4UserEventAction
{
  public:
    PENEventAction(PENRunAction* runaction);
   ~PENEventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
    void AddBulkEnergy(G4double);
	void SiPMTrue() { ifSiPM = true; }
	void BulkTrue() { ifBulk = true; }
    void DetectableTrue() { ifDetectable = true; }
    void SetAccelerate(G4bool ifacc) { ifAccelerate = ifacc; }
	void CountTotalSiPMPhoton(G4int ph) { TotalSiPMPhotonCount = TotalSiPMPhotonCount + ph; }
    void CountEscapedPhoton(G4int ph) { EscapedPhotonCount = EscapedPhotonCount + ph; }
    G4int GetTotalSiPMPhotonCnt() { return TotalSiPMPhotonCount; }

    G4int GetRowNb() { return RowNb; }
    G4int GetColumnNb() { return ColumnNb; }
    G4int GetContainerRowNb() { return RowNb; }
    G4int GetContainerColumnNb() { return ColumnNb; }
    G4int GetEscapedPhotonCnt() { return EscapedPhotonCount; }
    G4bool GetAcceletateStatus() { return ifAccelerate; }

    void AddToSiPM(G4int, G4int);
    void AddToSiPMSignal(G4int, G4int);
    G4int GetSiPMPhotonCount(G4int j, G4int k) { return SiPMPhotonCount[j][k]; }
    G4int GetSiPMSignalCount(G4int j, G4int k) { return SiPMSignalCount[j][k]; }

    void AddToContainerSiPM(G4int, G4int);
    void AddToContainerSiPMSignal(G4int, G4int);
    G4int GetContainerSiPMPhotonCount(G4int j, G4int k) { return ContainerSiPMPhotonCount[j][k]; }
    G4int GetContainerSiPMSignalCount(G4int j, G4int k) { return ContainerSiPMSignalCount[j][k]; }

  private:
    G4double edepBulk;
    G4int SiPMPhotonCount[100][4];
    G4int SiPMSignalCount[100][4];
    G4int ContainerSiPMPhotonCount[250][4];
    G4int ContainerSiPMSignalCount[250][4];
    G4int PhotonCut_0;
    G4int PhotonCut_1;
    G4int PhotonCut_2;
    G4int PhotonCut_3;
    G4int PhotonCut_4;
    G4int RowNb;
    G4int ColumnNb;
    G4int ContainerRowNb;
    G4int ContainerColumnNb;
    G4int SignalSiPMCount;
    G4int SignalSiPMCount_0;
    G4int SignalSiPMCount_1;
    G4int SignalSiPMCount_2;
    G4int SignalSiPMCount_3;
    G4int SignalSiPMCount_4;
    G4int ContainerSignalSiPMCount;
    G4int MinSignalSiPMCount;

    G4int Total;

    G4int ID;
	G4int TotalSiPMPhotonCount;
    G4int EscapedPhotonCount;
    G4double EnergyThreshold;
    
	G4bool ifSiPM;
	G4bool ifBulk;
    G4bool ifDetectable;
    G4bool ifAccelerate;

	PENRunAction* run;
    //TFile ResultFile;
    //TTree Distribution_Results;
};

#endif