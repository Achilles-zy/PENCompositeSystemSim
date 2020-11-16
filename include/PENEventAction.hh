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
    void AddToSiPM(G4int, G4int);
	void SiPMTrue() { ifSiPM = true; }
	void BulkTrue() { ifBulk = true; }
    void DetectableTrue() { ifDetectable = true; }
	void CountTotalSiPMPhoton(G4int ph) { TotalSiPMPhotonCount = TotalSiPMPhotonCount + ph; }
    void CountEscapedPhoton(G4int ph) { EscapedPhotonCount = EscapedPhotonCount + ph; }
    G4int GetSiPMPhotonCnt() { return TotalSiPMPhotonCount; }
    G4int GetEscapedPhotonCnt() { return EscapedPhotonCount; }

  private:
    G4double edepBulk;
    G4int SiPMPhotonCount[3][5];
    G4int PhotonCut_0;
    G4int PhotonCut_1;
    G4int PhotonCut_2;
    G4int PhotonCut_3;
    G4int SignalSiPMCount_0;
    G4int SignalSiPMCount_1;
    G4int SignalSiPMCount_2;
    G4int SignalSiPMCount_3;
    G4int MinSignalSiPMCount;

    G4int Total;

    G4int ID;
	G4int TotalSiPMPhotonCount;
    G4int EscapedPhotonCount;
    
	G4bool ifSiPM;
	G4bool ifBulk;
    G4bool ifDetectable;

	PENRunAction* run;
    //TFile ResultFile;
    //TTree Distribution_Results;
};

#endif