// PENSteppingAction.hh

#ifndef PENSteppingAction_h
#define PENSteppingAction_h 1

#include "G4UserSteppingAction.hh"

#include "G4Types.hh"

class PENDetectorConstruction;
class PENEventAction;
class PENRunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PENSteppingAction : public G4UserSteppingAction
{
  public:
      PENSteppingAction(PENEventAction*, PENRunAction*);
      ~PENSteppingAction() {};

    void UserSteppingAction(const G4Step*);
    inline G4double GetEfficiency(G4double wavelength);

  private:

    PENEventAction*      PENEvent;
    PENRunAction* PENRun;
    G4int SignalSiPMCount;
    G4int ContainerSignalSiPMCount;
    G4bool EnableAcc;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
