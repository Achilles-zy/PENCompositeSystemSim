#include "PENPrimaryGeneratorMessenger.hh"

#include "PENPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"
#include "G4SystemOfUnits.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PENPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;

PENPrimaryGeneratorMessenger::PENPrimaryGeneratorMessenger(PENPrimaryGeneratorAction* Gun)
  : G4UImessenger(),
    fAction(Gun),
    fSrcDir(0),
    cmdSetSrcType(0)
{
  fSrcDir = new G4UIdirectory("/PEN/src/");
  fSrcDir->SetGuidance("PrimaryGenerator control");

  cmdSetSrcType = new G4UIcmdWithAString("/PEN/src/type", this);
  cmdSetSrcType->SetGuidance("Choose the type of source");
  cmdSetSrcType->SetParameterName("SrcType",true);
  cmdSetSrcType->SetDefaultValue("PENShell");
  cmdSetSrcType->AvailableForStates(G4State_PreInit, G4State_Idle);

  cmdSetImprintID = new G4UIcmdWithAnInteger("/PEN/src/imprintid", this);
  cmdSetImprintID->SetGuidance("Choose a detector unit");
  cmdSetImprintID->SetParameterName("ImprintID", true);
  cmdSetImprintID->SetDefaultValue(1);
  cmdSetImprintID->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PENPrimaryGeneratorMessenger::~PENPrimaryGeneratorMessenger()
{
  delete cmdSetSrcType;
  delete cmdSetImprintID;
  delete fSrcDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PENPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == cmdSetSrcType) {
  	fAction->SetSrcType(newValue);
  }
  if (command == cmdSetImprintID) {
      fAction->SetImprintID(cmdSetImprintID->GetNewIntValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
