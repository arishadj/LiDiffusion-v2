#include "LiPrimaryGeneratorMessenger.hh"
#include "LiPrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LiPrimaryGeneratorMessenger::LiPrimaryGeneratorMessenger(
                                           LiPrimaryGeneratorAction* LiAct)
  :LiAction(LiAct),LiDir(0),ActDir(0)
{
  //Commands for user-specified beam characteristics
  LiDir = new G4UIdirectory("/LiDiffusionCode/");
  LiDir->SetGuidance("UI commands of this example");
 
  ActDir = new G4UIdirectory("/LiDiffusionCode/action/");
  ActDir->SetGuidance("primary generator control");

  PosxCmd = new G4UIcmdWithADoubleAndUnit("/LiDiffusionCode/action/setBeamPosx",this);	
  PosxCmd->SetGuidance("Set the x-value of the center of the beam, in mm");
  PosxCmd->SetParameterName("posx",true);
  PosxCmd->SetUnitCategory("Length");
  PosxCmd->SetRange("posx>=-10 && posx<=10");
  PosxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  PosyCmd = new G4UIcmdWithADoubleAndUnit("/LiDiffusionCode/action/setBeamPosy",this);	
  PosyCmd->SetGuidance("Set the y-value of the center of the beam, in mm");
  PosyCmd->SetParameterName("posy",true);
  PosyCmd->SetUnitCategory("Length");
  PosyCmd->SetRange("posy>=-10 && posy<=10");
  PosyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SigmaCmd = new G4UIcmdWithADoubleAndUnit("/LiDiffusionCode/action/setBeamSigma",this);	
  SigmaCmd->SetGuidance("Set the sigma value of the distribution of the beam, in mm");
  SigmaCmd->SetParameterName("sigma",true);
  SigmaCmd->SetUnitCategory("Length");
  SigmaCmd->SetRange("sigma>=0. && sigma<=10.");
  SigmaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LiPrimaryGeneratorMessenger::~LiPrimaryGeneratorMessenger()
{
  delete LiDir;
  delete ActDir;
  delete PosxCmd;
  delete PosyCmd;
  delete SigmaCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LiPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
 //Passing the user-specified values to the simulation
 if( command == PosxCmd )
   { LiAction->SetPosx(PosxCmd->GetNewDoubleValue(newValue));}

 if( command == PosyCmd )
   { LiAction->SetPosy(PosyCmd->GetNewDoubleValue(newValue));}

 if( command == SigmaCmd )
   { LiAction->SetSigma(SigmaCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
