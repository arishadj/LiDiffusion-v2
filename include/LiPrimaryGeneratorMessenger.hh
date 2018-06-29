#ifndef LiPrimaryGeneratorMessenger_h
#define LiPrimaryGeneratorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class LiPrimaryGeneratorAction;

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADouble ;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class LiPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    LiPrimaryGeneratorMessenger(LiPrimaryGeneratorAction* );
   ~LiPrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    LiPrimaryGeneratorAction*   LiAction;
    
    G4UIdirectory*              LiDir;
    G4UIdirectory*              ActDir;
    G4UIcmdWithADoubleAndUnit*  PosxCmd ;
    G4UIcmdWithADoubleAndUnit*  PosyCmd ;
    G4UIcmdWithADoubleAndUnit*  SigmaCmd ;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
