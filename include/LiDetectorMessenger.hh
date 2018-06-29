#ifndef LiDetectorMessenger_h
#define LiDetectorMessenger_h 1
#include "G4UImessenger.hh"
#include "globals.hh"

class LiDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class LiDetectorMessenger: public G4UImessenger
{
public:

  LiDetectorMessenger(LiDetectorConstruction* );
  ~LiDetectorMessenger();
     
  virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
   LiDetectorConstruction*      fDetector;
   
   G4UIdirectory*             LiDir;
   G4UIdirectory*             DetDir;
   G4UIcmdWithAString*        SampleMaterialCmd; 
   G4UIcmdWithAString*        CreateSampleMaterialCmd;
   G4UIcmdWithAString*        SubstrateMaterialCmd;
   G4UIcmdWithAString*        CreateSubstrateMaterialCmd;
   G4UIcmdWithADoubleAndUnit* SampleWidthCmd;
   G4UIcmdWithADoubleAndUnit* SampleDepthCmd;
   G4UIcmdWithADoubleAndUnit* SubstrateWidthCmd;
   G4UIcmdWithADoubleAndUnit* SubstrateDepthCmd;
   G4UIcmdWithADoubleAndUnit* DetWidthCmd;
   G4UIcmdWithADoubleAndUnit* DetDepthCmd;
   G4UIcmdWithADoubleAndUnit* DetRadiousCmd;
   G4UIcmdWithoutParameter*   UpdateCmd;

   G4UIcmdWithADoubleAndUnit* TimeCmd;
   G4UIcmdWithAnInteger*      DiagnosticsCmd;
   G4UIcmdWithAString*        FileNameCmd;
   G4UIcmdWithAString*        SRIMCmd;
   G4UIcmdWithADouble*        DiffRateCmd;
   G4UIcmdWithAString*        BoundaryCondCmd;
   G4UIcmdWithAnInteger*      NLiCmd;
   G4UIcmdWithADouble*        BeamOnCmd;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
