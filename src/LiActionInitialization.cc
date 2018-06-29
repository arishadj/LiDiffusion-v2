#include "LiActionInitialization.hh"
#include "LiPrimaryGeneratorAction.hh"
#include "LiRunAction.hh"
#include "LiEventAction.hh"
#include "LiSteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LiActionInitialization::LiActionInitialization(const LiDetectorConstruction* detector)
 : G4VUserActionInitialization(), fDetector(detector)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LiActionInitialization::~LiActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LiActionInitialization::BuildForMaster() const
{
  SetUserAction(new LiRunAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LiActionInitialization::Build() const
{
  // Actions
  //
  LiRunAction* runAction = new LiRunAction();  
  SetUserAction(runAction);
  
  SetUserAction(new LiPrimaryGeneratorAction(runAction));

  LiEventAction* eventAction = new LiEventAction(runAction);
  SetUserAction(eventAction);

  LiSteppingAction* steppingAction = new LiSteppingAction(eventAction);
  SetUserAction(steppingAction);
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
