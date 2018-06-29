#ifndef LiSteppingAction_h
#define LiSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class LiEventAction;

class G4LogicalVolume;

class LiSteppingAction : public G4UserSteppingAction
{
  public:
    LiSteppingAction(LiEventAction* eventAction);
    virtual ~LiSteppingAction();

    // method from the base class
    virtual void UserSteppingAction(const G4Step*);
    virtual void SetDiagnostics(G4int d) {diagnostics=d;};

  private:
    LiEventAction*   fEventAction;

    G4LogicalVolume* fScoringVolume;    
    G4LogicalVolume* fScoringVolumeBeta;
    G4LogicalVolume* fTargetVolume;
    G4LogicalVolume* fOut;

    G4int            diagnostics;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
