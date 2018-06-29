#ifndef LiRunAction_h
#define LiRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class LiRunAction : public G4UserRunAction
{
  public:
    LiRunAction();
    virtual ~LiRunAction();

    virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
 
  private:
    G4String         filename;
    G4int 	     diagnostics;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

