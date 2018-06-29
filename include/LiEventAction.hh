#ifndef LiEventAction_h
#define LiEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class LiRunAction;

/// Event action class
///
class LiEventAction : public G4UserEventAction
{
  public:
    LiEventAction(LiRunAction*);
    virtual ~LiEventAction();
    
    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    void AddEdepAlpha(G4double edep) { fEdepAlpha += edep; }
    void AddEdepAlphaFromAlphas(G4double edepA) { fEdepAlphaFromAlphas += edepA; }
    void AddEdepAlphaFromBetas(G4double edepB) { fEdepAlphaFromBetas += edepB; }
    void AddEdepBeta(G4double edep2) { fEdepBeta += edep2; }
    void SetDepth(G4double d) {depth=d ;}
    void AddAlphaEnergy(G4double e) {if (aEnergy1==0.) 
					{aEnergy1 = e;}
				     else{ aEnergy2=e;};}
    void AddBetaEnergy(G4double e) {bEnergy=e;}
    G4int GetDiagnostics(){return diagnostics;}

  private:
    LiRunAction*  runAct;

    G4double  fEdepAlpha;
    G4double  fEdepAlphaFromAlphas;
    G4double  fEdepAlphaFromBetas;
    G4double  fEdepBeta;
    G4double  depth;
    G4double  aEnergy1 ;
    G4double  aEnergy2 ;
    G4double  bEnergy ;
    G4int     diagnostics;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
