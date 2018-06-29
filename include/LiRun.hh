#ifndef LiRun_h
#define LiRun_h 1

#include "G4Run.hh"
#include "globals.hh"

class G4Event;

class LiRun : public G4Run
{
  public:
    LiRun();
    virtual ~LiRun();

    virtual void Merge(const G4Run*);
    
    void AddEdep (G4double edep); 
    G4double GetEdep()  const { return fEdep; }
    virtual G4int GetNumberOfEvents() {return N;};

  private:
    G4double  fEdep;
    G4int     N ;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

