#ifndef LiPrimaryGeneratorAction_h
#define LiPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class LiRunAction;
class G4ParticleGun;
class G4Event;

class LiPrimaryGeneratorMessenger;

class LiPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
   LiPrimaryGeneratorAction(LiRunAction* runA);    
   virtual ~LiPrimaryGeneratorAction();

   // method from the base class
   virtual void GeneratePrimaries(G4Event*);   

   void SetPosx(G4double x1){posx=x1;};
   void SetPosy(G4double y1){posy=y1;};
   void SetSigma(G4double s){sigma=s;};

   const G4ParticleGun* GetParticleGun() const { return fParticleGun; }
   void MakeAlphaDistribution();
   void MakeBetaDistribution();

   G4double getRandomAlphaEnergy() ;
   G4double getRandomBetaEnergy() ;
   G4double getDepth(G4int eventNb);

  private:
   LiPrimaryGeneratorMessenger*   primaryGeneratorMessenger;
   G4ParticleGun*  		  fParticleGun; // pointer a to G4 gun class
   LiRunAction*		  	  runAct;
  
   G4double 			  posx;
   G4double 			  posy;
   G4double 			  sigma;

   G4int                  	  nPointsAlpha;     
   G4int                  	  nPointsBeta;     
   G4int                  	  nPointsZ;     
   G4double 			  time;

   std::vector<G4double>  	  x;
   std::vector<G4double>  	  f;           //f(x)
   std::vector<G4double>  	  a;           //slopes
   std::vector<G4double>  	  Fc;          //cumulative of f
   G4double               	  fMax;        //max(f)

   std::vector<G4double>  	  xBeta;
   std::vector<G4double>  	  fBeta;           //f(x)
   std::vector<G4double>  	  aBeta;           //slopes
   std::vector<G4double>  	  FcBeta;          //cumulative of f
   G4double               	  fMaxBeta;        //max(f)
 
   std::vector<G4double>  	  xZ;
   std::vector<G4double>  	  fZ;           //f(x)
   std::vector<G4double>  	  nZ;
   std::vector<G4double>  	  FcZ;          //cumulative of f
   G4double                      NLi;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


