#include "LiPrimaryGeneratorAction.hh"
#include "LiPrimaryGeneratorMessenger.hh"
#include "LiDetectorConstruction.hh"
#include "LiRunAction.hh"
#include "LiRun.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4GenericIon.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include "G4DynamicParticle.hh"

#include "G4DecayTable.hh"
#include "G4VDecayChannel.hh"
#include "G4PhaseSpaceDecayChannel.hh"

#include "globals.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <iterator> 

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LiPrimaryGeneratorAction::LiPrimaryGeneratorAction(LiRunAction* runA)
: G4VUserPrimaryGeneratorAction(),fParticleGun(0),
  runAct(runA), posx(0.),posy(0.),sigma(2.) 
{
  primaryGeneratorMessenger = new LiPrimaryGeneratorMessenger(this);

  MakeAlphaDistribution();
  MakeBetaDistribution();
  
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LiPrimaryGeneratorAction::~LiPrimaryGeneratorAction()
{
  delete fParticleGun;
  delete primaryGeneratorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void LiPrimaryGeneratorAction::MakeBetaDistribution()
{
  //read the file with the (experimenetal) beta energy distribution, for the 8Li->8Be decay
  std::vector<G4double> tempBeta;

  G4double temp=0;
  nPointsBeta = 0;
  G4String fnameBeta = "DiffusionProfiles/EnergyHistograms/beta.txt";
  std::ifstream fieldBeta;
  fieldBeta.open(fnameBeta) ;
  if (fieldBeta.is_open())
  {
    std::string lineBeta;
    std::istringstream inStreamBeta; 
    while(!fieldBeta.eof())
    {			
      getline(fieldBeta, lineBeta) ;
      inStreamBeta.clear();
      inStreamBeta.str(lineBeta); 
      inStreamBeta >> temp ;
      tempBeta.push_back(temp);
      nPointsBeta++;
    }
    nPointsBeta--;
  }
  fieldBeta.close();

  //copy arrays in std::vector and compute fMax
  //
  xBeta.resize(nPointsBeta); fBeta.resize(nPointsBeta);
  fMaxBeta = 0.;
  for (G4int j=0; j<nPointsBeta; j++) {
    xBeta[j] = j+1; 
    fBeta[j] = tempBeta[j];
    if (fBeta[j]<0) fBeta[j]=0;	
	
    if (fMaxBeta < fBeta[j]) fMaxBeta = fBeta[j];
  };
  tempBeta.clear();

  //compute slopes
  //
  aBeta.resize(nPointsBeta);
  for (G4int j=0; j<nPointsBeta-1; j++) { 
    aBeta[j] = (fBeta[j+1] - fBeta[j])/(xBeta[j+1] - xBeta[j]);
  };
  
  //compute cumulative function
  //
  FcBeta.resize(nPointsBeta);  
  FcBeta[0] = 0.;
  for (G4int j=1; j<nPointsBeta; j++) {
    FcBeta[j] = FcBeta[j-1] + 0.5*(fBeta[j] + fBeta[j-1])*(xBeta[j] - xBeta[j-1]);
  };
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LiPrimaryGeneratorAction::MakeAlphaDistribution()
{
  //read the file with the (experimenetal) alpha energy distribution, for the 8Be->2a decay
  std::vector<G4double> tempAlphaE;
  std::vector<G4double> tempAlphaN;


  G4double temp1=0;
  G4double temp2=0;
  nPointsAlpha = 0;
  G4String fname = "DiffusionProfiles/EnergyHistograms/AlphaEnergySpectrum.txt";

  std::ifstream field;
  field.open(fname) ;
  if (field.is_open())
  {
    std::string line;
    std::istringstream inStream; 
    getline(field, line) ;
    inStream.clear();
    inStream.str(line);
    while(!field.eof()) 
    {			
      getline(field, line) ;
      inStream.clear();
      inStream.str(line); 
      inStream >> temp1>>temp2;
      tempAlphaE.push_back(temp1);
      tempAlphaN.push_back(temp2);
      nPointsAlpha++;
    }
    nPointsAlpha--;
  }
  field.close();
  
  //copy arrays in std::vector and compute fMax
  //
  x.resize(nPointsAlpha); f.resize(nPointsAlpha);
  fMax = 0.;
  for (G4int j=0; j<nPointsAlpha; j++) {
    x[j] = tempAlphaE[j]; 
    f[j] = tempAlphaN[j];
    if (f[j]<0) f[j]=0;	

    if (fMax < f[j]) fMax = f[j];
  };
  tempAlphaE.clear();
  tempAlphaN.clear();

  //compute slopes
  //
  a.resize(nPointsAlpha);
  for (G4int j=0; j<nPointsAlpha-1; j++) { 
    a[j] = (f[j+1] - f[j])/(x[j+1] - x[j]);
  };
  
  //compute cumulative function
  //
  Fc.resize(nPointsAlpha);  
  Fc[0] = 0.;
  for (G4int j=1; j<nPointsAlpha; j++) {
    Fc[j] = Fc[j-1] + 0.5*(f[j] + f[j-1])*(x[j] - x[j-1]);
  };     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LiPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //generate the primary alphas and beta for each decay, with correct distribution of energies, momentums and decay position
  G4double x0 = G4RandGauss::shoot(posx,sigma)*mm;	
  G4double y0 = G4RandGauss::shoot(posy,sigma)*mm;
  G4double z0 = 0.1*nm; 

  // depth of target
  const LiDetectorConstruction* detectorConstruction
      = static_cast<const LiDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  if(detectorConstruction->GetTime()>0)
  {
    z0=(getDepth(anEvent->GetEventID())+G4UniformRand())*nm;
  }
  else
  {
    G4int initDepth=0;
    G4int finalDepth=initDepth+20001; 
    z0=(initDepth+(G4UniformRand()*(finalDepth-initDepth)))*nm; 
  }  

  //define high energy beta 
  G4String electronName;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* electron
    = particleTable->FindParticle(electronName="e-");
  fParticleGun->SetParticleDefinition(electron);
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

  G4double betaMin =  0*deg;      //beta in [0,pi]
  G4double betaMax = 180*deg;
  G4double fCosBetaMin = std::cos(betaMin);
  G4double fCosBetaMax = std::cos(betaMax);
  
  G4double fPsiBetaMin = 0*deg;       //psi in [0, 2*pi]
  G4double fPsiBetaMax = 360*deg;

  //direction uniform in solid angle
  //
  G4double cosBeta = fCosBetaMin-G4UniformRand()*(fCosBetaMin-fCosBetaMax);
  G4double sinBeta = std::sqrt(1. - cosBeta*cosBeta);
  G4double psiBeta = fPsiBetaMin + G4UniformRand()*(fPsiBetaMax - fPsiBetaMin);

  G4double uxBeta = sinBeta*std::cos(psiBeta),
           uyBeta = sinBeta*std::sin(psiBeta),
           uzBeta = cosBeta;

  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(uxBeta,uyBeta,uzBeta));

  G4double energyOfBeta=(getRandomBetaEnergy())*keV;

  fParticleGun->SetParticleEnergy(energyOfBeta);
  fParticleGun->GeneratePrimaryVertex(anEvent);

  //define high energy alphas 
  G4int Z = 2, A = 4;
  G4double excitEnergy = 0.*MeV;
    
  G4ParticleDefinition* ion
       = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);

  fParticleGun->SetParticleDefinition(ion);
    fParticleGun->SetParticleCharge(2*eplus);

  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

  G4double alphaMin =  0*deg;      //alpha in [0,pi]
  G4double alphaMax = 180*deg;
  G4double fCosAlphaMin = std::cos(alphaMin);
  G4double fCosAlphaMax = std::cos(alphaMax);
  
  G4double fPsiMin = 0*deg;       //psi in [0, 2*pi]
  G4double fPsiMax = 360*deg;

  //direction uniform in solid angle
  //
  G4double cosAlpha = fCosAlphaMin-G4UniformRand()*(fCosAlphaMin-fCosAlphaMax);
  G4double sinAlpha = std::sqrt(1. - cosAlpha*cosAlpha);
  G4double psi = fPsiMin + G4UniformRand()*(fPsiMax - fPsiMin);

  G4double ux = sinAlpha*std::cos(psi),
           uy = sinAlpha*std::sin(psi),
           uz = cosAlpha;

  G4double energyOfAlpha=(getRandomAlphaEnergy())*keV;

  fParticleGun->SetParticleEnergy(energyOfAlpha);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz)); 
  fParticleGun->GeneratePrimaryVertex(anEvent);
 
  fParticleGun->SetParticleEnergy(energyOfAlpha);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(-1*ux,-1*uy,-1*uz));
  fParticleGun->GeneratePrimaryVertex(anEvent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double LiPrimaryGeneratorAction::getRandomAlphaEnergy()
{
  // tabulated function 
  // f is assumed positive, linear per segment, continuous
  //  
  G4double x_rndm = 0., y_rndm = 0., f_inter = -1.;
  
  while (y_rndm > f_inter) {
    //choose a point randomly
    x_rndm = x[0] + G4UniformRand()*(x[nPointsAlpha-1] - x[0]);
    y_rndm = G4UniformRand()*fMax;
    //find bin
    G4int j = nPointsAlpha-2;
    while ((x[j] > x_rndm) && (j > 0)) j--;
    //compute f(x_rndm) by linear interpolation 
   f_inter = f[j] + a[j]*(x_rndm - x[j]);
  };
  return x_rndm;
}

G4double LiPrimaryGeneratorAction::getRandomBetaEnergy()
{
  // tabulated function 
  // f is assumed positive, linear per segment, continuous
  //  
  G4double x_rndm = 0., y_rndm = 0., f_inter = -1.;
  while (y_rndm > f_inter) {
    //choose a point randomly
    x_rndm = xBeta[0] + G4UniformRand()*(xBeta[nPointsBeta-1] - xBeta[0]);
    y_rndm = G4UniformRand()*fMaxBeta;
    //find bin
    G4int j = nPointsBeta-2;
    while ((xBeta[j] > x_rndm) && (j > 0)) j--;
    //compute f(x_rndm) by linear interpolation
    f_inter = fBeta[j] + aBeta[j]*(x_rndm - xBeta[j]);
  };

  return x_rndm;
}

G4double LiPrimaryGeneratorAction::getDepth(G4int eventNb)
{
  // tabulated function 
  // f is assumed positive, linear per segment, continuous
  //  
  const LiDetectorConstruction* detectorConstruction
      = static_cast<const LiDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  const LiRun* run 
    = static_cast<const LiRun*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  nPointsZ=detectorConstruction->GetNpointsZ();
  xZ=detectorConstruction->GetXz();
  fZ=detectorConstruction->GetFz();   
  NLi=run->GetNumberOfEventToBeProcessed() ; 

  G4double xOut = 0.;

  nZ.resize(nPointsZ);
  FcZ.resize(nPointsZ);
  std::fill(nZ.begin(), nZ.end(), 0);


  for(G4int i=0;i<nPointsZ;i++)
  {
     double temp=static_cast<double>(fZ[i]*NLi);
     nZ[i]= static_cast<G4double>(std::round(temp));	

  }
  FcZ[0]=nZ[0];

  for(G4int i=1;i<nPointsZ;i++)
  {
     FcZ[i]=FcZ[i-1]+nZ[i];
  }

  bool FoundX=false;

  G4int j=0;
  while(!FoundX)
  {
    if(eventNb<=FcZ[j])
    {
       FoundX=true;
    }
    else
    {
       j++;
    }
  }

  xOut=xZ[j];

  return xOut;
}
