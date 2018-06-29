#include "LiEventAction.hh"
#include "LiRun.hh"
#include "LiRunAction.hh"
#include "Analysis.hh"
#include "LiDetectorConstruction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Threading.hh" 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LiEventAction::LiEventAction(LiRunAction* runA)
:G4UserEventAction(),runAct(runA)
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LiEventAction::~LiEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LiEventAction::BeginOfEventAction(const G4Event*)
{   
  //resetting parameters at the start of each event 
  fEdepAlpha   		 = 0.;
  fEdepAlphaFromAlphas   = 0.;
  fEdepAlphaFromBetas    = 0.;
  fEdepBeta    		 = 0.;
  depth        		 = 0.;
  aEnergy1     		 = 0.;
  aEnergy2     		 = 0.;
  bEnergy      		 = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LiEventAction::EndOfEventAction(const G4Event*)
{   
  // accumulate statistics in LiRun
  LiRun* run 
    = static_cast<LiRun*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->AddEdep(fEdepAlpha);

  //get the detector class, to pass information from it
  const LiDetectorConstruction* detectorConstruction
      = static_cast<const LiDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  diagnostics=detectorConstruction->GetDiagnostics();
  G4double time =detectorConstruction->GetTime()/s;

  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Fill the output file with the data of the current event
  analysisManager->FillNtupleDColumn(0, time); 
  analysisManager->FillNtupleDColumn(1, fEdepAlpha); 
  analysisManager->FillNtupleDColumn(2, fEdepBeta); 
  analysisManager->FillNtupleDColumn(3, depth); 
  analysisManager->AddNtupleRow(); 

  analysisManager->FillH1(1, depth); 

  // Optional information
  if(diagnostics==1)
  {
    analysisManager->FillH1(2, aEnergy1);
    analysisManager->FillH1(2, aEnergy2);
    analysisManager->FillH1(3, bEnergy);
  }

  analysisManager->FillH1(4, fEdepAlpha); 
  analysisManager->FillH1(5, fEdepAlphaFromAlphas); 
  analysisManager->FillH1(6, fEdepAlphaFromBetas); 
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
