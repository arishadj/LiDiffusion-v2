#include "LiRunAction.hh"
#include "LiEventAction.hh"
#include "LiPrimaryGeneratorAction.hh"
#include "LiDetectorConstruction.hh"
#include "LiRun.hh"
#include "Analysis.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <iterator>  
#include <cmath> 
#include <stdio.h>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LiRunAction::LiRunAction()
  : G4UserRunAction(),filename("outputFile.root"),diagnostics(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LiRunAction::~LiRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* LiRunAction::GenerateRun()
{
  return new LiRun; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LiRunAction::BeginOfRunAction(const G4Run* aRun)
{ 

  // Default settings for creating output ROOT file and its data structures 
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;
 
  analysisManager->SetNtupleMerging(true);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFileName("temp.root"); 
  analysisManager->SetFirstHistoId(1);

  analysisManager->CreateH1("depth","Depth distribution in sample", 1000, 0.,1000.);
  analysisManager->CreateH1("Ealpha","Energy of alphas at vertex", 100, 0., 10000);
  analysisManager->CreateH1("Ebeta","Energy of beta at vertex", 150, 0., 15000);
  analysisManager->CreateH1("EnergyDeposition","Energy deposition in the alpha ring detector", 800, 0., 8000);
  analysisManager->CreateH1("EnergyDepositionFromAlphas","Energy deposition in the alpha ring detector by alphas", 800, 0., 8000);
  analysisManager->CreateH1("EnergyDepositionFromBetas","Energy deposition in the alpha ring detector by betas", 800, 0., 8000);

  // Creating ntuple
  //
  analysisManager->CreateNtuple("EnergiesAndDepth","Edep and depth");
  analysisManager->CreateNtupleDColumn("time");
  analysisManager->CreateNtupleDColumn("EdepAlpha");
  //analysisManager->CreateNtupleDColumn("EdepAlphaFromAlphas");
  //analysisManager->CreateNtupleDColumn("EdepAlphaFromBetas");
  analysisManager->CreateNtupleDColumn("EdepBeta");
  analysisManager->CreateNtupleDColumn("depth");
  analysisManager->FinishNtuple();

  analysisManager->OpenFile("temp.root"); 

 G4int evts_to_process = aRun->GetNumberOfEventToBeProcessed();
    G4RunManager::GetRunManager()->SetPrintProgress((evts_to_process > 100)
                                                    ? evts_to_process/100
                                                    : 1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LiRunAction::EndOfRunAction(const G4Run*)
{

  // closing and renaming ROOT file at the end of run
  //
  auto  analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile(); 

  if(isMaster)
  {
  int result;

  const LiDetectorConstruction* detectorConstruction
      = static_cast<const LiDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  filename= detectorConstruction->GetOutputFileName();

  std::string actualName=filename;
  std::replace(actualName.begin(), actualName.end(), '.', '_');
  actualName+=".root";
  char oldname[] ="temp.root";
  result= rename( oldname , actualName.c_str());
  if ( result == 0 )
    {
      puts ( "File successfully renamed to " );
      puts(actualName.c_str());
    }		
  else
    perror( "Error renaming file" );
  }
  delete G4AnalysisManager::Instance(); 
}


