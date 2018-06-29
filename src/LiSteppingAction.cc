#include "LiSteppingAction.hh"
#include "LiEventAction.hh"
#include "LiDetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4TransportationManager.hh" 
#include "G4Navigator.hh" 
#include "G4ParticleDefinition.hh" 
#include "G4DecayTable.hh" 
#include "G4VDecayChannel.hh" 
#include "G4PhaseSpaceDecayChannel.hh" 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LiSteppingAction::LiSteppingAction(LiEventAction* eventAction)
: G4UserSteppingAction(),fEventAction(eventAction),fScoringVolume(0),
  fScoringVolumeBeta(0),fTargetVolume(0),fOut(0),diagnostics(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LiSteppingAction::~LiSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LiSteppingAction::UserSteppingAction(const G4Step* step)
{
  diagnostics=fEventAction->GetDiagnostics();

    const LiDetectorConstruction* detectorConstruction
      = static_cast<const LiDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  fScoringVolume = detectorConstruction->GetScoringVolume();   
  fScoringVolumeBeta = detectorConstruction->GetScoringVolumeBeta(); 
  fTargetVolume = detectorConstruction->GetTargetVolume(); 

  // get volume of the current step
  G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

    // get the decay depth
    if((step->GetTrack()->GetParentID()==0)&&(step->GetTrack()->GetCurrentStepNumber()==1))
    {
       G4ThreeVector pos = step->GetTrack()->GetVertexPosition() ;
       G4double z = pos[2] ;
       fEventAction->SetDepth(z/nm);
    }

  // Optional diagnostical outputs
  if(diagnostics==1)
  {
    if(step->GetTrack()->GetCurrentStepNumber()==1)
    {
      if((step->GetTrack()->GetDefinition()->GetParticleName()=="alpha"))                 
      {
	fEventAction->AddAlphaEnergy(step->GetTrack()->GetVertexKineticEnergy()/keV ) ;
      }
	
      if((step->GetTrack()->GetDefinition()->GetParticleName()=="e-")&&(step->GetTrack()->GetParentID()==0))
      {
	fEventAction->AddBetaEnergy(step->GetTrack()->GetVertexKineticEnergy()/keV ) ;	
      }
    }
  }

  // check if we are in scoring volume
  if ((volume != fScoringVolume)&&(volume != fScoringVolumeBeta)) return;

  // Energy deposition at Alpha Detector
  if(volume == fScoringVolume)
  {
     //register the total energy deposited
     fEventAction->AddEdepAlpha(step->GetTotalEnergyDeposit()/keV); 	

     //(for diagnostics) register the energy deposited from alphas alone and then from everything else
      if((step->GetTrack()->GetParentID()==1)||((step->GetTrack()->GetDefinition()->GetParticleName()=="e-")&&(step->GetTrack()->GetParentID()==0)))
      {
       fEventAction->AddEdepAlphaFromBetas(step->GetTotalEnergyDeposit()/keV); 	
      }
      else
      {
        fEventAction->AddEdepAlphaFromAlphas(step->GetTotalEnergyDeposit()/keV); 
      }	
  } 

  // Energy deposition at Beta Detector
  if(volume == fScoringVolumeBeta)
  {
     fEventAction->AddEdepBeta(step->GetTotalEnergyDeposit()/keV) ; 		     
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

