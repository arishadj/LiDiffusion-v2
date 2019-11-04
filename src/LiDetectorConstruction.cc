#include "LiDetectorConstruction.hh"
#include "LiDetectorMessenger.hh"

#include "G4Material.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4NistMaterialBuilder.hh"
#include "G4NistElementBuilder.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <iterator>  
#include <cmath> 
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>

#define pi   3.14

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LiDetectorConstruction::LiDetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0),fScoringVolumeBeta(0),fOut(0),fTargetVolume(0),
  detector_mat(0),beta_detector_mat(0),world_mat(0),sample_mat(0),
  substrate_mat(0),Shielding_mat(0),logicWorld(0),logicEnv(0),
  logicShapeSample(0),logicShapeSubstrate(0),A_logicRing(0),
  S_logic(0),B_logicDet(0), solidWorld(0),physWorld(0),
  solidEnv(0),physEnv(0),sample(0),physSample(0),substrate(0),
  physSubstrate(0),A_solidRing(0),physAlphaDet(0),Shielding_solid(0),
  physShielding(0),BetaDet_solid(0),physBetaDet(0),
  SampleVisAtt(0),SubstrateVisAtt(0),DetectorVisAtt(0)
{ 
 LiDetMessenger = new LiDetectorMessenger(this);

 //dummy initialization. The actual values are defined by the user/Messenger//
 SampleWidth    =    7*mm;  
 SampleDepth    =  500*um;      
 SubstrateWidth =    7*mm;        
 SubstrateDepth =  500*um;     
 A_cryst_dZ     =    0.4*mm;
 A_cryst_dY     =  0.1*mm;
 A_ring_R1      =   1.1*5.2*mm;

 MakeMaterials();
 SetSampleMaterial("TiO2");
 SetSubstrateMaterial("Galactic");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LiDetectorConstruction::~LiDetectorConstruction()
{
  delete LiDetMessenger;
}



void LiDetectorConstruction::MakeMaterials() 
{
  // Get nist material manager
  G4NistManager* man = G4NistManager::Instance();
 
  //This function illustrates the possible ways to define materials

  G4String symbol;        //a=mass of a mole;
  G4double a, z;          //z=mean number of protons;
  G4double density;//, temperature, pressure;

  G4int ncomponents, natoms;
  G4double fractionmass;          //abundance,

  //
  // define Elements
  //

  //G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
  G4Element* Ti  = man->FindOrBuildElement("Ti");
  G4Element* O  = man->FindOrBuildElement("O");
  G4Element* Zn = man->FindOrBuildElement("Zn");
  G4Element* S = man->FindOrBuildElement("S");
  G4Element* Ag = man->FindOrBuildElement("Ag");
  //G4Element* Zn = new G4Element("Zinc",symbol="Zn" , z= 30., a= 65.38*g/mole);
  //G4Element* S = new G4Element("Sulfur",symbol="S" , z= 16., a= 32.06*g/mole);
  //G4Element* Ag = new G4Element("Silver",symbol="Ag" , z= 47., a= 107.87*g/mole);

  //
  // define sample materials
  //
  G4Material* Vacuum =
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= 0.000000000001*mg/cm3,
                           kStateGas, 3.e-18*pascal, 2.73*kelvin);

  G4Material* TiO2 =new G4Material("TiO2",density= 4.23*g/cm3, ncomponents=2);
  TiO2->AddElement(Ti, natoms=1);
  TiO2->AddElement(O , natoms=2);

  G4Material* Al2O3 = man->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");

  G4Material* ZnS_Ag =new G4Material("ZnS_Ag",density= 4.090*g/cm3, ncomponents=3);
  ZnS_Ag->AddElement(Zn, fractionmass=0.449);
  ZnS_Ag->AddElement(S, fractionmass=0.498);
  ZnS_Ag->AddElement(Ag, fractionmass=0.053);

  // end of own definitions

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  world_mat         =  Vacuum;  
  detector_mat      =  ZnS_Ag ;
 // sample_mat      =  TiO2 ;
  //substrate_mat     =  Vacuum;//Al2O3 ;
  Shielding_mat     =  man->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  beta_detector_mat =  ZnS_Ag ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LiDetectorConstruction::DestroyMaterials()
{
  // Destroy all allocated elements and materials
  size_t i;
  G4MaterialTable* matTable = (G4MaterialTable*)G4Material::GetMaterialTable();
  for(i=0;i<matTable->size();i++)
  { delete (*(matTable))[i]; }
  matTable->clear();
  G4ElementTable* elemTable = (G4ElementTable*)G4Element::GetElementTable();
  for(i=0;i<elemTable->size();i++)
  { delete (*(elemTable))[i]; }
  elemTable->clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LiDetectorConstruction::Update() 
{
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  Construct();
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  CheckDepthDistribution();
  G4RunManager::GetRunManager()->ReinitializeGeometry();

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4VPhysicalVolume* LiDetectorConstruction::Construct()
{  

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  G4double beta_depth    =  10*mm;
  G4double beta_radious  =  5*mm;
  G4double beta_R1       =  0*mm;
  G4double beta_pos      =  10*cm; 
 
  G4double env_sizeXY    =  2*(A_ring_R1+ A_cryst_dY+1*mm) ;
  G4double env_sizeZ     =  2*(beta_pos+beta_depth+1*cm);

  G4double world_sizeXY  =  env_sizeXY+1*mm;
  G4double world_sizeZ   =  env_sizeZ+1*mm;

   //     
  // World
  //

  solidWorld =    
    new G4Box("World",                                           //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                   
  //     
  // Envelope (made of vacuum for diagnostic reasons)
  //  


  solidEnv =    
    new G4Box("Envelope",                               //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        world_mat,           //its material
                        "Envelope");         //its name
   physEnv=            
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
  //     
  // Sample
  //  

  G4ThreeVector SamplePosition = G4ThreeVector(0., 0., SampleDepth/2);
              
  sample =    
   new G4Box("sample", SampleWidth/2,SampleWidth/2,SampleDepth/2);
                      
  logicShapeSample =                         
    new G4LogicalVolume(sample,              //its solid
                        sample_mat,          //its material
                        "sample");           //its name
    
  physSample=         
  new G4PVPlacement(0,                       //no rotation
                    SamplePosition,          //at position
                    logicShapeSample,        //its logical volume
                    "sample",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //     
  // Substrate (optional, by default made of vacuum)
  //

  G4ThreeVector SubstratePos = G4ThreeVector(0., 0., (SampleDepth+SubstrateDepth/2));

  substrate =    
    new G4Box("substrate", SubstrateWidth/2,SubstrateWidth/2,SubstrateDepth/2);
       
  logicShapeSubstrate =                         
    new G4LogicalVolume(substrate,           //its solid
                        substrate_mat,       //its material
                        "substrate");        //its name
  physSubstrate=             
  new G4PVPlacement(0,                       //no rotation
                    SubstratePos,            //at position
                    logicShapeSubstrate,     //its logical volume
                    "substrate",             //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //     
  // Alpha ring detector
  //
 
  G4double twopi = 2*pi;
  G4double A_ring_R2 = (A_ring_R1+ A_cryst_dY);
      
 A_solidRing =
   new G4Cons("Alpha_Ring",0.9*A_ring_R1,0.9*A_ring_R2, A_ring_R1, A_ring_R2, 0.5*A_cryst_dZ, 0, twopi);

  A_logicRing =                         
    new G4LogicalVolume(A_solidRing,             //its solid
                        detector_mat,            //its material
                        "Alpha_Ring");           //its name
                    
  physAlphaDet=
    new G4PVPlacement(0,                                        //no rotation
                    G4ThreeVector(0., 0., -A_cryst_dZ/2*mm),    //at position
                    A_logicRing,                                //its logical volume
                    "Alpha_detector",                           //its name
                    logicEnv ,                                  //its mother  volume
                    false,                                      //no boolean operation
                    0,                                          //copy number
                    checkOverlaps);                             //overlaps checking

       
  //
  // Shielding of beta detector
  //

  Shielding_solid =
    new G4Tubs("Beta_Shield", 0, beta_radious, 0.1*cm, 0., twopi);
      
  S_logic =                         
    new G4LogicalVolume(Shielding_solid,           //its solid
                        Shielding_mat,             //its material
                        "Beta_Shield");            //its name
                    

  physShielding=
    new G4PVPlacement(0,                                     //no rotation
		   G4ThreeVector(0., 0.,(beta_pos-1*cm) ),   //at position
                    S_logic,                                 //its logical volume
                    "Beta_Shield",                           //its name
                   logicEnv ,                                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //
  // Beta detector
  //

  BetaDet_solid =
    new G4Tubs("Beta_det", beta_R1, beta_radious, 0.5*beta_depth, 0., twopi);
      
  B_logicDet =                         
    new G4LogicalVolume(BetaDet_solid,           //its solid
                        beta_detector_mat,         //its material
                        "Beta_det");             //its name
                    

  physBetaDet=
    new G4PVPlacement(0,                             //no rotation
                    G4ThreeVector(0., 0.,beta_pos ), //at position
                    B_logicDet,                      //its logical volume
                    "Beta_detector",                 //its name
                   logicEnv ,                        //its mother  volume
                    false,                           //no boolean operation
                    0,                               //copy number
                    checkOverlaps);                  //overlaps checking

  //visualization attributes
  logicEnv->SetVisAttributes(G4VisAttributes::Invisible);

  SampleVisAtt    = new G4VisAttributes(G4Colour(0.,0.2,0.5));
  SampleVisAtt->SetVisibility(true);
  logicShapeSample->SetVisAttributes(SampleVisAtt);

  SubstrateVisAtt = new G4VisAttributes(G4Colour(0.5,0.5,1.0));
  SubstrateVisAtt->SetVisibility(true);
  logicShapeSubstrate->SetVisAttributes(SubstrateVisAtt);

  DetectorVisAtt  = new G4VisAttributes(G4Colour(0.5,0.3,0.8));
  DetectorVisAtt->SetVisibility(true);
  A_logicRing->SetVisAttributes(DetectorVisAtt);  
  B_logicDet->SetVisAttributes(DetectorVisAtt);  

  fScoringVolume =  A_logicRing; 
  fScoringVolumeBeta =  B_logicDet; 
  fTargetVolume = logicShapeSample;   
  fOut = logicWorld;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

void LiDetectorConstruction::SetSampleMaterial(G4String mat)
{
  // search the material by its name
  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(mat);

  if (material && material != sample_mat) {
    sample_mat = material;
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

void LiDetectorConstruction::CreateSampleMaterial(std::vector<std::string> atoms)
{
  //user-defined sample material
  G4double d =atof(atoms[0].c_str());
  G4int n = atoms.size()-1;

  G4Material* custom =new G4Material("customSampleMaterial",d*g/cm3, n/2);
  for(G4int i=1;i<n;i=i+2)
  { 
    G4int nAtom=atoi(atoms[i+1].c_str());
    G4String elemName = atoms[i];
    G4Element* elem  = G4NistManager::Instance()->FindOrBuildElement(elemName.c_str());
    custom->AddElement(elem, nAtom);
  }

  if (custom && custom != sample_mat) {
    sample_mat = custom;
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

void LiDetectorConstruction::CreateSubstrateMaterial(std::vector<std::string> atoms)
{
  //user-defined substrate material
  G4double d =atof(atoms[0].c_str());
  G4int n = atoms.size()-1;

  G4Material* custom =new G4Material("customSubstrateMaterial",d*g/cm3, n/2);
  for(G4int i=1;i<n;i=i+2)
  { 
    G4int nAtom=atoi(atoms[i+1].c_str());
    G4String elemName = atoms[i];
    G4Element* elem  = G4NistManager::Instance()->FindOrBuildElement(elemName.c_str());
    custom->AddElement(elem, nAtom);
  }

  if (custom && custom != substrate_mat) {
    substrate_mat = custom;
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

void LiDetectorConstruction::SetSubstrateMaterial(G4String mat)
{
  // search the material by its name
  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(mat);

  if (material && material != substrate_mat) {
    substrate_mat = material;
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

void LiDetectorConstruction::CheckDepthDistribution()
{
  //check if the depth profile with the defined parameters (diff. coeff.,time, SRIM profile etc) already exists
  std::vector<G4double>  	  C0;   
  maxDepth=SampleDepth/nm;
  G4String pathAndSRIMname="DiffusionProfiles/SRIM/"+SRIMfilename+".txt";

  std::ostringstream strs;

  strs <<filename<<"D"<<D<<"t"<<time/s<<condition;
  nameOut=strs.str();

  std::ostringstream buffer;
  buffer <<"DiffusionProfiles/"<<SRIMfilename<<"/D"<<D<<"t"<<time/s<<condition;
  G4String profileName=buffer.str();
  std::replace(profileName.begin(), profileName.end(), '.', '_');
  profileName+=".txt";


  G4String Dirname="DiffusionProfiles/"+SRIMfilename;
  G4bool fileExists=false;
  DIR* dir = opendir(Dirname.c_str());
  if (opendir(Dirname.c_str()))
  {
    /* Directory exists. */
    if (std::ifstream(profileName.c_str()))
    {
        /* file exists. */
      closedir(dir);
      G4cout<<"Diffusion profile exists already."<<G4endl;
      fileExists=true;
    }
  }
  else
  {
    G4String commandMkDir="mkdir "+Dirname;
    int mkdir=system(commandMkDir);

   if(mkdir==0) G4cout<<"Could not create dir: "<<Dirname<<G4endl;
  }

  if(fileExists) 
  {
    C0=ReadDiffusionProfile(profileName);
    NpointsZ=static_cast<G4int>(C0.size());
    SetCoefficients(C0);
  }
  //if it doesn't exist, calculate it
  else
  {
    MakeDepthDistribution(pathAndSRIMname, profileName);
  }

}

void LiDetectorConstruction::MakeDepthDistribution(G4String pathAndSRIMname, G4String profileName)
{
  //calculates the new depth profile, given the defined parameters of the Li motion
  std::vector<G4double>  C0;   

  G4cout<<"Diffusion profile does not exist already, new profile will be created."<<G4endl;

  //C0 has the initial profile, Ct the depth profile at time t
  C0=ReadSRIMProfile(pathAndSRIMname);

  G4double Total=0;
  for (G4int i=0 ; i< static_cast<G4int>(C0.size()) ; i++) 
  {
    Total+=C0[i];
  }
  for(G4int i=0; i<static_cast<G4int>(C0.size());i++)
  {
    C0[i]=C0[i]/Total;
  }

  G4double rate=D*1e14;
  G4double lamda = 0.8249;
  G4double Dt=0.01;

  std::vector<G4double>  	  Ct;
  std::vector<G4double>  	  N;

  N.push_back(Libeam*Dt);
  G4double R=0;

  for(G4int i=0;i<static_cast<G4int>(C0.size());i++)
  {
    Ct.push_back(C0[i]*Libeam*Dt);
  }

  G4int upperTime=static_cast<G4int>((time/s)/Dt);

  for(G4int t_i=0;t_i<upperTime;t_i++)
  {
    G4double t=static_cast<G4double>(t_i)*Dt; 

    std::vector<G4double>  	  tempProfilePlus; 
    std::vector<G4double>  	  tempProfileMinus; 
    G4double  	  tempProfileZero=0; 
    tempProfilePlus.push_back(0);
    tempProfileMinus.push_back(0);

    R=N[t_i]-N[t_i]*exp(-Dt*lamda);

    for(G4int i=0;i<static_cast<G4int>(Ct.size());i++)
    {
      //remove the decays of previous Dt period
      Ct[i]=Ct[i]*exp(-Dt*lamda);
    }  

    if(condition=="accumulative")
    {   
      tempProfileZero=Ct[0];
      G4int d_i=0;

      do
      {
        tempProfilePlus.push_back(0);
	tempProfileMinus.push_back(0);         
        for(G4int j=1;j<static_cast<G4int>(Ct.size());j++)
	{
	  G4double distPlus=abs(j-d_i);
	  G4double distMinus=abs(j+d_i);
	  tempProfilePlus[d_i]+=(Ct[j]/(sqrt(4*pi*rate*Dt)))*exp(-pow(distPlus,2)/(4*rate*Dt));
	  tempProfileMinus[d_i]+=(Ct[j]/(sqrt(4*pi*rate*Dt)))*exp(-pow(distMinus,2)/(4*rate*Dt));
	}    
	d_i++;
      }while(tempProfilePlus[d_i-1]>1e-7);
    
      G4int newSize=std::max(tempProfilePlus.size(),tempProfileMinus.size()-tempProfilePlus.size());
      newSize=std::min(newSize,maxDepth);
      Ct.clear();
      Ct.resize(newSize);

      Ct[0]=tempProfileZero;

      for(G4int i=1;i<static_cast<G4int>(Ct.size());i++)
      {
        Ct[i]=tempProfilePlus[i];
      }  

      for(G4int i=0;i<static_cast<G4int>(tempProfileMinus.size());i++)
      {
	Ct[0]+=tempProfileMinus[i];
      }  

      if(static_cast<G4int>(tempProfilePlus.size())>maxDepth)
      {
	for(G4int i=maxDepth;i<=static_cast<G4int>(tempProfilePlus.size());i++)
	{
	  Ct[maxDepth-1]=Ct[maxDepth-1]+tempProfilePlus[i];
	}    
      }    
    }
    else if(condition=="reflective")
    {
      G4int d_i=0;

      do
      { 
	 tempProfilePlus.push_back(0);
	 tempProfileMinus.push_back(0);
	 for(G4int j=0;j<static_cast<G4int>(Ct.size());j++)
	 {
	   G4double distPlus=abs(j-d_i);
	   G4double distMinus=abs(j+d_i);
	   tempProfilePlus[d_i]+=(Ct[j]/(sqrt(4*pi*rate*Dt)))*exp(-pow(distPlus,2)/(4*rate*Dt));
	   tempProfileMinus[d_i]+=(Ct[j]/(sqrt(4*pi*rate*Dt)))*exp(-pow(distMinus,2)/(4*rate*Dt));
	 }    
	 d_i++;
       }while(tempProfilePlus[d_i-1]>1e-7);

       G4int newSize=static_cast<G4int>(std::min(static_cast<int>(tempProfilePlus.size()),static_cast<int>(maxDepth)));

       Ct.clear();
       Ct.resize(newSize);

       Ct[0]=tempProfilePlus[0];

       for(G4int i=1;i<static_cast<G4int>(Ct.size());i++)
       {
         Ct[i]=tempProfilePlus[i];
       }  

       for(G4int i=1;i<static_cast<G4int>(tempProfileMinus.size());i++)
       {
	 Ct[i]+=tempProfileMinus[i];
       }  

       if(static_cast<G4int>(tempProfilePlus.size())>maxDepth)
       {
         for(G4int i=maxDepth;i<=static_cast<G4int>(tempProfilePlus.size());i++)
	 {
	   Ct[maxDepth-1-i]=Ct[maxDepth-1-i]+tempProfilePlus[i];
	 }    
       }    
    } 
    else
    {
      G4cout<<"ATTENTION!!!!!! BOUNDARY CONDITION NOT RECOGNIZED, RESULTS MAY NOT BE VALID"<<G4endl;
    }

    if(t<beamOn)
    {
      for(G4int j=0;j<static_cast<G4int>(C0.size());j++)
      {
        Ct[j]=Ct[j]+C0[j]*Libeam*Dt;
      }
    }
    G4double tot=0;
   
    if(t<beamOn)
    {
      N.push_back(N[t_i]+Libeam*Dt-R);
    }
    else
    {
      N.push_back(N[t_i]-R);
    }
 
    for(G4int j=0;j<static_cast<G4int>(Ct.size());j++)
    {
      tot+=Ct[j];
    }

    for(int j=0;j<static_cast<G4int>(Ct.size());j++)
    {
      Ct[j]=Ct[j]*N[t_i]/tot;
    }
  }   
	
  std::ofstream myfile ;
  G4bool lineNbIsZero=true;

  myfile.open(profileName.c_str()) ;

  if (myfile.is_open())
  {
    if(lineNbIsZero)
    {
      myfile<<time/s <<"\n";
      myfile<<"#Depth (nm)   N/nm"<<"\n";
      lineNbIsZero=false;
    }
    for(G4int j=0;j<static_cast<G4int>(Ct.size());j++)
    {
      myfile <<j+0.5<<"            "<<Ct[j]<<"\n";
    }
  }

  myfile.close();
  if (myfile.is_open()) G4cout<<"NOT CLOSED"<<G4endl;

  NpointsZ=static_cast<G4int>(Ct.size());
  SetCoefficients(Ct);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

std::vector<G4double> LiDetectorConstruction::ReadDiffusionProfile(G4String profileName)
{
  //read the (already existing) diffusion profile
  std::vector<G4double>  	  C0;   
  G4double dummy1=0;
  G4double dummy2=0;
  std::ifstream file;
  file.open(profileName) ;
  if (file.is_open())
  {
    std::string line;
    std::istringstream inStream; 
    getline(file, line) ;
    inStream.clear();
    inStream.str(line); 
    getline(file, line) ;
    inStream.clear();
    inStream.str(line); 
    while (!file.eof())
    {
      getline(file, line) ;
      inStream.clear();
      inStream.str(line); 
      inStream >> dummy1>>dummy2;
      if(!file.eof()) 
      {  
	C0.push_back(dummy2);
      }
    }
  }
  file.close();

  return C0;
}

std::vector<G4double> LiDetectorConstruction::ReadSRIMProfile(G4String profileName)
{
  //read the SRIM profile
  std::vector<G4double>  	  C0;   
  G4double dummy1=0;
  G4double dummy2=0;
  std::ifstream file;
  file.open(profileName) ;
  if (file.is_open())
  {
    std::string line;
    std::istringstream inStream; 
    getline(file, line) ;
    inStream.clear();
    inStream.str(line); 
    while (!file.eof())
    {
      getline(file, line) ;
      inStream.clear();
      inStream.str(line); 
      inStream >> dummy1>>dummy2;
      if(!file.eof()) 
      {  
	C0.push_back(dummy2);
      }
    }
  }
  file.close();

  return C0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

void LiDetectorConstruction::SetCoefficients(std::vector<G4double> Ct)
{
  //copy arrays in std::vector and compute fMax
  //
  xZ.resize(NpointsZ); fZ.resize(NpointsZ); 
  fZ=Ct;
  G4double sum=0;

  for (G4int j=0; j<NpointsZ; j++) {
    sum=sum+fZ[j];
  }

  for (G4int j=0; j<NpointsZ; j++) {
    fZ[j]=fZ[j]/sum;
  }   

  for (G4int j=0; j<NpointsZ; j++) {
    xZ[j] = j; 
    if (fZ[j]<0) fZ[j]=0;	
  };
}
