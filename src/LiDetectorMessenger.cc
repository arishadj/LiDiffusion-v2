#include "LiDetectorMessenger.hh"
#include "LiDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"

#include <string>
#include <sstream>
#include <vector>
#include <iterator>

using namespace std;

LiDetectorMessenger::LiDetectorMessenger(LiDetectorConstruction * Det)
  :G4UImessenger(),fDetector(Det),LiDir(0),DetDir(0),
   SampleMaterialCmd(0),CreateSampleMaterialCmd(0), 
   SubstrateMaterialCmd(0),CreateSubstrateMaterialCmd(0),
   SampleWidthCmd(0),SampleDepthCmd(0),SubstrateWidthCmd(0),
   SubstrateDepthCmd(0),DetWidthCmd(0),DetDepthCmd(0),
   DetRadiousCmd(0),UpdateCmd(0),TimeCmd(0),DiagnosticsCmd(0),
   FileNameCmd(0),SRIMCmd(0),DiffRateCmd(0),BoundaryCondCmd(0),
   NLiCmd(0), BeamOnCmd(0)
 { 

  //Defining the commands for user-specified input related to the simulation's geometry/materials
  LiDir = new G4UIdirectory("/LiDiffusionCode/");
  LiDir->SetGuidance("UI commands of this example");

  DetDir = new G4UIdirectory("/LiDiffusionCode/geometry/");
  DetDir->SetGuidance("geometry settings");
   
  SampleMaterialCmd = new G4UIcmdWithAString("/LiDiffusionCode/geometry/setSampleMaterial",this);
  SampleMaterialCmd->SetGuidance("Select material of the sample.");
  SampleMaterialCmd->SetParameterName("choice",false);
  SampleMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SampleMaterialCmd->SetToBeBroadcasted(false);

  CreateSampleMaterialCmd = new G4UIcmdWithAString("/LiDiffusionCode/geometry/createSampleMaterial",this);
  CreateSampleMaterialCmd->SetGuidance("Create material of the sample.");
  CreateSampleMaterialCmd->SetGuidance("Set its density (in g/cm3) followed by each atom of the unit cell with its stoichiometric number.");
  CreateSampleMaterialCmd->SetGuidance("For example, to create rutile TiO2 write:");
  CreateSampleMaterialCmd->SetGuidance("4.23 Ti 1 O 2");
  CreateSampleMaterialCmd->SetParameterName("choice",false);
  CreateSampleMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateSampleMaterialCmd->SetToBeBroadcasted(false);

  SubstrateMaterialCmd = new G4UIcmdWithAString("/LiDiffusionCode/geometry/setSubstrateMaterial",this);
  SubstrateMaterialCmd->SetGuidance("Select material of the substrate.");
  SubstrateMaterialCmd->SetParameterName("choice",false);
  SubstrateMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SubstrateMaterialCmd->SetToBeBroadcasted(false);

  CreateSubstrateMaterialCmd = new G4UIcmdWithAString("/LiDiffusionCode/geometry/createSubstrateMaterial",this);
  CreateSubstrateMaterialCmd->SetGuidance("Create material of the substrate.");
  CreateSubstrateMaterialCmd->SetGuidance("Set its density (in g/cm3) followed by each atom of the unit cell with its stoichiometric number.");
  CreateSubstrateMaterialCmd->SetGuidance("For example, to create rutile TiO2 write:");
  CreateSubstrateMaterialCmd->SetGuidance("4.23 Ti-1 O-2");
  CreateSubstrateMaterialCmd->SetParameterName("choice",false);
  CreateSubstrateMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateSubstrateMaterialCmd->SetToBeBroadcasted(false);

  SampleWidthCmd = new G4UIcmdWithADoubleAndUnit("/LiDiffusionCode/geometry/setSampleWidth",this);
  SampleWidthCmd->SetGuidance("Set the width of the sample in mm");
  SampleWidthCmd->SetParameterName("sampleWidth",false);
  SampleWidthCmd->SetUnitCategory("Length");
  SampleWidthCmd->SetDefaultValue(8.);
  SampleWidthCmd->SetRange("sampleWidth>=0. && sampleWidth<=15.");
  SampleWidthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SampleWidthCmd->SetToBeBroadcasted(false);

  SampleDepthCmd = new G4UIcmdWithADoubleAndUnit("/LiDiffusionCode/geometry/setSampleDepth",this);	
  SampleDepthCmd->SetGuidance("Set the depth of the sample in um");
  SampleDepthCmd->SetParameterName("sampleDepth",false);
  SampleDepthCmd->SetUnitCategory("Length");
  SampleDepthCmd->SetDefaultValue(500.);
  SampleDepthCmd->SetRange("sampleDepth>=0. && sampleDepth<=3000.");
  SampleDepthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SampleDepthCmd->SetToBeBroadcasted(false);

  SubstrateWidthCmd = new G4UIcmdWithADoubleAndUnit("/LiDiffusionCode/geometry/setSubstrateWidth",this);	
  SubstrateWidthCmd->SetGuidance("Set the width of the substrate in mm");
  SubstrateWidthCmd->SetParameterName("substrateWidth",false);
  SubstrateWidthCmd->SetUnitCategory("Length");
  SubstrateWidthCmd->SetDefaultValue(8.);
  SubstrateWidthCmd->SetRange("substrateWidth>=0. && substrateWidth<=15.");
  SubstrateWidthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SubstrateWidthCmd->SetToBeBroadcasted(false);

  SubstrateDepthCmd = new G4UIcmdWithADoubleAndUnit("/LiDiffusionCode/geometry/setSubstrateDepth",this);	 
  SubstrateDepthCmd->SetGuidance("Set the depth of the substrate in um");
  SubstrateDepthCmd->SetParameterName("substrateDepth",false);
  SubstrateDepthCmd->SetUnitCategory("Length");
  SubstrateDepthCmd->SetDefaultValue(500.);
  SubstrateDepthCmd->SetRange("substrateDepth>=0. && substrateDepth<=3000.");
  SubstrateDepthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SubstrateDepthCmd->SetToBeBroadcasted(false);

  DetWidthCmd = new G4UIcmdWithADoubleAndUnit("/LiDiffusionCode/geometry/setDetectorWidth",this);	
  DetWidthCmd->SetGuidance("Set the width of the alpha detector");
  DetWidthCmd->SetParameterName("detWidth",false);
  DetWidthCmd->SetUnitCategory("Length");
  DetWidthCmd->SetDefaultValue(1.);
  DetWidthCmd->SetRange("detWidth>=0. && detWidth<=1500.");
  DetWidthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  DetWidthCmd->SetToBeBroadcasted(false);

  DetDepthCmd = new G4UIcmdWithADoubleAndUnit("/LiDiffusionCode/geometry/setDetectorDepth",this);	
  DetDepthCmd->SetGuidance("Set the depth of the detector in mm");
  DetDepthCmd->SetParameterName("detDepth",false);
  DetDepthCmd->SetUnitCategory("Length");
  DetDepthCmd->SetDefaultValue(0.1);
  DetDepthCmd->SetRange("detDepth>=0. && detDepth<=30.");
  DetDepthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  DetDepthCmd->SetToBeBroadcasted(false);

  DetRadiousCmd = new G4UIcmdWithADoubleAndUnit("/LiDiffusionCode/geometry/setDetectorRadious",this);	
  DetRadiousCmd->SetGuidance("Set the radious of the detector in mm");
  DetRadiousCmd->SetParameterName("detR",false);
  DetRadiousCmd->SetUnitCategory("Length");
  DetRadiousCmd->SetDefaultValue(10.);
  DetRadiousCmd->SetRange("detR>=0. && detR<=200.");
  DetRadiousCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  DetRadiousCmd->SetToBeBroadcasted(false);

  UpdateCmd =new G4UIcmdWithoutParameter("/LiDiffusionCode/geometry/update",this);
  UpdateCmd->SetGuidance("Update geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
  UpdateCmd->SetToBeBroadcasted(false);

  //Miscellenius commands (not related to the geometry).
  TimeCmd = new G4UIcmdWithADoubleAndUnit("/LiDiffusionCode/geometry/setTime",this);	
  TimeCmd->SetGuidance("Set time of study after beam implantation");
  TimeCmd->SetParameterName("time",true);
  TimeCmd->SetUnitCategory("Time");
  TimeCmd->SetDefaultValue(0.1);
  TimeCmd->SetRange("time>=0. && time<=6.");
  TimeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  TimeCmd->SetToBeBroadcasted(false);

  DiagnosticsCmd = new G4UIcmdWithAnInteger("/LiDiffusionCode/geometry/diagnostics",this);
  DiagnosticsCmd->SetGuidance("0 for diagnostics off, 1 for diagnostics on");
  DiagnosticsCmd->SetParameterName("diagnostics",true);
  DiagnosticsCmd->SetDefaultValue(0);
  DiagnosticsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  DiagnosticsCmd->SetToBeBroadcasted(false);

  FileNameCmd = new G4UIcmdWithAString("/LiDiffusionCode/geometry/setFilename",this);	
  FileNameCmd->SetGuidance("Choose output file name");
  FileNameCmd->SetParameterName("filename",true);
  FileNameCmd->SetDefaultValue("sample");
  FileNameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  FileNameCmd->SetToBeBroadcasted(false);

  SRIMCmd = new G4UIcmdWithAString("/LiDiffusionCode/geometry/setSRIMfilename",this);	
  SRIMCmd->SetGuidance("Choose input SRIM file name");
  SRIMCmd->SetParameterName("SRIM",true);
  SRIMCmd->SetDefaultValue("rutile5keV");
  SRIMCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SRIMCmd->SetToBeBroadcasted(false);

  DiffRateCmd = new G4UIcmdWithADouble("/LiDiffusionCode/geometry/setDiffusionRate",this);	//time 
  DiffRateCmd->SetGuidance("Set the diffusion rate (e.g. 1e-11) in cm^2/s");
  DiffRateCmd->SetParameterName("diffR",true);
  DiffRateCmd->SetDefaultValue(1e-11);
  DiffRateCmd->SetRange("diffR>=0. && diffR<=1");
  DiffRateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  DiffRateCmd->SetToBeBroadcasted(false);

  BoundaryCondCmd = new G4UIcmdWithAString("/LiDiffusionCode/geometry/setBoundaryCondition",this);	
  BoundaryCondCmd->SetGuidance("Choose 'reflective' or 'accumulative' boundary condition in the sample");
  BoundaryCondCmd->SetParameterName("boundaryCond",true);
  BoundaryCondCmd->SetDefaultValue("reflective");
  BoundaryCondCmd->SetCandidates("reflective accumulative");
  BoundaryCondCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  BoundaryCondCmd->SetToBeBroadcasted(false);

  NLiCmd = new G4UIcmdWithAnInteger("/LiDiffusionCode/geometry/LiIonsPerPulse",this);
  NLiCmd->SetGuidance("(integer) Number of Li ions per beam pulse");
  NLiCmd->SetParameterName("NLi",true);
  NLiCmd->SetDefaultValue(1000000);
  NLiCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  NLiCmd->SetToBeBroadcasted(false);

  BeamOnCmd = new G4UIcmdWithADouble("/LiDiffusionCode/geometry/setBeamOnTime",this);
  BeamOnCmd->SetGuidance("Set the beam On time in seconds");
  BeamOnCmd->SetParameterName("beamOn",true);
  BeamOnCmd->SetDefaultValue(1.);
  BeamOnCmd->SetRange("beamOn>=0. && beamOn<=100.");
  BeamOnCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  BeamOnCmd->SetToBeBroadcasted(false);
 }

LiDetectorMessenger::~LiDetectorMessenger()
{
  delete LiDir;
  delete DetDir;
  delete SampleMaterialCmd;
  delete CreateSampleMaterialCmd;
  delete  SubstrateMaterialCmd;
  delete CreateSubstrateMaterialCmd;
  delete SampleWidthCmd;
  delete SampleDepthCmd;
  delete SubstrateWidthCmd;
  delete SubstrateDepthCmd;
  delete DetWidthCmd;
  delete DetDepthCmd;
  delete DetRadiousCmd;
  delete UpdateCmd;

  delete TimeCmd;
  delete FileNameCmd;
  delete DiagnosticsCmd; 
  delete DiffRateCmd;
  delete BoundaryCondCmd;
  delete SRIMCmd;
  delete NLiCmd;
  delete BeamOnCmd;
}
   

void LiDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  //passing the commands' values to the DetectorConstruction class
  if( command == SampleMaterialCmd )
  { fDetector->SetSampleMaterial(newValue);}

  if( command == SubstrateMaterialCmd )
  { fDetector->SetSubstrateMaterial(newValue);}

  if( command == CreateSampleMaterialCmd )
  { 
	std::vector<std::string> elems;
        std::stringstream ss;
        ss.str(newValue);
        std::string item;
        while (std::getline(ss, item,' ')) {
           elems.push_back(item);
        }
	fDetector->CreateSampleMaterial(elems);
  }

    if( command == CreateSubstrateMaterialCmd )
   { 
	std::vector<std::string> elems;
        std::stringstream ss;
        ss.str(newValue);
        std::string item;
        while (std::getline(ss, item,' ')) {
           elems.push_back(item);
        }
	fDetector->CreateSubstrateMaterial(elems);
}

  if(command == SampleWidthCmd)
  {fDetector->SetSampleWidth(SampleWidthCmd->GetNewDoubleValue(newValue)) ;}

  if(command == SampleDepthCmd)
    {fDetector->SetSampleDepth(SampleDepthCmd->GetNewDoubleValue(newValue));}

  if(command == SubstrateWidthCmd)
    {fDetector->SetSubstrateWidth(SubstrateWidthCmd->GetNewDoubleValue(newValue));}

  if(command == SubstrateDepthCmd)
    {fDetector->SetSubstrateDepth(SubstrateDepthCmd->GetNewDoubleValue(newValue));}

  if(command == DetWidthCmd)
    {fDetector->SetDetectorWidth(DetWidthCmd->GetNewDoubleValue(newValue));}

  if(command == DetDepthCmd)
    {fDetector->SetDetectorDepth(DetDepthCmd->GetNewDoubleValue(newValue));}

  if(command == DetRadiousCmd)
    {fDetector->SetDetectorRadious(DetRadiousCmd->GetNewDoubleValue(newValue));}

  if(command == UpdateCmd)
    { fDetector->Update();}

  if(command == TimeCmd)
    {fDetector->SetTime(TimeCmd->GetNewDoubleValue(newValue));}
       
  if( command == FileNameCmd )
    { fDetector->SetFilename(newValue);}

  if( command == DiagnosticsCmd )
    { fDetector->SetDiagnostics(DiagnosticsCmd->GetNewIntValue(newValue));}

  if(command == DiffRateCmd)
    {fDetector->SetDiffRate(DiffRateCmd->GetNewDoubleValue(newValue));}

  if(command == BoundaryCondCmd)
    {fDetector->SetBoundaryCond(newValue);}

  if(command == SRIMCmd)
    {fDetector->SetSRIM(newValue);}

  if( command == NLiCmd )
    { fDetector->SetNLi(NLiCmd->GetNewIntValue(newValue));}

  if(command == BeamOnCmd)
    {fDetector->SetBeamOn(BeamOnCmd->GetNewDoubleValue(newValue));}
}
