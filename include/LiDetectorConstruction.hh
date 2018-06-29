#ifndef LiDetectorConstruction_h
#define LiDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class LiDetectorMessenger ;
class G4Box;
class G4Tubs;
class G4CutTubs;
class G4Cons;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4NistMaterialBuilder;
class G4VisAttributes;

/// Detector construction class to define materials and geometry.
class LiDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    LiDetectorConstruction();
    virtual ~LiDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    virtual void Update() ;
    virtual void MakeMaterials();
    virtual void DestroyMaterials();

    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
    G4LogicalVolume* GetScoringVolumeBeta() const { return fScoringVolumeBeta; }
    G4LogicalVolume* GetTargetVolume() const { return fTargetVolume; }
    G4LogicalVolume* GetOutVolume() const { return fOut; }
    G4double GetSampleDepth() const { return SampleDepth; }

    virtual void SetSampleWidth(G4double w){SampleWidth =w;};
    virtual void SetSampleDepth(G4double d){ SampleDepth=d;};
    virtual void SetSubstrateWidth(G4double w){SubstrateWidth =w;};
    virtual void SetSubstrateDepth(G4double d){ SubstrateDepth=d;};
    virtual void SetDetectorWidth(G4double w){A_cryst_dZ =w;};
    virtual void SetDetectorDepth(G4double d){A_cryst_dY =d;};
    virtual void SetDetectorRadious(G4double r){A_ring_R1 =r;};

    virtual void SetSampleMaterial(G4String);
    virtual void CreateSampleMaterial(std::vector<std::string>);
    virtual void CreateSubstrateMaterial(std::vector<std::string>);
    virtual void SetSubstrateMaterial(G4String);
    virtual void CheckDepthDistribution();
    virtual void MakeDepthDistribution(G4String, G4String);
    virtual void SetCoefficients(std::vector<G4double> C0);
    virtual std::vector<G4double> ReadDiffusionProfile(G4String profileName);
    virtual std::vector<G4double> ReadSRIMProfile(G4String profileName);

    void SetFilename(G4String f) {filename=f ;}
    void SetSRIM(G4String f) {SRIMfilename=f ;}
    void SetTime(G4double t) {time=t ;}
    void SetDiagnostics(G4int d) {diagnostics=d ;}
    void SetDiffRate(G4double r) {D=r ;}
    void SetBoundaryCond(G4String b) {condition=b ;}
    void SetBeamOn(G4double t) {beamOn=t ;}
    void SetNLi(G4int d) {Libeam=d ;}
    G4String GetFilename() const {return filename ;}
    G4String GetSRIM() const {return SRIMfilename ;}
    G4double GetDiffRate() const {return D;}
    G4String GetBoundaryCond() const {return condition ;}
    G4int    GetDiagnostics() const {return diagnostics;}
    G4double GetTime() const {return time ;}
    G4String GetOutputFileName() const {return nameOut;}
    G4int GetNLi() {return Libeam;}

    void SetXz(std::vector<G4double> x) { xZ=x;}
    void SetFz(std::vector<G4double> x) { fZ=x;}

    std::vector<G4double> GetXz() const {return xZ;}
    std::vector<G4double> GetFz() const {return fZ;}
    G4int		  GetNpointsZ() const {return NpointsZ;}

  protected:
    LiDetectorMessenger* LiDetMessenger;
    G4LogicalVolume*   fScoringVolume;
    G4LogicalVolume*   fScoringVolumeBeta;
    G4LogicalVolume*   fOut;
    G4LogicalVolume*   fTargetVolume;
 
    G4Material*        detector_mat  ;
    G4Material*        beta_detector_mat ;
    G4Material*        world_mat;
    G4Material*        sample_mat ;
    G4Material*        substrate_mat ;
    G4Material*        Shielding_mat;

    G4LogicalVolume*   logicWorld;
    G4LogicalVolume*   logicEnv;
    G4LogicalVolume*   logicShapeSample;
    G4LogicalVolume*   logicShapeSubstrate;
    G4LogicalVolume*   A_logicRing ;
    G4LogicalVolume*   S_logic;
    G4LogicalVolume*   B_logicDet;

    G4Box*             solidWorld;
    G4VPhysicalVolume* physWorld;

    G4Box*             solidEnv;
    G4VPhysicalVolume* physEnv;

    G4Box*             sample;
    G4VPhysicalVolume* physSample;

    G4Box*             substrate;
    G4VPhysicalVolume* physSubstrate;
 
    G4Cons*            A_solidRing;
    G4VPhysicalVolume* physAlphaDet;

    G4Tubs*            Shielding_solid;
    G4VPhysicalVolume* physShielding;

    G4Tubs*            BetaDet_solid;
    G4VPhysicalVolume* physBetaDet;

    G4VisAttributes * SampleVisAtt;
    G4VisAttributes * SubstrateVisAtt;
    G4VisAttributes * DetectorVisAtt;

    G4String          sample_mat_name;
    G4String          substrate_mat_name;
    G4double          SampleWidth;  
    G4double          SampleDepth;      
    G4double          SubstrateWidth;        
    G4double          SubstrateDepth;     
    G4double          A_cryst_dZ;
    G4double          A_cryst_dY;
    G4double          A_ring_R1;

    G4int            diagnostics ;
    G4String         filename;
    G4String         SRIMfilename;
    G4double         time;
    G4double         D;
    G4String         condition;
    G4double         beamOn;
    G4int            Libeam;
    G4String         nameOut;                

    G4int			  maxDepth;
    G4int			  NpointsZ;
    std::vector<G4double>  	  xZ;
    std::vector<G4double>  	  fZ;           //f(x)
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

