#ifndef LiActionInitialization_h
#define LiActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class LiDetectorConstruction;

/// Action initialization class.
///
class LiActionInitialization : public G4VUserActionInitialization
{
  public:
    LiActionInitialization(const LiDetectorConstruction* detector);
    virtual ~LiActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;

  private:
    const LiDetectorConstruction* fDetector;
};

#endif
    
