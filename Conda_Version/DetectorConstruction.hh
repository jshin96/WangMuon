#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

  private:
    // We keep pointers to logical volumes in case we want to attach 
    // SensitiveDetectors to them later instead of using SteppingAction.
//    G4LogicalVolume* fLogicDetectorIn1;
//    G4LogicalVolume* fLogicDetectorIn2;
//    G4LogicalVolume* fLogicDetectorOut1;
//    G4LogicalVolume* fLogicDetectorOut2;
    G4LogicalVolume* fLogicDetectorInner;
    G4LogicalVolume* fLogicDetectorOuter;
};

#endif
