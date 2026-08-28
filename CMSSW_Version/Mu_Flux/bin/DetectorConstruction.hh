#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1
#include "G4VUserDetectorConstruction.hh"
// This class builds the little world: three GEM boards and the magnet.
class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  DetectorConstruction() = default;
  ~DetectorConstruction() override = default;
  // Geant4 calls this once before the muons start moving.
  G4VPhysicalVolume* Construct() override;
};
#endif
