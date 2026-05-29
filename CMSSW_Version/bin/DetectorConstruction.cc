// DetectorConstruction.cc
#include "DetectorConstruction.hh"
#include "G4NistManager.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Ellipsoid.hh"


DetectorConstruction::DetectorConstruction() {}
DetectorConstruction::~DetectorConstruction() {}
G4VPhysicalVolume* DetectorConstruction::Construct() {
    // 1. Grab the NIST database pointer
    G4NistManager* nistManager = G4NistManager::Instance();
    
    // 2. Define the materials (using the pointer we just created)
    G4Material* air = nistManager->FindOrBuildMaterial("G4_AIR");
    G4Material* dirt = nistManager->FindOrBuildMaterial("G4_SILICON_DIOXIDE"); 
    
    G4Material* rock = nistManager->FindOrBuildMaterial("G4_ROCK_DENSITY_SINK") 
                       ? nistManager->FindOrBuildMaterial("G4_ROCK_DENSITY_SINK") 
                       : nistManager->FindOrBuildMaterial("G4_CONCRETE");
                       
    G4Material* scintillator = nistManager->FindOrBuildMaterial("G4_PLASTIC_SCINTILLATOR");

    // =================================================================
    // 1. The World (50x50x50 meters)
    // =================================================================
    G4Box* solidWorld = new G4Box("SolidWorld", 50*m, 50*m, 25*m);
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, air, "LogicWorld");
    G4VPhysicalVolume* physWorld = new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), logicWorld, "PhysWorld", nullptr, false, 0, true);

    // =================================================================
    // 2. The Dirt Mound (Ellipsoid: r=26.5m, base at z=0)
    // =================================================================
    G4Ellipsoid* solidMound = new G4Ellipsoid("SolidMound", 26.5*m, 26.5*m, 12.7*m, 0., 12.7*m);
    G4LogicalVolume* logicMound = new G4LogicalVolume(solidMound, dirt, "LogicMound");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), logicMound, "PhysMound", logicWorld, false, 0, true);

    // =================================================================
    // 3. The Nested Room Setup (Mound -> Rock Wall -> Air Room)
    // =================================================================
    G4double wallThickness = 20.0 * cm;

    // Half-dimensions of the inner air room (Original: 6.5m x 4.2m x 2.1m total)
    G4double roomX = 3.25 * m;
    G4double roomY = 2.1  * m;
    G4double roomZ = 1.05 * m;

    // Create the outer Rock Wall Box (Room dimensions + thickness)
    G4Box* solidRockWall = new G4Box("SolidRockWall",
                                     roomX + wallThickness,
                                     roomY + wallThickness,
                                     roomZ + wallThickness);
    G4LogicalVolume* logicRockWall = new G4LogicalVolume(solidRockWall, rock, "LogicRockWall");

    // Create the inner Air Room Box
    G4Box* solidRoom = new G4Box("SolidRoom", roomX, roomY, roomZ);
    G4LogicalVolume* logicRoom = new G4LogicalVolume(solidRoom, air, "LogicRoom");

    // Place the Air Room directly in the CENTER of the Rock Wall box
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), logicRoom, "PhysRoom", logicRockWall, false, 0, true);

    // Place the complete Rock Wall assembly inside the Mound.
    // Shifted up slightly to ensure the 20cm rock floor rests perfectly on the z=0 ground plane.
    G4double wallCenterZ = roomZ + wallThickness; // 1.05m + 0.2m = 1.25m
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, wallCenterZ), logicRockWall, "PhysRockWall", logicMound, false, 0, true);

    // =================================================================
    // 4. The Ground
    // =================================================================
    G4Box* solidGround = new G4Box("SolidGround", 50*m, 50*m, 12.5*m);
    G4LogicalVolume* logicGround = new G4LogicalVolume(solidGround, dirt, "LogicGround");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -12.5*m), logicGround, "PhysGround", logicWorld, false, 0, true);

    // =================================================================
    // 5. The Muon Detector ($2\pi$ Cylindrical Shell)
    // =================================================================
    G4double moundRadius = 26.5 * m;
    G4double detectorDistance = 2.0 * m;
    G4double detectorThickness = 2.0 * cm; // Thin plastic scintillator panel
    G4double detectorHeight = 1.0 * m;

    G4double innerRadius = moundRadius + detectorDistance; // 28.5m
    G4double outerRadius = innerRadius + detectorThickness;

    G4Tubs* solidDetector = new G4Tubs("SolidDetector",
                                       innerRadius,
                                       outerRadius,
                                       detectorHeight / 2.0, // Geant4 expects half-height
                                       0.0 * deg,
                                       360.0 * deg); // Full 2-pi coverage

    G4LogicalVolume* logicDetector = new G4LogicalVolume(solidDetector, scintillator, "LogicDetector");

    // Place it in the World volume, sitting on the ground (z goes from 0 to 1m, so center is at 0.5m)
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, detectorHeight / 2.0), logicDetector, "PhysDetector", logicWorld, false, 0, true);


    return physWorld;
}
