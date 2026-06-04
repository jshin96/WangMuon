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




    //------------- Soil Definition---------------------
    G4Element* elO  = nistManager->FindOrBuildElement("O");
    G4Element* elSi = nistManager->FindOrBuildElement("Si");
    G4Element* elAl = nistManager->FindOrBuildElement("Al");
    G4Element* elFe = nistManager->FindOrBuildElement("Fe");
    G4Element* elCa = nistManager->FindOrBuildElement("Ca");
    G4Element* elH  = nistManager->FindOrBuildElement("H"); 
    G4Element* elC  = nistManager->FindOrBuildElement("C");

    // 2.a Define the Soil Material
    G4Material* soil = new G4Material("DrySoil", 1.6 * g/cm3, 7);

    // 2.b Add elements by fraction of mass
    soil->AddElement(elO,  0.51); // 51%
    soil->AddElement(elSi, 0.28); // 28%
    soil->AddElement(elAl, 0.07); // 7%
    soil->AddElement(elFe, 0.05); // 5%
    soil->AddElement(elCa, 0.03); // 3%
    soil->AddElement(elH,  0.04); // 4% (Acts as the neutron moderator)
    soil->AddElement(elC,  0.02); // 2%





    //------------- Granite Rock Definition---------------------
    // 3.a Pull the additional elements for Granite from the NIST database
    G4Element* elK  = nistManager->FindOrBuildElement("K");
    G4Element* elNa = nistManager->FindOrBuildElement("Na");
    G4Element* elMg = nistManager->FindOrBuildElement("Mg");

    // 3.b Construct the Granite material (Density: 2.7 g/cm^3, 8 components)
    G4Material* rock = new G4Material("Granite", 2.7 * g/cm3, 8);
    rock->AddElement(elO,  0.488); // 48.8% Oxygen
    rock->AddElement(elSi, 0.333); // 33.3% Silicon
    rock->AddElement(elAl, 0.077); // 7.7%  Aluminum
    rock->AddElement(elFe, 0.027); // 2.7%  Iron
    rock->AddElement(elK,  0.026); // 2.6%  Potassium
    rock->AddElement(elNa, 0.025); // 2.5%  Sodium
    rock->AddElement(elCa, 0.021); // 2.1%  Calcium
    rock->AddElement(elMg, 0.003); // 0.3%  Magnesium


    // 4. Define predefined materials (using the pointer we just created)
    G4Material* air = nistManager->FindOrBuildMaterial("G4_AIR");
    G4Material* scintillator = nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");









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
    G4LogicalVolume* logicMound = new G4LogicalVolume(solidMound, soil, "LogicMound");
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
    G4LogicalVolume* logicGround = new G4LogicalVolume(solidGround, soil, "LogicGround");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -12.5*m), logicGround, "PhysGround", logicWorld, false, 0, true);

    // =================================================================
    // 5. The Muon Detector ($2\pi$ Cylindrical Shell)
    // =================================================================
    G4double moundRadius = 26.5 * m;
    G4double detectorDistance = 2.0 * m;
    G4double detectorThickness = 2.0 * cm; // Thin plastic scintillator panel
    G4double detectorHeight = 20.0 * m;

    G4double innerRadius = moundRadius + detectorDistance; // 28.5m
    G4double outerRadius = innerRadius + detectorThickness;

    G4Tubs* solidDetector = new G4Tubs("SolidDetector",
                                       innerRadius,
                                       outerRadius,
                                       detectorHeight / 2.0, // Geant4 expects half-height
                                       0.0 * deg,
                                       360.0 * deg); // Full 2-pi coverage

    G4LogicalVolume* logicDetector = new G4LogicalVolume(solidDetector, scintillator, "LogicDetector");

    // Place it in the World volume, sitting on the ground (z goes from 0 to 10m, so center is at 5m)
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, detectorHeight / 2.0), logicDetector, "PhysDetector", logicWorld, false, 0, true);

    // Flat detector 
    G4Box* flatDetector = new G4Box("FlatDetector", 20.0*m, 2.0*cm, detectorHeight / 2.0);

    G4LogicalVolume* logicflatDetector = new G4LogicalVolume(flatDetector, scintillator, "LogicFlatDetector");

    // Place it in the World volume, sitting on the ground (z goes from 0 to 1m, so center is at 0.5m)
    new G4PVPlacement(nullptr, G4ThreeVector(0, innerRadius+(10.0*m)+(2.0*cm), detectorHeight / 2.0), logicflatDetector, "PhysFlatDetector", logicWorld, false, 0, true);




    return physWorld;
}
