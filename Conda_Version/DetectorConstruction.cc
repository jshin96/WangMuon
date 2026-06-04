#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Ellipsoid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fLogicDetectorIn1(nullptr), fLogicDetectorIn2(nullptr),
  fLogicDetectorOut1(nullptr), fLogicDetectorOut2(nullptr)
{ }

DetectorConstruction::~DetectorConstruction()
{ }

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    // -----------------------------------------------------
    // 1. Materials & Elements
    // -----------------------------------------------------
    G4NistManager* nistManager = G4NistManager::Instance();
    G4Material* air = nistManager->FindOrBuildMaterial("G4_AIR");
    G4Material* scintillator = nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

    G4Element* elO  = nistManager->FindOrBuildElement("O");
    G4Element* elSi = nistManager->FindOrBuildElement("Si");
    G4Element* elAl = nistManager->FindOrBuildElement("Al");
    G4Element* elFe = nistManager->FindOrBuildElement("Fe");
    G4Element* elCa = nistManager->FindOrBuildElement("Ca");
    G4Element* elH  = nistManager->FindOrBuildElement("H");
    G4Element* elC  = nistManager->FindOrBuildElement("C");
    G4Element* elK  = nistManager->FindOrBuildElement("K");
    G4Element* elNa = nistManager->FindOrBuildElement("Na");
    G4Element* elMg = nistManager->FindOrBuildElement("Mg");

    G4Material* soil = new G4Material("DrySoil", 1.6 * g/cm3, 7);
    soil->AddElement(elO,  0.51); 
    soil->AddElement(elSi, 0.28); 
    soil->AddElement(elAl, 0.07); 
    soil->AddElement(elFe, 0.05); 
    soil->AddElement(elCa, 0.03); 
    soil->AddElement(elH,  0.04); 
    soil->AddElement(elC,  0.02); 

    G4Material* rock = new G4Material("Granite", 2.7 * g/cm3, 8);
    rock->AddElement(elO,  0.488); 
    rock->AddElement(elSi, 0.333); 
    rock->AddElement(elAl, 0.077); 
    rock->AddElement(elFe, 0.027); 
    rock->AddElement(elK,  0.026); 
    rock->AddElement(elNa, 0.025); 
    rock->AddElement(elCa, 0.021); 
    rock->AddElement(elMg, 0.003); 

    // -----------------------------------------------------
    // 1.5 World Volume
    // -----------------------------------------------------
    G4double worldSizeXY = 150 * m; // Expanded to fit wider detectors
    G4double worldSizeZ  = 100 * m;
    
    G4Box* solidWorld = new G4Box("World", worldSizeXY/2, worldSizeXY/2, worldSizeZ/2);
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, air, "World");
    G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, true);

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
    G4double roomX = 3.25 * m;
    G4double roomY = 2.1  * m;
    G4double roomZ = 1.05 * m;

    G4Box* solidRockWall = new G4Box("SolidRockWall",
                                     roomX + wallThickness,
                                     roomY + wallThickness,
                                     roomZ + wallThickness);
    G4LogicalVolume* logicRockWall = new G4LogicalVolume(solidRockWall, rock, "LogicRockWall");

    G4Box* solidRoom = new G4Box("SolidRoom", roomX, roomY, roomZ);
    G4LogicalVolume* logicRoom = new G4LogicalVolume(solidRoom, air, "LogicRoom");

    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), logicRoom, "PhysRoom", logicRockWall, false, 0, true);

    G4double wallCenterZ = roomZ + wallThickness; 
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, wallCenterZ), logicRockWall, "PhysRockWall", logicMound, false, 0, true);

    // =================================================================
    // 4. The Ground
    // =================================================================
    G4Box* solidGround = new G4Box("SolidGround", worldSizeXY/2, worldSizeXY/2, worldSizeZ/4);
    G4LogicalVolume* logicGround = new G4LogicalVolume(solidGround, soil, "LogicGround");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, - worldSizeZ/4), logicGround, "PhysGround", logicWorld, false, 0, true);

    // =================================================================
    // 5. The Detectors (4 Layers on Y-Axis, flanking the 26.5m mound)
    // =================================================================
    G4double detSizeX = 10.0 * m;  // 10m total width 
    G4double detSizeZ = 10.0 * m;  // 10m total height
    G4double detThickness = 1.0 * cm; 
    
    G4Box* solidDetector = new G4Box("DetectorShape", detSizeX/2, detThickness, detSizeZ/2);

    // INCOMING FLANK (-Y side, beyond -26.5m)
    fLogicDetectorIn1 = new G4LogicalVolume(solidDetector, scintillator, "DetectorIn1");
    new G4PVPlacement(0, G4ThreeVector(0, -32.0*m, detSizeZ/2), fLogicDetectorIn1, "DetectorIn1", logicWorld, false, 0, true);

    fLogicDetectorIn2 = new G4LogicalVolume(solidDetector, scintillator, "DetectorIn2");
    new G4PVPlacement(0, G4ThreeVector(0, -30.0*m, detSizeZ/2), fLogicDetectorIn2, "DetectorIn2", logicWorld, false, 0, true);

    // OUTGOING FLANK (+Y side, beyond +26.5m)
    fLogicDetectorOut1 = new G4LogicalVolume(solidDetector, scintillator, "DetectorOut1");
    new G4PVPlacement(0, G4ThreeVector(0, 30.0*m, detSizeZ/2), fLogicDetectorOut1, "DetectorOut1", logicWorld, false, 0, true);

    fLogicDetectorOut2 = new G4LogicalVolume(solidDetector, scintillator, "DetectorOut2");
    new G4PVPlacement(0, G4ThreeVector(0, 32.0*m, detSizeZ/2), fLogicDetectorOut2, "DetectorOut2", logicWorld, false, 0, true);

    // -----------------------------------------------------
    // 6. Visual Attributes
    // -----------------------------------------------------
    logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

    G4VisAttributes* groundVis = new G4VisAttributes(G4Colour(0.4, 0.25, 0.1, 0.4)); // Darker soil
    groundVis->SetForceSolid(true);
    logicGround->SetVisAttributes(groundVis);

    G4VisAttributes* moundVis = new G4VisAttributes(G4Colour(0.6, 0.4, 0.2, 0.3)); // Lighter soil
    moundVis->SetForceSolid(true);
    logicMound->SetVisAttributes(moundVis);

    G4VisAttributes* rockVis = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.8)); // Solid grey rock
    rockVis->SetForceSolid(true);
    logicRockWall->SetVisAttributes(rockVis);

    G4VisAttributes* roomVis = new G4VisAttributes(G4Colour(0.0, 0.8, 1.0, 1.0)); // Air interior
    roomVis->SetForceSolid(true);
    logicRoom->SetVisAttributes(roomVis);

    G4VisAttributes* detVis = new G4VisAttributes(G4Colour(0.0, 0.8, 1.0, 0.5)); // Cyan detectors
    detVis->SetForceSolid(true);
    fLogicDetectorIn1->SetVisAttributes(detVis);
    fLogicDetectorIn2->SetVisAttributes(detVis);
    fLogicDetectorOut1->SetVisAttributes(detVis);
    fLogicDetectorOut2->SetVisAttributes(detVis);

    return physWorld;
}
