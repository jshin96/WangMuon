// DetectorConstruction.cc
#include "DetectorConstruction.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Ellipsoid.hh"


DetectorConstruction::DetectorConstruction() {}
DetectorConstruction::~DetectorConstruction() {}
G4VPhysicalVolume* DetectorConstruction::Construct() {
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
    G4Material* dirt = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE"); 

    // 1. The World (50x50x50 meters)
    G4Box* solidWorld = new G4Box("SolidWorld", 50*m, 50*m, 25*m);
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, air, "LogicWorld");
    G4VPhysicalVolume* physWorld = new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), logicWorld, "PhysWorld", nullptr, false, 0, true);

    // 2. The Dirt Mound (Ellipsoid: r=20m, height=10m. Cut bottom at z=0, top at z=10)
    G4Ellipsoid* solidMound = new G4Ellipsoid("SolidMound", 26.5*m, 26.5*m, 12.7*m, 0., 12.7*m);
    G4LogicalVolume* logicMound = new G4LogicalVolume(solidMound, dirt, "LogicMound");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), logicMound, "PhysMound", logicWorld, false, 0, true);

    // 3. The Empty Room (10x5x3 meters). Center of mound is at z=5m.
    G4Box* solidRoom = new G4Box("SolidRoom", 3.25*m, 2.1*m, 1.05*m);
    G4LogicalVolume* logicRoom = new G4LogicalVolume(solidRoom, air, "LogicRoom");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 1.05*m), logicRoom, "PhysRoom", logicMound, false, 0, true);

    // 4. The Ground 
    G4Box* solidGround = new G4Box("SolidGround", 50*m,50*m,12.5*m);
    G4LogicalVolume* logicGround = new G4LogicalVolume(solidGround, dirt, "LogicGround");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -12.5*m), logicGround, "PhysGround", logicWorld, false, 0, true);

    // --- VISUALIZATION ATTRIBUTES ---
    
    // 1. Make the 50m World box completely invisible
    logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

    // 2. Make the Dirt Mound semi-transparent brown (so we can see inside it)
    // G4Colour(Red, Green, Blue, Alpha/Transparency)
    G4VisAttributes* dirtVis = new G4VisAttributes(G4Colour(0.6, 0.4, 0.2, 0.2)); 
    dirtVis->SetForceSolid(true);
    logicMound->SetVisAttributes(dirtVis);

    // 3. Make the Empty Room solid blue so it stands out brightly in the center
    G4VisAttributes* roomVis = new G4VisAttributes(G4Colour(0.0, 0.5, 1.0, 0.8));
    roomVis->SetForceSolid(true);
    logicRoom->SetVisAttributes(roomVis);

    // 4. Make the Ground solid, dark brown
    G4VisAttributes* groundVis = new G4VisAttributes(G4Colour(0.3, 0.2, 0.1, 0.4));
    groundVis->SetForceSolid(true);
    logicGround->SetVisAttributes(groundVis);

    return physWorld;
}
