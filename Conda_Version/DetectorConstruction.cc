#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Exception.hh"
#include "G4ios.hh"
#include <cmath>
#include <cstdlib>
#include <string>

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fLogicDetectorInner(nullptr), fLogicDetectorOuter(nullptr)
{ }

DetectorConstruction::~DetectorConstruction()
{ }

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    // Positive inclination means that terrain descends as global +Y increases:
    // z_ground(y) = -tan(inclination) * y.  The local visualisation defaults
    // to the requested 15 degree site slope; set the environment variable to
    // zero to draw the former horizontal configuration.
    G4double groundInclineDeg = 15.0;
    if (const char* inclineSetting = std::getenv("MOUND_GROUND_INCLINE_DEG")) {
        try {
            groundInclineDeg = std::stod(inclineSetting);
        } catch (const std::exception&) {
            G4Exception("DetectorConstruction::Construct", "InvalidGroundIncline",
                        FatalException,
                        "MOUND_GROUND_INCLINE_DEG must be a finite number of degrees.");
        }
    }
    if (!std::isfinite(groundInclineDeg) || std::abs(groundInclineDeg) >= 30.0) {
        G4Exception("DetectorConstruction::Construct", "InvalidGroundIncline",
                    FatalException,
                    "MOUND_GROUND_INCLINE_DEG must be finite and between -30 and 30 degrees.");
    }
    const G4double groundIncline = groundInclineDeg * deg;
    const G4double groundSlope = std::tan(groundIncline);
    G4cout << "Ground inclination: " << groundInclineDeg
           << " deg (positive is downhill in +Y)" << G4endl;

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
    // 2. Vertical-axis mound resting on the tilted terrain
    // =================================================================
    // The mound height is always measured along global Z, not along the
    // ground normal: z_top = z_ground(y) + H sqrt(1-r^2/R^2).
    G4VSolid* solidMound = nullptr;
    constexpr G4double moundRadius = 26.5*m;
    constexpr G4double moundHeight = 12.7*m;
    if (std::abs(groundInclineDeg) < 1.0e-12) {
        solidMound = new G4Ellipsoid("SolidMound", moundRadius, moundRadius,
                                     moundHeight, 0., moundHeight);
    } else {
        constexpr int radialSegments = 32;
        constexpr int azimuthSegments = 96;
        constexpr G4double pi = 3.14159265358979323846;
        auto surfacePoint = [&](int radialIndex, int azimuthIndex, bool upper) {
            const G4double fraction = static_cast<G4double>(radialIndex) / radialSegments;
            const G4double radius = moundRadius * fraction;
            const G4double phi = 2.0 * pi * azimuthIndex / azimuthSegments;
            const G4double x = radius * std::cos(phi);
            const G4double y = radius * std::sin(phi);
            const G4double height = upper
                ? moundHeight * std::sqrt(1.0 - fraction * fraction) : 0.0;
            return G4ThreeVector(x, y, -groundSlope * y + height);
        };
        auto* tiltedMound = new G4TessellatedSolid("SolidMound");
        const G4ThreeVector lowerCentre(0., 0., 0.);
        const G4ThreeVector upperCentre(0., 0., moundHeight);
        for (int phi = 0; phi < azimuthSegments; ++phi) {
            const int next = (phi + 1) % azimuthSegments;
            tiltedMound->AddFacet(new G4TriangularFacet(
                lowerCentre, surfacePoint(1, next, false), surfacePoint(1, phi, false), ABSOLUTE));
            tiltedMound->AddFacet(new G4TriangularFacet(
                upperCentre, surfacePoint(1, phi, true), surfacePoint(1, next, true), ABSOLUTE));
        }
        for (int radial = 1; radial < radialSegments; ++radial) {
            for (int phi = 0; phi < azimuthSegments; ++phi) {
                const int next = (phi + 1) % azimuthSegments;
                const auto lowerInner = surfacePoint(radial, phi, false);
                const auto lowerOuter = surfacePoint(radial + 1, phi, false);
                const auto lowerOuterNext = surfacePoint(radial + 1, next, false);
                const auto lowerInnerNext = surfacePoint(radial, next, false);
                tiltedMound->AddFacet(new G4TriangularFacet(
                    lowerInner, lowerOuterNext, lowerOuter, ABSOLUTE));
                tiltedMound->AddFacet(new G4TriangularFacet(
                    lowerInner, lowerInnerNext, lowerOuterNext, ABSOLUTE));

                const auto upperInner = surfacePoint(radial, phi, true);
                const auto upperOuter = surfacePoint(radial + 1, phi, true);
                const auto upperOuterNext = surfacePoint(radial + 1, next, true);
                const auto upperInnerNext = surfacePoint(radial, next, true);
                tiltedMound->AddFacet(new G4TriangularFacet(
                    upperInner, upperOuter, upperOuterNext, ABSOLUTE));
                tiltedMound->AddFacet(new G4TriangularFacet(
                    upperInner, upperOuterNext, upperInnerNext, ABSOLUTE));
            }
        }
        tiltedMound->SetSolidClosed(true);
        solidMound = tiltedMound;
    }
    G4LogicalVolume* logicMound = new G4LogicalVolume(solidMound, soil, "LogicMound");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), logicMound, "PhysMound", logicWorld, false, 0, true);

    // =================================================================
    // 3. The Nested Room Setup (Mound -> Rock Wall -> Air Room)
    // =================================================================
    G4double wallThickness = 50.0 * cm;
    G4double roomX = 8.0 * m;
    G4double roomY = 8.0 * m;
    G4double roomZ = 4.0 * m;

    G4Box* solidRockWall = new G4Box("SolidRockWall",
                                     roomX/2 + wallThickness,
                                     roomY/2 + wallThickness,
                                     roomZ/2 + wallThickness);
    G4LogicalVolume* logicRockWall = new G4LogicalVolume(solidRockWall, rock, "LogicRockWall");

    G4Box* solidRoom = new G4Box("SolidRoom", roomX/2, roomY/2, roomZ/2);
    G4LogicalVolume* logicRoom = new G4LogicalVolume(solidRoom, air, "LogicRoom");

    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), logicRoom, "PhysRoom", logicRockWall, false, 0, true);

    // Keep the room level, but lower it until the outer wall is just above
    // the highest inclined-ground point beneath its footprint. A level room
    // cannot touch the sloped ground everywhere without intersecting it.
    const G4double wallHalfY = roomY/2 + wallThickness;
    const G4double wallHalfZ = roomZ/2 + wallThickness;
    const G4double highestGroundUnderWall =
        std::abs(groundSlope) * wallHalfY;
    constexpr G4double roomGroundClearance = 1.0*mm;
    const G4double wallCenterZ =
        highestGroundUnderWall + roomGroundClearance + wallHalfZ;
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, wallCenterZ),
                      logicRockWall, "PhysRockWall", logicMound,
                      false, 0, true);
    G4cout << "Room outer-wall bottom Z="
           << (wallCenterZ - wallHalfZ)/m
           << " m (lowest non-overlapping level placement)" << G4endl;

    // =================================================================
    // 4. The Ground
    // =================================================================
    G4VSolid* solidGround = nullptr;
    if (std::abs(groundInclineDeg) < 1.0e-12) {
        solidGround = new G4Box("SolidGround", worldSizeXY/2, worldSizeXY/2, worldSizeZ/4);
    } else {
        // Terrain is limited to the mound footprint so it remains clear of
        // the local detector planes at y = +/-30 and +/-32 m.
        constexpr int terrainSegments = 96;
        constexpr G4double pi = 3.14159265358979323846;
        const G4double terrainBottomZ = -45.0*m;
        auto topPoint = [&](int index) {
            const G4double phi = 2.0 * pi * index / terrainSegments;
            const G4double x = moundRadius * std::cos(phi);
            const G4double y = moundRadius * std::sin(phi);
            return G4ThreeVector(x, y, -groundSlope * y);
        };
        auto bottomPoint = [&](int index) {
            const G4double phi = 2.0 * pi * index / terrainSegments;
            return G4ThreeVector(moundRadius * std::cos(phi),
                                 moundRadius * std::sin(phi), terrainBottomZ);
        };
        auto* tiltedGround = new G4TessellatedSolid("SolidGround");
        const G4ThreeVector topCentre(0., 0., 0.);
        const G4ThreeVector bottomCentre(0., 0., terrainBottomZ);
        for (int index = 0; index < terrainSegments; ++index) {
            const int next = (index + 1) % terrainSegments;
            const auto top = topPoint(index);
            const auto topNext = topPoint(next);
            const auto bottom = bottomPoint(index);
            const auto bottomNext = bottomPoint(next);
            tiltedGround->AddFacet(new G4TriangularFacet(topCentre, top, topNext, ABSOLUTE));
            tiltedGround->AddFacet(new G4TriangularFacet(bottomCentre, bottomNext, bottom, ABSOLUTE));
            tiltedGround->AddFacet(new G4TriangularFacet(top, bottom, bottomNext, ABSOLUTE));
            tiltedGround->AddFacet(new G4TriangularFacet(top, bottomNext, topNext, ABSOLUTE));
        }
        tiltedGround->SetSolidClosed(true);
        solidGround = tiltedGround;
    }
    G4LogicalVolume* logicGround = new G4LogicalVolume(solidGround, soil, "LogicGround");
    new G4PVPlacement(nullptr, G4ThreeVector(), logicGround, "PhysGround", logicWorld, false, 0, true);

    // =================================================================
    // 5. Full-2pi rigid detector cylinders
    // =================================================================
    // Use two ordinary cylindrical shells with their common axis normal to
    // the inclined ground. Their lower caps lie in the ground plane and the
    // tracking layers are separated by 30 cm in radius.
    constexpr G4double detectorHeight = 13.0*m;
    constexpr G4double detectorThickness = 1.0*cm;
    constexpr G4double innerRadius = 28.0*m;
    constexpr G4double detectorLayerSpacing = 30.0*cm;
    constexpr G4double outerRadius = innerRadius + detectorLayerSpacing;

    fLogicDetectorInner = new G4LogicalVolume(
        new G4Tubs("TiltedDetectorInner", innerRadius,
                   innerRadius + detectorThickness, detectorHeight/2.0,
                   0.0, CLHEP::twopi),
        scintillator, "DetectorInner");
    fLogicDetectorOuter = new G4LogicalVolume(
        new G4Tubs("TiltedDetectorOuter", outerRadius,
                   outerRadius + detectorThickness, detectorHeight/2.0,
                   0.0, CLHEP::twopi),
        scintillator, "DetectorOuter");

    auto* detectorRotation = new G4RotationMatrix();
    detectorRotation->rotateX(-groundIncline);
    const G4ThreeVector detectorAxis =
        (*detectorRotation) * G4ThreeVector(0., 0., 1.);
    const G4ThreeVector detectorCentre = detectorAxis * (detectorHeight/2.0);
    const G4Transform3D detectorTransform(*detectorRotation, detectorCentre);
    new G4PVPlacement(detectorTransform, fLogicDetectorInner,
                      "DetectorInner", logicWorld, false, 0, true);
    new G4PVPlacement(detectorTransform, fLogicDetectorOuter,
                      "DetectorOuter", logicWorld, false, 0, true);
    G4cout << "Detector cylinders: inclination=" << groundInclineDeg
           << " deg, inner/outer layer radii=" << innerRadius/m << "/"
           << outerRadius/m << " m, radial layer spacing="
           << detectorLayerSpacing/cm << " cm" << G4endl;

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
    fLogicDetectorInner->SetVisAttributes(detVis);
    fLogicDetectorOuter->SetVisAttributes(detVis);

    return physWorld;
}
