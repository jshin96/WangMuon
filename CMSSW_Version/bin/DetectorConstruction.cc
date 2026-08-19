#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4IntersectionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"
#include <cstdlib>
#include <cmath>
#include <string>
#include <algorithm>
#include <cctype>
G4double ReadFiniteEnvironmentDouble(const char* name, G4double fallback) {
    const char* setting = std::getenv(name);
    if (setting == nullptr) return fallback;
    try {
        std::size_t consumed = 0;
        const G4double value = std::stod(setting, &consumed);
        if (consumed != std::string(setting).size() || !std::isfinite(value)) {
            throw std::invalid_argument("not a finite number");
        }
        return value;
    } catch (const std::exception&) {
        const std::string message = std::string(name)
            + " must be a finite number.";
        G4Exception("PrimaryGeneratorAction", "InvalidGeneratorSetting",
                    FatalException, message.c_str());
        return fallback;
    }
}
DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fLogicDetectorInner(nullptr), fLogicDetectorOuter(nullptr)
//  fLogicDetectorInner(nullptr), fLogicDetectorIn2(nullptr),
//  fLogicDetectorOuter(nullptr), fLogicDetectorOut2(nullptr)
{ }

DetectorConstruction::~DetectorConstruction()
{ }

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    // Keep the source and detector setup identical for the room/no-room
    // comparison.  Only the material geometry changes between the two runs.
    const char* roomSetting = std::getenv("MOUND_INCLUDE_ROOM");
    const bool includeRoom = !(roomSetting && std::string(roomSetting) == "0");
    // Select the room fill at runtime so air-filled and soil-filled samples
    // use the identical compiled geometry and detector/source setup.
    std::string roomMaterialName = "air";
    if (const char* materialSetting = std::getenv("MOUND_ROOM_MATERIAL")) {
        roomMaterialName = materialSetting;
    }
    std::transform(roomMaterialName.begin(), roomMaterialName.end(),
                   roomMaterialName.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (roomMaterialName != "air" && roomMaterialName != "soil") {
        G4Exception("DetectorConstruction::Construct", "InvalidRoomMaterial",
                    FatalException,
                    "MOUND_ROOM_MATERIAL must be either 'air' or 'soil'.");
    }
    // Positive inclination means that the terrain falls as global +Y increases:
    // z_ground(y) = -tan(inclination) * y.  Leave it at zero for the
    // established horizontal-geometry baseline.
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
    G4cout << "Mound room geometry: " << (includeRoom ? "enabled" : "disabled")
           << "; room/wall materials: "
           << (roomMaterialName == "air" ? "air/granite" : "soil/soil")
           << G4endl;
    G4cout << "Ground inclination: " << groundInclineDeg
           << " deg (positive is downhill in +Y)" << G4endl;
    // -----------------------------------------------------
    // 1. Materials & Elements
    // -----------------------------------------------------
    G4NistManager* nistManager = G4NistManager::Instance();
    G4Material* air = nistManager->FindOrBuildMaterial("G4_AIR");
    G4Material* tungsten = nistManager->FindOrBuildMaterial("G4_W");
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
    // 2. The Dirt Mound
    // =================================================================
    // Preserve the established 12.7 m hemi-ellipsoid for horizontal terrain.
    // On inclined terrain, use one vertical analytic ellipsoid whose base is
    // at the lowest ground elevation and whose apex is moundHeight above the
    // highest ground elevation.  Intersect it with the half-space above the terrain so
    // the mound and ground touch but never overlap as sibling volumes.
    G4VSolid* solidMound = nullptr;
    constexpr G4double moundRadius = 19.5*m;
    constexpr G4double moundHeight = 8.0*m;
    constexpr G4double moundGroundSeparation = 1.0*um;
    G4double moundPlacementZ = 0.0;
    // The two points where the mound/ground contact curve reaches its
    // highest and lowest global elevations.  They also define the room's
    // centre along Y.
    // For horizontal terrain the entire base rim is a contact curve, so its
    // unambiguous centre is y=0.  On an incline both values below are solved
    // from the actual ellipse/plane intersection.
    G4double highestContactY = 0.0;
    G4double lowestContactY = 0.0;
    if (std::abs(groundInclineDeg) < 1.0e-12) {
        solidMound = new G4Ellipsoid("SolidMound", moundRadius, moundRadius,
                                     moundHeight, 0., moundHeight);
    } else {
        const G4double maximumGroundElevation =
            std::abs(groundSlope) * moundRadius;
        moundPlacementZ = -maximumGroundElevation;
        // Use the declared mound height as the clearance above the highest
        // part of the terrain; do not introduce a second height parameter.
        const G4double ellipsoidHeight =
            moundHeight + 2.0 * maximumGroundElevation;
        auto* moundEnvelope =
            new G4Ellipsoid("SolidMoundEnvelope", moundRadius, moundRadius,
                            ellipsoidHeight, 0., ellipsoidHeight);

        // In mound-local coordinates the ground plane is
        // z = -groundSlope*y - moundPlacementZ.  A rotated box supplies the
        // finite representation of the half-space above that plane.
        const G4double halfSpaceHalfXY = worldSizeXY;
        const G4double halfSpaceHalfZ = worldSizeZ;
        auto* aboveGround =
            new G4Box("SolidAboveGround", halfSpaceHalfXY, halfSpaceHalfXY,
                      halfSpaceHalfZ);
        auto* halfSpaceRotation = new G4RotationMatrix();
        halfSpaceRotation->rotateX(-groundIncline);
        // A one-micrometre numerical separation prevents Geant4's overlap
        // sampler from treating the coincident mound/ground faces as overlap.
        const G4ThreeVector planePoint(
            0., 0., -moundPlacementZ + moundGroundSeparation);
        const G4ThreeVector halfSpaceCentre =
            planePoint
            + (*halfSpaceRotation) * G4ThreeVector(0., 0., halfSpaceHalfZ);
        // Use the explicit object transform overload; the pointer overload
        // instead interprets its rotation/translation as a frame transform.
        const G4Transform3D halfSpaceTransform(*halfSpaceRotation,
                                               halfSpaceCentre);
        solidMound =
            new G4IntersectionSolid("SolidMound", moundEnvelope, aboveGround,
                                    halfSpaceTransform);

        // Solve the ellipse/ground-plane intersection at x=0.  The result
        // remains valid for any mound radius and nonzero ground inclination.
        const G4double heightRatio = maximumGroundElevation / ellipsoidHeight;
        const G4double contactDenominator =
            heightRatio*heightRatio + 1.0;
        const G4double uphillContactFactor =
            (heightRatio*heightRatio - 1.0) / contactDenominator;
        const G4double downhillContactFactor =
            (heightRatio*heightRatio + 1.0) / contactDenominator;
        const G4double downhillYDirection =
            std::copysign(1.0, groundSlope);
        highestContactY =
            downhillYDirection * moundRadius * uphillContactFactor;
        lowestContactY =
            downhillYDirection * moundRadius * downhillContactFactor;

        G4cout << "Inclined mound analytic envelope: base Z="
               << moundPlacementZ/m << " m, apex Z="
               << (moundPlacementZ + ellipsoidHeight)/m
               << " m (" << moundHeight/m
               << " m above highest ground)" << G4endl;
    }
    G4LogicalVolume* logicMound = new G4LogicalVolume(solidMound, soil, "LogicMound");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, moundPlacementZ), logicMound,
                      "PhysMound", logicWorld, false, 0, true);

    // =================================================================
    // 3. The Nested Room Setup (Mound -> Rock Wall -> Air Room)
    // =================================================================
    G4double wallThickness = 50.0 * cm;
    G4double roomX = 4.0 * m;
    G4double roomY = 4.0 * m;
    G4double roomZ = 4.0 * m;

    G4Box* solidRockWall = new G4Box("SolidRockWall",
                                     roomX/2 + wallThickness,
                                     roomY/2 + wallThickness,
                                     roomZ/2 + wallThickness);
    // The air-filled target includes its physical granite wall.  The
    // soil-filled control replaces both wall and cavity with mound soil, so
    // it has no material contrast or artificial internal boundary.
    G4Material* wallMaterial = roomMaterialName == "air" ? rock : soil;
    G4LogicalVolume* logicRockWall = new G4LogicalVolume(solidRockWall, wallMaterial,
                                                          "LogicRockWall");

    G4Box* solidRoom = new G4Box("SolidRoom", roomX/2, roomY/2, roomZ/2);
    G4Material* roomMaterial = roomMaterialName == "air" ? air : soil;
    G4LogicalVolume* logicRoom = new G4LogicalVolume(solidRoom, roomMaterial, "LogicRoom");

    // The room and its wall remain vertical in global Z.  Centre the room
    // using the requested separation of the two terrain-contact points.
    // Place the outer wall so its -Y, -Z edge lies on the inclined terrain.
    const G4double roomCentreY =
        0.5 * (lowestContactY + highestContactY);
    const G4double wallHalfY = roomY/2 + wallThickness;
    const G4double wallHalfZ = roomZ/2 + wallThickness;
    const G4double wallMinusYEdgeGroundZ =
        -groundSlope * (roomCentreY - wallHalfY);
    const G4double wallCenterZ = wallMinusYEdgeGroundZ + wallHalfZ;
    if (includeRoom) {
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), logicRoom,
                          "PhysRoom", logicRockWall, false, 0, true);
        new G4PVPlacement(nullptr,
                          G4ThreeVector(0, roomCentreY,
                                        wallCenterZ - moundPlacementZ),
                          logicRockWall,
                          "PhysRockWall", logicMound, false, 0, true);
        G4cout << "Room centre: y=" << roomCentreY/m
               << " m; outer-wall -Y,-Z edge terrain z="
               << wallMinusYEdgeGroundZ/m << " m"
               << G4endl;
    }

    // =================================================================
    // 4. The Ground
    // =================================================================
    G4VSolid* solidGround = nullptr;
    G4double groundPlacementZ = 0.0;
    if (std::abs(groundInclineDeg) < 1.0e-12) {
        solidGround = new G4Box("SolidGround", worldSizeXY/2, worldSizeXY/2,
                                worldSizeZ/4);
        // Put the upper face at global z=0 instead of centring the ground box
        // on the mound, which would overlap the entire horizontal mound.
        groundPlacementZ = -worldSizeZ/4;
    } else {
        // The terrain only needs to support the mound footprint.  A finite
        // sloped terrain disk avoids intersecting the detector cylinders at
        // 28--29 m radius while retaining the exact surface beneath the mound.
        constexpr int terrainSegments = 96;
        constexpr G4double pi = 3.14159265358979323846;
        const G4double terrainBottomZ = -25.0*m;
        auto terrainTopPoint = [&](int index) {
            const G4double phi = 2.0 * pi * index / terrainSegments;
            const G4double x = moundRadius * std::cos(phi);
            const G4double y = moundRadius * std::sin(phi);
            return G4ThreeVector(x, y, -groundSlope * y);
        };
        auto terrainBottomPoint = [&](int index) {
            const G4double phi = 2.0 * pi * index / terrainSegments;
            return G4ThreeVector(moundRadius * std::cos(phi),
                                 moundRadius * std::sin(phi), terrainBottomZ);
        };

        auto* tiltedGround = new G4TessellatedSolid("SolidGround");
        const G4ThreeVector terrainTopCentre(0., 0., 0.);
        const G4ThreeVector terrainBottomCentre(0., 0., terrainBottomZ);
        for (int index = 0; index < terrainSegments; ++index) {
            const int nextIndex = (index + 1) % terrainSegments;
            const auto top = terrainTopPoint(index);
            const auto topNext = terrainTopPoint(nextIndex);
            const auto bottom = terrainBottomPoint(index);
            const auto bottomNext = terrainBottomPoint(nextIndex);
            // Upper terrain surface, horizontal lower surface, and vertical rim.
            tiltedGround->AddFacet(new G4TriangularFacet(
                terrainTopCentre, top, topNext, ABSOLUTE));
            tiltedGround->AddFacet(new G4TriangularFacet(
                terrainBottomCentre, bottomNext, bottom, ABSOLUTE));
            tiltedGround->AddFacet(new G4TriangularFacet(
                top, bottom, bottomNext, ABSOLUTE));
            tiltedGround->AddFacet(new G4TriangularFacet(
                top, bottomNext, topNext, ABSOLUTE));
        }
        tiltedGround->SetSolidClosed(true);
        solidGround = tiltedGround;
    }
    G4LogicalVolume* logicGround = new G4LogicalVolume(solidGround, soil, "LogicGround");
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., groundPlacementZ),
                      logicGround, "PhysGround", logicWorld, false, 0, true);

    // =================================================================
    // 5. The Detectors (4 Layers on Y-Axis, flanking the 26.5m mound)
    // =================================================================



    /*
    G4double detSizeX = 40.0 * m;  
    G4double detSizeZ = 15.0 * m;  
    G4double detThickness = 1.0 * cm; 
    
    G4Box* InsolidDetector = new G4Box("DetectorShape", detSizeX/2, detThickness, detSizeZ/2);


    // INCOMING FLANK (-Y side, beyond -26.5m)
    fLogicDetectorIn1 = new G4LogicalVolume(InsolidDetector, scintillator, "DetectorIn1");
    new G4PVPlacement(0, G4ThreeVector(0, -29.0*m, detSizeZ/2), fLogicDetectorIn1, "DetectorIn1", logicWorld, false, 0, true);

    fLogicDetectorIn2 = new G4LogicalVolume(InsolidDetector, scintillator, "DetectorIn2");
    new G4PVPlacement(0, G4ThreeVector(0, -28.0*m, detSizeZ/2), fLogicDetectorIn2, "DetectorIn2", logicWorld, false, 0, true);

    // Flat Detector in +Y side
    G4Box* OutsolidDetector = new G4Box("DetectorShape", 30*detSizeX/2, detThickness, detSizeZ/2);
    // OUTGOING FLANK (+Y side, beyond +26.5m)
    fLogicDetectorOut1 = new G4LogicalVolume(OutsolidDetector, scintillator, "DetectorOut1");
    new G4PVPlacement(0, G4ThreeVector(0, 28.0*m, detSizeZ/2), fLogicDetectorOut1, "DetectorOut1", logicWorld, false, 0, true);

    fLogicDetectorOut2 = new G4LogicalVolume(OutsolidDetector, scintillator, "DetectorOut2");
    new G4PVPlacement(0, G4ThreeVector(0, 29.0*m, detSizeZ/2), fLogicDetectorOut2, "DetectorOut2", logicWorld, false, 0, true);
    */





    /*
    // Cylinderical out detector
    G4double wallHeight = 5.0 * m;
    G4double InnerRadius = 28.0 * m;
    G4double OuterRadius = 29.0 * m;
    
    // Define the arc to wrap around the +Y side
    G4double OutstartAngle = 0.0;   // Starts at the +X axis
    G4double spanAngle  = CLHEP::pi; // Sweeps through +Y axis and ends at -X axis
    G4double InstartAngle = CLHEP::pi;   // Starts at the -X axis


    // G4Tubs arguments: (Name, InnerRadius, OuterRadius, HalfHeight, StartAngle, SpanAngle)
    G4Tubs* InsolidDetector1 = new G4Tubs("CylInWall1", InnerRadius, InnerRadius + detThickness, wallHeight/2.0, InstartAngle, spanAngle);
    G4Tubs* InsolidDetector2 = new G4Tubs("CylInWall2", OuterRadius, OuterRadius + detThickness, wallHeight/2.0, InstartAngle, spanAngle);
    G4Tubs* OutsolidDetector1 = new G4Tubs("CylOutWall1", InnerRadius, InnerRadius + detThickness, wallHeight/2.0, OutstartAngle, spanAngle);
    G4Tubs* OutsolidDetector2 = new G4Tubs("CylOutWall2", OuterRadius, OuterRadius + detThickness, wallHeight/2.0, OutstartAngle, spanAngle);

    // Place the cylinders centered at X=0, Y=0 so the radius perfectly flanks the mound
    fLogicDetectorIn1 = new G4LogicalVolume(InsolidDetector1, scintillator, "DetectorIn1");
    new G4PVPlacement(0, G4ThreeVector(0, 0, wallHeight/2.0), fLogicDetectorIn1, "DetectorIn1", logicWorld, false, 0, true);

    fLogicDetectorIn2 = new G4LogicalVolume(InsolidDetector2, scintillator, "DetectorIn2");
    new G4PVPlacement(0, G4ThreeVector(0, 0, wallHeight/2.0), fLogicDetectorIn2, "DetectorIn2", logicWorld, false, 0, true);

    fLogicDetectorOut1 = new G4LogicalVolume(OutsolidDetector1, scintillator, "DetectorOut1");
    new G4PVPlacement(0, G4ThreeVector(0, 0, wallHeight/2.0), fLogicDetectorOut1, "DetectorOut1", logicWorld, false, 0, true);

    fLogicDetectorOut2 = new G4LogicalVolume(OutsolidDetector2, scintillator, "DetectorOut2");
    new G4PVPlacement(0, G4ThreeVector(0, 0, wallHeight/2.0), fLogicDetectorOut2, "DetectorOut2", logicWorld, false, 0, true);

    */









    // -----------------------------------------------------
    // 5. Full-2pi rigid detector cylinders
    // -----------------------------------------------------
    // Use two ordinary cylindrical shells with their common axis normal to
    // the inclined ground.  Their lower caps lie in the ground plane and the
    // layers are separated by 30 cm in radius.
//    constexpr G4double detectorHeight = 13.0*m;
    constexpr G4double detectorThickness = 1.0*cm;
    constexpr G4double innerRadius = 23.0*m;
    constexpr G4double detectorLayerSpacing = 400.0*cm;
    constexpr G4double outerRadius = innerRadius + detectorLayerSpacing;
    auto* detectorRotation = new G4RotationMatrix();
    detectorRotation->rotateX(-groundIncline);
    const G4ThreeVector detectorAxis =
        (*detectorRotation) * G4ThreeVector(0., 0., 1.);
    const G4double detectorHeight = ReadFiniteEnvironmentDouble(
        "MOUND_DETECTOR_WINDOW_HEIGHT_M", 1.0)*m;
    const G4double detectorElevationM = ReadFiniteEnvironmentDouble(
        "MOUND_DETECTOR_WINDOW_ELEVATION_M", 1.0)*m;
    if (detectorHeight <= 0.0) {
        G4Exception("PrimaryGeneratorAction", "InvalidDetectorWindowHeight",
                    FatalException,
                    "MOUND_DETECTOR_WINDOW_HEIGHT_M must be positive.");
    }
    const G4ThreeVector detectorCentre = detectorAxis * ((detectorHeight/2.0)+(detectorElevationM));
    const G4Transform3D detectorTransform(*detectorRotation, detectorCentre);
    G4double DetDegMin = -1.1;
    G4double DetDegMax = 1.1;
    if (const char* DetDegMinSetting = std::getenv("MOUND_DETECTOR_WINDOW_MIN_DEG")) {
        try {
            DetDegMin = std::stod(DetDegMinSetting);
        } catch (const std::exception&) {
            G4Exception("DetectorConstruction::Construct", "InvalidMOUND_DETECTOR_WINDOW_MIN_DEG",
                        FatalException,
                        "MOUND_DETECTOR_WINDOW_MIN_DEG must be a finite number of degrees.");
        }
    }
    if (const char* DetDegMaxSetting = std::getenv("MOUND_DETECTOR_WINDOW_MAX_DEG")) {
        try {
            DetDegMax = std::stod(DetDegMaxSetting);
        } catch (const std::exception&) {
            G4Exception("DetectorConstruction::Construct", "InvalidMOUND_DETECTOR_WINDOW_MAX_DEG",
                        FatalException,
                        "MOUND_DETECTOR_WINDOW_MAX_DEG must be a finite number of degrees.");
        }
    }
    if (!std::isfinite(DetDegMin) || !std::isfinite(DetDegMax)
        || DetDegMax <= DetDegMin || DetDegMax - DetDegMin > 360.0) {
        G4Exception("DetectorConstruction::Construct", "InvalidDetectorWindow",
                    FatalException,
                    "Detector-window bounds must be finite, ordered, and span at most 360 degrees.");
    }
    G4double sectorSpan = (DetDegMax-DetDegMin)*deg;

    const G4double sectorStarts[] = {
        CLHEP::pi/2.0 + (DetDegMin*deg),
        3.0*CLHEP::pi/2.0 + (DetDegMin*deg)
    };
    G4VisAttributes* detVis = new G4VisAttributes(G4Colour(0.0, 0.8, 1.0, 0.5));
    detVis->SetForceSolid(true);
    for (G4int copyNo = 0; copyNo < 2; ++copyNo) {
        auto* innerLogic = new G4LogicalVolume(
            new G4Tubs("TiltedDetectorInnerSector", innerRadius,
                     innerRadius + detectorThickness, detectorHeight/2.0,
                     sectorStarts[copyNo], sectorSpan),
            scintillator, "DetectorInner");

        auto* outerLogic = new G4LogicalVolume(
            new G4Tubs("TiltedDetectorOuterSector", outerRadius,
                     outerRadius + detectorThickness, detectorHeight/2.0,
                     sectorStarts[copyNo], sectorSpan),
            scintillator, "DetectorOuter");

        new G4PVPlacement(detectorTransform, innerLogic,
                          "DetectorInner", logicWorld, false, copyNo, true);
        new G4PVPlacement(detectorTransform, outerLogic,
                          "DetectorOuter", logicWorld, false, copyNo, true);
	if (copyNo == 0) {
		fLogicDetectorInner = innerLogic;
		fLogicDetectorOuter = outerLogic;
	}
	innerLogic->SetVisAttributes(detVis);
	outerLogic->SetVisAttributes(detVis);
    }



    // A 1 m x 1 m tungsten face immediately outside the inner +Y detector
    // sector.  Define and place it in the detector's tilted local frame so
    // its radial face remains parallel to that detector layer.
    constexpr G4double plateFaceSize = 1.0*m;
    constexpr G4double plateRadialThickness = 30.0*cm;
    constexpr G4double plateClearance = 1.0*um;
    auto* tungstenPlate = new G4Box("TungstenPlate",
                                    plateFaceSize/2.0,
                                    plateRadialThickness/2.0,
                                    plateFaceSize/2.0);
    auto* logicTungstenPlate = new G4LogicalVolume(tungstenPlate, tungsten,
                                                   "TungstenPlate");
    const G4double plateCentreRadius =
        innerRadius + detectorThickness + plateClearance
        + plateRadialThickness/2.0;
    // detectorTransform translates the cylinder origin by detectorCentre so
    // its lower cap is on the ground plane.  Apply that same translation to
    // the plate; rotation alone would incorrectly leave it at the lower cap.
    const G4ThreeVector plateCentre = detectorCentre
        + (*detectorRotation) * G4ThreeVector(0., plateCentreRadius+(1.0*cm), 0.);
    const G4Transform3D plateTransform(*detectorRotation, plateCentre);
    new G4PVPlacement(plateTransform, logicTungstenPlate,
                      "TungstenPlate", logicWorld, false, 0, true);







    // -----------------------------------------------------
    // 6. Visual Attributes
    // -----------------------------------------------------
    logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

    G4VisAttributes* groundVis = new G4VisAttributes(G4Colour(0.4, 0.25, 0.1, 0.4)); // Darker soil
    groundVis->SetForceSolid(true);
    logicGround->SetVisAttributes(groundVis);

    G4VisAttributes* moundVis = new G4VisAttributes(G4Colour(0.6, 0.4, 0.2, 0.6)); // Lighter soil
    moundVis->SetForceSolid(true);
    logicMound->SetVisAttributes(moundVis);

    G4VisAttributes* rockVis = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.8)); // Solid grey rock
    rockVis->SetForceSolid(true);
    logicRockWall->SetVisAttributes(rockVis);

    G4VisAttributes* roomVis = new G4VisAttributes(G4Colour(0.0, 0.8, 1.0, 1.0)); // Air interior
    roomVis->SetForceSolid(true);
    logicRoom->SetVisAttributes(roomVis);

//    G4VisAttributes* detVis = new G4VisAttributes(G4Colour(0.0, 0.8, 1.0, 0.5)); // Cyan detectors
//    detVis->SetForceSolid(true);
//    fLogicDetectorInner->SetVisAttributes(detVis);
//    fLogicDetectorOuter->SetVisAttributes(detVis);

    G4VisAttributes* TungstenPlateVis = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.8)); // Solid grey rock
    TungstenPlateVis->SetForceSolid(true);
    logicTungstenPlate->SetVisAttributes(TungstenPlateVis);
    return physWorld;
}
