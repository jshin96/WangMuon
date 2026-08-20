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

// Resolve convenient metal names while still permitting any Geant4/NIST
// material identifier (for example, G4_Al or G4_Au).
G4Material* ReadScatteringPlateMaterial(G4NistManager* nistManager) {
    std::string materialName = "tungsten";
    if (const char* setting = std::getenv("MOUND_SCATTERING_PLATE_MATERIAL")) {
        materialName = setting;
    }

    std::string normalizedName = materialName;
    std::transform(normalizedName.begin(), normalizedName.end(),
                   normalizedName.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    std::string nistName;
    if (normalizedName == "tungsten" || normalizedName == "w") {
        nistName = "G4_W";
    } else if (normalizedName == "lead" || normalizedName == "pb") {
        nistName = "G4_Pb";
    } else if (normalizedName == "copper" || normalizedName == "cu") {
        nistName = "G4_Cu";
    } else if (normalizedName == "iron" || normalizedName == "fe") {
        nistName = "G4_Fe";
    } else if (normalizedName.rfind("g4_", 0) == 0) {
        // Preserve the supplied NIST spelling, whose element symbols are case-sensitive.
        nistName = materialName;
    } else {
        G4Exception("DetectorConstruction::Construct", "InvalidPlateMaterial",
                    FatalException,
                    "MOUND_SCATTERING_PLATE_MATERIAL must be tungsten, lead, copper, iron, or a valid G4_* NIST material name.");
        return nullptr;
    }

    G4Material* material = nistManager->FindOrBuildMaterial(nistName, false);
    if (material == nullptr) {
        const std::string message = "Unknown Geant4/NIST plate material: " + nistName;
        G4Exception("DetectorConstruction::Construct", "InvalidPlateMaterial",
                    FatalException, message.c_str());
    }
    return material;
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
    G4Material* scatteringPlateMaterial = ReadScatteringPlateMaterial(nistManager);
    // A parameterised triple-GEM response is applied in EventAction.  The
    // active volume is nevertheless gas, so Geant4 provides the correct
    // energy deposition and passive-material scattering baseline.
    G4Material* argon = nistManager->FindOrBuildMaterial("G4_Ar");
    G4Material* carbonDioxide = nistManager->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
    G4Material* gemGas = G4Material::GetMaterial("GEMGasArCO2", false);
    if (!gemGas) {
        // 70/30 by volume Ar/CO2, expressed as approximate mass fractions.
        gemGas = new G4Material("GEMGasArCO2", 1.84*mg/cm3, 2,
                                kStateGas, 293.15*kelvin, 1.0*atmosphere);
        gemGas->AddMaterial(argon, 0.620);
        gemGas->AddMaterial(carbonDioxide, 0.380);
    }

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
    const G4double wallThickness = ReadFiniteEnvironmentDouble(
        "MOUND_ROOM_WALL_THICKNESS_CM", 50.0)*cm;
    const G4double roomX = ReadFiniteEnvironmentDouble(
        "MOUND_ROOM_SIZE_X_M", 4.0)*m;
    const G4double roomY = ReadFiniteEnvironmentDouble(
        "MOUND_ROOM_SIZE_Y_M", 4.0)*m;
    const G4double roomZ = ReadFiniteEnvironmentDouble(
        "MOUND_ROOM_SIZE_Z_M", 4.0)*m;
    if (roomX <= 0.0 || roomY <= 0.0 || roomZ <= 0.0
        || wallThickness <= 0.0) {
        G4Exception("DetectorConstruction::Construct", "InvalidRoomSize",
                    FatalException,
                    "MOUND_ROOM_SIZE_X_M, MOUND_ROOM_SIZE_Y_M, MOUND_ROOM_SIZE_Z_M, and MOUND_ROOM_WALL_THICKNESS_CM must be positive.");
    }

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
               << " m; interior size (x,y,z)=(" << roomX/m << ", "
               << roomY/m << ", " << roomZ/m << ") m"
               << "; wall thickness=" << wallThickness/cm << " cm"
               << "; outer-wall -Y,-Z edge terrain z="
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
    // 5. Flat GEM detector planes
    // -----------------------------------------------------
    // Planes are vertical, with their thin dimension along global Y.  Each
    // plane has an independent X width, Z height, Y centre, and elevation
    // above the inclined ground at that Y coordinate.
    // This makes the measured five-plane layout explicit instead of deriving
    // it from two angular sectors of cylindrical shells.
    constexpr G4double detectorThickness = 1.0*cm;
    struct FlatGEMConfig {
        G4String label;
        G4double width;
        G4double height;
        G4double y;
        G4double z;
    };
    const auto readGEMSetting = [](const G4String& label, const char* suffix,
                                   G4double fallback) {
        const std::string name = "MOUND_GEM_" + std::string(label) + suffix;
        return ReadFiniteEnvironmentDouble(name.c_str(), fallback);
    };
    const auto makeGEMConfig = [&readGEMSetting, groundSlope](const char* label,
                                                 G4double defaultY) {
        const G4String gemLabel(label);
        const G4double y = readGEMSetting(gemLabel, "_Y_M", defaultY)*m;
        const G4double elevation = readGEMSetting(gemLabel, "_Z_M", 1.0)*m;
        return FlatGEMConfig{
            gemLabel,
            readGEMSetting(gemLabel, "_WIDTH_M", 1.0)*m,
            readGEMSetting(gemLabel, "_HEIGHT_M", 1.0)*m,
            y,
            elevation - groundSlope*y};
    };
    const FlatGEMConfig gemIn1 = makeGEMConfig("IN1", -27.0);
    const FlatGEMConfig gemIn2 = makeGEMConfig("IN2", -23.0);
    const FlatGEMConfig gemOut1 = makeGEMConfig("OUT1", 23.0);
    const FlatGEMConfig gemOut2 = makeGEMConfig("OUT2", 27.0);
    const FlatGEMConfig gemOut3 = makeGEMConfig("OUT3", 28.04);
    const FlatGEMConfig gemConfigs[] = {
        gemIn1, gemIn2, gemOut1, gemOut2, gemOut3};
    for (const auto& config : gemConfigs) {
        if (config.width <= 0.0 || config.height <= 0.0) {
            G4Exception("DetectorConstruction::Construct", "InvalidFlatGEMSize",
                        FatalException,
                        "Every MOUND_GEM_<PLANE>_WIDTH_M and _HEIGHT_M must be positive.");
        }
    }

    G4VisAttributes* detVis = new G4VisAttributes(G4Colour(0.0, 0.8, 1.0, 0.5));
    detVis->SetForceSolid(true);
    for (const auto& config : gemConfigs) {
        const G4String volumeName = "Detector" + config.label;
        auto* logic = new G4LogicalVolume(
            new G4Box("Flat" + volumeName, config.width/2.0,
                      detectorThickness/2.0, config.height/2.0),
            gemGas, volumeName);
        new G4PVPlacement(nullptr, G4ThreeVector(0., config.y, config.z), logic,
                          volumeName, logicWorld, false, 0, true);
        logic->SetVisAttributes(detVis);
        if (config.label == "IN1") fLogicDetectorInner = logic;
        if (config.label == "OUT2") fLogicDetectorOuter = logic;
    }

    // The scattering plate is centred halfway between OUT2 and OUT3.  Its
    // position therefore follows the independently configured detector planes.
    const G4double plateRadialThickness = ReadFiniteEnvironmentDouble(
        "MOUND_TUNGSTEN_THICKNESS_CM", 3.0)*cm;
    if (plateRadialThickness <= 0.0) {
        G4Exception("DetectorConstruction::Construct", "InvalidPostTungstenGeometry",
                    FatalException,
                    "MOUND_TUNGSTEN_THICKNESS_CM must be positive.");
    }
    auto* scatteringPlate = new G4Box("ScatteringPlate",
                                      std::min(gemOut2.width, gemOut3.width)/2.0,
                                      plateRadialThickness/2.0,
                                      std::min(gemOut2.height, gemOut3.height)/2.0);
    auto* logicScatteringPlate = new G4LogicalVolume(scatteringPlate,
                                                      scatteringPlateMaterial,
                                                      "ScatteringPlate");
    const G4ThreeVector plateCentre(0., 0.5*(gemOut2.y + gemOut3.y),
                                    0.5*(gemOut2.z + gemOut3.z));
    new G4PVPlacement(nullptr, plateCentre, logicScatteringPlate,
                      "ScatteringPlate", logicWorld, false, 0, true);







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

    G4VisAttributes* scatteringPlateVis = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.8));
    scatteringPlateVis->SetForceSolid(true);
    logicScatteringPlate->SetVisAttributes(scatteringPlateVis);
    return physWorld;
}
