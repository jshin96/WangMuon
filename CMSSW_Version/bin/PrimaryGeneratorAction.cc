#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4Exception.hh"
#include "G4ParticleTable.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include <array>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>

namespace {
constexpr G4double kInnerDetectorRadius = 23.0*m;
constexpr G4double kDetectorLayerSpacing = 400.0*cm;
constexpr G4double kOuterDetectorRadius =
    kInnerDetectorRadius + kDetectorLayerSpacing;
constexpr G4double kDetectorThickness = 1.0*cm;
constexpr G4double kPhysicalDetectorHeight = 13.0*m;
constexpr G4double kMoundRadius = 20.0*m;
constexpr G4double kRadialSafetyMargin = 2.0*m;
constexpr G4double kSourceRadius =
    kOuterDetectorRadius + kDetectorThickness + kRadialSafetyMargin;
constexpr G4double kPi = CLHEP::pi;

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

G4long ReadMaximumAttempts() {
    constexpr G4long kDefaultMaximumAttempts = 5000;
    const char* setting = std::getenv("MOUND_MAX_GENERATION_ATTEMPTS");
    if (setting == nullptr) return kDefaultMaximumAttempts;
    try {
        std::size_t consumed = 0;
        const long value = std::stol(setting, &consumed);
        if (consumed != std::string(setting).size() || value <= 0) {
            throw std::invalid_argument("not positive");
        }
        return static_cast<G4long>(value);
    } catch (const std::exception&) {
        G4Exception("PrimaryGeneratorAction", "InvalidMaximumAttempts",
                    FatalException,
                    "MOUND_MAX_GENERATION_ATTEMPTS must be a positive integer.");
        return kDefaultMaximumAttempts;
    }
}

G4double FoldOpposingPairAngleDeg(G4double angleDeg) {
    while (angleDeg >= 90.0) angleDeg -= 180.0;
    while (angleDeg < -90.0) angleDeg += 180.0;
    return angleDeg;
}

struct ShellCrossing {
    G4double entryParameter = 0.0;
    G4double exitParameter = 0.0;
    G4double entryX = 0.0;
    G4double entryY = 0.0;
    G4double exitX = 0.0;
    G4double exitY = 0.0;
    G4double entryAxial = 0.0;
    G4double exitAxial = 0.0;
    G4double entryPairAngleDeg = 0.0;
    G4double exitPairAngleDeg = 0.0;
};

bool IntersectDetectorShell(G4double radius,
                            G4double posX, G4double posY,
                            G4double posAxial,
                            G4double dirX, G4double dirY,
                            G4double dirAxial,
                            ShellCrossing& crossing) {
    const G4double a = dirX*dirX + dirY*dirY;
    if (a <= 1.e-18) return false;
    const G4double b = 2.0*(posX*dirX + posY*dirY);
    const G4double c = posX*posX + posY*posY - radius*radius;
    const G4double discriminant = b*b - 4.0*a*c;
    if (discriminant <= 0.0) return false;

    const G4double root = std::sqrt(discriminant);
    crossing.entryParameter = (-b - root)/(2.0*a);
    crossing.exitParameter = (-b + root)/(2.0*a);
    if (crossing.entryParameter <= 0.0 || crossing.exitParameter <= 0.0) {
        return false;
    }

    crossing.entryX = posX + crossing.entryParameter*dirX;
    crossing.entryY = posY + crossing.entryParameter*dirY;
    crossing.exitX = posX + crossing.exitParameter*dirX;
    crossing.exitY = posY + crossing.exitParameter*dirY;
    crossing.entryAxial =
        posAxial + crossing.entryParameter*dirAxial;
    crossing.exitAxial =
        posAxial + crossing.exitParameter*dirAxial;
    crossing.entryPairAngleDeg = FoldOpposingPairAngleDeg(
        std::atan2(crossing.entryX, crossing.entryY)*180.0/kPi);
    crossing.exitPairAngleDeg = FoldOpposingPairAngleDeg(
        std::atan2(crossing.exitX, crossing.exitY)*180.0/kPi);
    return true;
}
} // namespace

PrimaryGeneratorAction::PrimaryGeneratorAction(RunAction* runAction)
    : fParticleGun(new G4ParticleGun(1)),
      fEcoMug(new EcoMug()),
      fRunAction(runAction),
      fGroundSlope(0.0),
      fInclineCos(1.0),
      fInclineSin(0.0),
      fWindowMinDeg(-30.0),
      fWindowMaxDeg(30.0),
      fAxialMin(0.0),
      fAxialMax(kPhysicalDetectorHeight),
      fMaximumAttempts(ReadMaximumAttempts()) {
    const G4double groundInclineDeg = ReadFiniteEnvironmentDouble(
        "MOUND_GROUND_INCLINE_DEG", 0.0);
    if (std::abs(groundInclineDeg) >= 30.0) {
        G4Exception("PrimaryGeneratorAction", "InvalidGroundIncline",
                    FatalException,
                    "MOUND_GROUND_INCLINE_DEG must be between -30 and 30 degrees.");
    }
    const G4double groundIncline = groundInclineDeg*kPi/180.0;
    fGroundSlope = std::tan(groundIncline);
    fInclineCos = std::cos(groundIncline);
    fInclineSin = std::sin(groundIncline);

    fWindowMinDeg = ReadFiniteEnvironmentDouble(
        "MOUND_DETECTOR_WINDOW_MIN_DEG", -30.0);
    fWindowMaxDeg = ReadFiniteEnvironmentDouble(
        "MOUND_DETECTOR_WINDOW_MAX_DEG", 30.0);
    if (fWindowMinDeg <= -90.0 || fWindowMaxDeg >= 90.0
        || fWindowMinDeg >= fWindowMaxDeg) {
        G4Exception("PrimaryGeneratorAction", "InvalidDetectorWindow",
                    FatalException,
                    "Detector-window angles must satisfy -90 < min < max < 90 degrees.");
    }

    const G4double windowHeightM = ReadFiniteEnvironmentDouble(
        "MOUND_DETECTOR_WINDOW_HEIGHT_M", 13.0);
    const G4double windowElevationM = ReadFiniteEnvironmentDouble(
        "MOUND_DETECTOR_WINDOW_ELEVATION_M", 6.5);
    if (windowHeightM <= 0.0) {
        G4Exception("PrimaryGeneratorAction", "InvalidDetectorWindowHeight",
                    FatalException,
                    "MOUND_DETECTOR_WINDOW_HEIGHT_M must be positive.");
    }
    fAxialMin = (windowElevationM - 0.5*windowHeightM)*m;
    fAxialMax = (windowElevationM + 0.5*windowHeightM)*m;
    if (fAxialMin < 0.0 || fAxialMax > kPhysicalDetectorHeight
        || fAxialMin >= fAxialMax) {
        G4Exception("PrimaryGeneratorAction", "DetectorWindowOutsideShell",
                    FatalException,
                    "The detector-window elevation/height must lie inside the 0--13 m physical cylinder.");
    }

    const G4double minimumMomentumGeV = ReadFiniteEnvironmentDouble(
        "MOUND_MUON_MIN_MOMENTUM_GEV", 10.0);
    const G4double maximumMomentumGeV = ReadFiniteEnvironmentDouble(
        "MOUND_MUON_MAX_MOMENTUM_GEV", 100.0);
    if (minimumMomentumGeV <= 0.0
        || minimumMomentumGeV >= maximumMomentumGeV) {
        G4Exception("PrimaryGeneratorAction", "InvalidMomentumWindow",
                    FatalException,
                    "Muon momenta must satisfy 0 < minimum < maximum.");
    }

    fParticleGun->SetParticleDefinition(
        G4ParticleTable::GetParticleTable()->FindParticle("mu-"));

    fEcoMug->SetUseCylinder();
    fEcoMug->SetCylinderRadius(kSourceRadius);
    // EcoMug samples vertical height above z_ground. In the tilted-cylinder
    // frame, axial height is cos(incline) times that vertical height.
    const G4double sourceVerticalMin = fAxialMin/fInclineCos;
    const G4double sourceVerticalMax = fAxialMax/fInclineCos;
    fEcoMug->SetCylinderHeight(sourceVerticalMax - sourceVerticalMin);
    fEcoMug->SetCylinderCenterPosition(
        {{0.0, 0.0, 0.5*(sourceVerticalMin + sourceVerticalMax)}});

    // EcoMug position phi is measured counter-clockwise from global +X.
    // A downhill detector angle alpha measured from +Y toward +X has an
    // uphill source point at phi = 270 deg - alpha.
    const G4double sourcePhiMinDeg = 270.0 - fWindowMaxDeg;
    const G4double sourcePhiMaxDeg = 270.0 - fWindowMinDeg;
    fEcoMug->SetCylinderMinPositionPhi(sourcePhiMinDeg*kPi/180.0);
    fEcoMug->SetCylinderMaxPositionPhi(sourcePhiMaxDeg*kPi/180.0);

    fEcoMug->SetMinimumMomentum(minimumMomentumGeV);
    fEcoMug->SetMaximumMomentum(maximumMomentumGeV);
    fEcoMug->SetMinimumTheta(70.0*kPi/180.0);
    fEcoMug->SetMaximumTheta(80.0*kPi/180.0);

    G4cout << "CONDITIONAL_PRIMARY_CONFIGURATION downhill_pair_window_deg=["
           << fWindowMinDeg << "," << fWindowMaxDeg << "]"
           << " uphill_source_phi_deg=[" << sourcePhiMinDeg << ","
           << sourcePhiMaxDeg << "]"
           << " axial_window_m=[" << fAxialMin/m << "," << fAxialMax/m
           << "] momentum_GeV=[" << minimumMomentumGeV << ","
           << maximumMomentumGeV << "] max_attempts=" << fMaximumAttempts
           << G4endl;
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fParticleGun;
    delete fEcoMug;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    G4long trials = 0;
    bool accepted = false;
    G4ThreeVector acceptedPosition;
    G4ThreeVector acceptedDirection;
    G4double acceptedMomentumGeV = 0.0;
    ShellCrossing acceptedOuterCrossing;

    const auto insideWindow = [&](const ShellCrossing& crossing) {
        const bool correctSides = crossing.entryY < 0.0
                               && crossing.exitY > 0.0;
        const bool axialAccepted =
            crossing.entryAxial >= fAxialMin
            && crossing.entryAxial <= fAxialMax
            && crossing.exitAxial >= fAxialMin
            && crossing.exitAxial <= fAxialMax;
        const bool azimuthAccepted =
            crossing.entryPairAngleDeg >= fWindowMinDeg
            && crossing.entryPairAngleDeg <= fWindowMaxDeg
            && crossing.exitPairAngleDeg >= fWindowMinDeg
            && crossing.exitPairAngleDeg <= fWindowMaxDeg;
        return correctSides && axialAccepted && azimuthAccepted;
    };

    while (trials < fMaximumAttempts && !accepted) {
        fEcoMug->Generate();
        ++trials;

        std::array<double, 3> position = fEcoMug->GetGenerationPosition();
        // Interpret EcoMug Z as height above the local inclined ground.
        position[2] += -fGroundSlope*position[1];

        const G4double theta = fEcoMug->GetGenerationTheta();
        const G4double phi = fEcoMug->GetGenerationPhi();
        const G4double dirX = std::sin(theta)*std::cos(phi);
        const G4double dirY = std::sin(theta)*std::sin(phi);
        const G4double dirZ = -std::abs(std::cos(theta));
        if (dirY <= 0.0) continue;

        // Transform the candidate into the rigid tilted-cylinder frame.
        const G4double localPosX = position[0];
        const G4double localPosY =
            fInclineCos*position[1] - fInclineSin*position[2];
        const G4double localPosAxial =
            fInclineSin*position[1] + fInclineCos*position[2];
        const G4double localDirX = dirX;
        const G4double localDirY =
            fInclineCos*dirY - fInclineSin*dirZ;
        const G4double localDirAxial =
            fInclineSin*dirY + fInclineCos*dirZ;

        ShellCrossing outerCrossing;
        ShellCrossing innerCrossing;
        if (!IntersectDetectorShell(
                kOuterDetectorRadius, localPosX, localPosY,
                localPosAxial, localDirX, localDirY, localDirAxial,
                outerCrossing)
            || !IntersectDetectorShell(
                kInnerDetectorRadius, localPosX, localPosY,
                localPosAxial, localDirX, localDirY, localDirAxial,
                innerCrossing)
            || !insideWindow(outerCrossing)
            || !insideWindow(innerCrossing)) {
            continue;
        }

        // Preserve the existing requirement that both mound-footprint
        // crossings are downstream and remain within the configured axial
        // detector band. This rejects grazing rays before Geant4 transport.
        const G4double a = dirX*dirX + dirY*dirY;
        if (a <= 1.e-18) continue;
        const G4double b =
            2.0*(position[0]*dirX + position[1]*dirY);
        const G4double c =
            position[0]*position[0] + position[1]*position[1]
            - kMoundRadius*kMoundRadius;
        const G4double discriminant = b*b - 4.0*a*c;
        if (discriminant <= 0.0) continue;
        const G4double root = std::sqrt(discriminant);
        const G4double entry = (-b - root)/(2.0*a);
        const G4double exit = (-b + root)/(2.0*a);
        if (entry <= 0.0 || exit <= 0.0) continue;
        const G4double entryAxial =
            localPosAxial + entry*localDirAxial;
        const G4double exitAxial =
            localPosAxial + exit*localDirAxial;
        if (entryAxial < fAxialMin || entryAxial > fAxialMax
            || exitAxial < fAxialMin || exitAxial > fAxialMax) {
            continue;
        }

        acceptedPosition =
            G4ThreeVector(position[0], position[1], position[2]);
        acceptedDirection = G4ThreeVector(dirX, dirY, dirZ);
        acceptedMomentumGeV = fEcoMug->GetGenerationMomentum();
        acceptedOuterCrossing = outerCrossing;
        accepted = true;
    }

    if (!accepted) {
        fRunAction->RecordPrimaryGeneration(event->GetEventID(), trials, false);
        G4cout << "PRIMARY_GENERATION_LIMIT event=" << event->GetEventID()
               << " cap=" << fMaximumAttempts << G4endl;
        G4RunManager::GetRunManager()->AbortEvent();
        return;
    }

    fParticleGun->SetParticlePosition(acceptedPosition);
    fParticleGun->SetParticleMomentumDirection(acceptedDirection);
    // Retain the established project convention: EcoMug momentum in GeV/c
    // is passed as particle-gun kinetic energy in GeV.
    fParticleGun->SetParticleEnergy(acceptedMomentumGeV*GeV);

    if (std::getenv("MOUND_PRINT_PRIMARY_DIRECTIONS") != nullptr) {
        G4cout << "CONDITIONAL_PRIMARY event=" << event->GetEventID()
               << " position_m=(" << acceptedPosition.x()/m << ","
               << acceptedPosition.y()/m << ","
               << acceptedPosition.z()/m << ")"
               << " direction=(" << acceptedDirection.x() << ","
               << acceptedDirection.y() << ","
               << acceptedDirection.z() << ")"
               << " entry_angle_deg="
               << acceptedOuterCrossing.entryPairAngleDeg
               << " exit_angle_deg="
               << acceptedOuterCrossing.exitPairAngleDeg
               << " entry_axial_m="
               << acceptedOuterCrossing.entryAxial/m
               << " exit_axial_m="
               << acceptedOuterCrossing.exitAxial/m
               << " momentum_GeV=" << acceptedMomentumGeV
               << " trials=" << trials << G4endl;
    }

    fParticleGun->GeneratePrimaryVertex(event);
    fRunAction->RecordPrimaryGeneration(event->GetEventID(), trials, true);
    G4AnalysisManager::Instance()->FillNtupleDColumn(0, trials);
}
