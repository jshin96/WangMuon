#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"

#include "Randomize.hh"
#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include <array>
#include <cmath>
#include <cstdlib>

namespace {
constexpr G4double kSourceGroundClearance = 0.5*m;
constexpr G4double kMaximumSourceHeight = 22.0*m;
constexpr G4double kDetectorActiveHeight = 13.0*m;
constexpr G4double kInnerDetectorRadius = 28.0*m;
constexpr G4double kDetectorLayerSpacing = 30.0*cm;
constexpr G4double kOuterDetectorRadius =
    kInnerDetectorRadius + kDetectorLayerSpacing;
constexpr G4double kDetectorThickness = 1.0*cm;
constexpr G4double kMoundRadius = 26.5*m;

G4double GroundSlopeFromEnvironment() {
    const char* inclineSetting = std::getenv("MOUND_GROUND_INCLINE_DEG");
    const G4double inclineDeg =
        inclineSetting ? std::atof(inclineSetting) : 15.0;
    return std::tan(inclineDeg*CLHEP::pi/180.0);
}

G4long MaxGenerationAttemptsFromEnvironment() {
    constexpr G4long kDefaultMaxAttempts = 500;
    const char* attemptSetting =
        std::getenv("MOUND_MAX_GENERATION_ATTEMPTS");
    if (attemptSetting == nullptr) return kDefaultMaxAttempts;

    const long requested = std::strtol(attemptSetting, nullptr, 10);
    return requested > 0 ? static_cast<G4long>(requested)
                         : kDefaultMaxAttempts;
}
}

PrimaryGeneratorAction::PrimaryGeneratorAction(RunAction* runAction)
    : fRunAction(runAction) {
    fParticleGun = new G4ParticleGun(1);
    fParticleGun->SetParticleDefinition(
        G4ParticleTable::GetParticleTable()->FindParticle("mu-"));

    fEcoMug = new EcoMug();
    fEcoMug->SetUseCylinder();

    // Generate just outside the 28.3 m detector layer. EcoMug's cylinder Z
    // is treated as height above the local inclined ground below.
    constexpr G4double kRadialSafetyMargin = 2.0*m;
    constexpr G4double kSourceRadius =
        kOuterDetectorRadius + kDetectorThickness + kRadialSafetyMargin;
    fEcoMug->SetCylinderRadius(kSourceRadius);
    fEcoMug->SetCylinderHeight(
        kMaximumSourceHeight - kSourceGroundClearance);
    fEcoMug->SetCylinderCenterPosition(
        {{0., 0.,
          0.5*(kMaximumSourceHeight + kSourceGroundClearance)}});

    fEcoMug->SetMinimumMomentum(10.0);
    fEcoMug->SetMaximumMomentum(100.0);
    fEcoMug->SetMinimumTheta(54.0*CLHEP::pi/180.0);
    fEcoMug->SetMaximumTheta(90.0*CLHEP::pi/180.0);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fParticleGun;
    delete fEcoMug;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    const G4double groundSlope = GroundSlopeFromEnvironment();
    const G4double inclineCos =
        1.0/std::sqrt(1.0 + groundSlope*groundSlope);
    const G4double inclineSin = groundSlope*inclineCos;
    const G4long maximumAttempts = MaxGenerationAttemptsFromEnvironment();
    G4long trials = 0;

    bool accepted = false;
    G4ThreeVector acceptedPosition;
    G4ThreeVector acceptedDirection;
    G4double acceptedMomentumGeV = 0.0;

    while (!accepted) {
        fEcoMug->Generate();
        ++trials;

        if (trials > maximumAttempts) {
            fRunAction->RecordPrimaryGeneration(
                event->GetEventID(), trials, false);
            G4cout << "PRIMARY_GENERATION_LIMIT event="
                   << event->GetEventID()
                   << " cap=" << maximumAttempts << G4endl;
            G4RunManager::GetRunManager()->AbortEvent();
            return;
        }

        std::array<double, 3> position = fEcoMug->GetGenerationPosition();
        position[2] += -groundSlope*position[1];

        const G4double theta = fEcoMug->GetGenerationTheta();
        const G4double phi = fEcoMug->GetGenerationPhi();
        const G4double dirX = std::sin(theta)*std::cos(phi);
        const G4double dirY = std::sin(theta)*std::sin(phi);
        const G4double dirZ = -std::abs(std::cos(theta));

        // +Y is downhill. Downward-going uphill tracks cannot remain inside
        // both ends of the finite tilted detector cylinder.
        if (dirY < 0.0) continue;

        // Convert the candidate ray into the rigid detector-cylinder frame.
        const G4double localPosX = position[0];
        const G4double localPosY =
            inclineCos*position[1] - inclineSin*position[2];
        const G4double localPosAxial =
            inclineSin*position[1] + inclineCos*position[2];
        const G4double localDirX = dirX;
        const G4double localDirY =
            inclineCos*dirY - inclineSin*dirZ;
        const G4double localDirAxial =
            inclineSin*dirY + inclineCos*dirZ;

        const auto crossesDetectorLayer = [&](G4double radius) {
            const G4double a =
                localDirX*localDirX + localDirY*localDirY;
            if (a <= 1.e-18) return false;
            const G4double b =
                2.0*(localPosX*localDirX + localPosY*localDirY);
            const G4double c =
                localPosX*localPosX + localPosY*localPosY
                - radius*radius;
            const G4double discriminant = b*b - 4.0*a*c;
            if (discriminant <= 0.0) return false;

            const G4double root = std::sqrt(discriminant);
            const G4double entry = (-b - root)/(2.0*a);
            const G4double exit = (-b + root)/(2.0*a);
            if (entry <= 0.0 || exit <= 0.0) return false;

            const G4double entryHeight =
                localPosAxial + entry*localDirAxial;
            const G4double exitHeight =
                localPosAxial + exit*localDirAxial;
            return entryHeight > 0.0
                && entryHeight < kDetectorActiveHeight
                && exitHeight > 0.0
                && exitHeight < kDetectorActiveHeight;
        };

        if (!crossesDetectorLayer(kOuterDetectorRadius)
            || !crossesDetectorLayer(kInnerDetectorRadius)) {
            continue;
        }

        // Require the same ray to cross both sides of the mound footprint
        // inside the detector's active axial band before Geant4 transport.
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

        const G4double entryHeight =
            localPosAxial + entry*localDirAxial;
        const G4double exitHeight =
            localPosAxial + exit*localDirAxial;
        if (entryHeight <= 0.0
            || entryHeight >= kDetectorActiveHeight
            || exitHeight <= 0.0
            || exitHeight >= kDetectorActiveHeight) {
            continue;
        }

        acceptedPosition =
            G4ThreeVector(position[0], position[1], position[2]);
        acceptedDirection = G4ThreeVector(dirX, dirY, dirZ);
        acceptedMomentumGeV = fEcoMug->GetGenerationMomentum();
        accepted = true;
    }

    fParticleGun->SetParticlePosition(acceptedPosition);
    fParticleGun->SetParticleMomentumDirection(acceptedDirection);
    fParticleGun->SetParticleEnergy(acceptedMomentumGeV*GeV);

    if (std::getenv("MOUND_PRINT_PRIMARY_DIRECTIONS") != nullptr) {
        G4cout << "COSMIC_PRIMARY event=" << event->GetEventID()
               << " position_m=(" << acceptedPosition.x()/m << ","
               << acceptedPosition.y()/m << ","
               << acceptedPosition.z()/m << ")"
               << " direction=(" << acceptedDirection.x() << ","
               << acceptedDirection.y() << ","
               << acceptedDirection.z() << ")"
               << " momentum_GeV=" << acceptedMomentumGeV
               << " trials=" << trials << G4endl;
    }

    fParticleGun->GeneratePrimaryVertex(event);
    fRunAction->RecordPrimaryGeneration(event->GetEventID(), trials, true);
    G4AnalysisManager::Instance()->FillNtupleDColumn(0, trials);
}
