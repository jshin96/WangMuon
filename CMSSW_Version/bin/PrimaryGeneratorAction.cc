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
constexpr G4double kMoundRadius = 20.0*m;
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
        const std::string message = std::string(name) + " must be a finite number.";
        G4Exception("PrimaryGeneratorAction", "InvalidGeneratorSetting", FatalException,
                    message.c_str());
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
        G4Exception("PrimaryGeneratorAction", "InvalidMaximumAttempts", FatalException,
                    "MOUND_MAX_GENERATION_ATTEMPTS must be a positive integer.");
        return kDefaultMaximumAttempts;
    }
}

G4long ReadRateIntegrationPoints() {
    constexpr G4long kDefaultRateIntegrationPoints = 100000;
    const char* setting = std::getenv("MOUND_ECOMUG_RATE_INTEGRATION_POINTS");
    if (setting == nullptr) return kDefaultRateIntegrationPoints;
    try {
        std::size_t consumed = 0;
        const long value = std::stol(setting, &consumed);
        if (consumed != std::string(setting).size() || value <= 1) {
            throw std::invalid_argument("not greater than one");
        }
        return static_cast<G4long>(value);
    } catch (const std::exception&) {
        G4Exception("PrimaryGeneratorAction", "InvalidRateIntegrationPoints", FatalException,
                    "MOUND_ECOMUG_RATE_INTEGRATION_POINTS must be an integer greater than one.");
        return kDefaultRateIntegrationPoints;
    }
}
} // namespace

PrimaryGeneratorAction::PrimaryGeneratorAction(RunAction* runAction)
    : fParticleGun(new G4ParticleGun(1)), fEcoMug(new EcoMug()),
      fRunAction(runAction), fGroundSlope(0.0), fIn1{}, fIn2{},
      fMaximumAttempts(ReadMaximumAttempts()) {
    const G4double groundInclineDeg = ReadFiniteEnvironmentDouble(
        "MOUND_GROUND_INCLINE_DEG", 15.0);
    if (std::abs(groundInclineDeg) >= 30.0) {
        G4Exception("PrimaryGeneratorAction", "InvalidGroundIncline", FatalException,
                    "MOUND_GROUND_INCLINE_DEG must be between -30 and 30 degrees.");
    }
    const G4double groundSlope = std::tan(groundInclineDeg*kPi/180.0);
    fGroundSlope = groundSlope;
    const auto readPlane = [groundSlope](const char* plane) {
        const std::string prefix = std::string("MOUND_GEM_") + plane;
        const G4double y = ReadFiniteEnvironmentDouble(
            (prefix + "_Y_M").c_str(), 0.0)*m;
        const G4double elevation = ReadFiniteEnvironmentDouble(
            (prefix + "_Z_M").c_str(), 1.0)*m;
        const G4double height = ReadFiniteEnvironmentDouble(
            (prefix + "_HEIGHT_M").c_str(), 1.0)*m;
        return FlatPlane{
            ReadFiniteEnvironmentDouble((prefix + "_WIDTH_M").c_str(), 1.0)*m,
            height,
            y, elevation - groundSlope*y + height/2.0};
    };
    fIn1 = readPlane("IN1");
    fIn2 = readPlane("IN2");
    if (fIn1.width <= 0.0 || fIn1.height <= 0.0 || fIn2.width <= 0.0
        || fIn2.height <= 0.0 || fIn1.y >= 0.0 || fIn2.y <= fIn1.y) {
        G4Exception("PrimaryGeneratorAction", "InvalidIncomingGEMLayout", FatalException,
                    "IN1/IN2 widths and heights must be positive; IN1 must be on -Y and upstream of IN2.");
    }

    const G4double minimumMomentumGeV = ReadFiniteEnvironmentDouble(
        "MOUND_MUON_MIN_MOMENTUM_GEV", 10.0);
    const G4double maximumMomentumGeV = ReadFiniteEnvironmentDouble(
        "MOUND_MUON_MAX_MOMENTUM_GEV", 100.0);
    if (minimumMomentumGeV <= 0.0 || minimumMomentumGeV >= maximumMomentumGeV) {
        G4Exception("PrimaryGeneratorAction", "InvalidMomentumWindow", FatalException,
                    "Muon momenta must satisfy 0 < minimum < maximum.");
    }

    const G4double sourceClearance = ReadFiniteEnvironmentDouble(
        "MOUND_SOURCE_CLEARANCE_M", 2.0)*m;
    const G4double sourceElevationMin = ReadFiniteEnvironmentDouble(
        "MOUND_SOURCE_Z_MIN_M", 1.0)*m;
    const G4double sourceElevationMax = ReadFiniteEnvironmentDouble(
        "MOUND_SOURCE_Z_MAX_M", 5.0)*m;
    const G4double sourceArcHalfWidthDeg = ReadFiniteEnvironmentDouble(
        "MOUND_SOURCE_UPSTREAM_ARC_HALF_WIDTH_DEG", 20.0);
    const G4double sourceDirectionHalfWidthDeg = ReadFiniteEnvironmentDouble(
        "MOUND_SOURCE_DIRECTION_PHI_HALF_WIDTH_DEG", 30.0);
    if (sourceClearance <= 0.0 || sourceElevationMax <= sourceElevationMin
        || sourceArcHalfWidthDeg <= 0.0 || sourceArcHalfWidthDeg >= 90.0
        || sourceDirectionHalfWidthDeg <= 0.0 || sourceDirectionHalfWidthDeg >= 90.0) {
        G4Exception("PrimaryGeneratorAction", "InvalidFlatSource", FatalException,
                    "Source clearance must be positive, Z max must exceed Z min, and source angular half-widths must be between 0 and 90 degrees.");
    }

    fParticleGun->SetParticleDefinition(
        G4ParticleTable::GetParticleTable()->FindParticle("mu-"));
    fEcoMug->SetUseCylinder();
    fEcoMug->SetCylinderRadius(std::abs(fIn2.y) + sourceClearance);
    // EcoMug samples source Z as elevation above local ground.  Each sampled
    // point is later translated by its own inclined-ground height.
    fEcoMug->SetCylinderHeight(sourceElevationMax - sourceElevationMin);
    fEcoMug->SetCylinderCenterPosition(
        {{0.0, 0.0, 0.5*(sourceElevationMin + sourceElevationMax)}});

    // EcoMug position phi is counter-clockwise from +X; -Y is 270 degrees.
    const G4double sourcePhiMinDeg = 270.0 - sourceArcHalfWidthDeg;
    const G4double sourcePhiMaxDeg = 270.0 + sourceArcHalfWidthDeg;
    fEcoMug->SetCylinderMinPositionPhi(sourcePhiMinDeg*kPi/180.0);
    fEcoMug->SetCylinderMaxPositionPhi(sourcePhiMaxDeg*kPi/180.0);
    // EcoMug's cylinder direction azimuth is local to the sampled source
    // point: Generate() adds the source-position phi before returning it.
    // A source point on the upstream (-Y) arc is near 270 degrees, so its
    // local direction must be near 180 degrees to produce a world direction
    // near +Y (90 degrees) and cross IN1 then IN2.
    fEcoMug->SetMinimumPhi((180.0 - sourceDirectionHalfWidthDeg)*kPi/180.0);
    fEcoMug->SetMaximumPhi((180.0 + sourceDirectionHalfWidthDeg)*kPi/180.0);
    fEcoMug->SetMinimumMomentum(minimumMomentumGeV);
    fEcoMug->SetMaximumMomentum(maximumMomentumGeV);
    fEcoMug->SetMinimumTheta(70.0*kPi/180.0);
    fEcoMug->SetMaximumTheta(80.0*kPi/180.0);

    const G4double horizontalRateHzM2 = ReadFiniteEnvironmentDouble(
        "MOUND_HORIZONTAL_MUON_RATE_HZ_M2", 129.0);
    if (horizontalRateHzM2 <= 0.0) {
        G4Exception("PrimaryGeneratorAction", "InvalidHorizontalMuonRate", FatalException,
                    "MOUND_HORIZONTAL_MUON_RATE_HZ_M2 must be positive.");
    }
    const G4long rateIntegrationPoints = ReadRateIntegrationPoints();
    fEcoMug->SetHorizontalRate(horizontalRateHzM2*EMUnits::hertz/EMUnits::m2);
    double sourceRatePerAreaHzM2 = 0.0;
    double sourceRatePerAreaErrorHzM2 = 0.0;
    // EcoMug generates inward cylinder rays around local phi=pi.  Its rate
    // integrator expects the equivalent outward-normal convention around
    // local phi=0 and otherwise clips cos(phi)<0 to zero.  Shift only this
    // copied integrator by pi; the physical generated tracks are unchanged.
    EcoMug rateCalculator(*fEcoMug);
    rateCalculator.SetMinimumPhi(fEcoMug->GetMinimumPhi() - kPi);
    rateCalculator.SetMaximumPhi(fEcoMug->GetMaximumPhi() - kPi);
    rateCalculator.GetAverageGenRateAndError(sourceRatePerAreaHzM2,
                                              sourceRatePerAreaErrorHzM2,
                                              rateIntegrationPoints);
    fRunAction->SetSourceFluxMetadata(sourceRatePerAreaHzM2/(EMUnits::hertz/EMUnits::m2),
                                      sourceRatePerAreaErrorHzM2/(EMUnits::hertz/EMUnits::m2),
                                      fEcoMug->GetGenSurfaceArea()/(m*m),
                                      rateIntegrationPoints);

    G4cout << "CONDITIONAL_PRIMARY_CONFIGURATION source_upstream_phi_deg=["
           << sourcePhiMinDeg << "," << sourcePhiMaxDeg << "] source_elevation_m=["
           << sourceElevationMin/m << "," << sourceElevationMax/m << "] in2_yz_m=["
           << fIn2.y/m << "," << fIn2.z/m << "] momentum_GeV=["
           << minimumMomentumGeV << "," << maximumMomentumGeV
           << "] max_attempts=" << fMaximumAttempts << G4endl;
    G4cout << "SOURCE_FLUX_NORMALIZATION horizontal_rate_hz_m2="
           << horizontalRateHzM2 << " source_rate_hz="
           << sourceRatePerAreaHzM2*fEcoMug->GetGenSurfaceArea()/(m*m)
           << " integration_points=" << rateIntegrationPoints << G4endl;
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fParticleGun;
    delete fEcoMug;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    G4long trials = 0;
    bool accepted = false;
    G4ThreeVector acceptedPosition, acceptedDirection, acceptedIn1Hit, acceptedIn2Hit;
    G4double acceptedMomentumGeV = 0.0;

    const auto intersectsPlane = [](const FlatPlane& plane,
                                    const G4ThreeVector& position,
                                    const G4ThreeVector& direction,
                                    G4ThreeVector& hit) {
        if (direction.y() <= 1.e-12) return false;
        const G4double parameter = (plane.y - position.y())/direction.y();
        if (parameter <= 0.0) return false;
        hit = position + parameter*direction;
        return std::abs(hit.x()) <= plane.width/2.0
            && std::abs(hit.z() - plane.z) <= plane.height/2.0;
    };

    while (trials < fMaximumAttempts && !accepted) {
        fEcoMug->Generate();
        ++trials;
        std::array<double, 3> position = fEcoMug->GetGenerationPosition();
        // Convert local source elevation to global Z using the ground height
        // at this exact sampled Y coordinate: z_ground(y)=-slope*y.
        position[2] += -fGroundSlope*position[1];
        const G4double theta = fEcoMug->GetGenerationTheta();
        const G4double phi = fEcoMug->GetGenerationPhi();
        const G4ThreeVector direction(std::sin(theta)*std::cos(phi),
                                      std::sin(theta)*std::sin(phi),
                                      -std::abs(std::cos(theta)));
        const G4ThreeVector start(position[0], position[1], position[2]);
        G4ThreeVector in1Hit, in2Hit;
        //if (!intersectsPlane(fIn1, start, direction, in1Hit)
        //    || !intersectsPlane(fIn2, start, direction, in2Hit)) continue;
        if (!intersectsPlane(fIn2, start, direction, in2Hit)) continue;

        // Do not test downstream planes: that would preferentially remove
        // the large-scattering tracks which carry the tomography signal.
        const G4double a = direction.x()*direction.x() + direction.y()*direction.y();
        const G4double b = 2.0*(start.x()*direction.x() + start.y()*direction.y());
        const G4double c = start.x()*start.x() + start.y()*start.y()
            - kMoundRadius*kMoundRadius;
        const G4double discriminant = b*b - 4.0*a*c;
        if (a <= 1.e-18 || discriminant <= 0.0) continue;
        const G4double entry = (-b - std::sqrt(discriminant))/(2.0*a);
        const G4double exit = (-b + std::sqrt(discriminant))/(2.0*a);
        if (entry <= 0.0 || exit <= 0.0) continue;

        acceptedPosition = start;
        acceptedDirection = direction;
        acceptedIn1Hit = in1Hit;
        acceptedIn2Hit = in2Hit;
        acceptedMomentumGeV = fEcoMug->GetGenerationMomentum();
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
    fParticleGun->SetParticleEnergy(acceptedMomentumGeV*GeV);
    if (std::getenv("MOUND_PRINT_PRIMARY_DIRECTIONS") != nullptr) {
        G4cout << "CONDITIONAL_PRIMARY event=" << event->GetEventID()
               << " position_m=(" << acceptedPosition.x()/m << ","
               << acceptedPosition.y()/m << "," << acceptedPosition.z()/m
               << ") direction=(" << acceptedDirection.x() << ","
               << acceptedDirection.y() << "," << acceptedDirection.z()
               << ") in1_hit_m=(" << acceptedIn1Hit.x()/m << ","
               << acceptedIn1Hit.y()/m << "," << acceptedIn1Hit.z()/m
               << ") in2_hit_m=(" << acceptedIn2Hit.x()/m << ","
               << acceptedIn2Hit.y()/m << "," << acceptedIn2Hit.z()/m
               << ") momentum_GeV=" << acceptedMomentumGeV
               << " trials=" << trials << G4endl;
    }

    fParticleGun->GeneratePrimaryVertex(event);
    fRunAction->RecordPrimaryGeneration(event->GetEventID(), trials, true);
    G4AnalysisManager::Instance()->FillNtupleDColumn(0, trials);
}
