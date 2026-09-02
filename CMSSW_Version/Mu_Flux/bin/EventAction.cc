#include "EventAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
namespace {
// The detector's small "wobble" and a helper for measuring turns.
G4double Resolution() {
  const char* text = std::getenv("GEM_INTRINSIC_POSITION_UM");
  if (!text) return 1000.0 * micrometer;
  char* end = nullptr;
  const auto value = std::strtod(text, &end);
  return (end != text && *end == '\0' && value > 0.0 && std::isfinite(value))
      ? value * micrometer : 1000.0 * micrometer;
}

G4double Angle(const G4ThreeVector& first, const G4ThreeVector& second) {
  if (first.mag2() == 0.0 || second.mag2() == 0.0) {
    return std::numeric_limits<G4double>::quiet_NaN();
  }
  return std::atan2(first.cross(second).mag(), first.dot(second));
}

// The detector construction uses this same job setting for a uniform +Y field.
// The sign is retained here so a future signed field setting remains valid;
// momentum magnitude only depends on |B|.
G4double MagneticFieldY() {
  const char* text = std::getenv("MUON_HELMHOLTZ_FIELD_T");
  if (!text) return 1.5 * tesla;
  char* end = nullptr;
  const auto value = std::strtod(text, &end);
  return (end != text && *end == '\0' && std::isfinite(value))
      ? value * tesla : std::numeric_limits<G4double>::quiet_NaN();
}

G4double MinimumBendSignificance() {
  const char* text = std::getenv("MUON_RECO_MIN_BEND_SIGMA");
  if (!text) return 3.0;
  char* end = nullptr;
  const auto value = std::strtod(text, &end);
  return (end != text && *end == '\0' && std::isfinite(value) && value >= 0.0)
      ? value : std::numeric_limits<G4double>::quiet_NaN();
}

G4double FieldAperture() {
  const char* text = std::getenv("MUON_HELMHOLTZ_APERTURE_CM");
  if (!text) return 30.0 * cm;
  char* end = nullptr;
  const auto value = std::strtod(text, &end);
  return (end != text && *end == '\0' && std::isfinite(value) && value > 0.0)
      ? value * cm : std::numeric_limits<G4double>::quiet_NaN();
}

G4double FieldHalfLength() {
  const char* text = std::getenv("MUON_HELMHOLTZ_LENGTH_M");
  if (!text) return 0.50 * m / 2.0;
  char* end = nullptr;
  const auto value = std::strtod(text, &end);
  return (end != text && *end == '\0' && std::isfinite(value) && value > 0.0)
      ? value * m / 2.0 : std::numeric_limits<G4double>::quiet_NaN();
}

// Find the interval where the unbent incoming ray lies in the finite field
// cylinder.  This is the field path used by the transport model below.
bool FieldSegment(const G4ThreeVector& position, const G4ThreeVector& direction,
                  G4double maximumPath, G4double& start, G4double& length) {
  const auto radius = FieldAperture();
  const auto halfLength = FieldHalfLength();
  if (!std::isfinite(radius) || !std::isfinite(halfLength)
      || maximumPath <= 0.0) return false;

  const auto transverseSpeed = direction.x() * direction.x()
      + direction.z() * direction.z();
  G4double radialLow = -std::numeric_limits<G4double>::infinity();
  G4double radialHigh = std::numeric_limits<G4double>::infinity();
  if (transverseSpeed <= 1.0e-16) {
    if (position.x() * position.x() + position.z() * position.z()
        > radius * radius) return false;
  } else {
    const auto linear = position.x() * direction.x()
        + position.z() * direction.z();
    const auto constant = position.x() * position.x() + position.z() * position.z()
        - radius * radius;
    const auto discriminant = linear * linear - transverseSpeed * constant;
    if (discriminant < 0.0) return false;
    const auto root = std::sqrt(discriminant);
    radialLow = (-linear - root) / transverseSpeed;
    radialHigh = (-linear + root) / transverseSpeed;
  }

  G4double axialLow = -std::numeric_limits<G4double>::infinity();
  G4double axialHigh = std::numeric_limits<G4double>::infinity();
  if (std::abs(direction.y()) <= 1.0e-12) {
    if (std::abs(position.y()) > halfLength) return false;
  } else {
    const auto first = (-halfLength - position.y()) / direction.y();
    const auto second = (halfLength - position.y()) / direction.y();
    axialLow = std::min(first, second);
    axialHigh = std::max(first, second);
  }

  start = std::max({0.0, radialLow, axialLow});
  const auto end = std::min({maximumPath, radialHigh, axialHigh});
  length = end - start;
  return length > 1.0e-9 * mm;
}

// Transport an incoming track through the known finite B_y cylinder.  The
// signed inverse momentum chooses the bend sign, so no truth charge is used.
G4double PredictedLateralBend(const G4ThreeVector& incomingDirection,
                              const G4ThreeVector& incomingPosition,
                              const G4ThreeVector& tangent,
                              const G4ThreeVector& bendingNormal,
                              G4double targetLongitudinal,
                              G4double signedInverseMomentum) {
  const auto nan = std::numeric_limits<G4double>::quiet_NaN();
  const auto fieldY = MagneticFieldY();
  if (!std::isfinite(fieldY) || std::abs(fieldY) <= 1.0e-12 * tesla) return nan;
  const auto longitudinalSpeed = incomingDirection.dot(tangent);
  if (longitudinalSpeed <= 1.0e-12) return nan;
  const auto maximumPath = targetLongitudinal / longitudinalSpeed;
  G4double fieldStart = 0.0;
  G4double fieldPath = 0.0;
  if (!FieldSegment(incomingPosition, incomingDirection, maximumPath,
                    fieldStart, fieldPath)) return nan;

  const G4ThreeVector fieldDirection(0.0, fieldY > 0.0 ? 1.0 : -1.0, 0.0);
  const auto entryPosition = incomingPosition + fieldStart * incomingDirection;
  const auto parallelDirection = incomingDirection.dot(fieldDirection) * fieldDirection;
  const auto perpendicularDirection = incomingDirection - parallelDirection;
  const auto turnDirection = perpendicularDirection.cross(fieldDirection);
  constexpr G4double kGeVPerTeslaMeter = 0.299792458;
  const auto angularRate = signedInverseMomentum * kGeVPerTeslaMeter
      * std::abs(fieldY / tesla) / m;

  G4ThreeVector exitPosition;
  G4ThreeVector exitDirection;
  if (std::abs(angularRate) <= 1.0e-16 / m) {
    exitPosition = entryPosition + fieldPath * incomingDirection;
    exitDirection = incomingDirection;
  } else {
    const auto turn = angularRate * fieldPath;
    exitPosition = entryPosition + parallelDirection * fieldPath
        + perpendicularDirection * (std::sin(turn) / angularRate)
        + turnDirection * ((1.0 - std::cos(turn)) / angularRate);
    exitDirection = parallelDirection + perpendicularDirection * std::cos(turn)
        + turnDirection * std::sin(turn);
  }

  const auto usedLongitudinal = (exitPosition - incomingPosition).dot(tangent);
  const auto exitLongitudinalSpeed = exitDirection.dot(tangent);
  if (exitLongitudinalSpeed <= 1.0e-12) return nan;
  const auto driftPath = (targetLongitudinal - usedLongitudinal)
      / exitLongitudinalSpeed;
  if (driftPath < 0.0) return nan;
  const auto predictedPosition = exitPosition + driftPath * exitDirection;
  return (predictedPosition - incomingPosition).dot(bendingNormal);
}

// Fit |q/p| to the GEMOut lateral residual after transport through the finite
// field.  The bend sign determines the charge sign; only |p| is stored.
G4double ReconstructedMomentum(const G4ThreeVector& incomingDirection,
                               const G4ThreeVector& upstreamPosition,
                               const G4ThreeVector& incomingPosition,
                               const G4ThreeVector& outgoingPosition) {
  const auto fieldY = MagneticFieldY();
  const auto nan = std::numeric_limits<G4double>::quiet_NaN();
  if (!std::isfinite(fieldY) || std::abs(fieldY) <= 1.0e-12 * tesla) return nan;

  const G4ThreeVector fieldDirection(0.0, fieldY > 0.0 ? 1.0 : -1.0, 0.0);
  const auto transverseDirection =
      incomingDirection - incomingDirection.dot(fieldDirection) * fieldDirection;
  const auto transverseMagnitude = transverseDirection.mag();
  if (transverseMagnitude <= 1.0e-12) return nan;
  const auto tangent = transverseDirection / transverseMagnitude;

  const auto displacement = outgoingPosition - incomingPosition;
  const auto chord = displacement - displacement.dot(fieldDirection) * fieldDirection;
  const auto transverseChord = chord.dot(tangent.cross(fieldDirection));
  const auto upstreamDisplacement = incomingPosition - upstreamPosition;
  const auto upstreamChord = upstreamDisplacement
      - upstreamDisplacement.dot(fieldDirection) * fieldDirection;
  const auto upstreamLeverArm = std::abs(upstreamChord.dot(tangent));
  if (upstreamLeverArm <= 1.0e-12 * mm) return nan;

  // The bend is a residual relative to a slope measured with GEMIn1/In2.
  // Propagate the independent position errors from all three GEMs into that
  // residual.  X is measured exactly in this model; for B || Y the relevant
  // detector uncertainty is its projection on the X-Z bending normal.
  const auto outgoingLeverArm = std::abs(chord.dot(tangent));
  const auto leverArmRatio = outgoingLeverArm / upstreamLeverArm;
  const auto bendingNormal = tangent.cross(fieldDirection);
  const auto lateralResolution = std::abs(bendingNormal.z()) * Resolution();
  const auto bendUncertainty = lateralResolution * std::sqrt(
      1.0 + (1.0 + leverArmRatio) * (1.0 + leverArmRatio)
          + leverArmRatio * leverArmRatio);
  const auto minimumSignificance = MinimumBendSignificance();
  if (!std::isfinite(minimumSignificance)
      || std::abs(transverseChord) < minimumSignificance * bendUncertainty) {
    return nan;
  }

  // Find the inverse momentum for which the finite-field prediction matches
  // the observed lateral bend.  The scan brackets the first physical root,
  // then bisection gives a stable one-dimensional fit.
  const auto targetLongitudinal = chord.dot(tangent);
  if (targetLongitudinal <= 1.0e-12 * mm) return nan;
  const auto observedMagnitude = std::abs(transverseChord);
  const auto bendSign = transverseChord >= 0.0 ? 1.0 : -1.0;
  auto residual = [&](G4double inverseMomentum) {
    const auto predicted = PredictedLateralBend(
        incomingDirection, incomingPosition, tangent, bendingNormal,
        targetLongitudinal, bendSign * inverseMomentum);
    return std::isfinite(predicted) ? bendSign * predicted - observedMagnitude : nan;
  };

  G4double lower = 0.0;
  G4double upper = 1.0e-6;
  G4double upperResidual = residual(upper);
  while (std::isfinite(upperResidual) && upperResidual < 0.0 && upper < 10.0) {
    lower = upper;
    upper *= 2.0;
    upperResidual = residual(upper);
  }
  if (!std::isfinite(upperResidual) || upperResidual < 0.0) return nan;
  for (int iteration = 0; iteration < 48; ++iteration) {
    const auto middle = 0.5 * (lower + upper);
    const auto middleResidual = residual(middle);
    if (!std::isfinite(middleResidual)) return nan;
    if (middleResidual < 0.0) lower = middle;
    else upper = middle;
  }
  return GeV / (0.5 * (lower + upper));
}
}  // namespace

void EventAction::BeginOfEventAction(const G4Event*) {
  // Start every muon with an empty notebook.
  const auto nan = std::numeric_limits<G4double>::quiet_NaN();
  for (auto& hit : fHit) {
    hit.pos = hit.mom = hit.reco = {nan, nan, nan};
    hit.hit = false;
  }
}

void EventAction::RecordGEMHit(int plane, const G4ThreeVector& position,
                               const G4ThreeVector& momentum) {
  // Keep the first crossing only.  Add random Y/Z fuzz to imitate a GEM.
  if (plane < 0 || plane > 2 || fHit[plane].hit) return;
  fHit[plane].pos = position;
  fHit[plane].mom = momentum;
  fHit[plane].reco = {position.x(),
                      position.y() + G4RandGauss::shoot(0.0, Resolution()),
                      position.z() + G4RandGauss::shoot(0.0, Resolution())};
  fHit[plane].hit = true;
}
void EventAction::EndOfEventAction(const G4Event* e) {
  // Save every primary muon that reaches GEMIn1.  Missing downstream hits
  // remain NaN, so the fixed tree schema also records partial trajectories.
  if (!fHit[0].hit) return;
  const auto nan = std::numeric_limits<G4double>::quiet_NaN();
  const bool hasTruthDeflection = fHit[1].hit && fHit[2].hit;
  const bool hasAllHits = fHit[0].hit && fHit[1].hit && fHit[2].hit;
  G4double truthDeflection = nan;
  G4double recoDeflection = nan;
  G4double recoMomentum = nan;
  G4double truthMomentum = nan;

  if (hasTruthDeflection) {
    truthDeflection = Angle(fHit[1].mom, fHit[2].mom);
  }
  if (fHit[1].hit) {
    truthMomentum = fHit[1].mom.mag() / GeV;
  }
  if (hasAllHits) {
    const auto recoIn = (fHit[1].reco - fHit[0].reco).unit();
    // With one post-field plane, the observable is the deflection chord from
    // the near incoming GEM to GEMOut.
    const auto recoOut = (fHit[2].reco - fHit[1].reco).unit();
    recoDeflection = Angle(recoIn, recoOut);
    recoMomentum = ReconstructedMomentum(
        recoIn, fHit[0].reco, fHit[1].reco, fHit[2].reco) / GeV;
  }

  // Put one complete or partial muon story into one row of the ROOT file.
  auto* a = G4AnalysisManager::Instance();
  int c = 0;
  a->FillNtupleIColumn(c++,e->GetEventID());
  for (const auto& hit : fHit) {
    a->FillNtupleDColumn(c++, hit.pos.x()/mm);
    a->FillNtupleDColumn(c++, hit.pos.y()/mm);
    a->FillNtupleDColumn(c++, hit.pos.z()/mm);
    a->FillNtupleDColumn(c++, hit.mom.x()/GeV);
    a->FillNtupleDColumn(c++, hit.mom.y()/GeV);
    a->FillNtupleDColumn(c++, hit.mom.z()/GeV);
    a->FillNtupleDColumn(c++, hit.reco.x()/mm);
    a->FillNtupleDColumn(c++, hit.reco.y()/mm);
    a->FillNtupleDColumn(c++, hit.reco.z()/mm);
  }
  a->FillNtupleDColumn(c++, truthDeflection);
  a->FillNtupleDColumn(c++, recoDeflection);
  a->FillNtupleDColumn(c++, recoMomentum);
  a->FillNtupleDColumn(c++, truthMomentum);
  a->FillNtupleDColumn(c++, hasAllHits ? truthMomentum - recoMomentum : nan);
  a->AddNtupleRow();
}
