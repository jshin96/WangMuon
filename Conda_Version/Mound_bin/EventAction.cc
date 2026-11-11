#include <limits>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <string>
#include "EventAction.hh"
#include "RunAction.hh"
#include "G4Event.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

namespace {
// GEM RESPONSE CONFIGURATION
//
// Every setting below can be overridden per job with `export NAME=value`.
// Invalid, missing, or non-positive values fall back to the stated default.
// Tune these values to laboratory/test-beam calibration for the actual chamber
// and front-end; they are intentionally detector-model inputs, not physics
// constants.
//
// GEM_EFFICIENCY (0.97): Probability that a geometrical truth crossing is
//   reconstructable before the charge threshold. It represents dead area,
//   dead channels, and non-charge-related reconstruction losses.
// GEM_W_VALUE_EV (26): Mean energy in eV required to make one primary
//   electron in the Ar/CO2 gas. It converts Geant4 Edep into mean ionisation.
// GEM_FANO_FACTOR (0.20): Primary-ionisation variance / mean. Larger values
//   broaden the number of primary electrons and therefore pulse charge.
// GEM_MEAN_GAIN (10000): Mean triple-GEM avalanche multiplication factor.
// GEM_RELATIVE_GAIN_RMS (0.50): RMS(gain)/mean(gain), used for the log-normal
//   Polya-like avalanche fluctuation. Larger values broaden pulse charge.
// GEM_ENC_ELECTRONS (1000): Equivalent noise charge, RMS, in electrons at the
//   front-end input. It is added as zero-mean Gaussian electronic noise.
// GEM_THRESHOLD_FC (3): Discriminator threshold in fC. A hit below it has
//   GEM*_Valid=0; increasing it reduces efficiency and worsens low-SNR timing.
// GEM_STRIP_PITCH_MM (0.40): Readout-strip/pad pitch. It sets the finite
//   segmentation term in the two in-plane position uncertainties.
// GEM_INTRINSIC_POSITION_UM (70): High-charge single-coordinate RMS limit
//   from diffusion, charge sharing, and centroid reconstruction.
// GEM_INTRINSIC_TIME_NS (4): High-charge RMS timing limit after calibration.
// GEM_THRESHOLD_TIME_JITTER_NS (20): Low-charge discriminator time-walk/jitter
//   scale. It is divided by charge/threshold and added in quadrature.
//
// Output: GEM*_Edep_keV is unsmeared Geant4 gas deposition; GEM*_Charge_fC is
// the noisy, amplified charge; GEM*_{X,Y,Z,Time_ns} are only populated when
// GEM*_Valid=1.  Position smearing is applied in the two detector-plane axes,
// never through the detector normal.
constexpr G4double kGEMEfficiency = 0.96;       // dead area / reconstruction loss
constexpr G4double kGasWValue = 26.0*eV;        // Ar/CO2 70/30 effective W value
constexpr G4double kFanoFactor = 0.20;
constexpr G4double kMeanGain = 1.0e4;
constexpr G4double kRelativeGainRms = 0.50;     // Polya-like avalanche width
constexpr G4double kElectronChargeFC = 1.602176634e-4; // fC per electron
constexpr G4double kENC = 1000.0;               // electrons RMS
constexpr G4double kThresholdFC = 0.0001;          // discriminator threshold
constexpr G4double kStripPitch = 0.40*mm;
constexpr G4double kIntrinsicPositionResolution = 1.000*mm;
constexpr G4double kIntrinsicTimeResolution = 4.0*ns;
constexpr G4double kThresholdTimeJitter = 20.0*ns;

G4double ReadResponseSetting(const char* name, G4double fallback) {
    const char* text = std::getenv(name);
    if (!text) return fallback;
    char* end = nullptr;
    const double value = std::strtod(text, &end);
    return (end != text && *end == '\0' && std::isfinite(value) && value > 0.0)
        ? value : fallback;
}

G4double ReadFiniteSetting(const char* name, G4double fallback) {
    const char* text = std::getenv(name);
    if (!text) return fallback;
    char* end = nullptr;
    const double value = std::strtod(text, &end);
    return (end != text && *end == '\0' && std::isfinite(value))
        ? value : fallback;
}

bool ReadBooleanSetting(const char* name, bool fallback) {
    const char* text = std::getenv(name);
    if (!text) return fallback;
    std::string value(text);
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (value == "1" || value == "true" || value == "yes" || value == "on") return true;
    if (value == "0" || value == "false" || value == "no" || value == "off") return false;
    return fallback;
}

bool ReadNeutronMode() {
    const char* text = std::getenv("MOUND_PRIMARY_MODE");
    if (!text) return false;
    std::string value(text);
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value == "neutron";
}

void FillTruth(G4AnalysisManager* manager, G4int& column,
               const EventAction::Hit& hit) {
    manager->FillNtupleDColumn(column++, hit.pos.x());
    manager->FillNtupleDColumn(column++, hit.pos.y());
    manager->FillNtupleDColumn(column++, hit.pos.z());
    manager->FillNtupleDColumn(column++, hit.mom.x());
    manager->FillNtupleDColumn(column++, hit.mom.y());
    manager->FillNtupleDColumn(column++, hit.mom.z());
}

G4double IncomingZenithAngleDegrees(const EventAction::Hit& hit) {
    const G4double momentumMagnitude = hit.mom.mag();
    if (momentumMagnitude <= 0.0) {
        return std::numeric_limits<G4double>::quiet_NaN();
    }
    // Cosmic-muon zenith is measured from the downward vertical, -Z.
    const G4double cosine = std::max(-1.0, std::min(1.0,
        -hit.mom.z()/momentumMagnitude));
    return std::acos(cosine)/deg;
}
}

EventAction::EventAction(RunAction* runAction)
    : G4UserEventAction(), fRunAction(runAction),
      fNeutronMode(ReadNeutronMode()),
      fTopNormalizationMode(ReadBooleanSetting("MOUND_TOP_NORMALIZATION_MODE", false)),
      fHitMound(false), fHitRoom(false), fHitRoomInterior(false) {
    // The side GEMs are perpendicular to the inclined terrain. Their local X
    // axis remains global X while their second measured axis follows the
    // terrain normal. The top monitor is horizontal and measures global X/Y.
    const G4double incline = ReadFiniteSetting("MOUND_GROUND_INCLINE_DEG", 15.0)*deg;
    fGEMTangentialAxis = G4ThreeVector(1., 0., 0.);
    fGEMAxialAxis = G4ThreeVector(0., std::sin(incline), std::cos(incline));
    fTopFirstAxis = G4ThreeVector(1., 0., 0.);
    fTopSecondAxis = G4ThreeVector(0., 1., 0.);
    Reset(fIn1); Reset(fIn2); Reset(fOut0); Reset(fOut1); Reset(fOut2); Reset(fOut3); Reset(fTop);
}

EventAction::~EventAction() {}

void EventAction::Reset(Hit& hit) {
    const G4double missing = std::numeric_limits<G4double>::quiet_NaN();
    hit.pos = hit.mom = hit.gemPos = G4ThreeVector(missing, missing, missing);
    hit.time = hit.kineticEnergy = hit.gemTime = hit.gemCharge = missing;
    hit.energy = 0.0;
    hit.hit = hit.gemValid = false;
}

void EventAction::SetTruth(Hit& hit, G4ThreeVector pos, G4ThreeVector mom,
                           G4double time, G4double kineticEnergy) {
    if (!hit.hit) {
        hit.pos = pos;
        hit.mom = mom;
        hit.time = time;
        hit.kineticEnergy = kineticEnergy;
        hit.hit = true;
    }
}

void EventAction::SetHitIn1(G4ThreeVector p, G4ThreeVector m, G4double t, G4double e) { SetTruth(fIn1, p, m, t, e); }
void EventAction::SetHitIn2(G4ThreeVector p, G4ThreeVector m, G4double t, G4double e) { SetTruth(fIn2, p, m, t, e); }
void EventAction::SetHitOut0(G4ThreeVector p, G4ThreeVector m, G4double t, G4double e) { SetTruth(fOut0, p, m, t, e); }
void EventAction::SetHitOut1(G4ThreeVector p, G4ThreeVector m, G4double t, G4double e) { SetTruth(fOut1, p, m, t, e); }
void EventAction::SetHitOut2(G4ThreeVector p, G4ThreeVector m, G4double t, G4double e) { SetTruth(fOut2, p, m, t, e); }
void EventAction::SetHitOut3(G4ThreeVector p, G4ThreeVector m, G4double t, G4double e) { SetTruth(fOut3, p, m, t, e); }
void EventAction::SetHitTop(G4ThreeVector p, G4ThreeVector m, G4double t, G4double e) { SetTruth(fTop, p, m, t, e); }
void EventAction::AddEnergyIn1(G4double e) { fIn1.energy += e; }
void EventAction::AddEnergyIn2(G4double e) { fIn2.energy += e; }
void EventAction::AddEnergyOut1(G4double e) { fOut1.energy += e; }
void EventAction::AddEnergyOut2(G4double e) { fOut2.energy += e; }
void EventAction::AddEnergyOut3(G4double e) { fOut3.energy += e; }
void EventAction::AddEnergyTop(G4double e) { fTop.energy += e; }

void EventAction::Digitize(Hit& hit, const G4ThreeVector& firstAxis,
                           const G4ThreeVector& secondAxis) {
    // Neutron mode is a penetration study, not a microscopic GEM-response
    // model. Treat any neutron plane crossing as an ideal detector hit.
    // The cosmic-muon response below remains byte-for-byte on its old path.
    if (fNeutronMode) {
        if (!hit.hit) return;
        hit.gemValid = true;
        hit.gemPos = hit.pos;
        hit.gemTime = hit.time;
        hit.gemCharge = 0.0;
        return;
    }
    const G4double efficiency = ReadResponseSetting("GEM_EFFICIENCY", kGEMEfficiency);
    const G4double wValue = ReadResponseSetting("GEM_W_VALUE_EV", kGasWValue/eV)*eV;
    const G4double fano = ReadResponseSetting("GEM_FANO_FACTOR", kFanoFactor);
    const G4double meanGain = ReadResponseSetting("GEM_MEAN_GAIN", kMeanGain);
    const G4double relativeGainRms = ReadResponseSetting(
        "GEM_RELATIVE_GAIN_RMS", kRelativeGainRms);
    const G4double enc = ReadResponseSetting("GEM_ENC_ELECTRONS", kENC);
    const G4double threshold = ReadResponseSetting("GEM_THRESHOLD_FC", kThresholdFC);
    const G4double stripPitch = ReadResponseSetting("GEM_STRIP_PITCH_MM",
        kStripPitch/mm)*mm;
    const G4double intrinsicPosition = ReadResponseSetting(
        "GEM_INTRINSIC_POSITION_UM", kIntrinsicPositionResolution/micrometer)*micrometer;
    const G4double intrinsicTime = ReadResponseSetting(
        "GEM_INTRINSIC_TIME_NS", kIntrinsicTimeResolution/ns)*ns;
    const G4double thresholdTimeJitter = ReadResponseSetting(
        "GEM_THRESHOLD_TIME_JITTER_NS", kThresholdTimeJitter/ns)*ns;
    if (!hit.hit || G4UniformRand() > efficiency) return;
    // Primary-ionisation statistics followed by a log-normal approximation
    // to the Polya avalanche-gain distribution.
    const G4double meanElectrons = hit.energy/wValue;
    const G4double primaryElectrons = std::max(0.0, G4RandGauss::shoot(
        meanElectrons, std::sqrt(fano*meanElectrons)));
    const G4double gainLogSigma = std::sqrt(std::log(1.0
        + relativeGainRms*relativeGainRms));
    const G4double gain = meanGain*std::exp(G4RandGauss::shoot(
        -0.5*gainLogSigma*gainLogSigma, gainLogSigma));
    const G4double chargeNoiseless = primaryElectrons*gain*kElectronChargeFC;
    hit.gemCharge = chargeNoiseless + G4RandGauss::shoot(0.0,
        enc*kElectronChargeFC);
    if (hit.gemCharge < threshold) return;

    hit.gemValid = true;
    const G4double signalToThreshold = hit.gemCharge/threshold;
    const G4double positionSigma = std::sqrt(
        intrinsicPosition*intrinsicPosition
        + std::pow(stripPitch/(std::sqrt(12.0)*signalToThreshold), 2));
    hit.gemPos = hit.pos
        + firstAxis*G4RandGauss::shoot(0.0, positionSigma)
        + secondAxis*G4RandGauss::shoot(0.0, positionSigma);
    const G4double timeSigma = std::sqrt(
        intrinsicTime*intrinsicTime
        + std::pow(thresholdTimeJitter/signalToThreshold, 2));
    hit.gemTime = G4RandGauss::shoot(hit.time, timeSigma);
}

void EventAction::BeginOfEventAction(const G4Event*) {
    Reset(fIn1); Reset(fIn2); Reset(fOut0); Reset(fOut1); Reset(fOut2); Reset(fOut3); Reset(fTop);
    fHitMound = false; fHitRoom = false; fHitRoomInterior = false;
}

void EventAction::EndOfEventAction(const G4Event* event) {
    // Side mode preserves the established IN2 truth trigger. Top mode records
    // only primaries crossing the horizontal normalization detector.
//    if (!(fIn2.hit && fOut1.hit && fOut2.hit && fOut3.hit)) return;
    if (fTopNormalizationMode ? !fTop.hit : !fIn2.hit) return;
    Digitize(fIn1, fGEMTangentialAxis, fGEMAxialAxis);
    Digitize(fIn2, fGEMTangentialAxis, fGEMAxialAxis);
    Digitize(fOut1, fGEMTangentialAxis, fGEMAxialAxis);
    Digitize(fOut2, fGEMTangentialAxis, fGEMAxialAxis);
    Digitize(fOut3, fGEMTangentialAxis, fGEMAxialAxis);
    Digitize(fTop, fTopFirstAxis, fTopSecondAxis);

    G4int trajectoryFlag = fHitRoom ? 2 : (fHitMound ? 1 : 0);
    auto* manager = G4AnalysisManager::Instance();
    G4int column = 1;
    manager->FillNtupleIColumn(column++, event->GetEventID());
    FillTruth(manager, column, fIn1); FillTruth(manager, column, fIn2);
    manager->FillNtupleDColumn(column++, IncomingZenithAngleDegrees(fIn2));
    // Out0 is truth-only: it neither receives GEM digitisation nor participates
    // in the event trigger or the five-GEM efficiency accounting.
    manager->FillNtupleDColumn(column++, fOut0.pos.x());
    manager->FillNtupleDColumn(column++, fOut0.pos.y());
    manager->FillNtupleDColumn(column++, fOut0.pos.z());
    FillTruth(manager, column, fOut1); FillTruth(manager, column, fOut2);
    FillTruth(manager, column, fOut3);
    manager->FillNtupleIColumn(column++, trajectoryFlag);
    manager->FillNtupleIColumn(column++, fOut1.hit ? 1 : 0);
    manager->FillNtupleIColumn(column++, fOut2.hit ? 1 : 0);
    manager->FillNtupleIColumn(column++, fOut3.hit ? 1 : 0);

    const Hit* hits[] = {&fIn1, &fIn2, &fOut1, &fOut2, &fOut3};
    for (const Hit* hit : hits) {
        manager->FillNtupleDColumn(column++, hit->gemPos.x());
        manager->FillNtupleDColumn(column++, hit->gemPos.y());
        manager->FillNtupleDColumn(column++, hit->gemPos.z());
        manager->FillNtupleDColumn(column++, hit->gemTime/ns);
        manager->FillNtupleDColumn(column++, hit->energy/keV);
        manager->FillNtupleDColumn(column++, hit->gemCharge);
        manager->FillNtupleIColumn(column++, hit->gemValid ? 1 : 0);
    }
    FillTruth(manager, column, fTop);
    manager->FillNtupleDColumn(column++, IncomingZenithAngleDegrees(fTop));
    manager->FillNtupleIColumn(column++, fTop.hit ? 1 : 0);
    manager->FillNtupleIColumn(column++, fHitRoomInterior ? 1 : 0);
    manager->FillNtupleDColumn(column++, fTop.gemPos.x());
    manager->FillNtupleDColumn(column++, fTop.gemPos.y());
    manager->FillNtupleDColumn(column++, fTop.gemPos.z());
    manager->FillNtupleDColumn(column++, fTop.gemTime/ns);
    manager->FillNtupleDColumn(column++, fTop.energy/keV);
    manager->FillNtupleDColumn(column++, fTop.gemCharge);
    manager->FillNtupleIColumn(column++, fTop.gemValid ? 1 : 0);
    manager->FillNtupleIColumn(column++, fNeutronMode ? 2112 : 13);
    manager->FillNtupleDColumn(column++, fIn2.kineticEnergy/GeV);
    manager->FillNtupleDColumn(column++, fOut1.kineticEnergy/GeV);
    //if (fIn1.gemValid && fIn2.gemValid && fOut1.gemValid
    //    && fOut2.gemValid && fOut3.gemValid) {
    if (fIn2.gemValid && fOut1.gemValid
        && fOut2.gemValid && fOut3.gemValid) {
        fRunAction->RecordFiveGEMEvent();
    }
    if (fTop.hit) {
        fRunAction->RecordTopTruthEvent();
        if (fTop.gemValid) fRunAction->RecordTopGEMEvent();
        if (fHitRoomInterior) fRunAction->RecordTopRoomEvent();
    }
    manager->AddNtupleRow();
}
