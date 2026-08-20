#include <limits>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include "EventAction.hh"
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

void FillTruth(G4AnalysisManager* manager, G4int& column,
               const EventAction::Hit& hit) {
    manager->FillNtupleDColumn(column++, hit.pos.x());
    manager->FillNtupleDColumn(column++, hit.pos.y());
    manager->FillNtupleDColumn(column++, hit.pos.z());
    manager->FillNtupleDColumn(column++, hit.mom.x());
    manager->FillNtupleDColumn(column++, hit.mom.y());
    manager->FillNtupleDColumn(column++, hit.mom.z());
}
}

EventAction::EventAction() : G4UserEventAction(), fHitMound(false), fHitRoom(false) {
    // Flat GEM planes have global-Y normals, so their measured in-plane axes
    // are global X and global Z regardless of terrain inclination.
    fGEMTangentialAxis = G4ThreeVector(1., 0., 0.);
    fGEMAxialAxis = G4ThreeVector(0., 0., 1.);
    Reset(fIn1); Reset(fIn2); Reset(fOut1); Reset(fOut2); Reset(fOut3);
}

EventAction::~EventAction() {}

void EventAction::Reset(Hit& hit) {
    const G4double missing = std::numeric_limits<G4double>::quiet_NaN();
    hit.pos = hit.mom = hit.gemPos = G4ThreeVector(missing, missing, missing);
    hit.time = hit.gemTime = hit.gemCharge = missing;
    hit.energy = 0.0;
    hit.hit = hit.gemValid = false;
}

void EventAction::SetTruth(Hit& hit, G4ThreeVector pos, G4ThreeVector mom,
                           G4double time) {
    if (!hit.hit) {
        hit.pos = pos; hit.mom = mom; hit.time = time; hit.hit = true;
    }
}

void EventAction::SetHitIn1(G4ThreeVector p, G4ThreeVector m, G4double t) { SetTruth(fIn1, p, m, t); }
void EventAction::SetHitIn2(G4ThreeVector p, G4ThreeVector m, G4double t) { SetTruth(fIn2, p, m, t); }
void EventAction::SetHitOut1(G4ThreeVector p, G4ThreeVector m, G4double t) { SetTruth(fOut1, p, m, t); }
void EventAction::SetHitOut2(G4ThreeVector p, G4ThreeVector m, G4double t) { SetTruth(fOut2, p, m, t); }
void EventAction::SetHitOut3(G4ThreeVector p, G4ThreeVector m, G4double t) { SetTruth(fOut3, p, m, t); }
void EventAction::AddEnergyIn1(G4double e) { fIn1.energy += e; }
void EventAction::AddEnergyIn2(G4double e) { fIn2.energy += e; }
void EventAction::AddEnergyOut1(G4double e) { fOut1.energy += e; }
void EventAction::AddEnergyOut2(G4double e) { fOut2.energy += e; }
void EventAction::AddEnergyOut3(G4double e) { fOut3.energy += e; }

void EventAction::Digitize(Hit& hit) {
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
        + fGEMTangentialAxis*G4RandGauss::shoot(0.0, positionSigma)
        + fGEMAxialAxis*G4RandGauss::shoot(0.0, positionSigma);
    const G4double timeSigma = std::sqrt(
        intrinsicTime*intrinsicTime
        + std::pow(thresholdTimeJitter/signalToThreshold, 2));
    hit.gemTime = G4RandGauss::shoot(hit.time, timeSigma);
}

void EventAction::BeginOfEventAction(const G4Event*) {
    Reset(fIn1); Reset(fIn2); Reset(fOut1); Reset(fOut2); Reset(fOut3);
    fHitMound = false; fHitRoom = false;
}

void EventAction::EndOfEventAction(const G4Event* event) {
    // Preserve the established incoming two-plane trigger; the GEM trigger
    // response is available separately through the GEM*_Valid branches.
    if (!(fIn1.hit && fIn2.hit)) return;
    Digitize(fIn1); Digitize(fIn2); Digitize(fOut1); Digitize(fOut2); Digitize(fOut3);

    G4int trajectoryFlag = fHitRoom ? 2 : (fHitMound ? 1 : 0);
    auto* manager = G4AnalysisManager::Instance();
    G4int column = 1;
    manager->FillNtupleIColumn(column++, event->GetEventID());
    FillTruth(manager, column, fIn1); FillTruth(manager, column, fIn2);
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
    manager->AddNtupleRow();
}
