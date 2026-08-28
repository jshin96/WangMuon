#include "PrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4Exception.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cmath>
#include <cstdlib>
#include <string>
constexpr G4double kPi = CLHEP::pi;
namespace {
// Read a safe number from the job settings.
G4double Env(const char* name, G4double fallback) {
  const char* text = std::getenv(name);
  if (!text) return fallback;
  try {
    std::size_t used = 0;
    const auto value = std::stod(text, &used);
    if (used != std::string(text).size() || !std::isfinite(value)) throw 0;
    return value;
  } catch (...) {
    G4Exception("PrimaryGeneratorAction", "InvalidBeam", FatalException,
                (std::string(name) + " must be finite.").c_str());
    return fallback;
  }
}
}  // namespace

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : fGun(new G4ParticleGun(1)),
      fSourceX(Env("MUON_BEAM_SOURCE_X_M", -1.0) * m),
      fHalfBeamSize(Env("MUON_BEAM_SIZE_CM", 8.0) * cm / 2),
      fMinEnergy(Env("MUON_BEAM_MIN_ENERGY_GEV", 1.0) * GeV),
      fMaxEnergy(Env("MUON_BEAM_MAX_ENERGY_GEV", 100.0) * GeV) {
      fBeamTheta(Env("MUON_BEAM_THETA_SPREAD_DEG", 5.0)),
      fBeamPhi(Env("MUON_BEAM_PHI_SPREAD_DEG", 5.0)) {
  if (fHalfBeamSize <= 0.0 || fMinEnergy <= 0.0 || fMaxEnergy < fMinEnergy) {
    G4Exception("PrimaryGeneratorAction", "InvalidBeam", FatalException,
                "Require positive beam size and 0 < min energy <= max energy.");
  }
  fGun->SetParticleDefinition(
      G4ParticleTable::GetParticleTable()->FindParticle("mu-"));
  std::float theta = fBeamTheta/180*G4UniformRand()*kPi;
  std::float phi = fBeamPhi/180*G4UniformRand()*kPi;
  fGun->SetParticleMomentumDirection({std::cos(phi)*std::cos(theta), std::sin(phi)*std::cos(theta), std::sin(theta)});
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  delete fGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event*e){
  // Pick a random little spot on the square beam face, but always shoot +X.
  const auto y = (2.0 * G4UniformRand() - 1.0) * fHalfBeamSize;
  const auto z = (2.0 * G4UniformRand() - 1.0) * fHalfBeamSize;
  fGun->SetParticlePosition({fSourceX, y, z});
  // Pick one energy between 1 and 100 GeV (unless the job script changes it).
  fGun->SetParticleEnergy(fMinEnergy + G4UniformRand() * (fMaxEnergy - fMinEnergy));
  fGun->GeneratePrimaryVertex(e);
}
