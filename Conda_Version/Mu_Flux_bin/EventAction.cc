#include "EventAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cmath>
#include <cstdlib>
#include <limits>
namespace {
// The detector's small "wobble" and a helper for measuring turns.
G4double Resolution() {
  const char* text = std::getenv("GEM_INTRINSIC_POSITION_UM");
  if (!text) return 100.0 * micrometer;
  char* end = nullptr;
  const auto value = std::strtod(text, &end);
  return (end != text && *end == '\0' && value > 0.0 && std::isfinite(value))
      ? value * micrometer : 100.0 * micrometer;
}

G4double Angle(const G4ThreeVector& first, const G4ThreeVector& second) {
  if (first.mag2() == 0.0 || second.mag2() == 0.0) {
    return std::numeric_limits<G4double>::quiet_NaN();
  }
  return std::atan2(first.cross(second).mag(), first.dot(second));
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
  // A useful answer needs all three boards to see the same muon.
  if(!(fHit[0].hit&&fHit[1].hit&&fHit[2].hit))return;
  const auto truthIn = fHit[1].mom.unit();
  const auto truthOut = fHit[2].mom.unit();
  const auto recoIn = (fHit[1].reco - fHit[0].reco).unit();
  // With one post-field plane, the observable is the deflection chord from
  // the near incoming GEM to GEMOut.  The truth angle uses local momenta.
  const auto recoOut = (fHit[2].reco - fHit[1].reco).unit();
  // Put one complete muon story into one row of the ROOT file.
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
  a->FillNtupleDColumn(c++, Angle(truthIn, truthOut));
  a->FillNtupleDColumn(c++, Angle(recoIn, recoOut));
  a->AddNtupleRow();
}
