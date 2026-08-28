#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4MuonMinus.hh"
#include "G4Step.hh"
#include "G4Track.hh"
SteppingAction::SteppingAction(EventAction* eventAction)
    : fEventAction(eventAction) {}

void SteppingAction::UserSteppingAction(const G4Step* s){
  // Ignore everything except the original negative muon when it enters a board.
  if (s->GetTrack()->GetDefinition() != G4MuonMinus::MuonMinusDefinition() ||
      s->GetTrack()->GetParentID() != 0 ||
      s->GetPreStepPoint()->GetStepStatus() != fGeomBoundary) return;
  auto* volume = s->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  if (!volume) return;
  int plane = -1;
  const auto& name = volume->GetName();
  if (name == "GEMIn1") plane = 0;
  else if (name == "GEMIn2") plane = 1;
  else if (name == "GEMOut") plane = 2;
  if (plane >= 0) {
    fEventAction->RecordGEMHit(plane, s->GetPreStepPoint()->GetPosition(),
                               s->GetPreStepPoint()->GetMomentum());
  }
}
