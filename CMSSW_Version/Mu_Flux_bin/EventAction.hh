#ifndef EventAction_h
#define EventAction_h 1
#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
// This class keeps the three detector notes for one muon, then saves them.
class EventAction : public G4UserEventAction {
public:
  EventAction() = default;
  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;
  // Remember where a muon first enters one GEM board.
  void RecordGEMHit(int plane,const G4ThreeVector&,const G4ThreeVector&);
private:
  // truth = the exact answer; reco = what a real, imperfect detector says.
  struct Hit {
    G4ThreeVector pos;
    G4ThreeVector mom;
    G4ThreeVector reco;
    G4bool hit = false;
  };
  Hit fHit[3];
};
#endif
