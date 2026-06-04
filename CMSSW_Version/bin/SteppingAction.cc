#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4Step.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

SteppingAction::SteppingAction(EventAction* eventAction)
: fEventAction(eventAction) {}

void SteppingAction::UserSteppingAction(const G4Step* step) {
    G4StepPoint* prePoint = step->GetPreStepPoint();
    G4StepPoint* postPoint = step->GetPostStepPoint();

    if (!prePoint || !postPoint) return;

    G4VPhysicalVolume* preVolume = prePoint->GetPhysicalVolume();
    G4VPhysicalVolume* postVolume = postPoint->GetPhysicalVolume();

    if (!preVolume || !postVolume) return;
    // ========================================================
    // 1. CYLINDER SENSOR
    // ========================================================
    if (preVolume->GetName() != "PhysDetector" && postVolume->GetName() == "PhysDetector") {
        G4ThreeVector pos = postPoint->GetPosition();
        double r = pos.perp();

        if (std::abs(r - 28.5 * m) < 1.0 * mm) {
            double h_arc = 28.5 * m * pos.phi();
            double energy = prePoint->GetKineticEnergy() / GeV;
            double pid = step->GetTrack()->GetDefinition()->GetPDGEncoding();

            // Send to memory bank
	    if (std::abs(pid) == 13) {
                fEventAction->SetInnerHit(h_arc / m, pos.z() / m, energy, pid);
	    }
        }
    }

    // ========================================================
    // 2. FLAT WALL SENSOR
    // ========================================================
    if (preVolume->GetName() != "PhysFlatDetector" && postVolume->GetName() == "PhysFlatDetector") {
        G4ThreeVector pos = postPoint->GetPosition();

        if (pos.y() > 38.5 * m && pos.y() < 38.52 * m) {
            double energy = prePoint->GetKineticEnergy() / GeV;
            double pid = step->GetTrack()->GetDefinition()->GetPDGEncoding();

            // Send to memory bank
	    if (std::abs(pid) == 13) {
                fEventAction->SetOuterHit(pos.x() / m, pos.z() / m, energy, pid);
	    }
            // Safely kill at the final backstop
            step->GetTrack()->SetTrackStatus(fStopAndKill);
        }
    }
}
