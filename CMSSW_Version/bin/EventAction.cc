#include "EventAction.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"

EventAction::EventAction()
: G4UserEventAction(),
  fPosIn(0.), fMomIn(0.), fPosOut(0.), fMomOut(0.),
  fHitIn(false), fHitOut(false)
{}

EventAction::~EventAction()
{}

void EventAction::BeginOfEventAction(const G4Event*)
{
    // Reset flags for the new muon
    fHitIn = false;
    fHitOut = false;
}

void EventAction::EndOfEventAction(const G4Event* event)
{
    // We only care about coincidence events (muons that didn't stop inside the mound)
    if (fHitIn && fHitOut) {
        
        // Get analysis manager
        auto analysisManager = G4AnalysisManager::Instance();

        // Fill the Ntuple (assuming ID 0)
        analysisManager->FillNtupleIColumn(0, event->GetEventID());
        
        // --- INCOMING DETECTOR (-Y Axis) ---
        analysisManager->FillNtupleDColumn(1, fPosIn.x());
        analysisManager->FillNtupleDColumn(2, fPosIn.y());
        analysisManager->FillNtupleDColumn(3, fPosIn.z());
        analysisManager->FillNtupleDColumn(4, fMomIn.x());
        analysisManager->FillNtupleDColumn(5, fMomIn.y());
        analysisManager->FillNtupleDColumn(6, fMomIn.z());

        // --- OUTGOING DETECTOR (+Y Axis) ---
        analysisManager->FillNtupleDColumn(7, fPosOut.x());
        analysisManager->FillNtupleDColumn(8, fPosOut.y());
        analysisManager->FillNtupleDColumn(9, fPosOut.z());
        analysisManager->FillNtupleDColumn(10, fMomOut.x());
        analysisManager->FillNtupleDColumn(11, fMomOut.y());
        analysisManager->FillNtupleDColumn(12, fMomOut.z());

        // Commit the row to the file
        analysisManager->AddNtupleRow();
    }
}
