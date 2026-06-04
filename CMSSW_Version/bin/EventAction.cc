#include "EventAction.hh"
#include "G4AnalysisManager.hh"

EventAction::EventAction() : G4UserEventAction() {}


void EventAction::BeginOfEventAction(const G4Event*) {
    // Reset EVERYTHING
    fPassedMound = false; fPassedWall = false; fPassedRoom = false;
    fHitInner = false; fHitOuter = false;
    fInnerH = 0; fInnerZ = 0; fInnerE = 0; fInnerPID = 0;
    fOuterH = 0; fOuterZ = 0; fOuterE = 0; fOuterPID = 0;
}

void EventAction::EndOfEventAction(const G4Event*) {
    // If NO detectors were hit, do not write anything to the file!
    if (!fHitInner && !fHitOuter) return; 

    // Determine the user's requested Trigger Category (1, 2, or 3)
    double triggerType = 0.0;
    if (fHitInner && !fHitOuter) triggerType = 1.0; // Inner Only
    if (!fHitInner && fHitOuter) triggerType = 2.0; // Outer Only
    if (fHitInner && fHitOuter)  triggerType = 3.0; // Coincidence!

    // Determine the Physics Category (Mound/Rock/Air)
    double physicsCat = (fPassedRoom || fPassedWall) ? 1.0 : (fPassedMound ? 2.0 : 3.0);

    // Write the unified row
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleDColumn(0, triggerType);
    analysisManager->FillNtupleDColumn(1, physicsCat);
    analysisManager->FillNtupleDColumn(2, fInnerH);
    analysisManager->FillNtupleDColumn(3, fInnerZ);
    analysisManager->FillNtupleDColumn(4, fInnerE);
    analysisManager->FillNtupleDColumn(5, fInnerPID);
    analysisManager->FillNtupleDColumn(6, fOuterH);
    analysisManager->FillNtupleDColumn(7, fOuterZ);
    analysisManager->FillNtupleDColumn(8, fOuterE);
    analysisManager->FillNtupleDColumn(9, fOuterPID);
    analysisManager->AddNtupleRow();
}
