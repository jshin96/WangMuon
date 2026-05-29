#include "EventAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4UImanager.hh"


EventAction::EventAction(RunAction* runAction) : fRunAction(runAction) {}

void EventAction::BeginOfEventAction(const G4Event*) {
    // Wipe the slate clean before the cosmic muon fires
    fRunAction->ClearEventData();
}

void EventAction::EndOfEventAction(const G4Event* event) {
    auto analysisManager = G4AnalysisManager::Instance();
    int eventID = event->GetEventID();   
    // --- NEW: AUTOMATIC SCREENSHOT ---
    G4UImanager* ui = G4UImanager::GetUIpointer();

    // Force the viewer to render the tracks
    ui->ApplyCommand("/vis/viewer/refresh");

    // 1. Explicitly lock the format to jpg
    ui->ApplyCommand("/vis/ogl/set/exportFormat jpg");
    
    // 2. Generate a custom filename (e.g., Event_0.png, Event_1.png)
    G4String exportCmd = "/vis/ogl/export Event_" + std::to_string(eventID);
    ui->ApplyCommand(exportCmd);
    G4String renameCmd = "mv Event_" + std::to_string(eventID) + "_*.jpg Event_" + std::to_string(eventID) + ".jpg 2>/dev/null";
    system(renameCmd.c_str());

    // ---------------------------------
    // Fill all the scalar columns manually by their column index
    analysisManager->FillNtupleDColumn(0, fRunAction->fEntryPosX);
    analysisManager->FillNtupleDColumn(1, fRunAction->fEntryPosY);
    analysisManager->FillNtupleDColumn(2, fRunAction->fEntryPosZ);
    analysisManager->FillNtupleDColumn(3, fRunAction->fEntryPx);
    analysisManager->FillNtupleDColumn(4, fRunAction->fEntryPy);
    analysisManager->FillNtupleDColumn(5, fRunAction->fEntryPz);
    analysisManager->FillNtupleDColumn(6, fRunAction->fEntryE);

    analysisManager->FillNtupleDColumn(7, fRunAction->fExitPosX);
    analysisManager->FillNtupleDColumn(8, fRunAction->fExitPosY);
    analysisManager->FillNtupleDColumn(9, fRunAction->fExitPosZ);
    analysisManager->FillNtupleDColumn(10, fRunAction->fExitPx);
    analysisManager->FillNtupleDColumn(11, fRunAction->fExitPy);
    analysisManager->FillNtupleDColumn(12, fRunAction->fExitPz);
    analysisManager->FillNtupleDColumn(13, fRunAction->fExitE);

    analysisManager->FillNtupleIColumn(14, fRunAction->fPassedRoom);
    analysisManager->FillNtupleDColumn(15, fRunAction->fDistRoom);
    analysisManager->FillNtupleDColumn(16, fRunAction->fDistDirt);

    // Commit the row (including the bound vectors) to the TTree
    analysisManager->AddNtupleRow();

}
