#include "RunAction.hh"
#include "G4AnalysisManager.hh"

RunAction::RunAction() {
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("root");
    analysisManager->SetNtupleMerging(true);
    
    analysisManager->CreateNtuple("MuonTree", "Cosmic Muon Data");
    
    // Columns 0-6: Entry Data (Scalars)
    analysisManager->CreateNtupleDColumn("EntryPosX");
    analysisManager->CreateNtupleDColumn("EntryPosY");
    analysisManager->CreateNtupleDColumn("EntryPosZ");
    analysisManager->CreateNtupleDColumn("EntryPx");
    analysisManager->CreateNtupleDColumn("EntryPy");
    analysisManager->CreateNtupleDColumn("EntryPz");
    analysisManager->CreateNtupleDColumn("EntryE");

    // Columns 7-13: Exit Data (Scalars)
    analysisManager->CreateNtupleDColumn("ExitPosX");
    analysisManager->CreateNtupleDColumn("ExitPosY");
    analysisManager->CreateNtupleDColumn("ExitPosZ");
    analysisManager->CreateNtupleDColumn("ExitPx");
    analysisManager->CreateNtupleDColumn("ExitPy");
    analysisManager->CreateNtupleDColumn("ExitPz");
    analysisManager->CreateNtupleDColumn("ExitE");

    // Columns 14-16: Room and Distance Data (Scalars)
    analysisManager->CreateNtupleIColumn("PassedRoom");
    analysisManager->CreateNtupleDColumn("DistRoom");
    analysisManager->CreateNtupleDColumn("DistDirt");

    // Columns 17-23: Vectors (Bound directly to memory)
    analysisManager->CreateNtupleDColumn("StepPosX", fStepPosX);
    analysisManager->CreateNtupleDColumn("StepPosY", fStepPosY);
    analysisManager->CreateNtupleDColumn("StepPosZ", fStepPosZ);
    analysisManager->CreateNtupleDColumn("StepPx", fStepPx);
    analysisManager->CreateNtupleDColumn("StepPy", fStepPy);
    analysisManager->CreateNtupleDColumn("StepPz", fStepPz);
    analysisManager->CreateNtupleDColumn("StepE", fStepE);
    
    analysisManager->FinishNtuple();
}

RunAction::~RunAction() {}

void RunAction::BeginOfRunAction(const G4Run*) {
    G4AnalysisManager::Instance()->OpenFile("wang_muon_data.root");
}

void RunAction::EndOfRunAction(const G4Run*) {
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
}

void RunAction::ClearEventData() {
    // Reset scalars to dummy values
    fEntryPosX = -999.0; fEntryPosY = -999.0; fEntryPosZ = -999.0;
    fEntryPx = -999.0;   fEntryPy = -999.0;   fEntryPz = -999.0; fEntryE = -999.0;

    fExitPosX = -999.0;  fExitPosY = -999.0;  fExitPosZ = -999.0;
    fExitPx = -999.0;    fExitPy = -999.0;    fExitPz = -999.0;  fExitE = -999.0;

    // Reset flags and distances
    fPassedRoom = 0; 
    fDistRoom = 0.0; 
    fDistDirt = 0.0;

    // Clear vectors
    fStepPosX.clear(); fStepPosY.clear(); fStepPosZ.clear();
    fStepPx.clear();   fStepPy.clear();   fStepPz.clear();   fStepE.clear();
}
