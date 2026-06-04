#include "RunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Run.hh"

RunAction::RunAction() {
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("csv");

    analysisManager->OpenFile("MoundScanData");
    analysisManager->CreateNtuple("MuonHits", "Tomography Data");
    analysisManager->CreateNtupleIColumn("EventID");
    analysisManager->CreateNtupleDColumn("X_in");
    analysisManager->CreateNtupleDColumn("Y_in");
    analysisManager->CreateNtupleDColumn("Z_in");
    analysisManager->CreateNtupleDColumn("PX_in");
    analysisManager->CreateNtupleDColumn("PY_in");
    analysisManager->CreateNtupleDColumn("PZ_in");
    analysisManager->CreateNtupleDColumn("X_out");
    analysisManager->CreateNtupleDColumn("Y_out");
    analysisManager->CreateNtupleDColumn("Z_out");
    analysisManager->CreateNtupleDColumn("PX_out");
    analysisManager->CreateNtupleDColumn("PY_out");
    analysisManager->CreateNtupleDColumn("PZ_out");

    analysisManager->FinishNtuple();
}

RunAction::~RunAction() {}

void RunAction::BeginOfRunAction(const G4Run*) {
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->OpenFile("wang_muon_data.root");
}

void RunAction::EndOfRunAction(const G4Run*) {
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
}
