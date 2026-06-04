#include "RunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Run.hh"

RunAction::RunAction() {
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("root");
    
    // Create the Data Tree
    // Full 2pi coverage detector
    analysisManager->CreateNtuple("DetectorHits", "Muon Hit Coordinates");
    analysisManager->CreateNtupleDColumn("TriggerType");     // 1=InnerOnly, 2=OuterOnly, 3=Coincidence
    analysisManager->CreateNtupleDColumn("PhysicsCategory"); // 1=Wall/Room, 2=Mound, 3=Air

    // Inner Detector Data
    analysisManager->CreateNtupleDColumn("InnerH");
    analysisManager->CreateNtupleDColumn("InnerZ");
    analysisManager->CreateNtupleDColumn("InnerEnergy");
    analysisManager->CreateNtupleDColumn("InnerPID");

    // Outer Detector Data
    analysisManager->CreateNtupleDColumn("OuterH");
    analysisManager->CreateNtupleDColumn("OuterZ");
    analysisManager->CreateNtupleDColumn("OuterEnergy");
    analysisManager->CreateNtupleDColumn("OuterPID");
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
