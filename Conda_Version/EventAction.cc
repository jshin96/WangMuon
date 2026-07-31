#include "EventAction.hh"
#include "G4Event.hh"
#include "G4AnalysisManager.hh"

EventAction::EventAction()
: G4UserEventAction(),
  fPosIn1(0.), fMomIn1(0.), fPosIn2(0.), fMomIn2(0.),
  fPosOut1(0.), fMomOut1(0.), fPosOut2(0.), fMomOut2(0.),
  fHitIn1(false), fHitIn2(false), fHitOut1(false), fHitOut2(false)
{}

EventAction::~EventAction() {}

void EventAction::BeginOfEventAction(const G4Event*)
{
    // Reset the 4-fold coincidence trigger
    fHitIn1 = false;
    fHitIn2 = false;
    fHitOut1 = false;
    fHitOut2 = false;
}

void EventAction::EndOfEventAction(const G4Event* event)
{
    // Hardware Trigger: Only record muons that survived the entire journey
    if (fHitIn1 && fHitIn2 && fHitOut1 && fHitOut2) {
        
        auto analysisManager = G4AnalysisManager::Instance();
        // Column 0 stores the EcoMug trial count from PrimaryGeneratorAction.
        int col = 1;

        analysisManager->FillNtupleIColumn(col++, event->GetEventID());
        
        // In 1
        analysisManager->FillNtupleDColumn(col++, fPosIn1.x());
        analysisManager->FillNtupleDColumn(col++, fPosIn1.y());
        analysisManager->FillNtupleDColumn(col++, fPosIn1.z());
        analysisManager->FillNtupleDColumn(col++, fMomIn1.x());
        analysisManager->FillNtupleDColumn(col++, fMomIn1.y());
        analysisManager->FillNtupleDColumn(col++, fMomIn1.z());

        // In 2
        analysisManager->FillNtupleDColumn(col++, fPosIn2.x());
        analysisManager->FillNtupleDColumn(col++, fPosIn2.y());
        analysisManager->FillNtupleDColumn(col++, fPosIn2.z());
        analysisManager->FillNtupleDColumn(col++, fMomIn2.x());
        analysisManager->FillNtupleDColumn(col++, fMomIn2.y());
        analysisManager->FillNtupleDColumn(col++, fMomIn2.z());

        // Out 1
        analysisManager->FillNtupleDColumn(col++, fPosOut1.x());
        analysisManager->FillNtupleDColumn(col++, fPosOut1.y());
        analysisManager->FillNtupleDColumn(col++, fPosOut1.z());
        analysisManager->FillNtupleDColumn(col++, fMomOut1.x());
        analysisManager->FillNtupleDColumn(col++, fMomOut1.y());
        analysisManager->FillNtupleDColumn(col++, fMomOut1.z());

        // Out 2
        analysisManager->FillNtupleDColumn(col++, fPosOut2.x());
        analysisManager->FillNtupleDColumn(col++, fPosOut2.y());
        analysisManager->FillNtupleDColumn(col++, fPosOut2.z());
        analysisManager->FillNtupleDColumn(col++, fMomOut2.x());
        analysisManager->FillNtupleDColumn(col++, fMomOut2.y());
        analysisManager->FillNtupleDColumn(col++, fMomOut2.z());

        analysisManager->AddNtupleRow();
    }
}
