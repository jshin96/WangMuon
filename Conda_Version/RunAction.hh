#ifndef RUN_ACTION_HH
#define RUN_ACTION_HH

#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include <vector>

class RunAction : public G4UserRunAction {
public:
    RunAction();
    ~RunAction() override;

    void BeginOfRunAction(const G4Run* run) override;
    void EndOfRunAction(const G4Run* run) override;

    // A public method to clear the vectors/scalars at the start of every new event
    void ClearEventData();

    // --- Scalars ---
    double fEntryPosX, fEntryPosY, fEntryPosZ;
    double fExitPosX, fExitPosY, fExitPosZ;
    double fEntryPx, fEntryPy, fEntryPz, fEntryE;
    double fExitPx, fExitPy, fExitPz, fExitE;
    
    int fPassedRoom; 
    double fDistRoom;
    double fDistDirt;

    // --- Vectors ---
    std::vector<double> fStepPosX, fStepPosY, fStepPosZ;
    std::vector<double> fStepPx, fStepPy, fStepPz, fStepE;
};
#endif
