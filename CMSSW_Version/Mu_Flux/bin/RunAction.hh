#ifndef RunAction_h
#define RunAction_h 1
#include "G4UserRunAction.hh"
// This class names the columns in the results file and opens/closes that file.
class RunAction : public G4UserRunAction { public: RunAction(); ~RunAction() override=default; void BeginOfRunAction(const G4Run*) override; void EndOfRunAction(const G4Run*) override; };
#endif
