#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

// This is the list of helper jobs Geant4 uses for every simulation event.
class ActionInitialization : public G4VUserActionInitialization {
  public:
    ActionInitialization();
    virtual ~ActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;
};

#endif
