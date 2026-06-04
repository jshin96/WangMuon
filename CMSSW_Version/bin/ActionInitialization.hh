#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class ActionInitialization : public G4VUserActionInitialization {
  public:
    ActionInitialization();
    virtual ~ActionInitialization(); // <-- This fixes the destructor error!

    virtual void BuildForMaster() const;
    virtual void Build() const;
};

#endif
