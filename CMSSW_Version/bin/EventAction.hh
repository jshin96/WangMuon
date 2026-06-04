#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class EventAction : public G4UserEventAction {
public:
    EventAction();
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

    // Setters for the SteppingAction to flip
    void SetPassedMound() { fPassedMound = true; }
    void SetPassedWall()  { fPassedWall = true; }
    void SetPassedRoom()  { fPassedRoom = true; }



    void SetInnerHit(double h, double z, double e, double pid) {
        if (fHitInner) return; // Only record the first particle to cross
        fHitInner = true; fInnerH = h; fInnerZ = z; fInnerE = e; fInnerPID = pid;
    }

    void SetOuterHit(double h, double z, double e, double pid) {
        if (fHitOuter) return; // Only record the first particle to cross
        fHitOuter = true; fOuterH = h; fOuterZ = z; fOuterE = e; fOuterPID = pid;
    }


    // Getters for the final logic check
    bool GetPassedMound() const { return fPassedMound; }
    bool GetPassedWall()  const { return fPassedWall; }
    bool GetPassedRoom()  const { return fPassedRoom; }

private:
    bool fPassedMound = false;
    bool fPassedWall = false;
    bool fPassedRoom = false;

    bool fHitInner = false;
    bool fHitOuter = false;

    // Inner Data
    double fInnerH = 0.0, fInnerZ = 0.0, fInnerE = 0.0, fInnerPID = 0.0;

    // Outer Data
    double fOuterH = 0.0, fOuterZ = 0.0, fOuterE = 0.0, fOuterPID = 0.0;
};

#endif
