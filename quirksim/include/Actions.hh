#ifndef QSIM_Actions_hh
#define QSIM_Actions_hh

// User actions: per-event quirk summaries, quirk hits in sensitive layers
// (one hit per continuous visit; each quirk keeps its own open hit, so a hit
// open when its quirk is suspended is parked, never merged with the partner's),
// detector energies, calorimeter towers, prescaled quirk steps.

#include "G4VUserActionInitialization.hh"
#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserStackingAction.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <chrono>
#include <unordered_map>
#include <vector>

// end reasons of a quirk track
enum QuirkEnd { kAlive = 0, kEnvelope = 1, kWorld = 2, kStopped = 3, kPartner = 4, kOther = 5, kTime = 6 };

struct OutBuffers {
  // events ntuple
  G4int evt, genEvt, nCross;
  G4double vx, vy, vz, cpu, wall, maxSep;
  std::vector<G4double> q_pdg, q_px0, q_py0, q_pz0, q_e0, q_x0, q_y0, q_z0, q_xe, q_ye, q_ze, q_te, q_pxe, q_pye,
      q_pze, q_ee, q_steps, q_susp, q_edep, q_path, q_reason, q_tExit, q_edepDet, eDet, eDetSens;
  std::vector<G4int> h_q, h_vol, h_det, h_layer, h_nsteps;
  std::vector<G4double> h_x0, h_y0, h_z0, h_t0, h_x1, h_y1, h_z1, h_t1, h_px, h_py, h_pz, h_dE, h_len;
  std::vector<G4int> tw_det, tw_ieta, tw_iphi;
  std::vector<G4double> tw_e, tw_eq;
  void clear();
};

class RunAction : public G4UserRunAction {
public:
  RunAction();
  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;
  OutBuffers& buffers() { return m_buf; }

private:
  OutBuffers m_buf;
};

class EventAction : public G4UserEventAction {
public:
  explicit EventAction(RunAction* run);
  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;

  struct Quirk {
    G4bool seen{false};
    G4int pdg{0};
    G4LorentzVector p0, pEnd;
    G4ThreeVector x0, xEnd;
    G4double tEnd{0}, edep{0}, path{0};
    long steps{0}, susp{0};
    G4int reason{kAlive};
    std::vector<G4double> tExit, edepDet;
    G4int openHit{-1}, openVol{-1};
  };
  struct Hit {
    G4int q, vol, nsteps{0};
    G4ThreeVector x0, x1, p;
    G4double t0, t1, dE{0}, len{0};
  };
  struct Tower {
    G4double e{0}, eq{0};
  };

  Quirk q[2];
  std::vector<Hit> hits;
  std::vector<G4double> eDet, eDetSens;  // all volumes, sensitive volumes only
  std::unordered_map<long, Tower> towers;
  G4ThreeVector lastPos[2];
  G4double maxSep{0};
  G4int nCross{0}, crossSign{0};
  long stepCount[2]{0, 0};

private:
  RunAction* m_run;
  std::chrono::steady_clock::time_point m_wall0;
  G4double m_cpu0{0};
};

class SteppingAction : public G4UserSteppingAction {
public:
  explicit SteppingAction(EventAction* ev) : m_ev(ev) {}
  void UserSteppingAction(const G4Step* step) override;

private:
  void quirkStep(const G4Step* step, int i);
  EventAction* m_ev;
};

class StackingAction : public G4UserStackingAction {
public:
  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* t) override;
};

class ActionInitialization : public G4VUserActionInitialization {
public:
  void BuildForMaster() const override;
  void Build() const override;
};

#endif
