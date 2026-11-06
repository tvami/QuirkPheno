#ifndef QUIRKCORE_QuirkPhysics_H
#define QUIRKCORE_QuirkPhysics_H

// Quirk physics: the quirk/antiquirk definitions and their processes
// (QuirkTransportation, hadron-like EM energy loss, QuirkWatcher).
// Standalone version of the CMSSW CMSQuirkPhysics, configured by QuirkConfig.

#include "quirkcore/QuirkConfig.h"

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class Quirk;

class QuirkPhysics : public G4VPhysicsConstructor {
public:
  explicit QuirkPhysics(const QuirkConfig& p);
  ~QuirkPhysics() override = default;

  void ConstructParticle() override;
  void ConstructProcess() override;

private:
  const G4int m_verbose;
  const G4int m_pdg;
  const G4double m_mass;
  const G4double m_charge;
  const G4double m_stringForce;
  const G4double m_firstStringLength;
  const G4double m_maxBoost;
  const G4double m_maxMergeT;
  const G4double m_maxMergeMag;
  const G4int m_looperTrials;
  const G4double m_rMax;
  const G4double m_zMax;
  const G4bool m_keepStopped;
  const G4bool m_moveTinySteps;
  const G4bool m_pairProduction;
  const G4double m_crossingLength;
  const G4double m_stopThreshold;
  Quirk* m_quirk{nullptr};
  Quirk* m_antiQuirk{nullptr};
};

#endif
