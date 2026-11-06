// Derived from atlas/athena Simulation/G4Extensions/Quirks (Copyright (C) 2002-2024 CERN for the
// benefit of the ATLAS collaboration) and cms-sw/cmssw, both Apache-2.0; modified.
#ifndef QUIRKCORE_QuirkStringStore_H
#define QUIRKCORE_QuirkStringStore_H

// Per-thread owner of the two coupled string states of the quirk pair.
// Index 0 is the quirk with positive PDG code, 1 the antiquirk.

#include "quirkcore/QuirkInfracolorForce.h"

class QuirkStringStore {
public:
  static QuirkStringStore& instance();

  QuirkInfracolorForce& stringFor(G4int pdg) { return (pdg > 0) ? m_strings[0] : m_strings[1]; }

  void configure(
      G4double stringForce, G4double firstStringLength, G4double maxBoost, G4double maxMergeT, G4double maxMergeMag);

  void clear() {
    m_strings[0].Clear();
    m_strings[1].Clear();
  }

  // true once per event: the first suspension waits for the urgent stack, later ones alternate
  G4bool firstSuspension(G4int eventID) {
    if (eventID == m_eventID)
      return false;
    m_eventID = eventID;
    return true;
  }

  QuirkStringStore(const QuirkStringStore&) = delete;
  QuirkStringStore& operator=(const QuirkStringStore&) = delete;

private:
  QuirkStringStore();

  QuirkInfracolorForce m_strings[2];
  G4int m_eventID{-1};
};

#endif
