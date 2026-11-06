// Derived from atlas/athena Simulation/G4Extensions/Quirks (Copyright (C) 2002-2024 CERN for the
// benefit of the ATLAS collaboration) and cms-sw/cmssw, both Apache-2.0; modified.
#ifndef QUIRKCORE_Quirk_H
#define QUIRKCORE_Quirk_H

// Quirk particle definition (lepton-like, infracolor charged).
// Carries no string state: the definition is shared between threads, the
// state lives in the thread-local QuirkStringStore.

#include "G4ParticleDefinition.hh"

class Quirk : public G4ParticleDefinition {
public:
  Quirk(const G4String& name, G4double mass, G4double charge, G4int pdg);
  ~Quirk() override = default;

  Quirk(const Quirk&) = delete;
  Quirk& operator=(const Quirk&) = delete;

  static bool isQuirk(const G4ParticleDefinition* p) { return nullptr != dynamic_cast<const Quirk*>(p); }
};

#endif
