#ifndef QUIRKCORE_QuirkConfig_H
#define QUIRKCORE_QuirkConfig_H

// Quirk parameters, set by the /quirk/ macro commands; values are in Geant4 internal units.

#include "globals.hh"
#include <CLHEP/Units/SystemOfUnits.h>

struct QuirkConfig {
  G4int verbose{0};
  G4int pdg{17};
  G4double mass{0.};                                // > 0 enables quirks
  G4double charge{-1. * CLHEP::eplus};              // of the positive PDG code
  G4double lambda{0.};                              // F = Lambda^2/(hbar c)
  G4double firstStringLength{1.e-6 * CLHEP::mm};
  G4double maxBoost{0.1};                           // max rapidity change per step (Athena: 0.01)
  G4double maxMerge{1.e-6 * CLHEP::mm};
  G4int looperTrials{10};
  G4double maxRadius{800. * CLHEP::cm};             // quirks killed outside (<= 0: never)
  G4double maxZ{1100. * CLHEP::cm};                 // (<= 0: no envelope kill)
  G4bool keepStopped{false};                        // keep quirks stopped by dE/dx while the partner moves
  G4double stopThreshold{10. * CLHEP::MeV};         // partner below it: pair stopped
  G4bool moveTinySteps{true};                       // move quirks on sub-tolerance steps (head-on crossings)
  G4double crossingLength{0.};                      // shorter incoming string: restart it (0: never, experimental)
  G4bool pairProduction{false};                     // hPairProd, no Geant4 element data for PDG 17 (em0033)
};

#endif
