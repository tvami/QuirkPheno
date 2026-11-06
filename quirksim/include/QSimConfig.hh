#ifndef QSIM_QSimConfig_hh
#define QSIM_QSimConfig_hh

// Application settings (/qsim/ commands); set in PreInit, read-only afterwards.

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <CLHEP/Units/SystemOfUnits.h>

#include <string>
#include <vector>

struct QSimConfig {
  // geometry and physics
  std::string card{"cards/cms_like.card"};
  std::string physList{"FTFP_BERT"};
  G4double prodCut{0.7 * CLHEP::mm};

  // generator: "gun" (fixed pair), "hepmc" (HepMC3 ascii) or "pythia" (built-in Pythia8)
  std::string genMode{"gun"};
  std::string hepmcFile;
  G4bool quirksOnly{false};  // pass only the quirks of the HepMC/Pythia event to Geant4
  G4double etaMax{5.5};      // CMSSW Generator eta cut
  G4int skipEvents{0};
  G4ThreeVector gunP1{100. * CLHEP::GeV, -60. * CLHEP::GeV, 1. * CLHEP::GeV};
  G4ThreeVector gunP2{50. * CLHEP::GeV, 60. * CLHEP::GeV, 1. * CLHEP::GeV};
  G4ThreeVector gunVertex{0., 0., 0.};
  G4double vtxSigmaXY{0.};  // Gaussian vertex smearing
  G4double vtxSigmaZ{0.};
  std::string pythiaCmnd;   // Pythia8 settings file
  std::vector<std::string> pythiaStrings;
  G4int pythiaSeed{12345};  // + thread id

  // tracking
  G4double maxTime{500. * CLHEP::ns};  // CMSSW MaxTrackTime, all tracks (<= 0: none)
  G4bool killSecondaries{false};       // quirks only, no showers

  // output
  std::string output{"quirksim.root"};
  G4int stepPrescale{0};  // write every N-th quirk step (0: none)
  G4bool towers{true};    // calorimeter towers
  G4int printEvery{1};
};

QSimConfig& qsimConfig();

#endif
