#ifndef QSIM_QuirkSetup_hh
#define QSIM_QuirkSetup_hh

// Shared quirk configuration and the Pythia8 setup used by quirksim and quirkgen.

#include "quirkcore/QuirkConfig.h"

#include <memory>
#include <string>
#include <vector>

QuirkConfig& quirkConfig();

namespace Pythia8 {
  class Pythia;
}

// Drell-Yan pair of the fourth-generation lepton (PDG 17) through gamma*/Z
// (process 141 with Zprime:gmZmode = 4), as the CMSSW quirk fragment;
// mass and seed are set here, the rest comes from the settings file.
#ifdef QSIM_WITH_PYTHIA8
std::unique_ptr<Pythia8::Pythia> makeQuirkPythia(int seedOffset);
#endif
std::vector<std::string> quirkPythiaMassSettings(double massGeV);

#endif
