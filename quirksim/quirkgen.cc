// quirkgen: quirk Drell-Yan pairs from Pythia8 into a HepMC3 ascii file,
// with the same settings as the built-in quirksim generator.
//
// usage: quirkgen <mass GeV> <nevents> <out.hepmc3> <settings.cmnd> [seed] [vtxSigmaZ_mm]

#include "QSimConfig.hh"
#include "QuirkSetup.hh"

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"

#include <CLHEP/Units/SystemOfUnits.h>
#include <cstdlib>
#include <iostream>
#include <random>

int main(int argc, char** argv) {
  if (argc < 5) {
    std::cerr << "usage: " << argv[0] << " <mass GeV> <nevents> <out.hepmc3> <settings.cmnd> [seed] [vtxSigmaZ_mm]\n";
    return 1;
  }
  quirkConfig().mass = std::atof(argv[1]) * CLHEP::GeV;
  const int n = std::atoi(argv[2]);
  qsimConfig().pythiaCmnd = argv[4];
  if (argc > 5)
    qsimConfig().pythiaSeed = std::atoi(argv[5]);
  const double sigz = argc > 6 ? std::atof(argv[6]) : 0.;
  auto pythia = makeQuirkPythia(0);
  Pythia8::Pythia8ToHepMC out(argv[3]);
  std::mt19937 rng(qsimConfig().pythiaSeed);
  std::normal_distribution<double> gaus(0., 1.);
  for (int i = 0; i < n;) {
    if (!pythia->next())
      continue;
    if (sigz > 0) {
      // Gaussian luminous region along z, as the CMS vertex smearing
      const double dz = sigz * gaus(rng);
      for (int k = 0; k < pythia->event.size(); ++k)
        pythia->event[k].vProdAdd(Pythia8::Vec4(0., 0., dz, 0.));
    }
    ++i;
    out.writeNextEvent(*pythia);
  }
  pythia->stat();
  std::cout << "quirkgen: wrote " << n << " events, sigma = " << pythia->info.sigmaGen() * 1e9 << " pb" << std::endl;
  return 0;
}
