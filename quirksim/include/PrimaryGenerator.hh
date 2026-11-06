#ifndef QSIM_PrimaryGenerator_hh
#define QSIM_PrimaryGenerator_hh

// Primaries from the fixed pair gun, a HepMC3 ascii file (shared reader,
// events handed out in file order) or a per-thread Pythia8.

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VUserEventInformation.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <memory>

class G4Event;
class G4ParticleDefinition;
namespace Pythia8 {
  class Pythia;
}

class EventInfo : public G4VUserEventInformation {
public:
  G4int genEvent{-1};  // HepMC event number, else -1
  G4ThreeVector vertex;
  void Print() const override {}
};

class PrimaryGenerator : public G4VUserPrimaryGeneratorAction {
public:
  PrimaryGenerator();
  ~PrimaryGenerator() override;
  void GeneratePrimaries(G4Event* ev) override;

private:
  void addParticle(G4Event* ev, G4int pdg, const G4ThreeVector& p, const G4ThreeVector& x, G4double t);
  void fromHepMC(G4Event* ev, EventInfo* info);
  void fromPythia(G4Event* ev, EventInfo* info);
  G4ThreeVector m_shift;
#ifdef QSIM_WITH_PYTHIA8
  std::unique_ptr<Pythia8::Pythia> m_pythia;
#endif
};

#endif
