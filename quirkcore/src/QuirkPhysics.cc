// Copyright (C) 2002-2025 CERN for the benefit of the ATLAS collaboration
#include "quirkcore/QuirkPhysics.h"
#include "quirkcore/Quirk.h"
#include "quirkcore/QuirkStringStore.h"
#include "quirkcore/QuirkTransportation.h"
#include "quirkcore/QuirkWatcher.h"

#include "G4ios.hh"
#include "G4Exception.hh"

#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4EventManager.hh"
#include "G4StackManager.hh"
#include "G4Threading.hh"
#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"
#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>

namespace {
  // string tension F = Lambda^2 / hbar c
  double stringForce(double lambda) { return lambda * lambda / CLHEP::hbarc; }
}  // namespace

QuirkPhysics::QuirkPhysics(const QuirkConfig& p)
    : G4VPhysicsConstructor("Quirk Physics"),
      m_verbose(p.verbose),
      m_pdg(std::abs(p.pdg)),
      m_mass(p.mass),
      m_charge(p.charge),
      m_stringForce(stringForce(p.lambda)),
      m_firstStringLength(p.firstStringLength),
      m_maxBoost(p.maxBoost),
      m_maxMergeT(p.maxMerge),
      m_maxMergeMag(m_maxMergeT),
      m_looperTrials(p.looperTrials),
      m_rMax(p.maxRadius),
      m_zMax(p.maxZ),
      m_keepStopped(p.keepStopped),
      m_moveTinySteps(p.moveTinySteps),
      // off by default: G4MuPairProductionModel has no element data for PDG 17 (em0033 above 8 m_Q)
      m_pairProduction(p.pairProduction),
      // must stay well below the first string, else the restart repeats at production
      m_crossingLength(std::min(p.crossingLength, 0.1 * m_firstStringLength)),
      m_stopThreshold(p.stopThreshold) {
  if (m_mass <= 0. || m_stringForce <= 0. || m_maxBoost <= 0. || m_maxBoost >= 1.) {
    G4ExceptionDescription ed;
    ed << "QuirkPhysics: invalid parameters QuirkMass= " << m_mass / CLHEP::GeV
       << " GeV, QuirkLambda= " << p.lambda / CLHEP::eV << " eV, QuirkMaxBoost= " << m_maxBoost;
    G4Exception("QuirkPhysics::QuirkPhysics", "QuirkConfiguration", FatalException, ed);
  }
  G4cout << "QuirkPhysics: PDG " << m_pdg << " mass " << m_mass / CLHEP::GeV << " GeV, charge "
         << m_charge / CLHEP::eplus << " e, Lambda " << std::sqrt(m_stringForce * CLHEP::hbarc) / CLHEP::eV
         << " eV (string force " << m_stringForce / (CLHEP::MeV / CLHEP::mm) << " MeV/mm), first string length "
         << m_firstStringLength / CLHEP::mm << " mm, max boost " << m_maxBoost << ", max merge " << m_maxMergeT / CLHEP::mm
         << " mm, killed at r > " << m_rMax / CLHEP::cm << " cm or |z| > " << m_zMax / CLHEP::cm
         << " cm (<= 0: never), pair production " << (m_pairProduction ? "on" : "off") << ", stopped quirks "
         << (m_keepStopped ? "kept while the partner moves" : "killed") << G4endl;
}

void QuirkPhysics::ConstructParticle() {
  // definitions are created once, on the master
  if (nullptr != m_quirk) {
    return;
  }
  G4ParticleTable* table = G4ParticleTable::GetParticleTable();
  for (int pdg : {m_pdg, -m_pdg}) {
    const G4ParticleDefinition* old = table->FindParticle(pdg);
    if (nullptr != old) {
      G4ExceptionDescription ed;
      ed << "QuirkPhysics: PDG code " << pdg << " is already defined as " << old->GetParticleName()
         << "; do not combine with other exotica physics";
      G4Exception("QuirkPhysics::ConstructParticle", "QuirkConfiguration", FatalException, ed);
    }
  }
  m_quirk = new Quirk("quirk", m_mass, m_charge, m_pdg);
  m_antiQuirk = new Quirk("antiquirk", m_mass, -m_charge, -m_pdg);
}

void QuirkPhysics::ConstructProcess() {
  // string state of this thread
  QuirkStringStore::instance().configure(m_stringForce, m_firstStringLength, m_maxBoost, m_maxMergeT, m_maxMergeMag);

  for (Quirk* q : {m_quirk, m_antiQuirk}) {
    G4ProcessManager* pmanager = q->GetProcessManager();
    if (nullptr == pmanager) {
      G4Exception("QuirkPhysics::ConstructProcess", "QuirkConfiguration", FatalException,
                  "QuirkPhysics: quirk without a process manager");
    }
    // replace whatever was attached by the generic constructors
    while (pmanager->GetProcessListLength() != 0) {
      pmanager->RemoveProcess(0);
    }
    auto transport = new QuirkTransportation(m_verbose, m_moveTinySteps, m_crossingLength);
    transport->SetThresholdTrials(m_looperTrials);
    pmanager->AddProcess(transport, -1, 0, 0);
    pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
    pmanager->AddProcess(new G4hIonisation, -1, 2, 2);
    pmanager->AddProcess(new G4hBremsstrahlung, -1, 3, 3);
    if (m_pairProduction) {
      pmanager->AddProcess(new G4hPairProduction, -1, 4, 4);
    }
    pmanager->AddProcess(new QuirkWatcher(m_rMax, m_zMax, m_keepStopped, m_stopThreshold), -1, -1, 5);
    if (m_verbose > 1) {
      pmanager->DumpInfo();
    }
  }

  // the two quirks alternate through one extra waiting stack; user stacking actions keep
  // these defaults for suspended tracks
  if (!G4Threading::IsMasterThread()) {
    G4StackManager* stack = G4EventManager::GetEventManager()->GetStackManager();
    stack->SetNumberOfAdditionalWaitingStacks(1);
    stack->SetDefaultClassification(fSuspendAndWait, fWaiting);
    stack->SetDefaultClassification(fSuspend, fWaiting_1);
  }
}
