#include "PrimaryGenerator.hh"
#include "QSimConfig.hh"
#include "QuirkSetup.hh"

#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4ParticleTable.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"
#include "G4Exception.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

#ifdef QSIM_WITH_HEPMC3
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#endif
#ifdef QSIM_WITH_PYTHIA8
#include "Pythia8/Pythia.h"
#endif

#include <map>
#include <set>

namespace {
  G4Mutex readerMutex = G4MUTEX_INITIALIZER;
  G4bool isQuirkPdg(G4int pdg) { return std::abs(pdg) == std::abs(quirkConfig().pdg); }

#ifdef QSIM_WITH_HEPMC3
  // one reader for all threads; Geant4 event k always gets file event k (after skipEvents),
  // so the per-event seeds and the output do not depend on the thread count
  std::map<G4int, HepMC3::GenEvent>& readAhead() {
    static std::map<G4int, HepMC3::GenEvent> m;
    return m;
  }
  G4int nRead = 0;
  HepMC3::ReaderAscii* reader() {
    static std::unique_ptr<HepMC3::ReaderAscii> r;
    if (!r) {
      const QSimConfig& c = qsimConfig();
      r = std::make_unique<HepMC3::ReaderAscii>(c.hepmcFile);
      if (r->failed())
        G4Exception("PrimaryGenerator", "QSimHepMC", FatalException, ("cannot read " + c.hepmcFile).c_str());
      HepMC3::GenEvent skip;
      for (int i = 0; i < c.skipEvents; ++i)
        r->read_event(skip);
    }
    return r.get();
  }
#endif
}  // namespace

PrimaryGenerator::PrimaryGenerator() = default;
PrimaryGenerator::~PrimaryGenerator() = default;

void PrimaryGenerator::addParticle(G4Event* ev, G4int pdg, const G4ThreeVector& p, const G4ThreeVector& x, G4double t) {
  const QSimConfig& c = qsimConfig();
  if (c.quirksOnly && !isQuirkPdg(pdg))
    return;
  if (p.perp() > 0 && std::abs(p.eta()) > c.etaMax)
    return;
  const G4ParticleDefinition* def = G4ParticleTable::GetParticleTable()->FindParticle(pdg);
  if (nullptr == def) {
    static G4ThreadLocal std::set<G4int>* warned = nullptr;
    if (!warned)
      warned = new std::set<G4int>;
    if (warned->insert(pdg).second)
      G4Exception("PrimaryGenerator", "QSimUnknownPDG", JustWarning, ("skipping PDG " + std::to_string(pdg)).c_str());
    return;
  }
  auto vtx = new G4PrimaryVertex(x + m_shift, t);
  vtx->SetPrimary(new G4PrimaryParticle(def, p.x(), p.y(), p.z()));
  ev->AddPrimaryVertex(vtx);
}

void PrimaryGenerator::GeneratePrimaries(G4Event* ev) {
  const QSimConfig& c = qsimConfig();
  auto info = new EventInfo;
  ev->SetUserInformation(info);
  m_shift = G4ThreeVector(c.vtxSigmaXY > 0 ? G4RandGauss::shoot(0., c.vtxSigmaXY) : 0.,
                          c.vtxSigmaXY > 0 ? G4RandGauss::shoot(0., c.vtxSigmaXY) : 0.,
                          c.vtxSigmaZ > 0 ? G4RandGauss::shoot(0., c.vtxSigmaZ) : 0.);
  info->vertex = m_shift;
  if (c.genMode == "gun") {
    const G4int pdg = std::abs(quirkConfig().pdg);
    addParticle(ev, pdg, c.gunP1, c.gunVertex, 0.);
    addParticle(ev, -pdg, c.gunP2, c.gunVertex, 0.);
    info->vertex += c.gunVertex;
  } else if (c.genMode == "hepmc") {
    fromHepMC(ev, info);
  } else if (c.genMode == "pythia") {
    fromPythia(ev, info);
  } else {
    G4Exception("PrimaryGenerator", "QSimGenMode", FatalException, ("unknown mode " + c.genMode).c_str());
  }
}

void PrimaryGenerator::fromHepMC(G4Event* ev, EventInfo* info) {
#ifdef QSIM_WITH_HEPMC3
  HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::MM);
  {
    G4AutoLock lock(&readerMutex);
    auto r = reader();
    auto& buf = readAhead();
    const G4int k = ev->GetEventID();
    while (buf.find(k) == buf.end()) {
      HepMC3::GenEvent next(HepMC3::Units::GEV, HepMC3::Units::MM);
      if (!r->read_event(next) || r->failed()) {
        G4Exception("PrimaryGenerator", "QSimHepMCEnd", RunMustBeAborted, "end of HepMC3 file");
        ev->SetEventAborted();
        return;
      }
      buf.emplace(nRead++, std::move(next));
    }
    evt = std::move(buf[k]);
    buf.erase(k);
  }
  evt.set_units(HepMC3::Units::GEV, HepMC3::Units::MM);
  info->genEvent = evt.event_number();
  for (const auto& p : evt.particles()) {
    if (p->status() != 1)
      continue;
    HepMC3::FourVector x(0, 0, 0, 0);
    if (auto v = p->production_vertex())
      x = v->position();
    const auto& m = p->momentum();
    addParticle(ev, p->pid(), G4ThreeVector(m.px(), m.py(), m.pz()) * GeV, G4ThreeVector(x.x(), x.y(), x.z()) * mm,
                x.t() * mm / c_light);
  }
#else
  (void)ev;
  (void)info;
  G4Exception("PrimaryGenerator", "QSimNoHepMC", FatalException, "built without HepMC3");
#endif
}

void PrimaryGenerator::fromPythia(G4Event* ev, EventInfo* info) {
#ifdef QSIM_WITH_PYTHIA8
  if (!m_pythia) {
    m_pythia = makeQuirkPythia(G4Threading::G4GetThreadId());
  }
  while (!m_pythia->next()) {
  }
  info->genEvent = -1;
  const Pythia8::Event& e = m_pythia->event;
  for (int i = 0; i < e.size(); ++i) {
    if (!e[i].isFinal())
      continue;
    addParticle(ev, e[i].id(), G4ThreeVector(e[i].px(), e[i].py(), e[i].pz()) * GeV,
                G4ThreeVector(e[i].xProd(), e[i].yProd(), e[i].zProd()) * mm, e[i].tProd() * mm / c_light);
  }
#else
  (void)ev;
  (void)info;
  G4Exception("PrimaryGenerator", "QSimNoPythia", FatalException, "built without Pythia8");
#endif
}
