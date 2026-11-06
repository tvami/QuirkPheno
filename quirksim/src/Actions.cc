#include "Actions.hh"
#include "Geometry.hh"
#include "PrimaryGenerator.hh"
#include "QSimConfig.hh"
#include "QuirkSetup.hh"

#include "quirkcore/Quirk.h"
#include "quirkcore/QuirkStringStore.h"
#include "quirkcore/QuirkTrackStatus.h"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StackManager.hh"
#include "G4EventManager.hh"
#include "G4Threading.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Version.hh"

#include <ctime>
#include <fstream>
#include <iomanip>

namespace {
  enum { kEvents = 0, kSteps = 1 };
  constexpr G4double kTowerEta = 0.087;
  constexpr G4int kTowerPhi = 72;

  G4double threadCpu() {
    timespec ts;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
  }
}  // namespace

void OutBuffers::clear() {
  for (auto* v : {&q_pdg, &q_px0, &q_py0, &q_pz0, &q_e0, &q_x0, &q_y0, &q_z0, &q_xe, &q_ye, &q_ze, &q_te, &q_pxe,
                  &q_pye, &q_pze, &q_ee, &q_steps, &q_susp, &q_edep, &q_path, &q_reason, &q_tExit, &q_edepDet, &eDet, &eDetSens,
                  &h_x0, &h_y0, &h_z0, &h_t0, &h_x1, &h_y1, &h_z1, &h_t1, &h_px, &h_py, &h_pz, &h_dE, &h_len, &tw_e,
                  &tw_eq})
    v->clear();
  for (auto* v : {&h_q, &h_vol, &h_det, &h_layer, &h_nsteps, &tw_det, &tw_ieta, &tw_iphi})
    v->clear();
}

// ---------------- run ----------------
RunAction::RunAction() {
  auto am = G4AnalysisManager::Instance();
  am->SetDefaultFileType("root");
  am->SetVerboseLevel(0);
  am->SetNtupleMerging(true);
  OutBuffers& b = m_buf;
  am->CreateNtuple("events", "quirk events (mm, ns, MeV)");
  am->CreateNtupleIColumn("evt");
  am->CreateNtupleIColumn("genEvt");
  am->CreateNtupleIColumn("nCross");
  for (const char* n : {"vx", "vy", "vz", "cpu", "wall", "maxSep"})
    am->CreateNtupleDColumn(n);
  const std::pair<const char*, std::vector<G4double>*> dv[] = {
      {"q_pdg", &b.q_pdg},     {"q_px0", &b.q_px0},   {"q_py0", &b.q_py0},       {"q_pz0", &b.q_pz0},
      {"q_e0", &b.q_e0},       {"q_x0", &b.q_x0},     {"q_y0", &b.q_y0},         {"q_z0", &b.q_z0},
      {"q_xe", &b.q_xe},       {"q_ye", &b.q_ye},     {"q_ze", &b.q_ze},         {"q_te", &b.q_te},
      {"q_pxe", &b.q_pxe},     {"q_pye", &b.q_pye},   {"q_pze", &b.q_pze},       {"q_ee", &b.q_ee},
      {"q_steps", &b.q_steps}, {"q_susp", &b.q_susp}, {"q_edep", &b.q_edep},     {"q_path", &b.q_path},
      {"q_reason", &b.q_reason}, {"q_tExit", &b.q_tExit}, {"q_edepDet", &b.q_edepDet}, {"eDet", &b.eDet}, {"eDetSens", &b.eDetSens},
      {"h_x0", &b.h_x0},       {"h_y0", &b.h_y0},     {"h_z0", &b.h_z0},         {"h_t0", &b.h_t0},
      {"h_x1", &b.h_x1},       {"h_y1", &b.h_y1},     {"h_z1", &b.h_z1},         {"h_t1", &b.h_t1},
      {"h_px", &b.h_px},       {"h_py", &b.h_py},     {"h_pz", &b.h_pz},         {"h_dE", &b.h_dE},
      {"h_len", &b.h_len},     {"tw_e", &b.tw_e},     {"tw_eq", &b.tw_eq}};
  for (auto& [n, v] : dv)
    am->CreateNtupleDColumn(n, *v);
  const std::pair<const char*, std::vector<G4int>*> iv[] = {{"h_q", &b.h_q},         {"h_vol", &b.h_vol},
                                                           {"h_det", &b.h_det},     {"h_layer", &b.h_layer},
                                                           {"h_nsteps", &b.h_nsteps}, {"tw_det", &b.tw_det},
                                                           {"tw_ieta", &b.tw_ieta}, {"tw_iphi", &b.tw_iphi}};
  for (auto& [n, v] : iv)
    am->CreateNtupleIColumn(n, *v);
  am->FinishNtuple();

  am->CreateNtuple("steps", "prescaled quirk steps (mm, ns, MeV)");
  for (const char* n : {"evt", "q", "step", "vol", "status"})
    am->CreateNtupleIColumn(n);
  for (const char* n : {"x", "y", "z", "t", "px", "py", "pz", "ekin", "edep"})
    am->CreateNtupleDColumn(n);
  am->FinishNtuple();
}

void RunAction::BeginOfRunAction(const G4Run*) {
  const QSimConfig& c = qsimConfig();
  G4AnalysisManager::Instance()->OpenFile(c.output);
  if (!G4Threading::IsMasterThread())
    return;
  // geometry table and settings next to the ROOT file
  const GeomInfo& g = geomInfo();
  const QuirkConfig& qc = quirkConfig();
  std::string side = c.output;
  if (auto p = side.rfind(".root"); p != std::string::npos)
    side.erase(p);
  std::ofstream js(side + "_meta.json");
  js << std::setprecision(10) << "{\n  \"geant4\": \"" << G4Version << "\",\n  \"physList\": \"" << c.physList
     << "\",\n  \"card\": \"" << c.card << "\",\n  \"mass_GeV\": " << qc.mass / GeV << ",\n  \"lambda_eV\": "
     << qc.lambda / eV << ",\n  \"maxBoost\": " << qc.maxBoost << ",\n  \"maxTime_ns\": " << c.maxTime / ns
     << ",\n  \"genMode\": \"" << c.genMode << "\",\n  \"hepmcFile\": \"" << c.hepmcFile
     << "\",\n  \"quirksOnly\": " << (c.quirksOnly ? "true" : "false") << ",\n  \"dets\": [";
  for (size_t i = 0; i < g.dets.size(); ++i)
    js << (i ? ", " : "") << '"' << g.dets[i] << '"';
  js << "],\n  \"boundaries\": [";
  for (size_t i = 0; i < g.boundaries.size(); ++i)
    js << (i ? ", " : "") << "{\"name\": \"" << g.boundaries[i].name << "\", \"r_mm\": " << g.boundaries[i].r
       << ", \"halfz_mm\": " << g.boundaries[i].halfz << "}";
  js << "],\n  \"vols\": [\n";
  for (size_t i = 0; i < g.vols.size(); ++i) {
    const VolInfo& v = g.vols[i];
    js << "    {\"name\": \"" << v.name << "\", \"det\": " << v.det << ", \"layer\": " << v.layer
       << ", \"sens\": " << v.sens << ", \"calo\": " << v.calo << ", \"rmin\": " << v.rmin << ", \"rmax\": " << v.rmax
       << ", \"zmin\": " << v.zmin << ", \"zmax\": " << v.zmax << "}" << (i + 1 < g.vols.size() ? ",\n" : "\n");
  }
  js << "  ]\n}\n";
}

void RunAction::EndOfRunAction(const G4Run*) {
  auto am = G4AnalysisManager::Instance();
  am->Write();
  am->CloseFile();
}

// ---------------- event ----------------
EventAction::EventAction(RunAction* run) : m_run(run) {}

void EventAction::BeginOfEventAction(const G4Event*) {
  const GeomInfo& g = geomInfo();
  for (auto& x : q) {
    x = Quirk();
    x.tExit.assign(g.boundaries.size(), -1.);
    x.edepDet.assign(g.dets.size(), 0.);
  }
  hits.clear();
  eDet.assign(g.dets.size(), 0.);
  eDetSens.assign(g.dets.size(), 0.);
  towers.clear();
  maxSep = 0;
  nCross = crossSign = 0;
  stepCount[0] = stepCount[1] = 0;
  m_wall0 = std::chrono::steady_clock::now();
  m_cpu0 = threadCpu();
}

void EventAction::EndOfEventAction(const G4Event* ev) {
  const QSimConfig& c = qsimConfig();
  const GeomInfo& g = geomInfo();
  auto am = G4AnalysisManager::Instance();
  OutBuffers& b = m_run->buffers();
  b.clear();
  const G4double cpu = threadCpu() - m_cpu0;
  const G4double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - m_wall0).count();
  auto info = dynamic_cast<const EventInfo*>(ev->GetUserInformation());
  const G4ThreeVector vtx = info ? info->vertex : G4ThreeVector();

  for (const auto& x : q) {
    b.q_pdg.push_back(x.pdg);
    b.q_px0.push_back(x.p0.px());
    b.q_py0.push_back(x.p0.py());
    b.q_pz0.push_back(x.p0.pz());
    b.q_e0.push_back(x.p0.e());
    b.q_x0.push_back(x.x0.x());
    b.q_y0.push_back(x.x0.y());
    b.q_z0.push_back(x.x0.z());
    b.q_xe.push_back(x.xEnd.x());
    b.q_ye.push_back(x.xEnd.y());
    b.q_ze.push_back(x.xEnd.z());
    b.q_te.push_back(x.tEnd);
    b.q_pxe.push_back(x.pEnd.px());
    b.q_pye.push_back(x.pEnd.py());
    b.q_pze.push_back(x.pEnd.pz());
    b.q_ee.push_back(x.pEnd.e());
    b.q_steps.push_back(x.steps);
    b.q_susp.push_back(x.susp);
    b.q_edep.push_back(x.edep);
    b.q_path.push_back(x.path);
    b.q_reason.push_back(x.reason);
    b.q_tExit.insert(b.q_tExit.end(), x.tExit.begin(), x.tExit.end());
    b.q_edepDet.insert(b.q_edepDet.end(), x.edepDet.begin(), x.edepDet.end());
  }
  b.eDet = eDet;
  b.eDetSens = eDetSens;
  for (const auto& h : hits) {
    b.h_q.push_back(h.q);
    b.h_vol.push_back(h.vol);
    b.h_det.push_back(g.vols[h.vol].det);
    b.h_layer.push_back(g.vols[h.vol].layer);
    b.h_nsteps.push_back(h.nsteps);
    b.h_x0.push_back(h.x0.x());
    b.h_y0.push_back(h.x0.y());
    b.h_z0.push_back(h.x0.z());
    b.h_t0.push_back(h.t0);
    b.h_x1.push_back(h.x1.x());
    b.h_y1.push_back(h.x1.y());
    b.h_z1.push_back(h.x1.z());
    b.h_t1.push_back(h.t1);
    b.h_px.push_back(h.p.x());
    b.h_py.push_back(h.p.y());
    b.h_pz.push_back(h.p.z());
    b.h_dE.push_back(h.dE);
    b.h_len.push_back(h.len);
  }
  for (const auto& [key, t] : towers) {
    b.tw_det.push_back(key / 1000000);
    b.tw_ieta.push_back((key / 1000) % 1000 - 500);
    b.tw_iphi.push_back(key % 1000);
    b.tw_e.push_back(t.e);
    b.tw_eq.push_back(t.eq);
  }
  am->FillNtupleIColumn(kEvents, 0, ev->GetEventID());
  am->FillNtupleIColumn(kEvents, 1, info ? info->genEvent : -1);
  am->FillNtupleIColumn(kEvents, 2, nCross);
  am->FillNtupleDColumn(kEvents, 3, vtx.x());
  am->FillNtupleDColumn(kEvents, 4, vtx.y());
  am->FillNtupleDColumn(kEvents, 5, vtx.z());
  am->FillNtupleDColumn(kEvents, 6, cpu);
  am->FillNtupleDColumn(kEvents, 7, wall);
  am->FillNtupleDColumn(kEvents, 8, maxSep);
  am->AddNtupleRow(kEvents);

  if (c.printEvery > 0 && ev->GetEventID() % c.printEvery == 0) {
    G4cout << std::setprecision(6) << "quirksim: event " << ev->GetEventID() << " gen " << (info ? info->genEvent : -1)
           << " steps " << q[0].steps << " " << q[1].steps << " suspensions " << q[0].susp << " " << q[1].susp
           << " max separation(mm) " << maxSep << " edep(MeV) " << q[0].edep << " " << q[1].edep << " time(s) " << cpu
           << " end " << q[0].reason << " " << q[1].reason << " hits " << hits.size() << G4endl;
  }
}

// ---------------- stepping ----------------
void SteppingAction::UserSteppingAction(const G4Step* step) {
  G4Track* t = step->GetTrack();
  const G4ParticleDefinition* def = t->GetParticleDefinition();
  if (Quirk::isQuirk(def)) {
    quirkStep(step, def->GetPDGEncoding() > 0 ? 0 : 1);
    return;
  }
  const G4double e = step->GetTotalEnergyDeposit();
  if (e > 0) {
    const G4StepPoint* pre = step->GetPreStepPoint();
    const VolInfo& vi = geomInfo().vols[pre->GetPhysicalVolume()->GetCopyNo()];
    m_ev->eDet[vi.det] += e;
    if (vi.sens)
      m_ev->eDetSens[vi.det] += e;
    if (vi.calo && qsimConfig().towers) {
      const G4ThreeVector x = 0.5 * (pre->GetPosition() + step->GetPostStepPoint()->GetPosition());
      const long ieta = std::lround(std::floor(x.eta() / kTowerEta)) + 500;
      const long iphi = std::lround(std::floor((x.phi() + CLHEP::pi) / (CLHEP::twopi / kTowerPhi))) % kTowerPhi;
      m_ev->towers[(vi.det * 1000L + ieta) * 1000L + iphi].e += e;
    }
  }
  const G4double tmax = qsimConfig().maxTime;
  if (tmax > 0 && step->GetPostStepPoint()->GetGlobalTime() > tmax)
    t->SetTrackStatus(fStopAndKill);
}

void SteppingAction::quirkStep(const G4Step* step, int i) {
  const QSimConfig& c = qsimConfig();
  const QuirkConfig& qc = quirkConfig();
  const GeomInfo& g = geomInfo();
  G4Track* t = step->GetTrack();
  EventAction::Quirk& Q = m_ev->q[i];
  const G4StepPoint* pre = step->GetPreStepPoint();
  const G4StepPoint* post = step->GetPostStepPoint();
  const G4int vol = pre->GetPhysicalVolume()->GetCopyNo();
  const VolInfo& vi = g.vols[vol];

  if (!Q.seen) {
    Q.seen = true;
    Q.pdg = t->GetParticleDefinition()->GetPDGEncoding();
    Q.p0 = G4LorentzVector(pre->GetMomentum(), pre->GetTotalEnergy());
    Q.x0 = pre->GetPosition();
  }
  ++Q.steps;
  ++m_ev->stepCount[i];
  G4TrackStatus status = t->GetTrackStatus();
  // CMSSW MaxTrackTime also applies to quirks; EndTracking frees the partner
  if (c.maxTime > 0 && post->GetGlobalTime() > c.maxTime && status != fStopAndKill) {
    t->SetTrackStatus(fStopAndKill);
    status = fStopAndKill;
    Q.reason = kTime;
  }
  if (status == fSuspend)
    ++Q.susp;  // as QuirkDebugWatcher
  const G4double e = step->GetTotalEnergyDeposit();
  Q.edep += e;
  Q.path += step->GetStepLength();
  Q.edepDet[vi.det] += e;
  m_ev->eDet[vi.det] += e;
  if (vi.sens)
    m_ev->eDetSens[vi.det] += e;
  const G4ThreeVector& x = post->GetPosition();
  if (e > 0 && vi.calo && c.towers) {
    const G4ThreeVector xm = 0.5 * (pre->GetPosition() + x);
    const long ieta = std::lround(std::floor(xm.eta() / kTowerEta)) + 500;
    const long iphi = std::lround(std::floor((xm.phi() + CLHEP::pi) / (CLHEP::twopi / kTowerPhi))) % kTowerPhi;
    auto& tw = m_ev->towers[(vi.det * 1000L + ieta) * 1000L + iphi];
    tw.e += e;
    tw.eq += e;
  }
  Q.xEnd = x;
  Q.tEnd = post->GetGlobalTime();
  Q.pEnd = G4LorentzVector(t->GetMomentum(), t->GetTotalEnergy());
  for (size_t b = 0; b < g.boundaries.size(); ++b) {
    if (Q.tExit[b] < 0 && (x.perp() > g.boundaries[b].r || std::abs(x.z()) > g.boundaries[b].halfz))
      Q.tExit[b] = Q.tEnd;
  }

  // hits: one per continuous visit of a sensitive volume, per quirk
  if (vi.sens) {
    if (Q.openVol != vol) {
      EventAction::Hit h;
      h.q = i;
      h.vol = vol;
      h.x0 = pre->GetPosition();
      h.t0 = pre->GetGlobalTime();
      h.p = pre->GetMomentum();
      m_ev->hits.push_back(h);
      Q.openHit = m_ev->hits.size() - 1;
      Q.openVol = vol;
    }
    EventAction::Hit& h = m_ev->hits[Q.openHit];
    ++h.nsteps;
    h.dE += e;
    h.len += step->GetStepLength();
    h.x1 = x;
    h.t1 = Q.tEnd;
    if (post->GetPhysicalVolume() != pre->GetPhysicalVolume() || status == fStopAndKill) {
      Q.openHit = Q.openVol = -1;
    }
  } else {
    Q.openHit = Q.openVol = -1;
  }

  // largest separation and crossings, from the last points of the two quirks (not at equal times)
  m_ev->lastPos[i] = x;
  if (m_ev->stepCount[0] > 0 && m_ev->stepCount[1] > 0) {
    const G4ThreeVector d = m_ev->lastPos[1] - m_ev->lastPos[0];
    m_ev->maxSep = std::max(m_ev->maxSep, d.mag());
    const G4ThreeVector axis = m_ev->q[1].p0.boostVector() - m_ev->q[0].p0.boostVector();
    const G4double proj = d.dot(axis);
    const G4int s = (proj > 0) - (proj < 0);
    if (s != 0) {
      if (m_ev->crossSign != 0 && s != m_ev->crossSign)
        ++m_ev->nCross;
      m_ev->crossSign = s;
    }
  }

  if (status == fStopAndKill && Q.reason == kAlive) {
    if (post->GetPhysicalVolume() == nullptr)
      Q.reason = kWorld;
    else if ((qc.maxRadius > 0 && x.perp() > qc.maxRadius) || (qc.maxZ > 0 && std::abs(x.z()) > qc.maxZ))
      Q.reason = kEnvelope;
    else if (t->GetKineticEnergy() <= 0.)
      Q.reason = kStopped;
    else if (!QuirkStringStore::instance().stringFor(Q.pdg).IsSourceAlive())
      Q.reason = kPartner;
    else
      Q.reason = kOther;
  }

  if (c.stepPrescale > 0 && (Q.steps % c.stepPrescale == 0 || status == fStopAndKill)) {
    auto am = G4AnalysisManager::Instance();
    am->FillNtupleIColumn(kSteps, 0, G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID());
    am->FillNtupleIColumn(kSteps, 1, i);
    am->FillNtupleIColumn(kSteps, 2, t->GetCurrentStepNumber());
    am->FillNtupleIColumn(kSteps, 3, vol);
    am->FillNtupleIColumn(kSteps, 4, status);
    am->FillNtupleDColumn(kSteps, 5, x.x());
    am->FillNtupleDColumn(kSteps, 6, x.y());
    am->FillNtupleDColumn(kSteps, 7, x.z());
    am->FillNtupleDColumn(kSteps, 8, Q.tEnd);
    am->FillNtupleDColumn(kSteps, 9, Q.pEnd.px());
    am->FillNtupleDColumn(kSteps, 10, Q.pEnd.py());
    am->FillNtupleDColumn(kSteps, 11, Q.pEnd.pz());
    am->FillNtupleDColumn(kSteps, 12, t->GetKineticEnergy());
    am->FillNtupleDColumn(kSteps, 13, e);
    am->AddNtupleRow(kSteps);
  }
}

// ---------------- stacking (only with killSecondaries) ----------------
G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(const G4Track* t) {
  // a resumed (suspended) quirk is not new: keep the default set by QuirkPhysics
  if (QuirkTrackStatus::isSuspended(t))
    return stackManager->GetDefaultClassification();
  if (t->GetParentID() > 0 && !Quirk::isQuirk(t->GetParticleDefinition()))
    return fKill;
  return stackManager->GetDefaultClassification();
}

// ---------------- initialization ----------------
void ActionInitialization::BuildForMaster() const { SetUserAction(new RunAction); }

void ActionInitialization::Build() const {
  auto run = new RunAction;
  SetUserAction(run);
  auto ev = new EventAction(run);
  SetUserAction(ev);
  SetUserAction(new SteppingAction(ev));
  SetUserAction(new PrimaryGenerator);
  if (qsimConfig().killSecondaries)
    SetUserAction(new StackingAction);
}
