// Derived from atlas/athena Simulation/G4Extensions/Quirks (Copyright (C) 2002-2024 CERN for the
// benefit of the ATLAS collaboration) and cms-sw/cmssw, both Apache-2.0; modified.
#ifndef QUIRKCORE_QuirkTrackStatus_H
#define QUIRKCORE_QuirkTrackStatus_H

// Track status helpers used by the quirk code.

#include "G4Track.hh"
#include "G4TrackStatus.hh"

struct QuirkTrackStatus {
  // suspended tracks are resumed later (quirk pairs), not new tracks
  static bool isSuspended(G4TrackStatus status) { return status == fSuspend || status == fSuspendAndWait; }
  static bool isSuspended(const G4Track* aTrack) { return isSuspended(aTrack->GetTrackStatus()); }
};

#endif
