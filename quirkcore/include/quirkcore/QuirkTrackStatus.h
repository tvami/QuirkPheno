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
