#ifndef CYLCOWWLIDS_H
#define CYLCOWWLIDS_H

class TrackerInfo;

// Cylindrical Cow with Lids
//
// Intended coverage: |eta| < 2.4 with D_z_beam_spot = +-3 cm (3 sigma)
// B-layer extends to 2.55.
// Layers 1 and 2 have somewhat longer barrels. It is assumed
// those will be needed / used for seed finding.
//
// Layers 3 - 9:
//   Barrel:     0.0 - 1.0
//   Transition: 1.0 - 1.4
//   Endcap:     1.4 - 2.4
//
// Run root test/CylCowWLids.C to get a plot and dumps of
// edge coordinates and etas.


// This could be translated into plugin

void Create_TrackerInfo(TrackerInfo& ti, bool verbose=false);

#endif
