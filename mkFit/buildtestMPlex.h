#ifndef _buildtest_mplex_
#define _buildtest_mplex_

#include "Track.h"

/*
void   generateTracks(std::vector<Track>& simtracks, int Ntracks);
void   make_validation_tree(const char         *fname,
                            std::vector<Track> &simtracks,
                            std::vector<Track> &rectracks);
*/

typedef std::pair<unsigned int,unsigned int> BinInfo;

double runBuildingTest(std::vector<Track>& simtracks/*, std::vector<Track>& rectracks*/);

void buildTestParallel(std::vector<Track>& evt_seeds,std::vector<Track>& evt_track_candidates,
		       std::vector<std::vector<Hit> >& evt_lay_hits,std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,
		       const int& nhits_per_seed,const unsigned int& maxCand,const float& chi2Cut,const float& nSigma,const float& minDPhi,
		       SMatrix36& projMatrix36,SMatrix63& projMatrix36T);
void processCandidates(std::pair<Track, TrackState>& cand,std::vector<std::pair<Track, TrackState> >& tmp_candidates,
		       unsigned int ilay,std::vector<std::vector<Hit> >& evt_lay_hits,std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,
		       const int& nhits_per_seed,const unsigned int& maxCand,const float& chi2Cut,const float& nSigma,const float& minDPhi,
		       SMatrix36& projMatrix36,SMatrix63& projMatrix36T);

double runBuildingTestPlex(std::vector<Track>& simtracks/*, std::vector<Track>& rectracks*/);
double runBuildingTestPlexBestHit(std::vector<Track>& simtracks/*, std::vector<Track>& rectracks*/);

#endif
