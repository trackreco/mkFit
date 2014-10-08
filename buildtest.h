#ifndef _buildtest_
#define _buildtest_

#include <map>
#include <string>
#include <vector>

#include "Track.h"

class TTree;
class TH1F;

typedef std::pair<unsigned int,unsigned int> BinInfo;

unsigned int getPhiPartition(float phi);
void runBuildingTest(bool saveTree, unsigned int nevts);
void runBuildingTestEvt(bool saveTree, TTree *tree,unsigned int& tk_nhits, float& tk_chi2, std::map<std::string,TH1F*>& validation_hists);
void buildTestSerial(std::vector<Track>& evt_seeds,std::vector<Track>& evt_track_candidates,
		     std::vector<std::vector<Hit> >& evt_lay_hits,std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,
		     const int& nhits_per_seed,const unsigned int& maxCand,const float& chi2Cut,const float& nSigma,const float& minDPhi,
		     SMatrix36& projMatrix36,SMatrix63& projMatrix36T,bool debug);
void buildTestParallel(std::vector<Track>& evt_seeds,std::vector<Track>& evt_track_candidates,
		       std::vector<std::vector<Hit> >& evt_lay_hits,std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,
		       const int& nhits_per_seed,const unsigned int& maxCand,const float& chi2Cut,const float& nSigma,const float& minDPhi,
		       SMatrix36& projMatrix36,SMatrix63& projMatrix36T,bool debug);
void processCandidates(std::pair<Track, TrackState>& cand,std::vector<std::pair<Track, TrackState> >& tmp_candidates,
		       unsigned int ilay,std::vector<std::vector<Hit> >& evt_lay_hits,std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,
		       const int& nhits_per_seed,const unsigned int& maxCand,const float& chi2Cut,const float& nSigma,const float& minDPhi,
		       SMatrix36& projMatrix36,SMatrix63& projMatrix36T,bool debug);
#endif

