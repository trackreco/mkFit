#ifndef _buildtest_
#define _buildtest_

#include <map>
#include <string>
#include <vector>

#include "Track.h"
#include "Geometry.h"

class TFile;
class TTree;
class TH1F;

typedef std::pair<unsigned int,unsigned int> BinInfo;

unsigned int getPhiPartition(float phi);
void runBuildingTest(bool saveTree, unsigned int nevts, Geometry* theGeom);
void runBuildingTestEvt(bool saveTree, TTree *tree,unsigned int& tk_nhits, float& tk_chi2, 
			TTree *tree_br, unsigned int& layer, unsigned int& branches, unsigned int& cands, 
			std::map<std::string,TH1F*>& validation_hists, Geometry* theGeom);
void buildTestSerial(std::vector<Track>& evt_seeds,std::vector<Track>& evt_track_candidates,
		     std::vector<std::vector<Hit> >& evt_lay_hits,std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,
		     const int& nhits_per_seed,const unsigned int& maxCand,const float& chi2Cut,const float& nSigma,const float& minDPhi,
		     SMatrix36& projMatrix36,SMatrix63& projMatrix36T,bool debug,
		     bool saveTree, TTree *tree_br, unsigned int& layer, unsigned int& branches, unsigned int& cands, Geometry* theGeom);
void buildTestParallel(std::vector<Track>& evt_seeds,std::vector<Track>& evt_track_candidates,
		       std::vector<std::vector<Hit> >& evt_lay_hits,std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,
		       const int& nhits_per_seed,const unsigned int& maxCand,const float& chi2Cut,const float& nSigma,const float& minDPhi,
		       SMatrix36& projMatrix36,SMatrix63& projMatrix36T,bool debug, Geometry* theGeom);
void processCandidates(std::pair<Track, TrackState>& cand,std::vector<std::pair<Track, TrackState> >& tmp_candidates,
		       unsigned int ilay,std::vector<std::vector<Hit> >& evt_lay_hits,std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,
		       const int& nhits_per_seed,const unsigned int& maxCand,const float& chi2Cut,const float& nSigma,const float& minDPhi,
		       SMatrix36& projMatrix36,SMatrix63& projMatrix36T,bool debug, Geometry* theGeom);
#endif
