#ifndef _buildtest_
#define _buildtest_

#include <map>
#include <string>

#include "Track.h"
#include "Matrix.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

typedef std::pair<unsigned int,unsigned int> BinInfo;

bool sortByPhi(Hit hit1,Hit hit2);
unsigned int getPhiPartition(float phi);
void runBuildingTest(bool saveTree, unsigned int nevts);
void runBuildingTest(bool saveTree, TTree *tree,unsigned int& tk_nhits, float& tk_chi2, std::map<std::string,TH1F*>& validation_hists);
void buildTestSerial(std::vector<Track>& evt_seeds,std::vector<Track>& evt_track_candidates,
		     std::vector<std::vector<Hit> >& evt_lay_hits,std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,
		     const int& nhits_per_seed,const unsigned int& maxCand,
		     SMatrix36& projMatrix36,SMatrix63& projMatrix36T,bool debug);

#endif
