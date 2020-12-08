#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <list>
#include <unordered_map>
#include "Event.h"
#include "LayerNumberConverter.h"

using namespace mkfit;

constexpr bool useMatched = false;

constexpr int cleanSimTrack_minSimHits = 3;
constexpr int cleanSimTrack_minRecHits = 2;

//check if this is the same as in the release
enum class HitType {
  Pixel = 0,
  Strip = 1,
  Glued = 2,
  Invalid = 3,
  Phase2OT = 4,
  Unknown = 99
};

/// track algorithm; partial copy from TrackBase.h
enum class TrackAlgorithm {
  undefAlgorithm = 0,
  ctf = 1, 
  duplicateMerge = 2,
  cosmics = 3,
  initialStep = 4,
  lowPtTripletStep = 5,
  pixelPairStep = 6,
  detachedTripletStep = 7,
  mixedTripletStep = 8,
  pixelLessStep = 9,
  tobTecStep = 10,
  jetCoreRegionalStep = 11,
  conversionStep = 12,
  muonSeededStepInOut = 13,
  muonSeededStepOutIn = 14,
  outInEcalSeededConv = 15,
  inOutEcalSeededConv = 16,
  nuclInter = 17,
  standAloneMuon = 18,
  globalMuon = 19,
  cosmicStandAloneMuon = 20,
  cosmicGlobalMuon = 21,
  // Phase1
  highPtTripletStep = 22,
  lowPtQuadStep = 23,
  detachedQuadStep = 24,
  reservedForUpgrades1 = 25,
  reservedForUpgrades2 = 26,
  bTagGhostTracks = 27,
  beamhalo = 28,
  gsf = 29,
  // HLT algo name
  hltPixel = 30,
  // steps used by PF
  hltIter0 = 31,
  hltIter1 = 32,
  hltIter2 = 33,
  hltIter3 = 34,
  hltIter4 = 35,
  // steps used by all other objects @HLT
  hltIterX = 36,
  // steps used by HI muon regional iterative tracking
  hiRegitMuInitialStep = 37,
  hiRegitMuLowPtTripletStep = 38,
  hiRegitMuPixelPairStep = 39,
  hiRegitMuDetachedTripletStep = 40,
  hiRegitMuMixedTripletStep = 41,
  hiRegitMuPixelLessStep = 42,
  hiRegitMuTobTecStep = 43,
  hiRegitMuMuonSeededStepInOut = 44,
  hiRegitMuMuonSeededStepOutIn = 45,
  algoSize = 46
};

typedef std::list<std::string> lStr_t;
typedef lStr_t::iterator       lStr_i;
void next_arg_or_die(lStr_t& args, lStr_i& i)
{
  lStr_i j = i;
  if (++j == args.end() || ((*j)[0] == '-' ))
    {
      std::cerr <<"Error: option "<< *i <<" requires an argument.\n";
      exit(1);
    }
  i = j;
}

bool next_arg_option(lStr_t& args, lStr_i& i)
{
  lStr_i j = i;
  if (++j == args.end() || ((*j)[0] == '-' ))
    {
      return false;
    }
  i = j;
  return true;
}

void printHelp(const char* av0){
  printf(
	 "Usage: %s [options]\n"
	 "Options:\n"
	 "  --input          <str>    input file\n"
	 "  --output         <str>    output file\n"
	 "  --verbosity      <num>    print details (0 quiet, 1 print counts, 2 print all; def: 0)\n"
	 "  --maxevt         <num>    maxevt events to write (-1 for everything in the file def: -1)\n"
	 "  --clean-sim-tracks        apply sim track cleaning (def: no cleaning)\n"
	 "  --write-all-events        write all events (def: skip events with 0 simtracks or seeds)\n"
	 "  --write-rec-tracks        write rec tracks (def: not written)\n"
	 "  --apply-ccc               apply cluster charge cut to strip hits (def: false)\n"
	 , av0);
}

int main(int argc, char *argv[])
{
  bool haveInput = false;
  std::string inputFileName;
  bool haveOutput = false;
  std::string outputFileName;

  bool cleanSimTracks = false;
  bool writeAllEvents = false;
  bool writeRecTracks = false;
  bool applyCCC       = false;

  int verbosity = 0;
  long long maxevt = -1;

  int cutValueCCC = 1620; //Nominal value (from first iteration of CMSSW) is 1620 

  lStr_t mArgs;
  for (int i = 1; i < argc; ++i)
    {
      mArgs.push_back(argv[i]);
    }

  lStr_i i  = mArgs.begin();
  while (i != mArgs.end())
    {
      lStr_i start = i;
      
      if (*i == "-h" || *i == "-help" || *i == "--help")
	{
	  printHelp(argv[0]);
	}
      else if (*i == "--input")
	{
	  next_arg_or_die(mArgs, i);
	  inputFileName = *i;
	  haveInput = true;
	}
      else if (*i == "--output")
	{
	  next_arg_or_die(mArgs, i);
	  outputFileName = *i;
	  haveOutput = true;
	}
      else if (*i == "--verbosity")
	{
	  next_arg_or_die(mArgs, i);
	  verbosity = std::atoi(i->c_str());
	}
      else if (*i == "--maxevt")
	{
	  next_arg_or_die(mArgs, i);
	  maxevt = std::atoi(i->c_str());
	}
      else if (*i == "--clean-sim-tracks")
	{
	  cleanSimTracks = true;
	}
      else if (*i == "--write-all-events")
	{
	  writeAllEvents = true;
	}
      else if (*i == "--write-rec-tracks")
	{
	  writeRecTracks = true;
	}
      else if (*i == "--apply-ccc")
	{
	  applyCCC = true;
	  if( next_arg_option(mArgs, i))
	    {
	      cutValueCCC = std::atoi(i->c_str());
	    }
	}
      else
	{
	  fprintf(stderr, "Error: Unknown option/argument '%s'.\n", i->c_str());
	  printHelp(argv[0]);
	  exit(1);
	}
      mArgs.erase(start, ++i);
    }//while arguments

  if (not haveOutput or not haveInput)
    {
      fprintf(stderr, "Error: both input and output are required\n");
      printHelp(argv[0]);
      exit(1);
    }

  using namespace std;

  LayerNumberConverter lnc(TkLayout::phase1);
  const unsigned int nTotalLayers = lnc.nLayers();
  Config::nTotalLayers = lnc.nLayers();

  vector<unordered_map<unsigned int, unsigned int>> module_shortId_hash(Config::nTotalLayers);

  int nstot = 0;
  std::vector<int> nhitstot(nTotalLayers, 0);

  TFile* f = TFile::Open(inputFileName.c_str());
  if (f == 0)
  {
    fprintf(stderr, "Failed opening input root file '%s'\n", inputFileName.c_str());
    exit(1);
  }
  
  TTree* t = (TTree*) f->Get("trackingNtuple/tree");

  unsigned long long event;
  t->SetBranchAddress("event",&event);
  
  //sim tracks
  std::vector<float>* sim_eta = 0;
  std::vector<float>* sim_px = 0;
  std::vector<float>* sim_py = 0;
  std::vector<float>* sim_pz = 0;
  std::vector<int>*   sim_parentVtxIdx = 0;
  std::vector<int>*   sim_q = 0;
  std::vector<int>*   sim_event = 0;
  std::vector<int>*   sim_bunchCrossing = 0;
  std::vector<int>*   sim_nValid = 0;//simHit count, actually
  t->SetBranchAddress("sim_eta",&sim_eta);
  t->SetBranchAddress("sim_px",&sim_px);
  t->SetBranchAddress("sim_py",&sim_py);
  t->SetBranchAddress("sim_pz",&sim_pz);
  t->SetBranchAddress("sim_parentVtxIdx",&sim_parentVtxIdx);
  t->SetBranchAddress("sim_q",&sim_q);
  t->SetBranchAddress("sim_event",&sim_event);
  t->SetBranchAddress("sim_bunchCrossing",&sim_bunchCrossing);
  t->SetBranchAddress("sim_nValid",&sim_nValid);

  std::vector<vector<int> >*   sim_trkIdx = 0;
  t->SetBranchAddress("sim_trkIdx", &sim_trkIdx);

  //simvtx
  std::vector<float>* simvtx_x = 0;
  std::vector<float>* simvtx_y = 0;
  std::vector<float>* simvtx_z = 0;
  t->SetBranchAddress("simvtx_x"       , &simvtx_x);
  t->SetBranchAddress("simvtx_y"       , &simvtx_y);
  t->SetBranchAddress("simvtx_z"       , &simvtx_z);


  //simhit
  std::vector<short>* simhit_process = 0;
  std::vector<int>* simhit_particle = 0;
  std::vector<int>* simhit_simTrkIdx = 0;
  std::vector<float>* simhit_x = 0;
  std::vector<float>* simhit_y = 0;
  std::vector<float>* simhit_z = 0;
  std::vector<float>* simhit_px = 0;
  std::vector<float>* simhit_py = 0;
  std::vector<float>* simhit_pz = 0;
  t->SetBranchAddress("simhit_process",   &simhit_process);
  t->SetBranchAddress("simhit_particle",   &simhit_particle);
  t->SetBranchAddress("simhit_simTrkIdx", &simhit_simTrkIdx);
  t->SetBranchAddress("simhit_x",        &simhit_x);
  t->SetBranchAddress("simhit_y",        &simhit_y);
  t->SetBranchAddress("simhit_z",        &simhit_z);
  t->SetBranchAddress("simhit_px",        &simhit_px);
  t->SetBranchAddress("simhit_py",        &simhit_py);
  t->SetBranchAddress("simhit_pz",        &simhit_pz);

  std::vector<std::vector<int> >* simhit_hitIdx = 0;
  t->SetBranchAddress("simhit_hitIdx", &simhit_hitIdx);
  std::vector<std::vector<int> >* simhit_hitType = 0;
  t->SetBranchAddress("simhit_hitType", &simhit_hitType);

  //rec tracks
  std::vector<int>*                trk_q = 0;
  std::vector<unsigned int>*       trk_nValid = 0;
  std::vector<int>*                trk_seedIdx = 0;
  std::vector<unsigned long long>* trk_algoMask = 0;
  std::vector<unsigned int>*       trk_algo = 0   ;
  std::vector<unsigned int>*       trk_originalAlgo = 0;
  std::vector<float>*              trk_nChi2 = 0;
  std::vector<float>*              trk_px = 0;
  std::vector<float>*              trk_py = 0;
  std::vector<float>*              trk_pz = 0;
  std::vector<float>*              trk_pt = 0;
  std::vector<float>*              trk_phi = 0;
  std::vector<float>*              trk_lambda = 0;
  std::vector<float>*              trk_refpoint_x = 0;
  std::vector<float>*              trk_refpoint_y = 0;
  std::vector<float>*              trk_refpoint_z = 0;
  std::vector<float>*              trk_dxyErr = 0;
  std::vector<float>*              trk_dzErr = 0;
  std::vector<float>*              trk_ptErr = 0;
  std::vector<float>*              trk_phiErr = 0;
  std::vector<float>*              trk_lambdaErr = 0;
  t->SetBranchAddress("trk_q",   &trk_q);
  t->SetBranchAddress("trk_nValid",   &trk_nValid);
  t->SetBranchAddress("trk_seedIdx",  &trk_seedIdx);
  t->SetBranchAddress("trk_algoMask", &trk_algoMask);
  t->SetBranchAddress("trk_algo", &trk_algo);
  t->SetBranchAddress("trk_originalAlgo", &trk_originalAlgo);
  t->SetBranchAddress("trk_nChi2", &trk_nChi2);
  t->SetBranchAddress("trk_px", &trk_px);
  t->SetBranchAddress("trk_py", &trk_py);
  t->SetBranchAddress("trk_pz", &trk_pz);
  t->SetBranchAddress("trk_pt", &trk_pt);
  t->SetBranchAddress("trk_phi", &trk_phi);
  t->SetBranchAddress("trk_lambda", &trk_lambda);
  t->SetBranchAddress("trk_refpoint_x", &trk_refpoint_x);
  t->SetBranchAddress("trk_refpoint_y", &trk_refpoint_y);
  t->SetBranchAddress("trk_refpoint_z", &trk_refpoint_z);
  t->SetBranchAddress("trk_dxyErr", &trk_dxyErr);
  t->SetBranchAddress("trk_dzErr", &trk_dzErr);
  t->SetBranchAddress("trk_ptErr", &trk_ptErr);
  t->SetBranchAddress("trk_phiErr", &trk_phiErr);
  t->SetBranchAddress("trk_lambdaErr", &trk_lambdaErr);

  std::vector<std::vector<int> >* trk_hitIdx = 0;
  t->SetBranchAddress("trk_hitIdx", &trk_hitIdx);
  std::vector<std::vector<int> >* trk_hitType = 0;
  t->SetBranchAddress("trk_hitType", &trk_hitType);

  //seeds
  std::vector<float>*   see_stateTrajGlbX = 0;
  std::vector<float>*   see_stateTrajGlbY = 0;
  std::vector<float>*   see_stateTrajGlbZ = 0;
  std::vector<float>*   see_stateTrajGlbPx = 0;
  std::vector<float>*   see_stateTrajGlbPy = 0;
  std::vector<float>*   see_stateTrajGlbPz = 0;
  std::vector<float>*   see_eta= 0;//PCA parameters
  std::vector<float>*   see_pt = 0;//PCA parameters
  std::vector<float>*   see_stateCcov00 = 0;
  std::vector<float>*   see_stateCcov01 = 0;
  std::vector<float>*   see_stateCcov02 = 0;
  std::vector<float>*   see_stateCcov03 = 0;
  std::vector<float>*   see_stateCcov04 = 0;
  std::vector<float>*   see_stateCcov05 = 0;
  std::vector<float>*   see_stateCcov11 = 0;
  std::vector<float>*   see_stateCcov12 = 0;
  std::vector<float>*   see_stateCcov13 = 0;
  std::vector<float>*   see_stateCcov14 = 0;
  std::vector<float>*   see_stateCcov15 = 0;
  std::vector<float>*   see_stateCcov22 = 0;
  std::vector<float>*   see_stateCcov23 = 0;
  std::vector<float>*   see_stateCcov24 = 0;
  std::vector<float>*   see_stateCcov25 = 0;
  std::vector<float>*   see_stateCcov33 = 0;
  std::vector<float>*   see_stateCcov34 = 0;
  std::vector<float>*   see_stateCcov35 = 0;
  std::vector<float>*   see_stateCcov44 = 0;
  std::vector<float>*   see_stateCcov45 = 0;
  std::vector<float>*   see_stateCcov55 = 0;
  std::vector<int>*     see_q = 0;
  std::vector<unsigned int>*     see_algo = 0;
  t->SetBranchAddress("see_stateTrajGlbX",&see_stateTrajGlbX);
  t->SetBranchAddress("see_stateTrajGlbY",&see_stateTrajGlbY);
  t->SetBranchAddress("see_stateTrajGlbZ",&see_stateTrajGlbZ);
  t->SetBranchAddress("see_stateTrajGlbPx",&see_stateTrajGlbPx);
  t->SetBranchAddress("see_stateTrajGlbPy",&see_stateTrajGlbPy);
  t->SetBranchAddress("see_stateTrajGlbPz",&see_stateTrajGlbPz);
  t->SetBranchAddress("see_eta",&see_eta);
  t->SetBranchAddress("see_pt",&see_pt);
  t->SetBranchAddress("see_stateCcov00",&see_stateCcov00);
  t->SetBranchAddress("see_stateCcov01",&see_stateCcov01);
  t->SetBranchAddress("see_stateCcov02",&see_stateCcov02);
  t->SetBranchAddress("see_stateCcov03",&see_stateCcov03);
  t->SetBranchAddress("see_stateCcov04",&see_stateCcov04);
  t->SetBranchAddress("see_stateCcov05",&see_stateCcov05);
  t->SetBranchAddress("see_stateCcov11",&see_stateCcov11);
  t->SetBranchAddress("see_stateCcov12",&see_stateCcov12);
  t->SetBranchAddress("see_stateCcov13",&see_stateCcov13);
  t->SetBranchAddress("see_stateCcov14",&see_stateCcov14);
  t->SetBranchAddress("see_stateCcov15",&see_stateCcov15);
  t->SetBranchAddress("see_stateCcov22",&see_stateCcov22);
  t->SetBranchAddress("see_stateCcov23",&see_stateCcov23);
  t->SetBranchAddress("see_stateCcov24",&see_stateCcov24);
  t->SetBranchAddress("see_stateCcov25",&see_stateCcov25);
  t->SetBranchAddress("see_stateCcov33",&see_stateCcov33);
  t->SetBranchAddress("see_stateCcov34",&see_stateCcov34);
  t->SetBranchAddress("see_stateCcov35",&see_stateCcov35);
  t->SetBranchAddress("see_stateCcov44",&see_stateCcov44);
  t->SetBranchAddress("see_stateCcov45",&see_stateCcov45);
  t->SetBranchAddress("see_stateCcov55",&see_stateCcov55);
  t->SetBranchAddress("see_q",&see_q);
  t->SetBranchAddress("see_algo",&see_algo);

  std::vector<std::vector<int> >* see_hitIdx = 0;
  t->SetBranchAddress("see_hitIdx", &see_hitIdx);
  std::vector<std::vector<int> >* see_hitType = 0;
  t->SetBranchAddress("see_hitType", &see_hitType);

  //pixel hits
  vector<unsigned short>*    pix_det = 0;
  vector<unsigned short>*    pix_lay = 0;
  vector<unsigned int>*      pix_detId = 0;
  vector<float>*  pix_x = 0;
  vector<float>*  pix_y = 0;
  vector<float>*  pix_z = 0;
  vector<float>*  pix_xx = 0;
  vector<float>*  pix_xy = 0;
  vector<float>*  pix_yy = 0;
  vector<float>*  pix_yz = 0;
  vector<float>*  pix_zz = 0;
  vector<float>*  pix_zx = 0;
  vector<int>*    pix_csize_col = 0;
  vector<int>*    pix_csize_row = 0;
  //these were renamed in CMSSW_9_1_0: auto-detect
  bool has910_det_lay = t->GetBranch("pix_det") == nullptr;
  if (has910_det_lay){
    t->SetBranchAddress("pix_subdet",&pix_det);
    t->SetBranchAddress("pix_layer",&pix_lay);
  } else {
    t->SetBranchAddress("pix_det",&pix_det);
    t->SetBranchAddress("pix_lay",&pix_lay);
  }
  t->SetBranchAddress("pix_detId",&pix_detId);
  t->SetBranchAddress("pix_x",&pix_x);
  t->SetBranchAddress("pix_y",&pix_y);
  t->SetBranchAddress("pix_z",&pix_z);
  t->SetBranchAddress("pix_xx",&pix_xx);
  t->SetBranchAddress("pix_xy",&pix_xy);
  t->SetBranchAddress("pix_yy",&pix_yy);
  t->SetBranchAddress("pix_yz",&pix_yz);
  t->SetBranchAddress("pix_zz",&pix_zz);
  t->SetBranchAddress("pix_zx",&pix_zx);
  t->SetBranchAddress("pix_clustSizeCol",&pix_csize_col);
  t->SetBranchAddress("pix_clustSizeRow",&pix_csize_row);

  vector<vector<int> >*    pix_simHitIdx = 0;
  t->SetBranchAddress("pix_simHitIdx", &pix_simHitIdx);
  vector<vector<float> >*    pix_chargeFraction = 0;
  t->SetBranchAddress("pix_chargeFraction", &pix_chargeFraction);

  //strip hits
  vector<short>*  glu_isBarrel = 0;
  vector<unsigned int>*    glu_det = 0;
  vector<unsigned int>*    glu_lay = 0;
  vector<unsigned int>*    glu_detId = 0;
  vector<int>*    glu_monoIdx = 0;
  vector<int>*    glu_stereoIdx = 0;
  vector<float>*  glu_x = 0;
  vector<float>*  glu_y = 0;
  vector<float>*  glu_z = 0;
  vector<float>*  glu_xx = 0;
  vector<float>*  glu_xy = 0;
  vector<float>*  glu_yy = 0;
  vector<float>*  glu_yz = 0;
  vector<float>*  glu_zz = 0;
  vector<float>*  glu_zx = 0;
  t->SetBranchAddress("glu_isBarrel",&glu_isBarrel);
  if (has910_det_lay){
    t->SetBranchAddress("glu_subdet",&glu_det);
    t->SetBranchAddress("glu_layer",&glu_lay);
  } else {
    t->SetBranchAddress("glu_det",&glu_det);
    t->SetBranchAddress("glu_lay",&glu_lay);
  }
  t->SetBranchAddress("glu_detId",&glu_detId);
  t->SetBranchAddress("glu_monoIdx",&glu_monoIdx);
  t->SetBranchAddress("glu_stereoIdx",&glu_stereoIdx);
  t->SetBranchAddress("glu_x",&glu_x);
  t->SetBranchAddress("glu_y",&glu_y);
  t->SetBranchAddress("glu_z",&glu_z);
  t->SetBranchAddress("glu_xx",&glu_xx);
  t->SetBranchAddress("glu_xy",&glu_xy);
  t->SetBranchAddress("glu_yy",&glu_yy);
  t->SetBranchAddress("glu_yz",&glu_yz);
  t->SetBranchAddress("glu_zz",&glu_zz);
  t->SetBranchAddress("glu_zx",&glu_zx);

  vector<short>*    str_isBarrel = 0;
  vector<short>*    str_isStereo = 0;
  vector<unsigned int>*    str_det = 0;
  vector<unsigned int>*    str_lay = 0;
  vector<unsigned int>*    str_detId = 0;
  vector<unsigned int>*    str_simType = 0;
  vector<float>*  str_x = 0;
  vector<float>*  str_y = 0;
  vector<float>*  str_z = 0;
  vector<float>*  str_xx = 0;
  vector<float>*  str_xy = 0;
  vector<float>*  str_yy = 0;
  vector<float>*  str_yz = 0;
  vector<float>*  str_zz = 0;
  vector<float>*  str_zx = 0;
  vector<float>*  str_chargePerCM = 0;
  vector<int>*    str_csize = 0;
  t->SetBranchAddress("str_isBarrel",&str_isBarrel);
  t->SetBranchAddress("str_isStereo",&str_isStereo);
  if (has910_det_lay){
    t->SetBranchAddress("str_subdet",&str_det);
    t->SetBranchAddress("str_layer",&str_lay);
  } else {
    t->SetBranchAddress("str_det",&str_det);
    t->SetBranchAddress("str_lay",&str_lay);
  }
  t->SetBranchAddress("str_detId",&str_detId);
  t->SetBranchAddress("str_simType",&str_simType);
  t->SetBranchAddress("str_x",&str_x);
  t->SetBranchAddress("str_y",&str_y);
  t->SetBranchAddress("str_z",&str_z);
  t->SetBranchAddress("str_xx",&str_xx);
  t->SetBranchAddress("str_xy",&str_xy);
  t->SetBranchAddress("str_yy",&str_yy);
  t->SetBranchAddress("str_yz",&str_yz);
  t->SetBranchAddress("str_zz",&str_zz);
  t->SetBranchAddress("str_zx",&str_zx);
  t->SetBranchAddress("str_chargePerCM",&str_chargePerCM);
  t->SetBranchAddress("str_clustSize", &str_csize);

  vector<vector<int> >*    str_simHitIdx = 0;
  t->SetBranchAddress("str_simHitIdx", &str_simHitIdx);
  vector<vector<float> >*    str_chargeFraction = 0;
  t->SetBranchAddress("str_chargeFraction", &str_chargeFraction);

  long long totentries = t->GetEntries();
  long long savedEvents = 0;

  DataFile data_file;
  int outOptions = DataFile::ES_Seeds;
  if (writeRecTracks) outOptions |= DataFile::ES_CmsswTracks;
  if (maxevt < 0) maxevt = totentries;
  data_file.OpenWrite(outputFileName, std::min(maxevt, totentries), outOptions);

  Event EE(0);

  int numFailCCC = 0;
  int numTotalStr = 0;
  // gDebug = 8;

  for (long long i = 0; savedEvents < maxevt && i<totentries && i<maxevt; ++i)
  {
    EE.Reset(i);

    cout << "process entry i=" << i << " out of " << totentries << ", saved so far " << savedEvents << ", with max=" << maxevt << endl;

    t->GetEntry(i);

    cout << "edm event=" << event << endl;

    for (unsigned int istr = 0; istr < str_lay->size(); ++istr) {
      if(str_chargePerCM->at(istr) < cutValueCCC) numFailCCC++;
      numTotalStr++;
    }

    auto nSims = sim_q->size();
    if (nSims==0) {
      cout << "branches not loaded" << endl; exit(1);
    }
    if (verbosity>0) std::cout<<__FILE__<<" "<<__LINE__
			      <<" nSims "<<nSims
			      <<" nSeeds "<<see_q->size()
			      <<" nRecT "<<trk_q->size()
			      <<std::endl;

    //find best matching tkIdx from a list of simhits indices
    auto bestTkIdx = [&](std::vector<int> const& shs, std::vector<float> const& shfs, int rhIdx, HitType rhType){
      //assume that all simhits are associated
      int ibest = -1;
      int shbest = -1;
      float hpbest = -1;
      float tpbest = -1;
      float hfbest = -1;
      
      float maxfrac = -1;
      int ish = -1;
      int nshs = shs.size();
      for (auto const sh : shs){
	ish++;
	auto tkidx = simhit_simTrkIdx->at(sh);
	//use only sh with available TP
	if (tkidx < 0) continue;
	
	auto hpx = simhit_px->at(sh);
	auto hpy = simhit_py->at(sh);
	auto hpz = simhit_pz->at(sh);
	auto hp = sqrt(hpx*hpx + hpy*hpy + hpz*hpz);

	//look only at hits with p> 50 MeV
	if (hp < 0.05f) continue;
	
	auto tpx = sim_px->at(tkidx);
	auto tpy = sim_py->at(tkidx);
	auto tpz = sim_pz->at(tkidx);
	auto tp = sqrt(tpx*tpx + tpy*tpy + tpz*tpz);

	//take only hits with hp> 0.5*tp
	if (hp < 0.5*tp) continue;

	//pick tkidx corresponding to max hp/tp; .. this is probably redundant
	if (maxfrac < hp/tp){
	  maxfrac = hp/tp;
	  ibest = tkidx;
	  shbest = sh;
	  hpbest = hp;
	  tpbest = tp;
	  hfbest = shfs[ish];
	}
	
      }

      //arbitration: a rechit with one matching sim is matched to sim if it's the first
      //FIXME: SOME BETTER SELECTION CAN BE DONE (it will require some more correlated knowledge)
      if (nshs == 1 && ibest >= 0){
	auto const& srhIdxV = simhit_hitIdx->at(shbest);
	auto const& srhTypeV = simhit_hitType->at(shbest);
	int ih = -1;
	for (auto itype : srhTypeV){
	  ih++;
	  if (HitType(itype) == rhType && srhIdxV[ih] != rhIdx){
	    ibest = -1;
	    break;
	  }
	}
      }
      
      if (ibest >= 0 && false){
	std::cout<<" best tkIdx "<<ibest<<" rh "<<rhIdx <<" for sh "<<shbest<<" out of "<<shs.size()
	<<" hp "<<hpbest
	<<" chF "<<hfbest
	<<" tp "<<tpbest
	<<" process "<<simhit_process->at(shbest)
	<<" particle "<<simhit_particle->at(shbest)			
	<<std::endl;
	if (rhType == HitType::Strip){
	  std::cout<<"    sh "<<simhit_x->at(shbest)<<", "<<simhit_y->at(shbest)<<", "<<simhit_z->at(shbest)
		   <<"  rh "<<str_x->at(rhIdx)<<", "<<str_y->at(rhIdx)<<", "<<str_z->at(rhIdx) <<std::endl;
	}
      }
      return ibest;
    };

    
    vector<Track> &simTracks_ = EE.simTracks_;
    vector<int> simTrackIdx_(sim_q->size(),-1);//keep track of original index in ntuple
    vector<int> seedSimIdx(see_q->size(),-1);
    for (unsigned int isim = 0; isim < sim_q->size(); ++isim) {

      //load sim production vertex data
      auto iVtx = sim_parentVtxIdx->at(isim);
      constexpr float largeValF = 9999.f;
      float sim_prodx = iVtx >= 0 ? simvtx_x->at(iVtx) : largeValF;
      float sim_prody = iVtx >= 0 ? simvtx_y->at(iVtx) : largeValF;
      float sim_prodz = iVtx >= 0 ? simvtx_z->at(iVtx) : largeValF;
      //if (fabs(sim_eta->at(isim))>0.8) continue;

      vector<int> const& trkIdxV = sim_trkIdx->at(isim);
      
      //if (trkIdx<0) continue;
      //FIXME: CHECK IF THE LOOP AND BEST SELECTION IS NEEDED.
      //Pick the first
      const int trkIdx = trkIdxV.empty() ? -1 : trkIdxV[0];
      
      int nlay = 0;
      if (trkIdx>=0) {	
	std::vector<int> hitlay(nTotalLayers, 0);
	auto const& hits = trk_hitIdx->at(trkIdx);
	auto const& hitTypes = trk_hitType->at(trkIdx);
	auto nHits = hits.size();
	for (auto ihit = 0U; ihit< nHits; ++ihit){
	  auto ihIdx = hits[ihit];
	  auto const ihType = HitType(hitTypes[ihit]);
	  
	  switch (ihType){
            case HitType::Pixel:{
              int ipix = ihIdx;
              if (ipix<0) continue;
              int cmsswlay = lnc.convertLayerNumber(pix_det->at(ipix),pix_lay->at(ipix),useMatched,-1,pix_z->at(ipix)>0);
              if (cmsswlay>=0 && cmsswlay<static_cast<int>(nTotalLayers)) hitlay[cmsswlay]++;	    
              break;
            }
            case HitType::Strip:{
              int istr = ihIdx;
              if (istr<0) continue;
              int cmsswlay = lnc.convertLayerNumber(str_det->at(istr),str_lay->at(istr),useMatched,str_isStereo->at(istr),str_z->at(istr)>0);
              if (cmsswlay>=0 && cmsswlay<static_cast<int>(nTotalLayers)) hitlay[cmsswlay]++;	    
              break;
            }
            case HitType::Glued:{
              if (useMatched) {
		int iglu = ihIdx;
		if (iglu<0) continue;
		int cmsswlay = lnc.convertLayerNumber(glu_det->at(iglu),glu_lay->at(iglu),useMatched,-1,glu_z->at(iglu)>0);
		if (cmsswlay>=0 && cmsswlay<static_cast<int>(nTotalLayers)) hitlay[cmsswlay]++;
              }	    
              break;
            }
            case HitType::Invalid: break;//FIXME. Skip, really?
            default: throw std::logic_error("Track type can not be handled");
	  }//hit type
	}//hits on track
	for (unsigned int i=0;i<nTotalLayers;i++) if (hitlay[i]>0) nlay++;
      }//count nlay layers on matching reco track

      //cout << Form("track q=%2i p=(%6.3f, %6.3f, %6.3f) x=(%6.3f, %6.3f, %6.3f) nlay=%i",sim_q->at(isim),sim_px->at(isim),sim_py->at(isim),sim_pz->at(isim),sim_prodx,sim_prody,sim_prodz,nlay) << endl;

      
      SVector3 pos(sim_prodx,sim_prody,sim_prodz);
      SVector3 mom(sim_px->at(isim),sim_py->at(isim),sim_pz->at(isim));
      SMatrixSym66 err;
      err.At(0,0) = sim_prodx*sim_prodx;
      err.At(1,1) = sim_prody*sim_prody;
      err.At(2,2) = sim_prodz*sim_prodz;
      err.At(3,3) = sim_px->at(isim)*sim_px->at(isim);
      err.At(4,4) = sim_py->at(isim)*sim_py->at(isim);
      err.At(5,5) = sim_pz->at(isim)*sim_pz->at(isim);
      TrackState state(sim_q->at(isim), pos, mom, err);
      state.convertFromCartesianToCCS();
      //create track: store number of reco hits in place of track chi2; fill hits later
      //              set label to be its own index in the output file
      Track track(state, float(nlay), simTracks_.size(), 0, nullptr);
      if (sim_bunchCrossing->at(isim) == 0){//in time
	if (sim_event->at(isim) == 0) track.setProdType(Track::ProdType::Signal);
	else track.setProdType(Track::ProdType::InTimePU);
      } else {
	track.setProdType(Track::ProdType::OutOfTimePU);
      }
      if (trkIdx>=0) {
	int seedIdx = trk_seedIdx->at(trkIdx);
	auto const& shTypes = see_hitType->at(seedIdx);
	if (std::count(shTypes.begin(), shTypes.end(), int(HitType::Pixel)) > 0) {
	  seedSimIdx[seedIdx] = simTracks_.size();
	}
      }
      if (cleanSimTracks){
	if (sim_nValid->at(isim) < cleanSimTrack_minSimHits) continue;
	if (cleanSimTrack_minRecHits > 0){
	  int nRecToSimHit = 0;
	  for (unsigned int ipix = 0; ipix < pix_lay->size() && nRecToSimHit < cleanSimTrack_minRecHits; ++ipix) {
	    int ilay = -1;
	    ilay = lnc.convertLayerNumber(pix_det->at(ipix),pix_lay->at(ipix),useMatched,-1,pix_z->at(ipix)>0);
	    if (ilay<0) continue;
	    int simTkIdxNt = bestTkIdx(pix_simHitIdx->at(ipix), pix_chargeFraction->at(ipix), ipix, HitType::Pixel);
	    if (simTkIdxNt >= 0) nRecToSimHit++;
	  }
	  if (useMatched) {
	    for (unsigned int iglu = 0; iglu < glu_lay->size() && nRecToSimHit < cleanSimTrack_minRecHits; ++iglu) {
	      if (glu_isBarrel->at(iglu)==0) continue;
	      int igluMono = glu_monoIdx->at(iglu);
	      int simTkIdxNt = bestTkIdx(str_simHitIdx->at(igluMono), str_chargeFraction->at(igluMono),       igluMono, HitType::Strip);
	      if (simTkIdxNt >= 0) nRecToSimHit++;
	    }
	  }
	  for (unsigned int istr = 0; istr < str_lay->size() && nRecToSimHit < cleanSimTrack_minRecHits; ++istr) {
	    int ilay = -1;
	    ilay = lnc.convertLayerNumber(str_det->at(istr),str_lay->at(istr),useMatched,str_isStereo->at(istr),str_z->at(istr)>0);
	    if (useMatched && str_isBarrel->at(istr)==1 && str_isStereo->at(istr)) continue;
	    if (ilay==-1) continue;
	    int simTkIdxNt = bestTkIdx(str_simHitIdx->at(istr), str_chargeFraction->at(istr), istr, HitType::Strip);
	    if (simTkIdxNt >= 0) nRecToSimHit++;
	  }
	  if (nRecToSimHit < cleanSimTrack_minRecHits) continue;
	}//count rec-to-sim hits
      }//cleanSimTracks

      simTrackIdx_[isim] = simTracks_.size();
      simTracks_.push_back(track);  
      
    }

    if (simTracks_.empty() and not writeAllEvents) continue;
    
    vector<Track> &seedTracks_ = EE.seedTracks_;
    vector<vector<int> > pixHitSeedIdx(pix_lay->size());
    for (unsigned int is = 0; is<see_q->size(); ++is) {
      auto isAlgo = TrackAlgorithm(see_algo->at(is));
      if (isAlgo != TrackAlgorithm::initialStep 
          && isAlgo != TrackAlgorithm::hltIter0 ) continue;//select seed in acceptance
      //if (see_pt->at(is)<0.5 || fabs(see_eta->at(is))>0.8) continue;//select seed in acceptance
      SVector3 pos = SVector3(see_stateTrajGlbX->at(is),see_stateTrajGlbY->at(is),see_stateTrajGlbZ->at(is));
      SVector3 mom = SVector3(see_stateTrajGlbPx->at(is),see_stateTrajGlbPy->at(is),see_stateTrajGlbPz->at(is));
      SMatrixSym66 err;
      err.At(0,0) = see_stateCcov00->at(is);
      err.At(0,1) = see_stateCcov01->at(is);
      err.At(0,2) = see_stateCcov02->at(is);
      err.At(0,3) = see_stateCcov03->at(is);
      err.At(0,4) = see_stateCcov04->at(is);
      err.At(0,5) = see_stateCcov05->at(is);
      err.At(1,1) = see_stateCcov11->at(is);
      err.At(1,2) = see_stateCcov12->at(is);
      err.At(1,3) = see_stateCcov13->at(is);
      err.At(1,4) = see_stateCcov14->at(is);
      err.At(1,5) = see_stateCcov15->at(is);
      err.At(2,2) = see_stateCcov22->at(is);
      err.At(2,3) = see_stateCcov23->at(is);
      err.At(2,4) = see_stateCcov24->at(is);
      err.At(2,5) = see_stateCcov25->at(is);
      err.At(3,3) = see_stateCcov33->at(is);
      err.At(3,4) = see_stateCcov34->at(is);
      err.At(3,5) = see_stateCcov35->at(is);
      err.At(4,4) = see_stateCcov44->at(is);
      err.At(4,5) = see_stateCcov45->at(is);
      err.At(5,5) = see_stateCcov55->at(is);
      TrackState state(see_q->at(is), pos, mom, err);
      state.convertFromCartesianToCCS();
      Track track(state, 0, seedSimIdx[is], 0, nullptr);
      auto const& shTypes = see_hitType->at(is);
      auto const& shIdxs = see_hitIdx->at(is);
      if (! ( (isAlgo == TrackAlgorithm::initialStep || isAlgo == TrackAlgorithm::hltIter0)
	     && std::count(shTypes.begin(), shTypes.end(), int(HitType::Pixel))>=3)) continue;//check algo and nhits
      for (unsigned int ip=0; ip<shTypes.size(); ip++) {
	unsigned int ipix = shIdxs[ip];
	//cout << "ipix=" << ipix << " seed=" << seedTracks_.size() << endl;
	pixHitSeedIdx[ipix].push_back(seedTracks_.size());
      }
      seedTracks_.push_back(track);
    }

    if (seedTracks_.empty() and not writeAllEvents) continue;

    vector<Track> &cmsswTracks_ = EE.cmsswTracks_;
    vector<vector<int> > pixHitRecIdx(pix_lay->size());
    vector<vector<int> > strHitRecIdx(str_lay->size());
    vector<vector<int> > gluHitRecIdx(glu_lay->size());
    for (unsigned int ir = 0; ir<trk_q->size(); ++ir) {
      //check the origin; redundant for initialStep ntuples
      if ((trk_algoMask->at(ir) & ( (1 << int(TrackAlgorithm::initialStep )) | (1 << int(TrackAlgorithm::hltIter0 )) )) == 0){
	if (verbosity > 1){
	  std::cout<<"track "<<ir<<" failed algo selection for "<< int(TrackAlgorithm::initialStep) <<": mask "<<trk_algoMask->at(ir)
		   <<" origAlgo "<<trk_originalAlgo->at(ir)<<" algo "<<trk_algo->at(ir)
		   <<std::endl;
	}
	continue;
      }
      //fill the state in CCS upfront
      SMatrixSym66 err;
      /*	
	vx = -dxy*sin(phi) - pt*cos(phi)/p*pz/p*dz;
	vy =  dxy*cos(phi) - pt*sin(phi)/p*pz/p*dz;
	vz = dz*pt*pt/p/p;
	//partial: ignores cross-terms
	c(vx,vx) = c(dxy,dxy)*sin(phi)*sin(phi) + c(dz,dz)*pow(pt*cos(phi)/p*pz/p ,2);
	c(vx,vy) = -c(dxy,dxy)*cos(phi)*sin(phi) + c(dz,dz)*cos(phi)*sin(phi)*pow(pt/p*pz/p, 2);
	c(vy,vy) = c(dxy,dxy)*cos(phi)*cos(phi) + c(dz,dz)*pow(pt*sin(phi)/p*pz/p ,2);
	c(vx,vz) = -c(dz,dz)*pt*pt/p/p*pt/p*pz/p*cos(phi);
	c(vy,vz) = -c(dz,dz)*pt*pt/p/p*pt/p*pz/p*sin(phi);
	c(vz,vz) = c(dz,dz)*pow(pt*pt/p/p, 2);
      */
      float pt = trk_pt->at(ir);
      float pz = trk_pz->at(ir);
      float p2 = pt*pt + pz*pz;
      float phi = trk_phi->at(ir);
      float sP = sin(phi);
      float cP = cos(phi);
      float dxyErr2 = trk_dxyErr->at(ir); dxyErr2 *= dxyErr2;
      float dzErr2 = trk_dzErr->at(ir); dzErr2 *= dzErr2;
      float dzErrF2 = trk_dzErr->at(ir)*(pt*pz/p2); dzErr2 *= dzErr2;
      err.At(0,0) =  dxyErr2 *sP*sP + dzErrF2 *cP*cP;
      err.At(0,1) = -dxyErr2 *cP*sP + dzErrF2 *cP*sP;
      err.At(1,1) =  dxyErr2 *cP*cP + dzErrF2 *sP*sP;
      err.At(0,2) = -dzErrF2*cP*pt/pz;
      err.At(1,2) = -dzErrF2*sP*pt/pz;
      err.At(2,2) =  dzErr2*std::pow((pt*pt/p2), 2);
      err.At(3,3) = std::pow(trk_ptErr->at(ir)/pt/pt, 2);
      err.At(4,4) = std::pow(trk_phiErr->at(ir), 2);
      err.At(5,5) = std::pow(trk_lambdaErr->at(ir), 2);
      SVector3 pos = SVector3(trk_refpoint_x->at(ir), trk_refpoint_y->at(ir), trk_refpoint_z->at(ir));
      SVector3 mom = SVector3(1.f/pt, phi, M_PI_2 - trk_lambda->at(ir));
      TrackState state(trk_q->at(ir), pos, mom, err);
      Track track(state, trk_nChi2->at(ir), trk_seedIdx->at(ir), 0, nullptr);//hits are filled later
      auto const& hTypes = trk_hitType->at(ir);
      auto const& hIdxs =  trk_hitIdx->at(ir);
      for (unsigned int ip=0; ip<hTypes.size(); ip++) {
	unsigned int hidx = hIdxs[ip];
	switch( HitType(hTypes[ip]) ) {
	  case HitType::Pixel:{
	    //cout << "pix=" << hidx << " track=" << cmsswTracks_.size() << endl;
	    pixHitRecIdx[hidx].push_back(cmsswTracks_.size());
	    break;
	  }
	  case HitType::Strip:{
	    //cout << "pix=" << hidx << " track=" << cmsswTracks_.size() << endl;
	    strHitRecIdx[hidx].push_back(cmsswTracks_.size());
	    break;
	  }
	  case HitType::Glued:{
	    if (not useMatched ) throw std::logic_error("Tracks have glued hits, but matchedHit load is not configured");
	    //cout << "pix=" << hidx << " track=" << cmsswTracks_.size() << endl;
	    gluHitRecIdx[hidx].push_back(cmsswTracks_.size());
	    break;
	  }
	  case HitType::Invalid: break;//FIXME. Skip, really?
	  default: throw std::logic_error("Track hit type can not be handled");
	}//switch( HitType
      }
      cmsswTracks_.push_back(track);
    }



    
    vector<vector<Hit> > &layerHits_   = EE.layerHits_;
    vector<MCHitInfo>    &simHitsInfo_ = EE.simHitsInfo_;
    int totHits = 0;
    layerHits_.resize(nTotalLayers);
    for (unsigned int ipix = 0; ipix < pix_lay->size(); ++ipix) {
      int ilay = -1;
      ilay = lnc.convertLayerNumber(pix_det->at(ipix),pix_lay->at(ipix),useMatched,-1,pix_z->at(ipix)>0);
      if (ilay<0) continue;

      unsigned int imoduleid;
      {
        auto ii = module_shortId_hash[ilay].emplace(pix_detId->at(ipix), (unsigned int) module_shortId_hash[ilay].size());
        imoduleid = ii.first->second;
      }

      int simTkIdxNt = bestTkIdx(pix_simHitIdx->at(ipix), pix_chargeFraction->at(ipix), ipix, HitType::Pixel); 
      int simTkIdx = simTkIdxNt >= 0 ? simTrackIdx_[simTkIdxNt] : -1; //switch to index in simTracks_

      //cout << Form("pix lay=%i det=%i x=(%6.3f, %6.3f, %6.3f)",ilay+1,pix_det->at(ipix),pix_x->at(ipix),pix_y->at(ipix),pix_z->at(ipix)) << endl;
      SVector3 pos(pix_x->at(ipix),pix_y->at(ipix),pix_z->at(ipix));
      SMatrixSym33 err;
      err.At(0,0) = pix_xx->at(ipix);
      err.At(1,1) = pix_yy->at(ipix);
      err.At(2,2) = pix_zz->at(ipix);
      err.At(0,1) = pix_xy->at(ipix);
      err.At(0,2) = pix_zx->at(ipix);
      err.At(1,2) = pix_yz->at(ipix);
      if (simTkIdx>=0){
	simTracks_[simTkIdx].addHitIdx(layerHits_[ilay].size(), ilay, 0);
      }
      for (unsigned int is=0;is<pixHitSeedIdx[ipix].size();is++) {
	//cout << "xxx ipix=" << ipix << " seed=" << pixHitSeedIdx[ipix][is] << endl;
      	seedTracks_[pixHitSeedIdx[ipix][is]].addHitIdx(layerHits_[ilay].size(), ilay, 0);//per-hit chi2 is not known
      }
      for (unsigned int ir=0;ir<pixHitRecIdx[ipix].size();ir++) {
	//cout << "xxx ipix=" << ipix << " recTrack=" << pixHitRecIdx[ipix][ir] << endl;
      	cmsswTracks_[pixHitRecIdx[ipix][ir]].addHitIdx(layerHits_[ilay].size(), ilay, 0);//per-hit chi2 is not known
      }
      Hit hit(pos, err, totHits);
      hit.setupAsPixel(imoduleid, pix_csize_row->at(ipix), pix_csize_col->at(ipix));
      layerHits_[ilay].push_back(hit);
      MCHitInfo hitInfo(simTkIdx, ilay, layerHits_[ilay].size()-1, totHits);
      simHitsInfo_.push_back(hitInfo);
      totHits++;
    }

    if (useMatched) {
      for (unsigned int iglu = 0; iglu < glu_lay->size(); ++iglu) {
	if (glu_isBarrel->at(iglu)==0) continue;
	int igluMono = glu_monoIdx->at(iglu);
	int simTkIdxNt = bestTkIdx(str_simHitIdx->at(igluMono), str_chargeFraction->at(igluMono),	igluMono, HitType::Strip);
	int simTkIdx = simTkIdxNt >= 0 ? simTrackIdx_[simTkIdxNt] : -1; //switch to index in simTracks_

	int ilay = lnc.convertLayerNumber(glu_det->at(iglu),glu_lay->at(iglu),useMatched,-1,glu_z->at(iglu)>0);
	// cout << Form("glu lay=%i det=%i bar=%i x=(%6.3f, %6.3f, %6.3f)",ilay+1,glu_det->at(iglu),glu_isBarrel->at(iglu),glu_x->at(iglu),glu_y->at(iglu),glu_z->at(iglu)) << endl;
	SVector3 pos(glu_x->at(iglu),glu_y->at(iglu),glu_z->at(iglu));
	SMatrixSym33 err;
	err.At(0,0) = glu_xx->at(iglu);
	err.At(1,1) = glu_yy->at(iglu);
	err.At(2,2) = glu_zz->at(iglu);
	err.At(0,1) = glu_xy->at(iglu);
	err.At(0,2) = glu_zx->at(iglu);
	err.At(1,2) = glu_yz->at(iglu);	
	if (simTkIdx>=0){
	  simTracks_[simTkIdx].addHitIdx(layerHits_[ilay].size(), ilay, 0);
	}
	for (unsigned int ir=0;ir<gluHitRecIdx[iglu].size();ir++) {
	  //cout << "xxx iglu=" << iglu << " recTrack=" << gluHitRecIdx[iglu][ir] << endl;
	  cmsswTracks_[gluHitRecIdx[iglu][ir]].addHitIdx(layerHits_[ilay].size(), ilay, 0);//per-hit chi2 is not known
	}

        // QQQQ module-id-in-layer, adc and phi/theta spans are not done for matched hits.
        // Will we ever use / need this?
        assert (false && "Implement module-ids, cluster adc and spans for matched hits!");

	Hit hit(pos, err, totHits);
	layerHits_[ilay].push_back(hit);
	MCHitInfo hitInfo(simTkIdx, ilay, layerHits_[ilay].size()-1, totHits);
	simHitsInfo_.push_back(hitInfo);
	totHits++;
      }
    }

    vector<int> strIdx;
    strIdx.resize(str_lay->size());
    for (unsigned int istr = 0; istr < str_lay->size(); ++istr) {
      int ilay = -1;
      ilay = lnc.convertLayerNumber(str_det->at(istr),str_lay->at(istr),useMatched,str_isStereo->at(istr),str_z->at(istr)>0);
      if (useMatched && str_isBarrel->at(istr)==1 && str_isStereo->at(istr)) continue;
      if (ilay==-1) continue;

      unsigned int imoduleid;
      {
        auto ii = module_shortId_hash[ilay].emplace(str_detId->at(istr), (unsigned int) module_shortId_hash[ilay].size() );
        imoduleid = ii.first->second;
      }

      int simTkIdxNt = bestTkIdx(str_simHitIdx->at(istr), str_chargeFraction->at(istr), istr, HitType::Strip);
      int simTkIdx = simTkIdxNt >= 0 ? simTrackIdx_[simTkIdxNt] : -1; //switch to index in simTracks_

      bool passCCC = applyCCC ? (str_chargePerCM->at(istr) > cutValueCCC) : true;

      //if (str_onTrack->at(istr)==0) continue;//do not consider hits that are not on track!
      SVector3 pos(str_x->at(istr),str_y->at(istr),str_z->at(istr));
      SMatrixSym33 err;
      err.At(0,0) = str_xx->at(istr);
      err.At(1,1) = str_yy->at(istr);
      err.At(2,2) = str_zz->at(istr);
      err.At(0,1) = str_xy->at(istr);
      err.At(0,2) = str_zx->at(istr);
      err.At(1,2) = str_yz->at(istr);
      if (simTkIdx>=0){
        if(passCCC) simTracks_[simTkIdx].addHitIdx(layerHits_[ilay].size(), ilay, 0);
        else simTracks_[simTkIdx].addHitIdx( -9, ilay,0);
      }
      for (unsigned int ir=0;ir<strHitRecIdx[istr].size();ir++) {
	//cout << "xxx istr=" << istr << " recTrack=" << strHitRecIdx[istr][ir] << endl;
	if(passCCC) cmsswTracks_[strHitRecIdx[istr][ir]].addHitIdx(layerHits_[ilay].size(), ilay, 0);//per-hit chi2 is not known
	else cmsswTracks_[strHitRecIdx[istr][ir]].addHitIdx(-9,ilay,0); 
      }
      if(passCCC)
	{
          Hit hit(pos, err, totHits);
          hit.setupAsStrip(imoduleid, str_chargePerCM->at(istr), str_csize->at(istr));
          layerHits_[ilay].push_back(hit);
	  MCHitInfo hitInfo(simTkIdx, ilay, layerHits_[ilay].size()-1, totHits);
	  simHitsInfo_.push_back(hitInfo);
	  totHits++;
	}
    }
    // Loop over inactive detIDs and add dummy hits as placemarkers
    int numInactive = 497;
    float x_val[] = {-11.1437, -11.0704, -11.144, -11.0692, -11.1434, -11.0706, -5.72967, -5.71181, -5.72804, -5.71023, -5.72943, -5.71092, 0.0195421, -0.0177409, 0.0179494, -0.0185438, 0.0177897, -0.0167593, 5.76417, 5.67638, 5.76343, 5.67444, 5.7644, 5.67457, 11.1754, 11.0365, 11.176, 11.0387, 11.1762, 11.0366, 20.0668, 20.2125, 20.0671, 20.212, 20.0682, 20.2128, 20.2337, 27.0626, 27.064, 27.0622, -23.5754, -27.063, -20.2125, -5.62838, -5.64126, -5.62909, -5.64083, -5.62706, -5.64011, 11.0857, 10.9591, 11.0848, 26.838, 21.8219, 21.6986, 21.8247, 21.6962, -31.1809, 29.6338, 32.2734, 32.0516, 32.2747, 32.0481, 32.2739, 32.051, 31.1725, 30.9641, 31.1713, 30.9636, 31.1706, 30.9615, 29.0085, 28.8211, 29.0091, 28.8222, 29.0085, 28.8194, 25.8568, 25.6959, 25.8557, 25.6943, 25.857, 25.6968, 21.8229, 21.6977, 21.823, 21.6966, 21.8243, 21.6986, 5.93158, 21.6695, 21.8501, 21.67, -11.4746, -35.2654, -35.4932, -11.5463, 0.0178074, 11.4727, 11.4717, 11.4724, 35.3361, 35.335, 35.337, -35.3361, -8.55623, -8.55597, -8.55558, 8.5613, 8.56005, 8.56115, 38.5695, 38.5694, 38.5691, -38.5662, 27.4117, 35.3359, 35.3354, 35.3352, -2.8744, -2.87438, -8.56007, -8.5605, -14.0719, 42.5298, 42.5288, 42.529, 38.5678, 38.5677, 38.5699, 35.4882, 35.4894, 35.4876, 31.7446, 31.7486, 31.7478, 27.4143, 27.414, 27.4144, -5.91253, -5.9127, -11.7173, -11.7173, 11.7167, 11.7179, 11.7169, 43.3362, 34.0909, 29.7316, 29.7317, 29.7347, 8.68899, 8.68892, 8.68857, -37.9509, -37.9533, -37.9518, -34.0914, -34.0897, -34.0908, -2.90874, 19.7873, 19.788, 19.7887, 34.0939, 34.0926, 34.0932, 51.3061, 51.3058, 51.3086, 47.4772, 47.4774, 47.477, -34.2421, -34.2439, -34.2439, -24.8584, -24.8583, -24.8581, 2.87931, 2.88112, 2.87787, 24.8546, 24.8554, 24.8581, 29.7341, 29.7364, 29.7352, 34.24, 34.2398, 34.2386, 38.3154, 38.3148, 38.3167, 41.9098, 41.906, 41.9086, 47.4754, 47.4736, 47.4752, 49.3788, 49.3776, 49.3771, 50.6617, 50.6606, 50.6619, -41.2585, -2.87928, -2.87982, -2.87925, -8.60563, -8.6056, -14.2237, -50.6578, -50.6618, -50.6613, 51.3069, 51.3064, 51.3062, -28.6725, -28.6729, -24.8345, -24.8348, -14.338, -14.3371, 35.4581, 7.42116, -7.4253, -20.2763, -27.6924, -27.6984, -27.6915, -24.8286, -24.833, -28.6711, -28.6714, 35.4581, 35.4568, 25.9563, 25.9573, 35.4569, 35.4566, -31.7917, 6.99932, 24.8307, 14.3372, 0.00152803, 27.6954, 20.2758, 27.6962, -25.9558, -44.1905, -39.8686, 35.4599, 41.596, 41.2957, -54.8982, -56.214, -34.5951, -34.3605, -34.5969, -34.3619, -34.596, -34.3605, 19.742, 19.6074, 31.7205, 31.9236, 31.2137, 31.4173, 31.7232, 31.9254, 31.214, 31.4167, 31.7214, 31.9228, 31.2134, 31.4176, 45.3948, 41.3913, 40.7652, 41.1417, -69.7359, -54.2529, 66.4609, 67.0093, 66.4629, 67.0094, 66.4624, 67.0074, -31.0145, -31.0139, -31.0128, -37.6902, -38.0134, -37.688, -38.0142, -37.6889, -38.0126, -48.0652, -48.0656, -48.0649, -22.3924, -22.2031, -22.3925, -22.2018, -22.3906, -22.2025, 15.9295, 16.0481, 15.9284, 16.0484, 15.9279, 16.0468, 6.54018, 6.58916, 6.54207, 6.58934, 6.53923, 6.59019, -2.4289, -2.44734, -2.42713, -2.44529, -2.42896, -2.44629, -11.1966, -11.282, -11.1956, -11.2818, -11.1966, -11.2815, -20.6803, -20.832, -20.6793, -20.8339, -20.6809, -20.8332, -28.4443, -28.6618, -28.4428, -28.6604, -28.4436, -28.6635, -38.0288, -38.3073, -38.027, -38.3086, -38.0273, -38.3088, -44.4479, -44.7887, -44.4477, -44.7899, -44.4489, -44.791, -60.092, -50.7717, -51.1216, -50.771, -51.1214, -50.7724, -51.1225, 74.6744, 74.6754, 74.6748, -43.7157, -98.3338, 74.1814, -78.8052, -78.8037, -78.8041, -55.0393, -55.0381, -55.3674, -55.0392, 10.6238, 53.9803, -10.6236, -10.6232, -14.0301, -14.0482, -13.8291, -20.3143, -15.5547, -21.5021, -27.178, -97.8437, -98.7579, 33.9081, 33.9143, 14.0649, 37.911, 22.3328, 22.3467, 29.1128, 29.1228, -5.28616, -22.3461, 9.40876, 84.4641, 88.2672, 91.5238, 94.2145, 96.3277, -65.8071, -65.8179, -62.5264, -62.5367, -57.7016, -57.7158, 43.9798, 43.9587, 51.4871, -47.3737, -0.000144759, -53.9809, -71.6644, 56.0852, 81.7523, 7.77351, 27.1749, -27.1772, -21.4996, -15.5514, -22.6651, 22.6671, 57.7272, 57.7147, 62.5452, 62.5363, 65.8243, 65.8209, 84.4663, 88.2685, 91.5255, 94.219, 96.3251, -4.77591, 5.28351, -67.4809, 27.178, 44.1926, 61.0504, 80.144, -7.77224, -14.0457, -56.0892, -51.4543, -37.9119, -33.9061, -33.9159, 33.9207, -33.9123, 67.4825, 67.482, 67.4757, 67.4824, 97.844, 98.7595, 99.062, 98.7572, 97.8436, 20.3118, 13.8259, 27.1747, 21.4976, 15.5492, 39.6082, 31.3345, 22.6666, -57.7284, -57.7167, -62.548, -62.5356, -65.8249, -65.8176, 43.9607, -98.7578, -97.8404, -62.5391, -42.5548, -39.8664, -53.9792, -51.9037, -49.1696, -78.6792, -75.6478, -71.6667, -9.40948, -45.816, -49.1661, -51.9027, -53.9829, -71.6622, -75.6488, -78.6794, -80.1505};
    float y_val[] = {-21.2725, -21.0531, -21.2715, -21.0532, -21.2703, -21.0551, -23.3197, -23.0904, -23.3186, -23.0924, -23.3197, -23.0892, -24.0148, -23.785, -24.0142, -23.7868, -24.0147, -23.786, -23.3112, -23.0993, -23.3101, -23.1013, -23.3103, -23.1, -21.2556, -21.071, -21.2556, -21.0715, -21.2543, -21.0696, 18.0446, 18.222, 18.0435, 18.223, 18.0454, 18.2203, -18.1964, 2.86258, 2.86417, 2.86248, 13.5902, -2.86266, -18.2219, -26.3927, -26.6237, -26.3937, -26.6238, -26.3919, -26.6211, -24.8547, -24.6591, -24.8531, -2.83899, 23.9658, 23.7731, 23.9659, 23.7736, 8.85309, 19.3395, 3.00909, 2.95217, 3.00765, 2.95159, 3.00812, 2.95341, 8.88737, 8.79049, 8.88832, 8.79181, 8.88563, 8.79022, 14.4647, 14.3307, 14.4631, 14.3303, 14.4622, 14.3303, 19.5499, 19.3828, 19.5465, 19.3832, 19.5487, 19.3848, 23.9659, 23.7747, 23.9669, 23.7739, 23.9649, 23.7735, 31.6348, -23.7991, -23.9413, -23.7983, 33.4742, 2.93869, 2.92307, -33.6903, -35.6133, -33.4744, -33.4748, -33.4758, 19.298, 19.298, 19.2988, -19.2971, -39.3416, -39.3428, -39.3438, -39.3427, -39.3416, -39.3412, 19.984, 19.9825, 19.9826, 19.9861, -33.6955, 19.2964, 19.2988, 19.2973, 40.1602, 40.1592, 39.3412, 39.3411, 37.7226, 8.83617, 8.83643, 8.83588, 19.9835, 19.9816, 19.9822, 25.0494, 25.0469, 25.05, 29.6463, 29.6466, 29.6471, 33.6936, 33.6939, 33.6955, 43.0329, 43.0333, 41.8269, 41.8269, -41.8273, -41.8277, -41.8276, -2.96715, 34.093, 37.9527, 37.9548, 37.9527, 47.4232, 47.423, 47.4225, 29.7307, 29.7328, 29.7311, -34.0941, -34.0924, -34.0931, -48.1236, -43.9653, -43.964, -43.9634, -34.0892, -34.09, -34.0911, 2.8797, 2.87883, 2.87892, 19.6634, 19.6628, 19.6615, -38.3157, -38.3151, -38.3157, -44.9729, -44.9752, -44.9757, -51.3075, -51.3059, -51.3054, -44.9771, -44.9767, -44.9769, -41.9105, -41.91, -41.9115, -38.3173, -38.3187, -38.3181, -34.2422, -34.2425, -34.2447, -29.7386, -29.7382, -29.7381, -19.6664, -19.6695, -19.6675, -14.2278, -14.2285, -14.2286, -8.61095, -8.61121, -8.60904, 24.9411, 51.3064, 51.3066, 51.3064, 50.6615, 50.6612, 49.3792, 8.60907, 8.60985, 8.60812, -2.88336, -2.8816, -2.88325, 0.00231905, 0.00164029, -14.3352, -14.3381, -24.8325, -24.8274, 9.50117, 27.6955, 27.7004, 20.272, 7.42389, 7.42185, -7.42046, 14.3351, 14.34, 0.000224108, 0.00163443, 9.50071, 9.50062, -25.9562, -25.9569, -9.50064, -9.50028, 18.3523, -44.1928, 14.3385, 24.8284, 28.6724, 7.42112, -20.2752, -7.42073, 25.9586, -6.99647, -20.307, -9.49495, -43.1813, -42.8702, 30.342, 20.8582, -52.8401, -52.4799, -52.8422, -52.4781, -52.8397, -52.4779, -59.9953, -59.5842, 60.0639, 60.4472, 59.1058, 59.4874, 60.0642, 60.4442, 59.1038, 59.4862, 60.0658, 60.4455, 59.1058, 59.488, -54.7564, 58.3729, 57.4878, 58.0196, 6.51583, -44.9789, 43.212, 43.5659, 43.2136, 43.5676, 43.2121, 43.5672, 72.9564, 72.9558, 72.9564, 66.081, 66.6467, 66.0839, 66.6485, 66.0814, 66.6509, -63.8573, -63.8565, -63.8572, -73.3848, -72.7624, -73.3844, -72.7629, -73.3865, -72.7624, 86.622, 87.2637, 86.6201, 87.2652, 86.6218, 87.2628, 84.6206, 85.2719, 84.6219, 85.2727, 84.6223, 85.272, 88.0416, 88.6947, 88.0412, 88.6916, 88.0401, 88.6923, 84.1302, 84.779, 84.131, 84.7778, 84.1341, 84.7774, 85.6107, 86.2449, 85.611, 86.246, 85.613, 86.2456, 79.9676, 80.5795, 79.9667, 80.5815, 79.9661, 80.5807, 79.4416, 80.0292, 79.4426, 80.0268, 79.4397, 80.0282, 72.3059, 72.8587, 72.3055, 72.8598, 72.3058, 72.8603, -77.1276, -79.7898, -80.3404, -79.7891, -80.3405, -79.7909, -80.3396, -64.1194, -64.1206, -64.1187, 83.8649, 4.24153, -63.6943, 71.9736, 71.9715, 71.9738, -94.4026, -94.4006, -94.9632, -94.3997, 25.6543, 15.5542, -25.6518, -25.653, -33.9201, -33.9122, -42.5559, -39.8635, -53.98, -51.9003, -49.1658, -15.4991, -7.7704, 14.0644, 14.048, -33.9071, -91.5204, 29.1327, 29.1218, 22.3626, 22.3466, 67.4823, 29.1242, -55.374, 51.7611, 44.9756, 37.9103, 30.6095, 23.1255, 15.826, 15.8017, 25.9271, 25.9037, 35.3896, 35.3682, 51.4525, 51.4701, 43.9398, 66.7683, -44.7417, 15.5519, 39.609, 3.15081, -4.59066, 98.7567, 49.1702, 49.1728, 51.9029, 53.9787, 78.6719, -78.6789, -35.3427, -35.3657, -25.8792, -25.9032, -15.7747, -15.7999, -51.7614, -44.9729, -37.9114, -30.6127, -23.1268, 36.3965, 67.4833, -5.33807, -49.1669, 7.00213, 54.5615, 58.227, 98.7594, 33.915, -3.14899, -43.9811, 91.521, 14.0648, 14.046, 14.0281, -14.0491, -5.28342, -5.3091, 5.33891, 5.30946, -15.4941, -7.77258, 0.000605578, 7.77053, 15.4951, 39.8682, 42.5523, 49.1669, 51.901, 53.9814, 71.6644, 75.6496, 78.6804, 35.3452, 35.3679, 25.8801, 25.9024, 15.7763, 15.8009, 51.4704, -7.77202, -15.498, 25.9004, -13.8265, -20.3172, -15.5524, -21.4974, -27.1754, -22.6707, -31.335, -39.608, 55.3778, -32.5089, 27.1748, 21.5013, 15.5516, 39.6043, 31.3311, 22.6677, -58.2338};
    float z_val[] = {-3.69155, -3.69222, -24.2702, -24.2706, -44.7355, -44.7341, -3.69251, -3.69422, -24.2704, -24.2704, -44.7359, -44.7362, -3.69459, -3.69321, -24.2704, -24.2736, -44.7336, -44.7339, -3.69262, -3.69517, -24.2717, -24.27, -44.7361, -44.7355, -3.69264, -3.69299, -24.2713, -24.2702, -44.7354, -44.7359, -16.0189, -16.0198, -39.4868, -39.4873, -60.5457, -60.5443, -39.488, 7.57632, 31.0889, 54.4767, 54.4751, 54.4752, 54.4757, 7.57665, 7.5784, 31.0904, 31.0922, 54.4755, 54.4762, 31.0896, 54.4757, 54.4756, 31.0917, -18.0208, -18.0227, -38.9252, -38.9272, -59.7644, -54.6371, 2.93783, 2.93463, 23.8884, 23.888, 44.7723, 44.7743, 2.93706, 2.93523, 23.8874, 23.889, 44.7732, 44.7715, 2.93476, 2.93504, 23.8882, 23.8874, 44.7725, 44.7744, 2.93652, 2.93546, 23.8887, 23.8857, 44.7744, 44.7745, 2.93779, 2.93808, 23.8893, 23.8874, 44.772, 44.7747, 2.93403, 2.93574, 23.8866, 44.7735, 38.0178, 38.0181, 38.0182, 60.5463, 14.8634, 14.8638, 38.018, 60.5449, -3.60883, -24.6281, -45.6184, -24.6263, -3.60923, -24.6262, -45.6189, -3.60959, -24.626, -45.6194, -15.2393, -37.9121, -60.5459, -60.548, -15.241, 17.4205, 38.4261, 59.3937, 17.4202, 38.4254, 38.4263, 59.3946, 17.417, 7.45891, 30.1437, 52.7903, 7.45753, 30.144, 52.7889, 7.4578, 30.1428, 52.791, 7.45681, 30.1452, 52.7897, 7.4601, 30.1441, 52.7912, 30.144, 52.7903, 30.1445, 52.7912, 7.46011, 30.1445, 52.7913, 7.45816, -60.5439, -17.8213, -39.1976, -60.5434, -17.8212, -39.195, -60.5457, -17.8208, -39.1982, -60.5448, -17.8203, -39.1972, -60.5444, -60.5444, -17.82, -39.1994, -60.5476, -17.8218, -39.1977, -60.5451, -7.60125, -30.4055, -53.1809, -7.60271, -30.4036, -53.1809, -7.60026, -30.4031, -53.1805, -7.59977, -30.4051, -53.1786, -7.60133, -30.405, -53.1775, -7.60057, -30.4046, -53.1781, -7.59986, -30.4031, -53.1772, -7.60301, -30.4042, -53.1798, -7.60118, -30.4049, -53.1776, -7.60024, -30.4051, -53.1789, -7.60163, -30.4036, -53.1807, -7.60062, -30.4037, -53.1785, -7.60081, -30.4041, -53.1787, 46.1123, 14.9868, 37.7807, 60.545, 14.9866, 37.778, 60.5465, 14.9878, 37.781, 60.5454, 14.9875, 37.7805, 60.5453, -77.9162, -77.6707, -77.913, -77.6688, -77.911, -77.6667, -84.0185, -93.2305, -93.2303, -93.2298, -93.2275, -93.4737, -93.2281, -90.8638, -90.6162, -90.8635, -90.6202, -96.7286, -96.9736, -96.7286, -96.9733, -96.7287, -96.9736, -94.3632, -89.7741, 77.6687, 77.6679, 77.6723, 93.2289, 93.2286, 93.2284, 96.7257, 89.7749, 89.7781, 109.924, -98.9812, -98.9814, 9.23208, 9.23148, 9.2303, 9.23261, 45.9754, 45.9784, 82.2417, 82.2428, 9.23105, 9.23093, -9.23047, -9.23002, -27.3133, -27.3155, -45.977, -45.9751, -63.3833, -63.3831, -82.2425, -82.2425, -98.9808, -98.9809, -82.2424, 63.3841, 80.4689, 98.9801, 9.22942, 80.4714, -9.23022, -27.3158, -45.0665, -63.3826, -80.4702, -98.9817, -9.23199, -45.0675, -80.4704, -9.23001, -27.3151, -45.0674, -63.383, -80.4694, -98.982, -27.3144, -63.3836, -98.9817, 9.23285, 27.3162, 45.9763, 63.3816, 82.2429, 98.9819, 9.22911, 27.3154, 45.066, 63.3821, 80.4685, 98.9804, 9.23176, 27.3162, 45.068, 63.384, 80.4709, 98.9815, 9.23047, 27.3157, 45.0644, 63.3817, 80.4677, 98.9805, 9.23359, 27.3152, 45.0648, 63.3812, 80.4691, 98.9804, 9.23187, 27.3155, 45.0669, 63.3842, 80.4709, 98.9798, 9.23124, 27.3141, 45.0638, 63.3827, 80.47, 98.9796, 9.2296, 27.315, 45.0656, 63.3859, 80.4697, 98.9815, 9.23241, 27.3147, 45.0672, 63.3817, 80.4695, 98.9811, -9.23154, -9.23086, -27.3144, -45.0636, -63.3823, -80.4712, -98.9799, -27.3155, -63.3836, -98.9809, 27.3159, 45.9749, 98.9817, -9.23034, -45.9775, -82.2405, -27.3155, -63.3818, -82.2422, -98.9798, -134.613, -151.335, -148.284, -148.613, -151.36, -151.03, -148.308, -148.809, -151.335, -150.836, -151.338, -155.716, -155.218, -193.36, -193.03, -210.862, -208.382, -204.4, -204.071, -203.698, -203.368, -200.644, -203.369, -203.168, -227.379, -226.884, -227.379, -226.882, -227.379, -226.069, -226.381, -226.807, -227.117, -226.069, -226.381, -219.648, -219.959, -218.91, -222.587, -220.135, -250.33, -250.265, -243.166, -242.6, -240.214, -271.836, -271.823, -271.327, -271.838, -271.754, 137.266, 133.574, 133.889, 134.307, 134.624, 133.577, 133.888, 134.879, 134.382, 134.878, 134.381, 134.878, 130.197, 127.144, 126.406, 151.335, 155.141, 158.099, 155.712, 155.213, 179.033, 171.67, 169.143, 190.882, 193.363, 193.032, 210.861, 210.532, 200.642, 200.958, 199.91, 200.223, 201.213, 200.716, 201.214, 200.715, 201.212, 226.808, 227.308, 229.835, 229.335, 229.836, 229.766, 229.267, 229.763, 226.068, 226.38, 226.803, 227.116, 226.068, 226.382, 219.22, 219.718, 220.215, 247.622, 247.307, 247.807, 250.339, 249.838, 250.337, 250.265, 249.767, 250.262, 243.17, 242.667, 271.836, 271.337, 271.837, 271.767, 271.265, 271.764, 262.21};
    
    float det_val[] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}; 
    float lay_val[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9};
    float isStereo_val[] = {1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float phiSpan_val[] = {0.29759, 0.254004, 0.297598, 0.254011, 0.297614, 0.253986, 0.297604, 0.254001, 0.297624, 0.253987, 0.297605, 0.254017, 0.29759, 0.254018, 0.297593, 0.253998, 0.297588, 0.254009, 0.297605, 0.254004, 0.29762, 0.253985, 0.297615, 0.253999, 0.29759, 0.254001, 0.297588, 0.253988, 0.297603, 0.254014, 0.265137, 0.222321, 0.265141, 0.22232, 0.265121, 0.222327, 0.222332, 0.222324, 0.222313, 0.222327, 0.222337, 6.06086, 0.222322, 0.265142, 0.222313, 0.265131, 0.222313, 0.265154, 0.222336, 0.222313, 0.265157, 0.222327, 0.265127, 0.221166, 0.188142, 0.22115, 0.188151, 0.221157, 0.202669, 0.221157, 0.188141, 0.221147, 0.188161, 0.221153, 0.188145, 0.221147, 0.188137, 0.221155, 0.18814, 0.221162, 0.188153, 0.221147, 0.18814, 0.221146, 0.188134, 0.221152, 0.18815, 0.221142, 0.188146, 0.22116, 0.188153, 0.221145, 0.188135, 0.221159, 0.18814, 0.221154, 0.188148, 0.221156, 0.18814, 0.188148, 0.188144, 0.221159, 0.188147, 0.20267, 6.08052, 6.11303, 0.170151, 0.170154, 0.202671, 0.202669, 0.202663, 0.15057, 0.150574, 0.150566, 0.150571, 0.150574, 0.15057, 0.150567, 0.150566, 0.15057, 0.150572, 0.13959, 0.139593, 0.139593, 0.139596, 0.139597, 0.150573, 0.150572, 0.150574, 0.150584, 0.150574, 0.150576, 0.150578, 0.15057, 0.139594, 0.139598, 0.139597, 0.139597, 0.139599, 0.139592, 0.139594, 0.139594, 0.139595, 0.139604, 0.139592, 0.139594, 0.139596, 0.139597, 0.139592, 0.139598, 0.139613, 0.139599, 0.139596, 0.139596, 0.139594, 0.139595, 0.139595, 0.125805, 0.12581, 0.125805, 0.125803, 0.125807, 0.125808, 0.125809, 0.125814, 0.125806, 0.125812, 0.125803, 0.125809, 0.125805, 0.125811, 0.125806, 0.125809, 0.12581, 0.125808, 0.125808, 0.125805, 0.118052, 0.118052, 0.118046, 0.118049, 0.118049, 0.11805, 0.118052, 0.11805, 0.118049, 0.118054, 0.11805, 0.118049, 0.118049, 0.118053, 0.118053, 0.11805, 0.11805, 0.118047, 0.118052, 0.11805, 0.118049, 0.118052, 0.11805, 0.118053, 0.118052, 0.118053, 0.118046, 0.118047, 0.118054, 0.11805, 0.11805, 0.11805, 0.118049, 0.11805, 0.118052, 0.118052, 0.118048, 0.11805, 0.118049, 0.12581, 0.118051, 0.118056, 0.118052, 0.118059, 0.118051, 0.11806, 0.118058, 0.118049, 0.118051, 0.11805, 0.118052, 0.118051, 5.97924, 6.01706, 0.303916, 0.26607, 0.30394, 0.266165, 0.266066, 0.303971, 0.303947, 0.304018, 0.303953, 0.266079, 0.304024, 0.304, 0.2661, 5.9792, 6.01705, 0.28861, 0.266081, 0.288621, 0.266069, 0.288623, 0.266079, 0.288626, 0.160762, 0.266125, 0.266151, 0.266125, 0.303967, 0.303952, 0.303956, 0.288621, 0.160782, 0.160798, 0.266105, 0.186006, 0.157156, 0.149155, 0.186014, 0.176639, 0.149148, 0.176632, 0.14915, 0.176639, 0.149152, 0.176633, 0.149149, 0.164294, 0.136897, 0.166946, 0.139099, 0.16429, 0.136901, 0.16695, 0.139102, 0.164289, 0.136901, 0.166947, 0.139098, 0.156929, 0.130793, 0.1328, 0.15693, 0.159354, 0.132802, 0.118091, 0.117129, 0.118087, 0.117127, 0.118089, 0.11713, 0.118089, 0.11809, 0.11809, 0.123047, 0.122005, 0.123044, 0.122001, 0.123048, 0.121999, 0.117131, 0.117131, 0.117131, 0.122005, 0.123046, 0.122006, 0.123046, 0.122003, 0.123046, 0.106311, 0.10553, 0.106313, 0.105529, 0.106312, 0.105532, 0.110314, 0.109473, 0.110312, 0.109472, 0.110312, 0.109473, 0.10631, 0.105529, 0.106311, 0.105532, 0.106312, 0.105531, 0.110315, 0.109473, 0.110315, 0.109474, 0.11031, 0.109475, 0.106313, 0.105533, 0.106313, 0.105531, 0.10631, 0.105531, 0.110311, 0.109475, 0.110312, 0.109473, 0.110313, 0.109472, 0.106311, 0.105533, 0.106311, 0.105535, 0.106314, 0.105533, 0.110312, 0.109476, 0.110312, 0.109474, 0.110311, 0.109472, 0.0957792, 0.0990157, 0.0983386, 0.0990169, 0.0983386, 0.0990143, 0.0983386, 0.0951459, 0.0951443, 0.095146, 0.0990145, 6.18804, 0.0957787, 0.0877542, 0.0877562, 0.0877547, 0.0857089, 0.0857105, 0.0852032, 0.0857108, 0.266134, 0.114874, 0.295504, 0.266148, 0.288731, 0.266088, 0.160757, 0.160781, 0.114877, 0.114876, 0.114879, 0.0804391, 0.0804393, 0.288724, 0.26607, 0.288721, 0.0804403, 0.288717, 0.266078, 0.288707, 0.26607, 0.180699, 0.266072, 0.114891, 0.0804381, 0.0804359, 0.0804349, 0.0804412, 0.0804353, 0.180729, 0.15963, 0.180697, 0.159617, 0.180682, 0.15962, 0.180709, 0.159623, 0.180706, 0.114861, 0.160771, 0.114869, 0.114837, 0.114873, 0.114839, 0.0804383, 0.11486, 0.11486, 0.114863, 0.114874, 0.114865, 0.114842, 0.180712, 0.159621, 0.180704, 0.159618, 0.180708, 0.159618, 0.080436, 0.0804361, 0.0804342, 0.0804323, 0.0804396, 0.288705, 0.180703, 6.10249, 0.114874, 0.16077, 0.114844, 0.0804386, 0.0804358, 0.266084, 6.16832, 0.180702, 0.0804392, 0.288731, 0.266075, 0.288745, 0.266089, 0.180706, 0.159617, 0.180715, 0.159617, 0.0804398, 0.0804351, 0.0804373, 0.0804409, 0.0804404, 0.160764, 0.160764, 0.114869, 0.114866, 0.114873, 0.114838, 0.114837, 0.11484, 0.180693, 0.159612, 0.180691, 0.159629, 0.180693, 0.159635, 0.15962, 0.0804372, 0.0804396, 0.159624, 0.160758, 0.160772, 0.114871, 0.114866, 0.114861, 0.114845, 0.114839, 0.114835, 0.114876, 0.114864, 0.114872, 0.114867, 0.114863, 0.114844, 0.114844, 0.114841, 0.0804579};

    double spacing = 0.02;

    for(int i =0; i < numInactive; i++){
      int ilay = -1;
      bool posEndcap = (z_val[i] > 0);
      ilay = lnc.convertLayerNumber(det_val[i], lay_val[i], useMatched, isStereo_val[i], posEndcap);  

      int num_hits = (phiSpan_val[i]/spacing) + 1;
      double phi_center = std::atan(y_val[i] / x_val[i]);
      if((x_val[i] < 0) && (y_val[i] < 0) )  phi_center = (phi_center-3.14159);
      else if ((x_val[i] < 0) && (y_val[i] > 0) )  phi_center = (phi_center+3.14159);

      double r_val = std::sqrt(x_val[i]*x_val[i]+y_val[i]*y_val[i]);
      for(int k=0; k < num_hits; k++){
	double phi_val = phi_center - 0.5 * phiSpan_val[i] + 0.5 * spacing + k * spacing;  
	if(phi_val > phi_center + 0.5*phiSpan_val[i]){
	  phi_val = phi_center + 0.5*phiSpan_val[i] - 0.4*spacing;
	}
	double x_ = r_val * std::cos(phi_val);
	double y_ = r_val * std::sin(phi_val);

	SVector3 pos(x_, y_, z_val[i]);
	SMatrixSym33 err;
	float dummyErr = 40.0;
	float dummyIdx = -7;
	err.At(0,0) = dummyErr;
	err.At(1,1) = dummyErr;
	err.At(2,2) = dummyErr;
	err.At(0,1) = dummyErr;
	err.At(0,2) = dummyErr;
	err.At(1,2) = dummyErr;
	Hit hit(pos, err, dummyIdx);
	layerHits_[ilay].push_back(hit);
	//	std::cout << "Dummy hit " << i << " added on layer " << ilay << ", " << x_ << ", " << y_ << ", " << z_val[i] << "(" << x_val[i] << ", " << y_val[i] << ", " << z_val[i] << ")"<< std::endl;
      }
    }
    //End of dummy hit treatment 

    // Seed % hit statistics
    nstot += seedTracks_.size();
    for (unsigned int il = 0; il<layerHits_.size(); ++il) {
      int nh = layerHits_[il].size();
      nhitstot[il]+=nh;
    }

    if (verbosity>0)
    {
      int nt = simTracks_.size();

      int nl = layerHits_.size();
    
      int nm = simHitsInfo_.size();

      int ns = seedTracks_.size();

      int nr = cmsswTracks_.size();

      printf("number of simTracks %i\n",nt);
      printf("number of layerHits %i\n",nl);
      printf("number of simHitsInfo %i\n",nm);
      printf("number of seedTracks %i\n",ns);
      printf("number of recTracks %i\n",nr);

      if (verbosity>1) {
	printf("\n");
	for (int il = 0; il<nl; ++il) {
	  int nh = layerHits_[il].size();
	  for (int ih=0; ih<nh; ++ih ) {
	    printf("lay=%i idx=%i mcid=%i x=(%6.3f, %6.3f, %6.3f) r=%6.3f\n",il+1,ih,layerHits_[il][ih].mcHitID(),layerHits_[il][ih].x(),layerHits_[il][ih].y(),layerHits_[il][ih].z(),sqrt(pow(layerHits_[il][ih].x(),2)+pow(layerHits_[il][ih].y(),2)));
	  }
	}

	for (int i=0;i<nt;++i) {
	  float spt = sqrt(pow(simTracks_[i].px(),2)+pow(simTracks_[i].py(),2));
	  printf("sim track id=%i q=%2i p=(%6.3f, %6.3f, %6.3f) x=(%6.3f, %6.3f, %6.3f) pT=%7.4f nTotal=%i nFound=%i \n",i,
		 simTracks_[i].charge(),simTracks_[i].px(),simTracks_[i].py(),simTracks_[i].pz(),simTracks_[i].x(),simTracks_[i].y(),simTracks_[i].z(),spt,
		 simTracks_[i].nTotalHits(),simTracks_[i].nFoundHits());
	  int nh = simTracks_[i].nTotalHits();
	  for (int ih=0;ih<nh;++ih){
	    int hidx = simTracks_[i].getHitIdx(ih);
	    int hlay = simTracks_[i].getHitLyr(ih);
	    float hx = layerHits_[hlay][hidx].x();
	    float hy = layerHits_[hlay][hidx].y();
	    float hz = layerHits_[hlay][hidx].z();
	    printf("track #%4i hit #%2i idx=%4i lay=%2i x=(% 8.3f, % 8.3f, % 8.3f) r=%8.3f\n",
		   i,ih,hidx,hlay,hx,hy,hz, sqrt(hx*hx+hy*hy));
	  }
	}
	
	for (int i=0;i<ns;++i) {
	  printf("seed id=%i label=%i q=%2i pT=%6.3f p=(%6.3f, %6.3f, %6.3f) x=(%6.3f, %6.3f, %6.3f)\n",i,
		 seedTracks_[i].label(),seedTracks_[i].charge(),seedTracks_[i].pT(),seedTracks_[i].px(),seedTracks_[i].py(),seedTracks_[i].pz(),seedTracks_[i].x(),seedTracks_[i].y(),seedTracks_[i].z());
	  int nh = seedTracks_[i].nTotalHits();
	  for (int ih=0;ih<nh;++ih) printf("seed #%i hit #%i idx=%i\n",i,ih,seedTracks_[i].getHitIdx(ih));
	}

	if (writeRecTracks){
	  for (int i=0;i<nr;++i) {
	    float spt = sqrt(pow(cmsswTracks_[i].px(),2)+pow(cmsswTracks_[i].py(),2));
	    printf("rec track id=%i label%i chi2=%6.3f q=%2i p=(%6.3f, %6.3f, %6.3f) x=(%6.3f, %6.3f, %6.3f) pT=%7.4f nTotal=%i nFound=%i \n",
		   i, cmsswTracks_[i].label(), cmsswTracks_[i].chi2(),
		   cmsswTracks_[i].charge(),cmsswTracks_[i].px(),cmsswTracks_[i].py(),cmsswTracks_[i].pz(),cmsswTracks_[i].x(),cmsswTracks_[i].y(),cmsswTracks_[i].z(),spt,
		   cmsswTracks_[i].nTotalHits(),cmsswTracks_[i].nFoundHits());
	    int nh = cmsswTracks_[i].nTotalHits();
	    for (int ih=0;ih<nh;++ih){
	      int hidx = cmsswTracks_[i].getHitIdx(ih);
	      int hlay = cmsswTracks_[i].getHitLyr(ih);
	      float hx = layerHits_[hlay][hidx].x();
	      float hy = layerHits_[hlay][hidx].y();
	      float hz = layerHits_[hlay][hidx].z();
	      printf("track #%4i hit #%2i idx=%4i lay=%2i x=(% 8.3f, % 8.3f, % 8.3f) r=%8.3f\n",
		     i,ih,hidx,hlay,hx,hy,hz, sqrt(hx*hx+hy*hy));
	    }
	  }
	}//if (writeRecTracks){
	
      }//verbosity>1
    }//verbosity>0
    EE.write_out(data_file);
      
    savedEvents++;
    printf("end of event %lli\n",savedEvents);
  }

  data_file.CloseWrite(savedEvents);
  printf("\nSaved %lli events\n\n",savedEvents);

  printf("Average number of seeds per event %f\n",float(nstot)/float(savedEvents));
  for (unsigned int il=0;il<nhitstot.size();++il)
    printf("Average number of hits in layer %3i = %7.2f\n", il, float(nhitstot[il])/float(savedEvents)); //Includes those that failed the cluster charge cut

  printf("Out of %i hits, %i failed the cut",numTotalStr,numFailCCC);

  //========================================================================

  printf("\n\n================================================================\n");
  printf("=== Max module id for %d layers\n", Config::nTotalLayers);
  printf("================================================================\n");
  for (int ii = 0; ii < Config::nTotalLayers; ++ii)
  {
    printf("Layer%2d : %d\n", ii, (int) module_shortId_hash[ii].size());
  }
}
