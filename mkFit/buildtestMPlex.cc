#include "buildtestMPlex.h"

#include "MkFitter.h"
#include "Matrix.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "Simulation.h"
#include "Config.h"
#include "BinInfoUtils.h"

#include <omp.h>

#if defined(USE_VTUNE_PAUSE)
#include "ittnotify.h"
#endif

//#define PRINTOUTS_FOR_PLOTS

inline bool sortByHitsChi2(const std::pair<Track, TrackState>& cand1, const std::pair<Track, TrackState>& cand2)
{
  if (cand1.first.nFoundHits()==cand2.first.nFoundHits()) return cand1.first.chi2()<cand2.first.chi2();
  return cand1.first.nFoundHits()>cand2.first.nFoundHits();
}

inline bool sortCandByHitsChi2(const Track& cand1, const Track& cand2)
{
  if (cand1.nFoundHits()==cand2.nFoundHits()) return cand1.chi2()<cand2.chi2();
  return cand1.nFoundHits()>cand2.nFoundHits();
}

inline bool sortByPhi(const Hit& hit1, const Hit& hit2)
{
  return std::atan2(hit1.y(),hit1.x())<std::atan2(hit2.y(),hit2.x());
}

inline bool sortByEta(const Hit& hit1, const Hit& hit2){
  return hit1.eta()<hit2.eta();
}
inline bool sortTracksByEta(const Track& track1, const Track& track2){
  return track1.momEta()<track2.momEta();
}
inline bool sortTracksByPhi(const Track& track1, const Track& track2){
  return track1.momPhi()<track2.momPhi();
}
struct sortTracksByPhiStruct {
  std::vector<std::vector<Track> >* track_candidates;
  sortTracksByPhiStruct(std::vector<std::vector<Track> >* track_candidates_) { track_candidates=track_candidates_; }
  bool operator() (const std::pair<int,int>& track1, const std::pair<int,int>& track2) {
    return (*track_candidates)[track1.first][track1.second].posPhi()<(*track_candidates)[track2.first][track2.second].posPhi();
  }
};

// within a layer with a "reasonable" geometry, ordering by Z is the same as eta
static bool sortByZ(const Hit& hit1, const Hit& hit2){
  return hit1.z()<hit2.z();
}


//==============================================================================
// Common functions
//==============================================================================

class MkBuilder
{
private:
  Event         *m_event;
  EventOfHits    m_event_of_hits;
  MkFitter      *m_mkfp_arr[NUM_THREADS];

  std::vector<Track> m_recseeds;

  int m_cnt=0, m_cnt1=0, m_cnt2=0, m_cnt_8=0, m_cnt1_8=0, m_cnt2_8=0, m_cnt_nomc=0;

public:
  MkBuilder();
  ~MkBuilder();

  // --------

  void begin_event(Event& ev, const char* build_type);

  void fit_seeds();

  void end_event();

  // --------

  void quality_reset();
  void quality_process(Track& tkcand);
  void quality_print();

  // --------

  void FindTracksBestHit(EventOfCandidates& event_of_cands);

  void FindTracks(EventOfCombCandidates& event_of_comb_cands);
};

//------------------------------------------------------------------------------
// MkBuilder - constructor and destructor
//------------------------------------------------------------------------------

MkBuilder::MkBuilder() :
  m_event(0),
  m_event_of_hits(Config::nLayers)
{
  for (int i = 0; i < NUM_THREADS; ++i)
  {
    m_mkfp_arr[i] = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(Config::nlayers_per_seed);
  }
}

MkBuilder::~MkBuilder()
{
   for (int i = 0; i < NUM_THREADS; ++i)
   {
     _mm_free(m_mkfp_arr[i]);
   }
}

//------------------------------------------------------------------------------
// MkBuilder - common functions
//------------------------------------------------------------------------------

void MkBuilder::begin_event(Event& ev, const char* build_type)
{
  m_event = &ev;

  std::vector<Track>& simtracks = ev.simTracks_;

  std::cout << "Building tracks with '" << build_type << "', total simtracks=" << simtracks.size() << std::endl;
#ifdef DEBUG
  //unit test for eta partitioning
  for (int i=0;i<60; ++i)
    {
      float eta = -1.5 + 0.05*i;
      int b1, b2;
      int cnt = getBothEtaBins(eta, b1, b2);
      std::cout << "eta=" << eta << " bin=" << getEtaBin(eta) << " hb1=" << b1 << " hb2=" << b2 << std::endl;
    }
  //dump sim tracks
  for (int itrack=0;itrack<simtracks.size();++itrack) {
    Track track = simtracks[itrack];
    std::cout << "MX - simtrack with nHits=" << track.nFoundHits() << " chi2=" << track.chi2()  << " pT=" << sqrt(track.momentum()[0]*track.momentum()[0]+track.momentum()[1]*track.momentum()[1]) <<" phi="<< track.momPhi() <<" eta=" << track.momEta() << std::endl;
  }
#endif

#ifdef PRINTOUTS_FOR_PLOTS
  for (int itrack=0;itrack<simtracks.size();++itrack) {
    Track track = simtracks[itrack];
    std::cout << "MX - simtrack with nHits=" << track.nFoundHits() << " chi2=" << track.chi2()  << " pT=" << sqrt(track.momentum()[0]*track.momentum()[0]+track.momentum()[1]*track.momentum()[1]) <<" phi="<< track.momPhi() <<" eta=" << track.momEta() << std::endl;
  }
#endif

  m_event_of_hits.Reset();

  for (int itrack=0; itrack < simtracks.size(); ++itrack)
  {
    if (simtracks[itrack].label() != itrack)
    {
      printf("Bad label for simtrack %d -- %d\n", itrack, simtracks[itrack].label());
    }

    //fill vector of hits in each layer (assuming there is one hit per layer in hits vector)
    for (int ilay = 0; ilay < simtracks[itrack].nTotalHits(); ++ilay)
    {
      m_event_of_hits.InsertHit(simtracks[itrack].hitsVector(ev.layerHits_)[ilay], ilay);
    }
  }

  m_event_of_hits.SortByPhiBuildPhiBins();

  m_recseeds.resize(simtracks.size());
}

void MkBuilder::fit_seeds()
{
  std::vector<Track>& simtracks = m_event->simTracks_;

  int theEnd = simtracks.size();

#pragma omp parallel for num_threads(NUM_THREADS)
  for (int itrack = 0; itrack < theEnd; itrack += NN)
  {
    int end = std::min(itrack + NN, theEnd);

    MkFitter *mkfp = m_mkfp_arr[omp_get_thread_num()];

    mkfp->SetNhits(3);//just to be sure (is this needed?)

    mkfp->InputTracksAndHits(simtracks, m_event->layerHits_, itrack, end);

    mkfp->FitTracks();

    mkfp->OutputFittedTracksAndHitIdx(m_recseeds, itrack, end);
  }

  //ok now, we should have all seeds fitted in recseeds
#ifdef DEBUG
  std::cout << "found total seeds=" << m_recseeds.size() << std::endl;
  for (int iseed = 0; iseed < m_recseeds.size(); ++iseed)
  {
    Track& seed = m_recseeds[iseed];
    std::cout << "MX - found seed with nHits=" << seed.nFoundHits() << " chi2=" << seed.chi2() << " posEta=" << seed.posEta() << " posPhi=" << seed.posPhi() << " posR=" << seed.radius() << " pT=" << seed.pT() << std::endl;
  }
#endif

#ifdef PRINTOUTS_FOR_PLOTS
  std::cout << "found total seeds=" << m_recseeds.size() << std::endl;
  for (int iseed = 0; iseed < m_recseeds.size(); ++iseed)
  {
    Track& seed = m_recseeds[iseed];
    std::cout << "MX - found seed with nHits=" << seed.nFoundHits() << " chi2=" << seed.chi2() << " posEta=" << seed.posEta() << " posPhi=" << seed.posPhi() << " posR=" << seed.radius() << " pT=" << seed.pt() << std::endl;
  }
#endif
}

void MkBuilder::end_event()
{

  m_event = 0;
}

//------------------------------------------------------------------------------

void MkBuilder::quality_reset()
{
  m_cnt = m_cnt1 = m_cnt2 = m_cnt_8 = m_cnt1_8 = m_cnt2_8 = m_cnt_nomc = 0;
}

void MkBuilder::quality_process(Track &tkcand)
{
  int mctrk = tkcand.label();
  if (mctrk < 0 || mctrk >= Config::nTracks)
  {
    ++m_cnt_nomc;
    // std::cout << "XX bad track idx " << mctrk << "\n";
    return;
  }
  float pt    = tkcand.pT();
  float ptmc  = m_event->simTracks_[mctrk].pT() ;
  float pr    = pt / ptmc;

  ++m_cnt;
  if (pr > 0.9 && pr < 1.1) ++m_cnt1;
  if (pr > 0.8 && pr < 1.2) ++m_cnt2;

  if (tkcand.nFoundHits() >= 8)
  {
    ++m_cnt_8;
    if (pr > 0.9 && pr < 1.1) ++m_cnt1_8;
    if (pr > 0.8 && pr < 1.2) ++m_cnt2_8;
  }

#ifdef DEBUG
  std::cout << "MXBH - found track with nFoundHits=" << tkcand.nFoundHits() << " chi2=" << tkcand.chi2() << " pT=" << pt <<" pTmc="<< ptmc << std::endl;
#endif

#ifdef PRINTOUTS_FOR_PLOTS
  std::cout << "MX - found track with nFoundHits=" << tkcand.nFoundHits() << " chi2=" << tkcand.chi2() << " pT=" << pt <<" pTmc="<< ptmc << std::endl;
#endif
}

void MkBuilder::quality_print()
{
  std::cout << "found tracks=" << m_cnt   << "  in pT 10%=" << m_cnt1   << "  in pT 20%=" << m_cnt2   << "     no_mc_assoc="<< m_cnt_nomc <<std::endl;
  std::cout << "  nH >= 8   =" << m_cnt_8 << "  in pT 10%=" << m_cnt1_8 << "  in pT 20%=" << m_cnt2_8 << std::endl;
}

//------------------------------------------------------------------------------
// FindTracksBestHit
//------------------------------------------------------------------------------

void MkBuilder::FindTracksBestHit(EventOfCandidates& event_of_cands)
{
  // partition recseeds into eta bins
  for (int iseed = 0; iseed < m_recseeds.size(); ++iseed)
  {
    if (m_recseeds[iseed].label() != iseed)
    {
      printf("Bad label for recseed %d -- %d\n", iseed, m_recseeds[iseed].label());
    }

    event_of_cands.InsertCandidate(m_recseeds[iseed]);
  }

  //dump seeds
#ifdef DEBUG
  for (int ebin = 0; ebin < Config::nEtaBin; ++ebin)
  {
    EtaBinOfCandidates &etabin_of_candidates = event_of_cands.m_etabins_of_candidates[ebin]; 
    for (int iseed = 0; iseed < etabin_of_candidates.m_fill_index; iseed++)
    {
      Track& seed = etabin_of_candidates.m_candidates[iseed];
      std::cout << "MX - found seed with nFoundHits=" << seed.nFoundHits() << " chi2=" << seed.chi2() 
                << " x=" << seed.position()[0] << " y=" << seed.position()[1] << " z=" << seed.position()[2] 
                << " px=" << seed.momentum()[0] << " py=" << seed.momentum()[1] << " pz=" << seed.momentum()[2] 
                << " pT=" << sqrt(seed.momentum()[0]*seed.momentum()[0]+seed.momentum()[1]*seed.momentum()[1]) 
                << std::endl;
    }
  }
#endif

  //parallel section over seeds; num_threads can of course be smaller
  int nseeds = m_recseeds.size();
  //#pragma omp parallel num_threads(Config::nEtaBin)
#pragma omp parallel num_threads(1)//fixme: set to one for debugging (to be revisited anyway - what if there are more threads than eta bins?)
   for (int ebin = 0; ebin < Config::nEtaBin; ++ebin)
   {
     // XXXX Could have nested paralellism, like NUM_THREADS/nEtaBins (but rounding sucks here).
     // XXXX So one should really have TBB, for this and for the above.
     // vectorized loop
     EtaBinOfCandidates &etabin_of_candidates = event_of_cands.m_etabins_of_candidates[ebin];

     for (int itrack = 0; itrack < etabin_of_candidates.m_fill_index; itrack += NN)
       {
	 int end = std::min(itrack + NN, etabin_of_candidates.m_fill_index);
	 
#ifdef DEBUG
         std::cout << std::endl;
	 std::cout << "processing track=" << itrack << " etabin=" << ebin << " findex=" << etabin_of_candidates.m_fill_index << " thn=" << omp_get_thread_num() << std::endl;
#endif

	 MkFitter *mkfp = m_mkfp_arr[omp_get_thread_num()];

	 mkfp->SetNhits(3);//just to be sure (is this needed?)

	 mkfp->InputTracksAndHitIdx(etabin_of_candidates.m_candidates, itrack, end);

	 //ok now we start looping over layers
	 //loop over layers, starting from after the seed
	 //consider inverting loop order and make layer outer, need to trade off hit prefetching with copy-out of candidates
	 for (int ilay = Config::nlayers_per_seed; ilay < m_event_of_hits.m_n_layers; ++ilay)
	   {
	     BunchOfHits &bunch_of_hits = m_event_of_hits.m_layers_of_hits[ilay].m_bunches_of_hits[ebin];	     

             // XXX This should actually be done in some other thread for the next layer while
             // this thread is crunching the current one.
             // For now it's done in MkFitter::AddBestHit(), two loops before the data is needed.
             // for (int i = 0; i < bunch_of_hits.m_fill_index; ++i)
             // {
             //   _mm_prefetch((char*) & bunch_of_hits.m_hits[i], _MM_HINT_T1);
             // }

             //propagate to layer
#ifdef DEBUG
	     std::cout << "propagate to lay=" << ilay+1 << " start from x=" << mkfp->getPar(0, 0, 0) << " y=" << mkfp->getPar(0, 0, 1) << " z=" << mkfp->getPar(0, 0, 2)<< " r=" << std::sqrt(getRad2(mkfp->getPar(0, 0, 0), mkfp->getPar(0, 0, 1)))
		       << " px=" << mkfp->getPar(0, 0, 3) << " py=" << mkfp->getPar(0, 0, 4) << " pz=" << mkfp->getPar(0, 0, 5) << " pT=" << std::sqrt(getRad2(mkfp->getPar(0, 0, 3), mkfp->getPar(0, 0, 4))) << std::endl;
#endif
	     mkfp->PropagateTracksToR(4.*(ilay+1), end - itrack);
	 
#ifdef DEBUG
	     std::cout << "propagate to lay=" << ilay+1 << " arrive at x=" << mkfp->getPar(0, 1, 0) << " y=" << mkfp->getPar(0, 1, 1) << " z=" << mkfp->getPar(0, 1, 2)<< " r=" << std::sqrt(getRad2(mkfp->getPar(0, 1, 0), mkfp->getPar(0, 1, 1))) << std::endl;
	     std::cout << "now get hit range" << std::endl;
#endif
	   
	     mkfp->SelectHitRanges(bunch_of_hits, end - itrack);
	     
// #ifdef PRINTOUTS_FOR_PLOTS
// 	     std::cout << "MX number of hits in window in layer " << ilay << " is " <<  mkfp->getXHitEnd(0, 0, 0)-mkfp->getXHitBegin(0, 0, 0) << std::endl;
// #endif
	     //make candidates with best hit
#ifdef DEBUG
	     std::cout << "make new candidates" << std::endl;
#endif
	     mkfp->AddBestHit(bunch_of_hits);
	     
	     mkfp->SetNhits(ilay + 1);  //here again assuming one hit per layer (is this needed?)
	 
	   }//end of layer loop
	 
	 mkfp->OutputFittedTracksAndHitIdx(etabin_of_candidates.m_candidates, itrack, end);	 
       }//end of seed loop
     
   }//end of parallel section over seeds
}

//------------------------------------------------------------------------------
// FindTracks
//------------------------------------------------------------------------------

void MkBuilder::FindTracks(EventOfCombCandidates& event_of_comb_cands)
{
  // MT: partition recseeds into eta bins
  for (int iseed = 0; iseed < m_recseeds.size(); ++iseed)
  {
    if (m_recseeds[iseed].label() != iseed)
    {
      printf("Bad label for recseed %d -- %d\n", iseed, m_recseeds[iseed].label());
    }

    event_of_comb_cands.InsertSeed(m_recseeds[iseed]);
  }

  //dump seeds
#ifdef DEBUG
  for (int ebin = 0; ebin < Config::nEtaBin; ++ebin)
  {
    EtaBinOfCombCandidates &etabin_of_comb_candidates = event_of_comb_cands.m_etabins_of_comb_candidates[ebin]; 
    for (int iseed = 0; iseed < etabin_of_comb_candidates.m_fill_index; iseed++)
    {
      Track& seed = etabin_of_comb_candidates.m_candidates[iseed].front();
      std::cout << "MX - found seed with nFoundHits=" << seed.nFoundHits() << " chi2=" << seed.chi2() 
                << " x=" << seed.position()[0] << " y=" << seed.position()[1] << " z=" << seed.position()[2] 
                << " px=" << seed.momentum()[0] << " py=" << seed.momentum()[1] << " pz=" << seed.momentum()[2] 
                << " pT=" << sqrt(seed.momentum()[0]*seed.momentum()[0]+seed.momentum()[1]*seed.momentum()[1]) 
                << std::endl;
    }
  }
#endif

#ifdef DEBUG
  omp_lock_t writelock;

  omp_init_lock(&writelock);
#endif

  //the logic below is as follows:
  //- threads can be either over eta bins (a) or over seeds in one eta bin (b)
  //- for (a) we need the same number of eta bins in each thread
  //- for (b) we need the same number of threads in each eta bin
  assert( (Config::nEtaBin % NUM_THREADS == 0) || (NUM_THREADS % Config::nEtaBin == 0) );

  //parallel section over seeds
  //int nseeds=recseeds.size();
#pragma omp parallel num_threads(NUM_THREADS)
   {
     int thread_num = omp_get_thread_num();
     int num_threads = omp_get_num_threads();
     std::cout << "Thread " << thread_num << "/" << num_threads << std::endl;
     int n_th_per_eta_bin = num_threads/Config::nEtaBin;
     int n_eta_bin_per_th = Config::nEtaBin/num_threads;

     int th_start_ebin=-1,th_end_ebin=-1;

     if (n_th_per_eta_bin>=1) {

       // case (b): there is only one eta bin per thread (i.e. >1 thread per eta bin), we'll split seeds in different threads below
       th_start_ebin = thread_num/n_th_per_eta_bin;
       th_end_ebin = th_start_ebin+1;

     } else {

       //case (a): define first and last eta bin for this thread
       int ebin_idx_in_th = thread_num * n_eta_bin_per_th;
       th_start_ebin = ebin_idx_in_th;
       th_end_ebin = th_start_ebin + n_eta_bin_per_th;       

     }

#ifdef DEBUG
     omp_set_lock(&writelock);
     std::cout << "th_start_ebin-a=" << thread_num * n_eta_bin_per_th << " th_end_ebin-a=" << thread_num * n_eta_bin_per_th + n_eta_bin_per_th;
     if (n_th_per_eta_bin>=1) {
       std::cout << " th_start_ebin-b=" << thread_num/n_th_per_eta_bin << " th_end_ebin-b=" << thread_num/n_th_per_eta_bin+1;
     }
     std::cout << std::endl;
     omp_unset_lock(&writelock);
#endif

     //loop over eta bins
     for (int ebin = th_start_ebin; ebin < th_end_ebin; ++ebin)
     {

       EtaBinOfCombCandidates &etabin_of_comb_candidates = event_of_comb_cands.m_etabins_of_comb_candidates[ebin];

       int th_start_seed=-1,th_end_seed=-1;
       int nseeds_ebin = event_of_comb_cands.m_etabins_of_comb_candidates[ebin].m_fill_index;
       if (th_end_ebin==th_start_ebin+1) {
	 // case (b): define first and last seed in this eta bin for this thread
	 int th_idx_in_ebin = thread_num % n_th_per_eta_bin;       
	 th_start_seed = th_idx_in_ebin * nseeds_ebin / n_th_per_eta_bin;
	 th_end_seed = std::min( (th_idx_in_ebin + 1) * nseeds_ebin / n_th_per_eta_bin, nseeds_ebin );
       } else {
	 // case (a): we process >=1 full eta bins in this thread, se we need to loop over all seeds in each eta bin
	 th_start_seed=0;
	 th_end_seed= etabin_of_comb_candidates.m_fill_index;
       }

#ifdef DEBUG
       omp_set_lock(&writelock);
       std::cout << "th_start_seed-a=" << 0 << " th_end_seed-a=" << etabin_of_comb_candidates.m_fill_index;
       if (n_th_per_eta_bin>=1) {
         std::cout << " th_start_seed-b=" << (thread_num % n_th_per_eta_bin) * nseeds_ebin / n_th_per_eta_bin << " th_end_seed-b=" << std::min( ( (thread_num % n_th_per_eta_bin)+ 1) * nseeds_ebin / n_th_per_eta_bin, nseeds_ebin );
       }
       std::cout << std::endl;
       omp_unset_lock(&writelock);
#endif

       //ok now we start looping over layers
       //loop over layers, starting from after the seeD
       for (int ilay = Config::nlayers_per_seed; ilay < m_event_of_hits.m_n_layers; ++ilay)
	 {
	   BunchOfHits &bunch_of_hits = m_event_of_hits.m_layers_of_hits[ilay].m_bunches_of_hits[ebin];	     
	   
#ifdef DEBUG
	   std::cout << "processing lay=" << ilay+1 << std::endl;
#endif
	   
	   //prepare unrolled vector to loop over
	   std::vector<std::pair<int,int> > seed_cand_idx;
	   for (int iseed = th_start_seed; iseed != th_end_seed; ++iseed) 
	     {
	       for (int ic = 0; ic<etabin_of_comb_candidates.m_candidates[iseed].size(); ++ic)
		 {
		   seed_cand_idx.push_back(std::pair<int,int>(iseed,ic));
		 }
	     }
	   int theEndCand = seed_cand_idx.size();     
	   
	   std::vector<std::vector<Track> > tmp_candidates(th_end_seed-th_start_seed);     
	   for (int iseed=0;iseed<tmp_candidates.size();++iseed)
           {
             // XXXX MT: Tried adding 25 to reserve below as I was seeing some
             // time spent in push_back ... but it didn't really help.
             // We need to optimize this by throwing away and replacing the worst
             // candidate once a better one arrives. This will also avoid sorting.
	     tmp_candidates[iseed].reserve(2*Config::maxCand);//factor 2 seems reasonable to start with
	   }
	   
	   //vectorized loop
	   for (int itrack = 0; itrack < theEndCand; itrack += NN)
	     {
	       
	       int end = std::min(itrack + NN, theEndCand);
	       
#ifdef DEBUG
	       std::cout << "processing track=" << itrack << std::endl;
#endif
	       
	       MkFitter *mkfp = m_mkfp_arr[omp_get_thread_num()];
	       
	       mkfp->SetNhits(ilay);//here again assuming one hit per layer
	       
	       //fixme find a way to deal only with the candidates needed in this thread
	       mkfp->InputTracksAndHitIdx(etabin_of_comb_candidates.m_candidates, seed_cand_idx, itrack, end);
	       
	       //propagate to layer
#ifdef DEBUG
	       std::cout << "propagate to lay=" << ilay+1 << " start from x=" << mkfp->getPar(0, 0, 0) << " y=" << mkfp->getPar(0, 0, 1) << " z=" << mkfp->getPar(0, 0, 2)<< " r=" << std::sqrt(getRad2(mkfp->getPar(0, 0, 0), mkfp->getPar(0, 0, 1)))
			 << " px=" << mkfp->getPar(0, 0, 3) << " py=" << mkfp->getPar(0, 0, 4) << " pz=" << mkfp->getPar(0, 0, 5) << " pT=" << std::sqrt(getRad2(mkfp->getPar(0, 0, 3), mkfp->getPar(0, 0, 4))) << std::endl;
#endif
    		 mkfp->PropagateTracksToR(4.*(ilay+1), end - itrack);
	       
#ifdef DEBUG
	       std::cout << "propagate to lay=" << ilay+1 << " arrive at x=" << mkfp->getPar(0, 1, 0) << " y=" << mkfp->getPar(0, 1, 1) << " z=" << mkfp->getPar(0, 1, 2)<< " r=" << std::sqrt(getRad2(mkfp->getPar(0, 1, 0), mkfp->getPar(0, 1, 1))) << std::endl;
	       std::cout << "now get hit range" << std::endl;
#endif
	       
	       mkfp->SelectHitRanges(bunch_of_hits, end - itrack);
	       
//#ifdef PRINTOUTS_FOR_PLOTS
//std::cout << "MX number of hits in window in layer " << ilay << " is " <<  mkfp->getXHitEnd(0, 0, 0)-mkfp->getXHitBegin(0, 0, 0) << std::endl;
//#endif

	       //make candidates with best hit
#ifdef DEBUG
	       std::cout << "make new candidates" << std::endl;
#endif
	       mkfp->FindCandidates(bunch_of_hits,tmp_candidates,th_start_seed);
	       
	     }//end of vectorized loop
	   
	   //clean exceeding candidates per seed
	   for (int is=0;is<tmp_candidates.size();++is)
	     {
	       if (tmp_candidates[is].size()>Config::maxCand)
		 {
#ifdef DEBUG
		   std::cout << "erase extra candidates" << std::endl;
#endif	     
		   std::sort(tmp_candidates[is].begin(), tmp_candidates[is].end(), sortCandByHitsChi2);
		   tmp_candidates[is].erase(tmp_candidates[is].begin()+Config::maxCand,tmp_candidates[is].end());
		 }
	     } 
	   //now swap with input candidates
	   for (int is=0;is<tmp_candidates.size();++is)
	     {
	       if (tmp_candidates[is].size()>0)
		 {
		   etabin_of_comb_candidates.m_candidates[th_start_seed+is].swap(tmp_candidates[is]);
		   tmp_candidates[is].clear();
		 }
	       else 
		 {
		   //we do nothing in the SM version here, I think we should put these in the output and avoid keeping looping over them
		 }
	     }
	   
	 }//end of layer loop

       //final sorting
       int nCandsBeforeEnd = 0;
       for (int iseed=th_start_seed;iseed<th_end_seed;++iseed) 
	 {
	   std::vector<Track>& finalcands = etabin_of_comb_candidates.m_candidates[iseed];
	   if (finalcands.size()==0) continue;
	   std::sort(finalcands.begin(), finalcands.end(), sortCandByHitsChi2);
	 }

     }//end of loop over eta bins

   }//end of parallel section over seeds
}


//==============================================================================
// runBuildingTestPlexBestHit
//==============================================================================

double runBuildingTestPlexBestHit(Event& ev, EventOfCandidates& event_of_cands)
{
  MkBuilder builder;

  builder.begin_event(ev, __func__);

  double time = dtime();

#ifdef USE_VTUNE_PAUSE
  __itt_resume();
#endif

  builder.fit_seeds();

  //EventOfCandidates event_of_cands;
  event_of_cands.Reset();

  builder.FindTracksBestHit(event_of_cands);

#ifdef USE_VTUNE_PAUSE
  __itt_pause();
#endif

   time = dtime() - time;

   builder.quality_reset();

   for (int ebin = 0; ebin < Config::nEtaBin; ++ebin)
   {
     EtaBinOfCandidates &etabin_of_candidates = event_of_cands.m_etabins_of_candidates[ebin]; 

     for (int itrack = 0; itrack < etabin_of_candidates.m_fill_index; itrack++)
     {
       builder.quality_process(etabin_of_candidates.m_candidates[itrack]);
     }
   }

   builder.quality_print();

   builder.end_event();

   return time;
}


//==============================================================================
// runBuildingTestPlex
//==============================================================================

double runBuildingTestPlex(Event& ev, EventOfCombCandidates& event_of_comb_cands)
{
  MkBuilder builder;

  builder.begin_event(ev, __func__);

  double time = dtime();

#ifdef USE_VTUNE_PAUSE
  __itt_resume();
#endif

  builder.fit_seeds();

  //EventOfCombCandidates event_of_comb_cands;
  event_of_comb_cands.Reset();

  builder.FindTracks(event_of_comb_cands);

#ifdef USE_VTUNE_PAUSE
  __itt_pause();
#endif

   time = dtime() - time;

   //dump tracks
   //std::cout << "found total tracks=" << recseeds.size() << std::endl;
   builder.quality_reset();
   {
     int cnt=0, cnt1=0, cnt2=0, cnt_8=0, cnt1_8=0, cnt2_8=0, cnt_nomc=0;
     for (int ebin = 0; ebin < Config::nEtaBin; ++ebin)
     {
       EtaBinOfCombCandidates &etabin_of_comb_candidates = event_of_comb_cands.m_etabins_of_comb_candidates[ebin]; 

       for (int iseed = 0; iseed < etabin_of_comb_candidates.m_fill_index; iseed++)
       {
	 //take the first one!
         builder.quality_process(etabin_of_comb_candidates.m_candidates[iseed].front());
       }
     }
   }

   builder.quality_print();

   builder.end_event();

   return time;
}
