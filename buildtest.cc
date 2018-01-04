#include "buildtest.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "Event.h"
//#define DEBUG
#include "Debug.h"

#include <cmath>
#include <iostream>

//typedef std::pair<Track, TrackState> cand_t;
typedef Track cand_t;
typedef TrackVec::const_iterator TrkIter;

#ifndef TBB
typedef std::vector<cand_t> candvec;
#else
#include "tbb/tbb.h"
#include <mutex>
// concurrent_vector is only needed if we parallelize the candidate loops;
// not needed if we only parallelize over seeds.
//typedef tbb::concurrent_vector<cand_t> candvec;
typedef std::vector<cand_t> candvec;
static std::mutex evtlock;
#endif
typedef candvec::const_iterator canditer;

void extendCandidate(const BinInfoMap & segmentMap, const Event& ev, const cand_t& cand, candvec& tmp_candidates, int ilay, int seedID, bool debug);

void processCandidates(const BinInfoMap & segmentMap, Event& ev, candvec& candidates, int ilay, int seedID, const bool debug)
{
  auto& evt_track_candidates(ev.candidateTracks_);

  candvec tmp_candidates;
  tmp_candidates.reserve(3*candidates.size()/2);

  if (candidates.size() > 0)
  {
    //loop over running candidates
    for (auto&& cand : candidates) {
      extendCandidate(segmentMap, ev, cand, tmp_candidates, ilay, seedID, debug);
    }
    if (tmp_candidates.size()>Config::maxCandsPerSeed) {
      dprint("huge size=" << tmp_candidates.size() << " keeping best "<< Config::maxCandsPerSeed << " only");
      std::partial_sort(tmp_candidates.begin(),tmp_candidates.begin()+(Config::maxCandsPerSeed-1),tmp_candidates.end(),sortByHitsChi2);
      tmp_candidates.resize(Config::maxCandsPerSeed); // thread local, so ok not thread safe
    } else if (tmp_candidates.size()==0) {
      // save the best candidate from the previous iteration and then swap in
      // the empty new candidate list; seed will be skipped on future iterations
      auto best = std::min_element(candidates.begin(),candidates.end(),sortByHitsChi2);
#ifdef TBB
      std::lock_guard<std::mutex> evtguard(evtlock); // should be rare
#endif
      best->setLabel(evt_track_candidates.size());
      dprint("no more candidates, saving track from seed " << seedID << " label " << best->label() << " hits " 
                                       << best->nFoundHits() << " parameters " << best->parameters() << std::endl);
      evt_track_candidates.push_back(*best);
      ev.candidateTracksExtra_.emplace_back(seedID);
    }
    dprint("swapping with size=" << tmp_candidates.size());
    candidates.swap(tmp_candidates);
    tmp_candidates.clear();
  }
}

void buildTracksBySeeds(const BinInfoMap & segmentMap, Event& ev)
{
  auto& evt_track_candidates(ev.candidateTracks_);
  const auto& evt_lay_hits(ev.layerHits_);
  const auto& evt_seeds(ev.seedTracks_);
  const auto& evt_seeds_extra(ev.seedTracksExtra_);
  bool debug(true);

  std::vector<candvec> track_candidates(evt_seeds.size());
  for (auto iseed = 0; iseed < evt_seeds.size(); iseed++) {
    const auto& seed(evt_seeds[iseed]);
    track_candidates[iseed].push_back(seed);
  }

#ifdef TBB
  //loop over seeds
  parallel_for( tbb::blocked_range<size_t>(0, evt_seeds.size()), 
      [&](const tbb::blocked_range<size_t>& seediter) {
    for (auto iseed = seediter.begin(); iseed != seediter.end(); ++iseed) {
      const auto& seed(evt_seeds[iseed]);
#else
    for (auto iseed = 0; iseed != evt_seeds.size(); ++iseed) {
      const auto& seed(evt_seeds[iseed]);
#endif
      dprint("processing seed # " << iseed << " par=" << seed.parameters());
      dprint("from MC track par = " << ev.simTracks_[iseed].parameters()); // will be wrong with real seeding
      TrackState seed_state = seed.state();
      //seed_state.errors *= 0.01;//otherwise combinatorics explode!!!
      //should consider more than 1 candidate...
      auto&& candidates(track_candidates[iseed]);
      for (int ilay=Config::nlayers_per_seed;ilay<evt_lay_hits.size();++ilay) {//loop over layers, starting from after the seed
        dprint("going to layer #" << ilay << " with N cands=" << track_candidates.size());
        processCandidates(segmentMap, ev, candidates, ilay, evt_seeds_extra[iseed].seedID(), debug);
      }
      //end of layer loop
    }//end of process seeds loop
#ifdef TBB
  });
#endif
  for (auto iseed = 0; iseed < track_candidates.size(); ++iseed) {
    auto& cand = track_candidates[iseed];
    if (cand.size()>0) {
      // only save one track candidate per seed, one with lowest chi2
      //std::partial_sort(cand.begin(),cand.begin()+1,cand.end(),sortByHitsChi2);
      auto best = std::min_element(cand.begin(),cand.end(),sortByHitsChi2);
      best->setLabel(evt_track_candidates.size());
      dprint("Saving track from seed " << iseed << " label " << best->label() << " hits " 
                                       << best->nFoundHits() << " parameters " << best->parameters() << std::endl);
      evt_track_candidates.push_back(*best);
      ev.candidateTracksExtra_.emplace_back(evt_seeds_extra[iseed].seedID());
    }
  }
}		

void buildTracksByLayers(const BinInfoMap & segmentMap, Event& ev)
{
  auto& evt_track_candidates(ev.candidateTracks_);
  const auto& evt_lay_hits(ev.layerHits_);
  const auto& evt_seeds(ev.seedTracks_);
  const auto& evt_seeds_extra(ev.seedTracksExtra_);
  bool debug(true);

  std::vector<candvec> track_candidates(evt_seeds.size());
  for (auto iseed = 0; iseed < evt_seeds.size(); iseed++) {
    const auto& seed(evt_seeds[iseed]);
    track_candidates[iseed].push_back(seed);
  }

  //loop over layers, starting from after the seed
  for (auto ilay = Config::nlayers_per_seed; ilay < evt_lay_hits.size(); ++ilay) {
    dprint("going to layer #" << ilay << " with N cands = " << track_candidates.size());

#ifdef TBB
    //loop over seeds
    parallel_for( tbb::blocked_range<size_t>(0, evt_seeds.size()), 
        [&](const tbb::blocked_range<size_t>& seediter) {
      for (auto iseed = seediter.begin(); iseed != seediter.end(); ++iseed) {
        auto&& candidates(track_candidates[iseed]);
        processCandidates(segmentMap, ev, candidates, ilay, evt_seeds_extra[iseed].seedID(), debug);
      }
    }); //end of process seeds loop
#else
    //process seeds
    for (auto iseed = 0; iseed != evt_seeds.size(); ++iseed) {
      const auto& seed(evt_seeds[iseed]);
      auto&& candidates(track_candidates[iseed]);
      processCandidates(segmentMap, ev, candidates, ilay, evt_seeds_extra[iseed].seedID(), debug);
    }
#endif
  } //end of layer loop

  //std::lock_guard<std::mutex> evtguard(evtlock);
  for (auto iseed = 0; iseed < track_candidates.size(); ++iseed) {
    auto& cand = track_candidates[iseed];
    if (cand.size()>0) {
      // only save one track candidate per seed, one with lowest chi2
      //std::partial_sort(cand.begin(),cand.begin()+1,cand.end(),sortByHitsChi2);
      auto best = std::min_element(cand.begin(),cand.end(),sortByHitsChi2);
      best->setLabel(evt_track_candidates.size());
      evt_track_candidates.push_back(*best);
      ev.candidateTracksExtra_.emplace_back(evt_seeds_extra[iseed].seedID());
    }
  }
}

void extendCandidate(const BinInfoMap & segmentMap, const Event& ev, const cand_t& cand, candvec& tmp_candidates, int ilayer, int seedID, bool debug)
{
  std::vector<int> branch_hit_indices; // temp variable for validation... could be used for cand hit builder engine!
  const Track& tkcand = cand;
  const TrackState& updatedState = cand.state();
  const auto& evt_lay_hits(ev.layerHits_);
  const auto& segLayMap(segmentMap[ilayer]);

  const PropagationFlags pflags(PF_none);

  dprint("processing candidate with nHits=" << tkcand.nFoundHits());
#ifdef LINEARINTERP
  TrackState propState = propagateHelixToR(updatedState, ev.geom_.Radius(ilayer), pflags);
#else
#ifdef TBB
#error "Invalid combination of options (thread safety)"
#endif
  TrackState propState = propagateHelixToLayer(updatedState, ilayer,ev.geom_, pflags);
#endif // LINEARINTERP
#ifdef CHECKSTATEVALID
  if (!propState.valid) {
    return;
  }
#endif
  const float predx  = propState.parameters.At(0);
  const float predy  = propState.parameters.At(1);
  const float predz  = propState.parameters.At(2);

  const float eta  = getEta(std::sqrt(getRad2(predx,predy)),predz);
  const float deta = std::sqrt(std::abs(getEtaErr2(predx,predy,predz,propState.errors.At(0,0),propState.errors.At(1,1),propState.errors.At(2,2),propState.errors.At(0,1),propState.errors.At(0,2),propState.errors.At(1,2))));
  const float nSigmaDeta = std::min(std::max(Config::nSigma*deta,(float) Config::minDEta), (float) Config::maxDEta); // something to tune -- minDEta = 0.0
  const auto etaBinMinus = getEtaPartition(eta-nSigmaDeta);
  const auto etaBinPlus  = getEtaPartition(eta+nSigmaDeta);

  const float phi    = getPhi(predx,predy); //std::atan2(predy,predx); 
  const float dphi = std::sqrt(std::abs(getPhiErr2(predx,predy,propState.errors.At(0,0),propState.errors.At(1,1),propState.errors.At(0,1))));
  const float nSigmaDphi = std::min(std::max(Config::nSigma*dphi,(float) Config::minDPhi), (float) Config::maxDPhi);
  const auto phiBinMinus = getPhiPartition(normalizedPhi(phi-nSigmaDphi));
  const auto phiBinPlus  = getPhiPartition(normalizedPhi(phi+nSigmaDphi));

  dprint("propState at layer: " << ilayer << ": " << propState.parameters);
  dcall(dumpMatrix(propState.errors));
  // get candidate hits for this track candidate at this layer
  std::vector<int> cand_hit_indices = getCandHitIndices(etaBinMinus,etaBinPlus,phiBinMinus,phiBinPlus,segLayMap);

#ifdef LINEARINTERP
    const float minR = ev.geom_.Radius(ilayer);
    float maxR = minR;
    for (auto&& cand_hit_idx : cand_hit_indices){
      const float candR = evt_lay_hits[ilayer][cand_hit_idx].r();
      if (candR > maxR) maxR = candR;
    }
    const float deltaR = maxR - minR;
    dprint("min, max, delta: " << minR << ", " << maxR << ", " << deltaR);
    const TrackState propStateMin = propState;
    const TrackState propStateMax = propagateHelixToR(updatedState, maxR, pflags);
#ifdef CHECKSTATEVALID
    if (!propStateMax.valid) {
      return;
    }
#endif
#endif
  
    for (auto&& cand_hit_idx : cand_hit_indices){
      const Hit hitCand = evt_lay_hits[ilayer][cand_hit_idx];
      const MeasurementState hitMeas = hitCand.measurementState();

#ifdef LINEARINTERP
      if (hitCand.r() > minR && deltaR > 0.0f) {
        const float ratio = (hitCand.r() - minR)/deltaR;
        dprint("ratio " << ratio);
        propState.parameters = (1.0-ratio)*propStateMin.parameters + ratio*propStateMax.parameters;
      }
      dprint("interp propState " << propState.parameters << std::endl
                                 << propStateMax.parameters - propStateMin.parameters);
      dcall(dumpMatrix(propState.errors));
      dprint(hitMeas.parameters());
      dcall(dumpMatrix(hitMeas.errors()));
#endif
      dprint(propState.position() - hitMeas.parameters());
      const float chi2 = computeChi2(propState,hitMeas);
      dprint("found hit with index: " << cand_hit_idx << " from sim track " 
	     << ev.simHitsInfo_[evt_lay_hits[ilayer][cand_hit_idx].mcHitID()].mcTrackID()
	     << " chi2=" << chi2 << std::endl);
    
      if ((chi2<Config::chi2Cut)&&(chi2>0.)) {//fixme 
        const TrackState tmpUpdatedState = updateParameters(propState, hitMeas);
        Track tmpCand = tkcand.clone();
        tmpCand.addHitIdx(cand_hit_idx, ilayer, chi2);
        tmpCand.setState(tmpUpdatedState);
        tmp_candidates.push_back(tmpCand);
        branch_hit_indices.push_back(cand_hit_idx); // validation
      }
    }//end of consider hits on layer loop

  //add also the candidate for no hit found
  if (tkcand.nFoundHits()==ilayer) {//only if this is the first missing hit
    dprint("adding candidate with no hit");
    Track tmpCand = tkcand.clone();
    tmpCand.addHitIdx(-1, ilayer, 0.0f);
    tmp_candidates.push_back(tmpCand);  // fix this once moving to fix indices
    branch_hit_indices.push_back(Config::nTracks); // since tracks go from 0-Config::nTracks -1, the ghost index is just the one beyond
  }
}
