#include "CandCloner.h"

#include "HitStructures.h"

//#define DEBUG
#include "Debug.h"

namespace
{
inline bool sortCandListByHitsChi2(const mkfit::MkFinder::IdxChi2List& cand1,
                                   const mkfit::MkFinder::IdxChi2List& cand2)
{
  if (cand1.nhits == cand2.nhits) return cand1.chi2 < cand2.chi2;
  return cand1.nhits > cand2.nhits;
}

inline bool sortCandListByScore(const mkfit::MkFinder::IdxChi2List& cand1,
				const mkfit::MkFinder::IdxChi2List& cand2)
{


  float validHitBonus_=2.5;
  float missingHitPenalty_=20.0;
  int nfoundhits[2] = {cand1.nhits,cand2.nhits};
  int nmisshits[2] = {cand1.nholes,cand2.nholes};
  float chi2[2] = {cand1.chi2,cand2.chi2};
  float pt[2] = {cand1.pt,cand2.pt};
  float score[2] = {0.f,0.f};
  for(int c=0; c<2; ++c){
    score[c] = (validHitBonus_)*nfoundhits[c] - (missingHitPenalty_)*nmisshits[c] - chi2[c];
    if(pt[c]<0.9f) score[c] -= 0.5f*(validHitBonus_)*nfoundhits[c];
    else if(nfoundhits[c]>8) score[c] += (validHitBonus_)*nfoundhits[c];
  }
  return score[0]>score[1];

}
}

namespace mkfit {

//==============================================================================

void CandCloner::ProcessSeedRange(int is_beg, int is_end)
{
  // Process new hits for a range of seeds.

  // bool debug = true;

  dprintf("\nCandCloner::ProcessSeedRange is_beg=%d, is_end=%d\n", is_beg, is_end);

  //1) sort the candidates
  for (int is = is_beg; is < is_end; ++is)
  {
    std::vector<MkFinder::IdxChi2List>& hitsForSeed = m_hits_to_add[is];

    std::vector<CombCandidate> &cands = mp_event_of_comb_candidates->m_candidates;

#ifdef DEBUG
    int th_start_seed = m_start_seed;

    dprint("  seed n " << is << " with input candidates=" << hitsForSeed.size());
    for (int ih = 0; ih < (int) hitsForSeed.size(); ih++)
    {
      dprint("trkIdx=" << hitsForSeed[ih].trkIdx << " hitIdx=" << hitsForSeed[ih].hitIdx << " chi2=" <<  hitsForSeed[ih].chi2 << std::endl
             << "    "
             << "original pt=" << cands[th_start_seed+is][hitsForSeed[ih].trkIdx].pT() << " "
             << "nTotalHits="  << cands[th_start_seed+is][hitsForSeed[ih].trkIdx].nTotalHits() << " "
             << "nFoundHits="  << cands[th_start_seed+is][hitsForSeed[ih].trkIdx].nFoundHits() << " "
             << "chi2="        << cands[th_start_seed+is][hitsForSeed[ih].trkIdx].chi2());
    }
#endif

    if ( ! hitsForSeed.empty())
    {
      //sort the new hits
      //std::sort(hitsForSeed.begin(), hitsForSeed.end(), sortCandListByHitsChi2);
      std::sort(hitsForSeed.begin(), hitsForSeed.end(), sortCandListByScore);

      int num_hits = std::min((int) hitsForSeed.size(), Config::maxCandsPerSeed);

      // This is from buffer, we know it was cleared after last usage.
      std::vector<Track> &cv = t_cands_for_next_lay[is - is_beg];

      for (int ih = 0; ih < num_hits; ih++)
      {
        const MkFinder::IdxChi2List &h2a = hitsForSeed[ih];

        cv.push_back( cands[ m_start_seed + is ][ h2a.trkIdx ] );

        if (h2a.hitIdx != -4) // Could also skip storing of -3 hits.
        {
          cv.back().addHitIdx(h2a.hitIdx, m_layer, 0);
          cv.back().setChi2(h2a.chi2);

          if (h2a.hitIdx >= 0)
          {
            mp_kalman_update_list->push_back(std::pair<int,int>(m_start_seed + is, ih));
          }
        }
      }

      // Copy the best -2 cands back to the current list while there is room.
      if (num_hits < Config::maxCandsPerSeed)
      {
        const std::vector<Track> &ov = cands[m_start_seed + is];
        const int max_m2 = ov.size();

        int cur_m2 = 0;
        while (cur_m2 < max_m2 && ov[cur_m2].getLastHitIdx() != -2) ++cur_m2;
        while (cur_m2 < max_m2 && num_hits < Config::maxCandsPerSeed)
        {
          cv.push_back( ov[cur_m2++] );
          ++num_hits;
        }
      }

      //cands[ m_start_seed + is ].swap(cv); // segfaulting w/ backwards fit input tracks -- using loop below now
      cands[m_start_seed + is].resize(cv.size());
      for (size_t ii = 0; ii < cv.size(); ++ii)
      {
	memcpy( & cands[m_start_seed + is][ii], & cv[ii], sizeof(Track));
      }
      cv.clear();
    }
    else // hitsForSeed.empty()
    {
      if (cands[m_start_seed + is].m_state == CombCandidate::Finding)
      {
        cands[m_start_seed + is].clear();
      }

      // Cross-check for what is left once there are no more changes for a whole seed.
      //
      // for (auto &cand : cands[ m_start_seed + is ])
      // {
      //   if (cand.getLastHitIdx() < 0 && cand.getLastHitIdx() != -2)
      //   {
      //     printf("  We got cand with index %d\n", cand.getLastHitIdx());
      //     print("offending track (unknown index)", -666, cand, true);
      //   }
      // }
    }
  }
}

} // end namespace mkfit
