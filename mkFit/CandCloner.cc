#include "CandCloner.h"

#include "HitStructures.h"

//#define DEBUG
#include "Debug.h"

namespace
{
inline bool sortCandListByHitsChi2(const mkfit::IdxChi2List& cand1,
                                   const mkfit::IdxChi2List& cand2)
{
  if (cand1.nhits == cand2.nhits) return cand1.chi2 < cand2.chi2;
  return cand1.nhits > cand2.nhits;
}

inline bool sortCandListByScore(const mkfit::IdxChi2List& cand1,
				const mkfit::IdxChi2List& cand2)
{
  return mkfit::sortByScoreStruct(cand1,cand2);
}
}

namespace mkfit {

//==============================================================================

void CandCloner::ProcessSeedRange(int is_beg, int is_end)
{
  // Process new hits for a range of seeds.

  // bool debug = true;

  dprintf("\nCandCloner::ProcessSeedRange is_beg=%d, is_end=%d\n", is_beg, is_end);

  std::vector<CombCandidate> &cands = mp_event_of_comb_candidates->m_candidates;

  //1) sort the candidates
  for (int is = is_beg; is < is_end; ++is)
  {
    std::vector<IdxChi2List>& hitsForSeed = m_hits_to_add[is];

    CombCandidate          &ccand  = cands[m_start_seed + is];
    std::vector<TrackCand> &extras = (*mp_extra_cands)[is];
    auto extra_i = extras.begin();
    auto extra_e = extras.end();

    // Extras sorted now
    // std::sort(extras.begin(), extras.end(), sortByScoreCand);

#ifdef DEBUG
    dprint("  seed n " << is << " with input candidates=" << hitsForSeed.size());
    for (int ih = 0; ih < (int) hitsForSeed.size(); ih++)
    {
      dprint("trkIdx=" << hitsForSeed[ih].trkIdx << " hitIdx=" << hitsForSeed[ih].hitIdx << " chi2=" <<  hitsForSeed[ih].chi2 << std::endl
             << "    "
             << "original pt=" << ccand[hitsForSeed[ih].trkIdx].pT() << " "
             << "nTotalHits="  << cdand[hitsForSeed[ih].trkIdx].nTotalHits() << " "
             << "nFoundHits="  << ccand[hitsForSeed[ih].trkIdx].nFoundHits() << " "
             << "chi2="        << ccand[hitsForSeed[ih].trkIdx].chi2());
    }
#endif

    if ( ! hitsForSeed.empty())
    {
      //sort the new hits
      std::sort(hitsForSeed.begin(), hitsForSeed.end(), sortCandListByScore);

      int num_hits = std::min((int) hitsForSeed.size(), Config::maxCandsPerSeed);

      // This is from buffer, we know it was cleared after last usage.
      std::vector<TrackCand> &cv = t_cands_for_next_lay[is - is_beg];

      int n_pushed = 0;

      for (int ih = 0; ih < num_hits; ih++)
      {
        const IdxChi2List &h2a = hitsForSeed[ih];

        TrackCand cc( ccand[h2a.trkIdx] );
        cc.addHitIdx(h2a.hitIdx, m_layer, 0);
        cc.setChi2(h2a.chi2);
        cc.setCandScore(h2a.score);
        // h2a already carries correct score
        // cc.setCandScore( getScoreCand( cc );

        if (h2a.hitIdx == -2)
        {
          if (h2a.score > ccand.m_best_short_cand.getCandScore())
          {
            ccand.m_best_short_cand = cc;
          }
          continue;
        }

        // Could also skip storing of cands with last -3 hit.

        // Squeeze in extra tracks that are better than current one.
        while (extra_i != extra_e && sortByScoreTrackCand(*extra_i, cc) && n_pushed < Config::maxCandsPerSeed)
        {
          cv.emplace_back(*extra_i);
          ++n_pushed;
          ++extra_i;
        }

        if (n_pushed >= Config::maxCandsPerSeed)
          break;

        cv.emplace_back( cc );
        ++n_pushed;

        if (h2a.hitIdx >= 0)
        {
          mp_kalman_update_list->push_back(std::pair<int,int>(m_start_seed + is, n_pushed - 1));
        }
      }

      // Add remaining extras as long as there is still room for them.
      while (extra_i != extra_e && n_pushed < Config::maxCandsPerSeed)
      {
        cv.emplace_back(*extra_i);
        ++n_pushed;
        ++extra_i;
      }

      //ccand.swap(cv); // segfaulting w/ backwards fit input tracks -- using loop below now
      ccand.resize(cv.size());
      for (size_t ii = 0; ii < cv.size(); ++ii)
      {
        ccand[ii] = cv[ii];
      }
      cv.clear();
    }
    else // hitsForSeed.empty()
    {
      if (ccand.m_state == CombCandidate::Finding)
      {
        ccand.clear();

        while (extra_i != extra_e)
        {
          ccand.emplace_back(*extra_i);
          ++extra_i;
        }
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

    extras.clear();
  }
}

} // end namespace mkfit
