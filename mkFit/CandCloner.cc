#include "CandCloner.h"

namespace
{
bool sortCandListByHitsChi2(const MkFitter::IdxChi2List& cand1,
                            const MkFitter::IdxChi2List& cand2)
{
  if (cand1.nhits == cand2.nhits) return cand1.chi2 < cand2.chi2;
  return cand1.nhits > cand2.nhits;
}
}

//==============================================================================

void CandCloner::ProcessSeedRange(int is_beg, int is_end)
{
  // Like cut-n-pasta but work on a range of seeds.
  //
  // XXXX Blak, this is silly ... I can do 1 as the seeds become available and once number
  // of candidates to process is large enough step into that section (2).
  // 3 is just the extension of 2, assuming all candidates from a given seed get processed.
  // Furthermore ... once 3 is done, those seeds can proceed to the next layer.
  // Hi, funny idea ... how about having a thread per layer, passing candidtes in :)
  //
  // Talked with G ... the plan is to have num candidates 2 4 or 8 ... and thus a reasonable
  // LCM with the vector unit width.

  const int is_num = is_end - is_beg;

  // printf("CandCloner::ProcessSeedRange is_beg=%d, is_end=%d, is_num=%d\n", is_beg, is_end, is_num);

  //1) sort the candidates
  for (int is = is_beg; is < is_end; ++is)
  {
    std::vector<MkFitter::IdxChi2List>& hitsToAddForThisSeed = m_hits_to_add[is];
#ifdef DEBUG
    EtaBinOfCombCandidates &etabin_of_comb_candidates = *mp_etabin_of_comb_candidates;
    int th_start_seed = m_start_seed;

    std::cout << "dump seed n " << is << " with input candidates=" << hitsToAddForThisSeed.size() << std::endl;
    for (int ih = 0; ih<hitsToAddForThisSeed.size(); ih++)
    {
      std::cout << "trkIdx=" << hitsToAddForThisSeed[ih].trkIdx << " hitIdx=" << hitsToAddForThisSeed[ih].hitIdx << " chi2=" <<  hitsToAddForThisSeed[ih].chi2 << std::endl;
      std::cout << "original pt=" << etabin_of_comb_candidates.m_candidates[th_start_seed+is][hitsToAddForThisSeed[ih].trkIdx].pT() << " " 
                << "nTotalHits=" << etabin_of_comb_candidates.m_candidates[th_start_seed+is][hitsToAddForThisSeed[ih].trkIdx].nTotalHits() << " " 
                << "nFoundHits=" << etabin_of_comb_candidates.m_candidates[th_start_seed+is][hitsToAddForThisSeed[ih].trkIdx].nFoundHits() << " " 
                << "chi2=" << etabin_of_comb_candidates.m_candidates[th_start_seed+is][hitsToAddForThisSeed[ih].trkIdx].chi2() << " " 
                << std::endl;
    }
#endif

    //sort the new hits
    std::sort(hitsToAddForThisSeed.begin(), hitsToAddForThisSeed.end(), sortCandListByHitsChi2);

    const int num_seeds = std::min((int) hitsToAddForThisSeed.size(), Config::maxCandsPerSeed);

    if (num_seeds > 0)
    {
      std::vector<Track> &cv = t_cands_for_next_lay[is - is_beg];

      for (int ih = 0; ih < num_seeds; ih++)
      {
        MkFitter::IdxChi2List &h2a = hitsToAddForThisSeed[ih];
        cv.push_back( mp_etabin_of_comb_candidates->m_candidates[ m_start_seed + is ][ h2a.trkIdx ] );
        cv.back().addHitIdx(h2a.hitIdx, 0);
        cv.back().setChi2(h2a.chi2);
      }

      mp_etabin_of_comb_candidates->m_candidates[ m_start_seed + is ].swap(cv);
      cv.clear();
    }
    else
    {
      // XXXX MT This is just to get what we have to the end.
      // Would be more efficient to put the best candidate into some result list.
      // On the other hand, we might find more hits later on but could as well
      // screw up track parameters / chi2.
      for (auto &t : mp_etabin_of_comb_candidates->m_candidates[ m_start_seed + is ])
      {
        t.addHitIdx(-1, 0);
      }
    }
  }
}

void CandCloner::DoWorkInSideThread(CandClonerWork_t work)
{
  int beg     = work.first;
  int the_end = work.second;

  // printf("CandCloner::DoWorkInSideThread working on beg=%d to the_end=%d\n", beg, the_end);

  while (beg != the_end)
  {
    int end = std::min(beg + s_max_seed_range, the_end);

    // printf("CandCloner::DoWorkInSideThread processing %4d -> %4d\n", beg, end);

    ProcessSeedRange(beg, end);

    beg = end;
  }
}
