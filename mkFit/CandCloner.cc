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
    //sort the damn thing
    std::sort(hitsToAddForThisSeed.begin(), hitsToAddForThisSeed.end(), sortCandListByHitsChi2);
  }

  //2) now create the candidates for the best maxCandsPerSeed, we'll do it vectorized with MkFitter
  //2a) first we need to unroll the vectors

  // XXXX std::vector<std::pair<int,MkFitter::IdxChi2List> > seed_newcand_idx;

  for (int is = is_beg; is < is_end; ++is)
  {
    std::vector<MkFitter::IdxChi2List>& hitsToAddForThisSeed = m_hits_to_add[is];
    for (int ih = 0; ih < hitsToAddForThisSeed.size() && ih < Config::maxCandsPerSeed; ih++)
    {
      t_seed_newcand_idx.push_back(std::pair<int,MkFitter::IdxChi2List>(is, hitsToAddForThisSeed[ih]));
    }
  }

  int theEndNewCand = t_seed_newcand_idx.size();     

  //2b) vectorized loop
  for (int itrack = 0; itrack < theEndNewCand; itrack += NN)
  {
    int end = std::min(itrack + NN, theEndNewCand);

    ////m_fitter->SetNhits(ilay);//here again assuming one hit per layer

    if (Config::g_PropagateAtEnd) 
    {
      m_fitter->InputTracksAndHitIdx(mp_etabin_of_comb_candidates->m_candidates, t_seed_newcand_idx, itrack, end, true);
      
      //now we need to update with the hit in bunch_of_hits.m_hits[ hitsToAddForThisSeed[ih].hitIdx ]
      //also fill the temp vector of candidates
      m_fitter->UpdateWithHit(*mp_bunch_of_hits, t_seed_newcand_idx, itrack, end);
      
      if (m_layer!=9) 
      {
	//propagate to next layer
#ifdef DEBUG
	MkFitter *mkfp = m_fitter;
	int ilay = m_layer;
	
	std::cout << "propagate to lay=" << ilay+2 << " start from x=" << mkfp->getPar(0, 0, 0) << " y=" << mkfp->getPar(0, 0, 1) << " z=" << mkfp->getPar(0, 0, 2)<< " r=" << getHypot(mkfp->getPar(0, 0, 0), mkfp->getPar(0, 0, 1))
		  << " px=" << mkfp->getPar(0, 0, 3) << " py=" << mkfp->getPar(0, 0, 4) << " pz=" << mkfp->getPar(0, 0, 5) << " pT=" << getHypot(mkfp->getPar(0, 0, 3), mkfp->getPar(0, 0, 4)) << std::endl;
#endif
	//fixme: doesn't need itrack, end?
	// MT 2015-08-13: YES, shows up in valgrind ... added as N_proc argument,
	//    N-to-process, as we always run from mplex entry 0 to maximum.
	//    Maybe also needs to be done in FindCandidates, will investigate later.
	m_fitter->PropagateTracksToR(m_rad, end - itrack);
	  
	m_fitter->CopyOutClone(t_seed_newcand_idx, t_cands_for_next_lay,
			       m_start_seed + is_beg, itrack, end, true);
      }
      else
      {
	m_fitter->CopyOutClone(t_seed_newcand_idx, t_cands_for_next_lay,
			       m_start_seed + is_beg, itrack, end);	
      }

    } 
    else 
    {
      m_fitter->InputTracksAndHitIdx(mp_etabin_of_comb_candidates->m_candidates, t_seed_newcand_idx, itrack, end);

      //propagate to layer
#ifdef DEBUG
      MkFitter *mkfp = m_fitter;
      int ilay = m_layer;
      
      std::cout << "propagate to lay=" << ilay+1 << " start from x=" << mkfp->getPar(0, 0, 0) << " y=" << mkfp->getPar(0, 0, 1) << " z=" << mkfp->getPar(0, 0, 2)<< " r=" << getHypot(mkfp->getPar(0, 0, 0), mkfp->getPar(0, 0, 1))
		<< " px=" << mkfp->getPar(0, 0, 3) << " py=" << mkfp->getPar(0, 0, 4) << " pz=" << mkfp->getPar(0, 0, 5) << " pT=" << getHypot(mkfp->getPar(0, 0, 3), mkfp->getPar(0, 0, 4)) << std::endl;
#endif
      //fixme: doesn't need itrack, end?
      // MT 2015-08-13: YES, shows up in valgrind ... added as N_proc argument,
      //    N-to-process, as we always run from mplex entry 0 to maximum.
      //    Maybe also needs to be done in FindCandidates, will investigate later.
      m_fitter->PropagateTracksToR(m_rad, end - itrack);
      
      //now we need to update with the hit in bunch_of_hits.m_hits[ hitsToAddForThisSeed[ih].hitIdx ]
      //also fill the temp vector of candidates
      m_fitter->UpdateWithHit(*mp_bunch_of_hits, t_seed_newcand_idx, t_cands_for_next_lay,
			      m_start_seed + is_beg, itrack, end);
      
    }

  }

  //3) now swap the temp vector with the input candidates so that they are used as input for the next layer
  for (int is_rel = 0; is_rel < is_num; ++is_rel)
  {
    if (t_cands_for_next_lay[is_rel].size() > 0)
    {
      mp_etabin_of_comb_candidates->m_candidates[m_start_seed + is_beg + is_rel].swap(t_cands_for_next_lay[is_rel]);
      t_cands_for_next_lay[is_rel].clear();
    }
    else 
    {
      //we do nothing in the SM version here, I think we should put these in the output and avoid keeping looping over them
    }
  }

  // Cleanup member temporaries, t_cands_for_next_lay elements are already cleared after each swap.
  
  t_seed_newcand_idx.clear();

  //----------------------- END CLONE ENGINE -----------------------//
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
