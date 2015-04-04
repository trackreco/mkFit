#include "Hit.h"

// for each layer
//   Config::nEtaBin vectors of hits, resized to large enough N
//   filled with corresponding hits
//   sorted in phi
//   and with a corresponding BinInfoSomething
//
// My gut feeling is that we could have them all in one place and just store
// indices into small enough eta-phi bins (oh, and eta is a horrible separation variable)
// 10*50000*64/1024/1024
//    30.51757812500000000000
// Or make them smaller ... if we could use short indices that might make a difference.
//
// Need to add eta/phi ... or something significantly better to the Hit struct
//
// At least ... fast eta/phi estimators.
//
// Oh, and we should use radix sort. Could one vectorize this?

// Need a good "array of pods" class with aligned alloc and automatic growth.
// For now just implement the no-resize / no-destroy basics in the BoH.

namespace Config
{
  const int g_NEvents           = 10;
  const int g_NTracks           = 20000;
  const int g_MaxHitsPerBunch   = std::max(100, g_NTracks * 2 / Config::nEtaPart);

  const int g_MaxCandsPerSeed   = 6;
  const int g_MaxCandsPerEtaBin = std::max(100, g_MaxCandsPerSeed * g_NTracks / Config::nEtaPart);
  // Effective eta bin is one half of nEtaPart -- so the above is twice the "average".
  // Note that last and first bin are 3/4 nEtaPart ... but can be made 1/4 by adding
  // additional bins on each end.
}

typedef std::pair<int, int> PhiBinInfo_t;

//==============================================================================

inline bool sortHitsByPhiMT(const Hit& h1, const Hit& h2)
{
  return std::atan2(h1.position()[1],h1.position()[0])<std::atan2(h2.position()[1],h2.position()[0]);
}

inline bool sortTrksByPhiMT(const Track& t1, const Track& t2)
{
  return t1.momPhi() < t2.momPhi();
}

//==============================================================================

class BunchOfHits
{
public:
  // This eventually becomes a vector of pointers / refs / multi indices
  //std::vector<Hit>          m_hits;
  Hit                      *m_hits;
  std::vector<PhiBinInfo_t> m_phi_bin_infos;

  int     m_real_size;
  int     m_fill_index;

public:
  BunchOfHits() :
    //m_hits          (Config::g_MaxHitsPerBunch),
    m_phi_bin_infos (Config::nPhiPart),
    m_real_size     (Config::g_MaxHitsPerBunch),
    m_fill_index    (0)
  {
    m_hits = (Hit*) _mm_malloc(sizeof(Hit)*Config::g_MaxHitsPerBunch, 64);
  }

  ~BunchOfHits()
  {
    _mm_free(m_hits);
  }
  
  void Reset()
  {
    for (auto &bi : m_phi_bin_infos)
    {
      bi.first  = -1;
      bi.second =  0;
    }

    m_fill_index = 0;
  }

  void InsertHit(const Hit& hit)
  {
    assert (m_fill_index < m_real_size); // or something

    m_hits[m_fill_index] = hit;
    ++m_fill_index;
  }

  void SortByPhiBuildPhiBins()
  {
    // std::sort(m_hits.begin(), m_hits.begin() + m_fill_index, sortByPhiMT);
    std::sort(&m_hits[0], &m_hits[m_fill_index], sortHitsByPhiMT);

    int last_bin = -1;
    int idx      =  0;
    for (int i = 0; i < m_fill_index; ++i)
    {
      Hit &h = m_hits[i];

      int bin = getPhiPartition(h.phi());

      if (bin != last_bin) 
      {
        m_phi_bin_infos[bin].first  = idx;
        // PhiBinInfo.second set to 0 in Reset()
      }
      ++m_phi_bin_infos[bin].second;

      last_bin = bin;
      ++idx;
    }
  }
};

//==============================================================================

class LayerOfHits
{
public:
  std::vector<BunchOfHits> m_bunches_of_hits;

public:
  LayerOfHits() :
    m_bunches_of_hits(Config::nEtaBin)
  {}

  void Reset()
  {
    for (auto &i : m_bunches_of_hits)
    {
      i.Reset();
    }
  }

  void InsertHit(const Hit& hit)
  {
    int b1, b2;
    int cnt = Config::getBothEtaBins(hit.eta(), b1, b2);

    if (b1 != -1) m_bunches_of_hits[b1].InsertHit(hit);
    if (b2 != -1) m_bunches_of_hits[b2].InsertHit(hit);
  }

  void SortByPhiBuildPhiBins()
  {
    for (auto &i : m_bunches_of_hits)
    {
      i.SortByPhiBuildPhiBins();
    }
  }
};

//==============================================================================

class EventOfHits
{
public:
  std::vector<LayerOfHits> m_layers_of_hits;
  int                      m_n_layers;

public:
  EventOfHits(int n_layers) :
    m_layers_of_hits(n_layers),
    m_n_layers(n_layers)
  {}

  void Reset()
  {
    for (auto &i : m_layers_of_hits)
    {
      i.Reset();
    }
  }

  void InsertHit(const Hit& hit, int layer)
  {
    m_layers_of_hits[layer].InsertHit(hit);
  }

  void SortByPhiBuildPhiBins()
  {
    for (auto &i : m_layers_of_hits)
    {
      i.SortByPhiBuildPhiBins();
    }
  }
};


//==============================================================================
//==============================================================================

// This should actually be a BunchOfCandidates that share common hit vector.
// At the moment this is an EtaBin ...

class EtaBinOfCandidates
{
public:
  std::vector<Track> m_candidates;

  int     m_real_size;
  int     m_fill_index;

public:
  EtaBinOfCandidates() :
    m_candidates (Config::g_MaxCandsPerEtaBin),
    m_real_size  (Config::g_MaxCandsPerEtaBin),
    m_fill_index (0)
  {}

  void Reset()
  {
    m_fill_index = 0;
  }

  void InsertTrack(const Track& track)
  {
    assert (m_fill_index < m_real_size); // or something

    m_candidates[m_fill_index] = track;
    ++m_fill_index;
  }

  void SortByPhi()
  {
    std::sort(m_candidates.begin(), m_candidates.begin() + m_fill_index, sortTrksByPhiMT);
  }
};

class EventOfCandidates
{
public:
  std::vector<EtaBinOfCandidates> m_etabins_of_candidates;

public:
  EventOfCandidates() :
    m_etabins_of_candidates(Config::nEtaBin)
  {}

  void Reset()
  {
    for (auto &i : m_etabins_of_candidates)
    {
      i.Reset();
    }
  }

  void InsertCandidate(const Track& track)
  {
    // XXXX assuming vertex at origin.
    // XXXX the R condition is trying to get rid of bad seeds (as a quick hack)
    int bin = Config::getEtaBin(track.momEta());
    float r = track.posR();
    if (bin != -1 && r > 11.9 && r < 12.1)
    {
      m_etabins_of_candidates[bin].InsertTrack(track);
    }
  }

  void SortByPhi()
  {
    for (auto &i : m_etabins_of_candidates)
    {
      i.SortByPhi();
    }
  }
};



//-------------------------------------------------------
// for combinatorial version, switch to vector of vectors
//-------------------------------------------------------

class EtaBinOfCombCandidates
{
public:
  std::vector<std::vector<Track> > m_candidates;

  //these refer to seeds
  int     m_real_size;
  int     m_fill_index;

public:
  EtaBinOfCombCandidates() :
    m_candidates (Config::g_MaxCandsPerEtaBin),
    m_real_size  (Config::g_MaxCandsPerEtaBin),
    m_fill_index (0)
  {
    for (int s=0;s<m_real_size;++s) m_candidates[s].reserve(Config::g_MaxCandsPerSeed);//we should never exceed this
  }

  void Reset()
  {
    m_fill_index = 0;
  }

  void InsertSeed(const Track& seed)
  {
    assert (m_fill_index < m_real_size); // or something

    m_candidates[m_fill_index].push_back(seed);
    ++m_fill_index;
  }

  void InsertTrack(const Track& track, int seed_index)
  {
    assert (seed_index <= m_fill_index); // or something

    m_candidates[seed_index].push_back(track);
  }

  /* void SortByPhi() */
  /* { */
  /*   std::sort(m_candidates.begin(), m_candidates.begin() + m_fill_index, sortTrksByPhiMT); */
  /* } */
};

class EventOfCombCandidates
{
public:
  std::vector<EtaBinOfCombCandidates> m_etabins_of_comb_candidates;

public:
  EventOfCombCandidates() :
    m_etabins_of_comb_candidates(Config::nEtaBin)
  {}

  void Reset()
  {
    for (auto &i : m_etabins_of_comb_candidates)
    {
      i.Reset();
    }
  }

  void InsertSeed(const Track& seed)
  {
    // XXXX assuming vertex at origin.
    // XXXX the R condition is trying to get rid of bad seeds (as a quick hack)
    int bin = Config::getEtaBin(seed.momEta());
    float r = seed.posR();
    if (bin != -1 && r > 11.9 && r < 12.1)
      {
	m_etabins_of_comb_candidates[bin].InsertSeed(seed);
      } 
#ifdef DEBUG
    else std::cout << "excluding seed with r=" << r << " etaBin=" << bin << std::endl;
#endif
  }

  void InsertCandidate(const Track& track, int seed_index)
  {
    // XXXX assuming vertex at origin.
    // XXXX the R condition is trying to get rid of bad seeds (as a quick hack)
    int bin = Config::getEtaBin(track.momEta());
    m_etabins_of_comb_candidates[bin].InsertTrack(track,seed_index);
  }

};
