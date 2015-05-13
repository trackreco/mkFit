#include "BinInfoUtils.h"

hitIndices getCandHitIndices(const float eta, const float etaWindow, const float phi, const float phiWindow, const unsigned int ilayer, const BinInfoMap& segmentMap){    
  hitIndices cand_hit_idx;
#ifdef ETASEG
  const unsigned int etaBinMinus = getEtaPartition(normalizedEta(eta-etaWindow));
  const unsigned int etaBinPlus  = getEtaPartition(normalizedEta(eta+etaWindow));
#else
  const auto etaBinMinus = 0U;
  const auto etaBinPlus  = 0U;
#endif    
  const unsigned int phiBinMinus = getPhiPartition(normalizedPhi(phi-phiWindow));
  const unsigned int phiBinPlus  = getPhiPartition(normalizedPhi(phi+phiWindow));

  for (unsigned int ieta = etaBinMinus; ieta <= etaBinPlus; ++ieta){
    const BinInfo binInfoMinus = segmentMap[ilayer][ieta][int(phiBinMinus)];
    const BinInfo binInfoPlus  = segmentMap[ilayer][ieta][int(phiBinPlus)];
    
    // Branch here from wrapping
    if (phiBinMinus<=phiBinPlus){
      const auto firstIndex = binInfoMinus.first;
      const auto maxIndex   = binInfoPlus.first+binInfoPlus.second;
      for (auto ihit  = firstIndex; ihit < maxIndex; ++ihit){
	cand_hit_idx.push_back(ihit);
      }
    } 
    else { // loop wrap around end of array for phiBinMinus > phiBinPlus
      const auto firstIndex = binInfoMinus.first;
      const auto etaBinSize = segmentMap[ilayer][ieta][62].first+segmentMap[ilayer][ieta][62].second;
      for (auto ihit  = firstIndex; ihit < etaBinSize; ++ihit){
	cand_hit_idx.push_back(ihit);
      }

      const auto etaBinStart= segmentMap[ilayer][ieta][0].first;
      const auto maxIndex   = binInfoPlus.first+binInfoPlus.second;
      for (auto ihit  = etaBinStart; ihit < maxIndex; ++ihit){
	cand_hit_idx.push_back(ihit);
      }
    }
  }
  return cand_hit_idx;
}
