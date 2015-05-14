#include "BinInfoUtils.h"

hitIndices getCandHitIndices(const unsigned int & etaBinMinus, const unsigned int & etaBinPlus, const unsigned int & phiBinMinus, const unsigned int & phiBinPlus, const BinInfoLayerMap& segLayMap){    
  hitIndices cand_hit_idx;
  for (unsigned int ieta = etaBinMinus; ieta <= etaBinPlus; ++ieta){
    const BinInfo binInfoMinus = segLayMap[ieta][int(phiBinMinus)];
    const BinInfo binInfoPlus  = segLayMap[ieta][int(phiBinPlus)];
    
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
      const auto etaBinSize = segLayMap[ieta][Config::nPhiPart-1].first+segLayMap[ieta][Config::nPhiPart-1].second;
      for (auto ihit  = firstIndex; ihit < etaBinSize; ++ihit){
	cand_hit_idx.push_back(ihit);
      }

      const auto etaBinStart= segLayMap[ieta][0].first;
      const auto maxIndex   = binInfoPlus.first+binInfoPlus.second;
      for (auto ihit  = etaBinStart; ihit < maxIndex; ++ihit){
	cand_hit_idx.push_back(ihit);
      }
    }
  }
  return cand_hit_idx;
}
