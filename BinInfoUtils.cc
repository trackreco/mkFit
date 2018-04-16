#include "BinInfoUtils.h"

namespace mkfit {

std::vector<int> getCandHitIndices(const int & etaBinMinus, const int & etaBinPlus, const int & phiBinMinus, const int & phiBinPlus, const BinInfoLayerMap& segLayMap){    
  std::vector<int> cand_hit_indices;
  for (int ieta = etaBinMinus; ieta <= etaBinPlus; ++ieta){
    const BinInfo binInfoMinus = segLayMap[ieta][phiBinMinus];
    const BinInfo binInfoPlus  = segLayMap[ieta][phiBinPlus];
    
    // Branch here from wrapping
    if (phiBinMinus<=phiBinPlus){
      const auto firstIndex = binInfoMinus.first;
      const auto maxIndex   = binInfoPlus.first+binInfoPlus.second;
      for (auto ihit  = firstIndex; ihit < maxIndex; ++ihit){
	cand_hit_indices.push_back(ihit);
      }
    } 
    else { // loop wrap around end of array for phiBinMinus > phiBinPlus
      const auto firstIndex = binInfoMinus.first;
      const auto etaBinSize = segLayMap[ieta][Config::nPhiPart-1].first+segLayMap[ieta][Config::nPhiPart-1].second;
      for (auto ihit  = firstIndex; ihit < etaBinSize; ++ihit){
	cand_hit_indices.push_back(ihit);
      }

      const auto etaBinStart= segLayMap[ieta][0].first;
      const auto maxIndex   = binInfoPlus.first+binInfoPlus.second;
      for (auto ihit  = etaBinStart; ihit < maxIndex; ++ihit){
	cand_hit_indices.push_back(ihit);
      }
    }
  }
  return cand_hit_indices;
}

} // end namespace mkfit
