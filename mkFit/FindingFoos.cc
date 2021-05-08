#include "FindingFoos.h"
#include "MkBase.h"
#include "KalmanUtilsMPlex.h"

namespace {
  using namespace mkfit;
  const FindingFoos  s_fndfoos_brl( kalmanPropagateAndComputeChi2,       kalmanPropagateAndUpdate,       &MkBase::PropagateTracksToR );
  const FindingFoos  s_fndfoos_ec ( kalmanPropagateAndComputeChi2Endcap, kalmanPropagateAndUpdateEndcap, &MkBase::PropagateTracksToZ );
}

namespace mkfit {

const FindingFoos& FindingFoos::get_barrel_finding_foos() { return s_fndfoos_brl; }
const FindingFoos& FindingFoos::get_endcap_finding_foos() { return s_fndfoos_ec;  }

const FindingFoos& FindingFoos::get_finding_foos(bool is_barrel)
{
    return is_barrel ? s_fndfoos_brl : s_fndfoos_ec;
}

}
