//-------------------
// CMS 2017 geometry
//-------------------

#include "../Config.h"
#include "../TrackerInfo.h"

namespace
{
#include "CMS-2017.acc"
}

void* TrackerInfoCrator_ptr = (void*) Create_TrackerInfo;
