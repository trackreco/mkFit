#include "Hit.h"

std::atomic<int> MCHitInfo::mcHitIDCounter_(0);

void MCHitInfo::reset()
{
  mcHitIDCounter_ = 0;
}
