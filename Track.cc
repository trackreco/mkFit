#include "Track.h"

// find the simtrack that provided the most hits
std::pair<unsigned int,unsigned int> Track::SimTrackIDInfo() const
{
  std::vector<unsigned int> mctrack;
  for (auto&& ihit : hits_){
    mctrack.push_back(ihit.mcTrackID());
  }
  std::sort(mctrack.begin(), mctrack.end()); // ensures all elements are checked properly

  auto mtrk(mctrack[0]), mcount(0U), m(mctrack[0]), c(0U);

  for (auto i : mctrack) {
    if (i == m) { ++c; } else { c = 0; }
    if (c >= mcount) { mtrk = m; mcount = c; }
    m = i;
  }

  std::pair<unsigned int, unsigned int> simIDInfo (mtrk,mcount);
  return simIDInfo;
}

void Track::write_out(FILE *fp)
{
  Track t = clone_for_io();
  fwrite(&t, sizeof(Track), 1, fp);

  int nh = nHits();
  fwrite(&nh, sizeof(int), 1, fp);

  fwrite(&hits_[0], sizeof(Hit), nh, fp);
}

void Track::read_in(FILE *fp)
{
  fread(this, sizeof(Track), 1, fp);

  int nh = nHits();
  fread(&nh, sizeof(int), 1, fp);

  hits_.resize(nh);
  fread(&hits_[0], sizeof(Hit), nh, fp);
}
