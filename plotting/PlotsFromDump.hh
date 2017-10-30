#ifndef _PlotsFromDump_
#define _PlotsFromDump_

#include "plotting/Common.hh"

struct TestOpts
{
  TString arch;
  TString suffix;
  Color_t color;
  Marker_t marker;
};
typedef std::vector<TestOpts>  TOVec;

namespace
{
  TOVec tests;
  void setupTests()
  {
    tests.push_back({"SNB","NVU1_NTH1",kBlue,kOpenDiamond});
    tests.push_back({"SNB","NVU8int_NTH24",kBlack,kOpenCross});
    tests.push_back({"KNC","NVU1_NTH1",kRed,kOpenCircle});
    tests.push_back({"KNC","NVU16int_NTH240",kMagenta,kOpenSquare});
    tests.push_back({"KNL","NVU1_NTH1",kGreen+1,kOpenTriangleUp});
    tests.push_back({"KNL","NVU16int_NTH256",kOrange+1,kOpenTriangleDown});
  }
};

class PlotsFromDump
{
public:
  PlotsFromDump(const TString & sample, const TString & build);
  ~PlotsFromDump();
  void RunPlotsFromDump();

private:
  const TString sample;
  const TString build;

  TString label;
};

#endif
