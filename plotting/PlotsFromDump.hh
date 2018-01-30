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
    tests.push_back({"KNL","NVU1_NTH1",kGreen+1,kOpenTriangleUp});
    tests.push_back({"KNL","NVU16int_NTH256",kOrange+1,kOpenTriangleDown});
  }
};

struct PlotOpts
{
  TString name;
  TString xtitle;
  TString ytitle;
  TString outname;
};
typedef std::vector<PlotOpts> POVec;

namespace 
{
  POVec plots;
  void setupPlots()
  {
    plots.push_back({"h_MXNH","Number of Hits Found","Fraction of Tracks","nHits"});
    plots.push_back({"h_MXPT","p_{T}^{mkFit}","Fraction of Tracks","pt"});
    plots.push_back({"h_MXPHI","#phi^{mkFit}","Fraction of Tracks","phi"});
    plots.push_back({"h_MXETA","#eta^{mkFit}","Fraction of Tracks","eta"});

    plots.push_back({"h_DCNH","nHits^{mkFit}-nHits^{CMSSW}","Fraction of Tracks","dnHits"});
    plots.push_back({"h_DCPT","p_{T}^{mkFit}-p_{T}^{CMSSW}","Fraction of Tracks","dpt"});
    plots.push_back({"h_DCPHI","#phi^{mkFit}-#phi^{CMSSW}","Fraction of Tracks","dphi"});
    plots.push_back({"h_DCETA","#eta^{mkFit}-#eta^{CMSSW}","Fraction of Tracks","deta"});
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
