#ifndef _PlotsFromDump_
#define _PlotsFromDump_

#include "Common.hh"

struct TestOpts
{
  TestOpts() {}
  TestOpts(const TString & arch, const TString & suffix, const Color_t color, const Marker_t marker)
    : arch(arch), suffix(suffix), color(color), marker(marker) {}

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
    // N.B.: Consult ./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build.sh for info on which VU and TH tests were used for making text dumps

    tests.emplace_back("SNB","NVU1_NTH1",kBlue,kOpenDiamond);
    tests.emplace_back("SNB","NVU8int_NTH24",kBlack,kOpenCross);
    tests.emplace_back("KNL","NVU1_NTH1",kGreen+1,kOpenTriangleUp);
    tests.emplace_back("KNL","NVU16int_NTH256",kOrange+1,kOpenTriangleDown);
    tests.emplace_back("SKL-SP","NVU1_NTH1",kRed+1,kOpenTriangleUp);
    tests.emplace_back("SKL-SP","NVU16int_NTH64",kMagenta+1,kOpenTriangleDown);
  }
};

struct PlotOpts
{
  PlotOpts() {}
  PlotOpts(const TString & name, const TString & xtitle, const TString & ytitle, const TString & outname)
    : name(name), xtitle(xtitle), ytitle(ytitle), outname(outname) {}

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
    // N.B. Consult plotting/makePlotsFromDump.py for info on hist names

    plots.emplace_back("h_MXNH","Number of Hits Found","Fraction of Tracks","nHits");
    plots.emplace_back("h_MXPT","p_{T}^{mkFit}","Fraction of Tracks","pt");
    plots.emplace_back("h_MXPHI","#phi^{mkFit}","Fraction of Tracks","phi");
    plots.emplace_back("h_MXETA","#eta^{mkFit}","Fraction of Tracks","eta");

    plots.emplace_back("h_DCNH","nHits^{mkFit}-nHits^{CMSSW}","Fraction of Tracks","dnHits");
    plots.emplace_back("h_DCPT","p_{T}^{mkFit}-p_{T}^{CMSSW}","Fraction of Tracks","dpt");
    plots.emplace_back("h_DCPHI","#phi^{mkFit}-#phi^{CMSSW}","Fraction of Tracks","dphi");
    plots.emplace_back("h_DCETA","#eta^{mkFit}-#eta^{CMSSW}","Fraction of Tracks","deta");
  }
};

class PlotsFromDump
{
public:
  PlotsFromDump(const TString & sample, const TString & build, const TString & suite);
  ~PlotsFromDump();
  void RunPlotsFromDump();

private:
  const TString sample;
  const TString build;
  const TString suite;

  TString label;
};

#endif
