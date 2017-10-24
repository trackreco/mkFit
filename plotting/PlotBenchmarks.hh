#ifndef _PlotBenchmarks_
#define _PlotBenchmarks_

#include "plotting/Common.hh"

typedef std::vector<TGraphErrors*> TGEVec;

class PlotBenchmarks
{
public:
  PlotBenchmarks(const TString & arch, const TString & sample);
  ~PlotBenchmarks();
  void RunBenchmarkPlots();
  void MakeOverlay(const TString & text, const TString & title, const TString & xtitle,
		   const TString & ytitle, const Float_t xmax, const Float_t ymin, const Float_t ymax);
  void GetGraphs(TGEVec & graphs, const TString & text, const TString & title, const TString & xtitle,
		 const TString & ytitle, const Float_t xmax, const Float_t ymin, const Float_t ymax);

private:
  const TString arch;
  const TString sample;

  ArchEnum ARCH;
  TFile * file;
};

#endif
