#include "TString.h"
#include "TColor.h"
#include "TStyle.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TAxis.h"
#include "TH1F.h"

#include <iostream>
#include <vector>

struct BuildOpts
{
  TString name;
  Color_t color;
  TString label;
};
typedef std::vector<BuildOpts> BOVec;

namespace
{
  BOVec builds;
  void setupBuilds()
  {
    builds.push_back({"BH",kBlue,"Best Hit"});
    builds.push_back({"STD",kGreen+1,"Standard"});
    builds.push_back({"CE",kRed,"Clone Engine"});
  }
};

enum ArchEnum {SNB, KNC, KNL};
