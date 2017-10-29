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
#include "TLine.h"

#include <iostream>
#include <vector>

enum ArchEnum {SNB, KNC, KNL};

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

struct ArchOpts
{
  Int_t vumin;
  Int_t vumax;

  Int_t thmin;
  Int_t thmax;

  Double_t vutimemin;
  Double_t vutimemax;

  Double_t thtimemin;
  Double_t thtimemax;

  Double_t vuspeedupmin;
  Double_t vuspeedupmax;

  Double_t thspeedupmin;
  Double_t thspeedupmax;

  Double_t thmeiftimemin;
  Double_t thmeiftimemax;

  Double_t thmeifspeedupmin;
  Double_t thmeifspeedupmax;
};

namespace 
{
  ArchOpts arch_opt;
  void setupArch(ArchEnum ARCH)
  {
    if      (ARCH == SNB)
    {
      arch_opt.vumin = 1;
      arch_opt.vumax = 8;

      arch_opt.thmin = 1;
      arch_opt.thmax = 24;
      
      arch_opt.vutimemin = 0.;
      arch_opt.vutimemax = 1.;

      arch_opt.thtimemin = 0.005;
      arch_opt.thtimemax = 1.;

      arch_opt.vuspeedupmin = 0.;
      arch_opt.vuspeedupmax = arch_opt.vumax;

      arch_opt.thspeedupmin = 0.;
      arch_opt.thspeedupmax = arch_opt.thmax;

      arch_opt.thmeiftimemin = 0.02;
      arch_opt.thmeiftimemax = 0.5;

      arch_opt.thmeifspeedupmin = 0.;
      arch_opt.thmeifspeedupmax = arch_opt.thmax;
    }
    else if (ARCH == KNC)
    {
      arch_opt.vumin = 1;
      arch_opt.vumax = 16;

      arch_opt.thmin = 1;
      arch_opt.thmax = 240;

      arch_opt.vutimemin = 0.;
      arch_opt.vutimemax = 12.;

      arch_opt.thtimemin = 0.01;
      arch_opt.thtimemax = 2.;

      arch_opt.vuspeedupmin = 0.;
      arch_opt.vuspeedupmax = arch_opt.vumax;

      arch_opt.thspeedupmin = 0.;
      arch_opt.thspeedupmax = 40.;

      arch_opt.thmeiftimemin = 0.1;
      arch_opt.thmeiftimemax = arch_opt.thtimemax;

      arch_opt.thmeifspeedupmin = 0.;
      arch_opt.thmeifspeedupmax = 40;
    }
    else if (ARCH == KNL)
    {
      arch_opt.vumin = 1;
      arch_opt.vumax = 16;
      
      arch_opt.thmin = 1;
      arch_opt.thmax = 256;

      arch_opt.vutimemin = 0.;
      arch_opt.vutimemax = 2.;

      arch_opt.thtimemin = 0.005;
      arch_opt.thtimemax = 1.;

      arch_opt.vuspeedupmin = 0.;
      arch_opt.vuspeedupmax = arch_opt.vumax;

      arch_opt.thspeedupmin = 0.;
      arch_opt.thspeedupmax = 40.;

      arch_opt.thmeiftimemin = 0.03;
      arch_opt.thmeiftimemax = arch_opt.thtimemax;

      arch_opt.thmeifspeedupmin = 0.;
      arch_opt.thmeifspeedupmax = 40.;
    }
  }
};
