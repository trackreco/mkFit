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
#include "TF1.h"

#include <iostream>
#include <vector>

enum ArchEnum {SNB, KNL};

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
    builds.push_back({"FV",kMagenta,"Full Vector"});
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

      arch_opt.thtimemin = 0.001;
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
    else if (ARCH == KNL)
    {
      arch_opt.vumin = 1;
      arch_opt.vumax = 16;
      
      arch_opt.thmin = 1;
      arch_opt.thmax = 256;

      arch_opt.vutimemin = 0.;
      arch_opt.vutimemax = 2.;

      arch_opt.thtimemin = 0.001;
      arch_opt.thtimemax = 1.;

      arch_opt.vuspeedupmin = 0.;
      arch_opt.vuspeedupmax = arch_opt.vumax;

      arch_opt.thspeedupmin = 0.;
      arch_opt.thspeedupmax = 60.;

      arch_opt.thmeiftimemin = 0.01;
      arch_opt.thmeiftimemax = arch_opt.thtimemax;

      arch_opt.thmeifspeedupmin = 0.;
      arch_opt.thmeifspeedupmax = 40.;
    }
  }
};

void GetMinMaxHist(const TH1F * hist, Double_t & min, Double_t & max)
{
  for (Int_t ibin = 1; ibin <= hist->GetNbinsX(); ibin++)
  {
    const Double_t content = hist->GetBinContent(ibin);
    
    if (content < min && content != 0.0) min = content;
    if (content > max) max = content;
  }
}
 
void SetMinMaxHist(TH1F * hist, const Double_t min, const Double_t max, const Bool_t isLogy)
{
  hist->SetMinimum(isLogy?min/2.0:min/1.05);
  hist->SetMaximum(isLogy?max*2.0:max*1.05);
}
