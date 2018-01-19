#include "plotting/StackValidation.hh"

StackValidation::StackValidation(const TString & label, const TString & extra, const Bool_t cmsswComp) : label(label), extra(extra), cmsswComp(cmsswComp)
{
  gStyle->SetOptStat(0);

  // setup builds
  setupBuilds();

  files.resize(builds.size());
  for (UInt_t b = 0; b < files.size(); b++)
  {
    files[b] = TFile::Open("validation_"+label+"_"+builds[b].name+extra+"/plots.root");
  }

  // setup rates
  setupRates(cmsswComp);

  // setup ptcuts
  setupPtCuts();
}

StackValidation::~StackValidation()
{
  for (auto & file : files) delete file;
}

void StackValidation::MakeValidationStacks()
{
  StackValidation::MakeRatioStacks("build");

  if (cmsswComp)
  {    
    StackValidation::MakeRatioStacks("fit");
    StackValidation::MakeCMSSWKinematicDiffStacks("build");
    StackValidation::MakeCMSSWKinematicDiffStacks("fit");
  }
}

void StackValidation::MakeRatioStacks(const TString & trk)
{
  // Variables to plot
  std::vector<TString> vars = {"pt","phi","eta"};

  for (UInt_t i = 0; i < rates.size(); i++)
  {
    for (UInt_t j = 0; j < vars.size(); j++)
    {
      for (UInt_t p = 0; p < ptcuts.size(); p++)
      {
	TCanvas * canv = new TCanvas();
	canv->cd();

	TLegend * leg = new TLegend(0.85,0.80,1.0,1.0);
	
	std::vector<TGraphAsymmErrors*> graphs(builds.size());
	for (UInt_t b = 0; b < builds.size(); b++)
	{
	  graphs[b] = ((TEfficiency*)files[b]->Get(rates[i].dir+"/"+rates[i].rate+"_"+rates[i].sORr+"_"+vars[j]+"_"+trk+"_pt"+Form("%3.1f",ptcuts[p])))->CreateGraph();
	  graphs[b]->SetLineColor(builds[b].color);
	  graphs[b]->SetMarkerColor(builds[b].color);

	  graphs[b]->Draw(b>0?"PZ SAME":"APZ");

	  if (!rates[i].rate.Contains("ineff",TString::kExact)) graphs[b]->GetYaxis()->SetRangeUser(0.0,1.05);
	  else graphs[b]->GetYaxis()->SetRangeUser(0.0,0.25);
	  
	  leg->AddEntry(graphs[b],builds[b].label.Data(),"LEP");
	}
	
	leg->Draw("SAME");
	canv->SaveAs(label+"_"+rates[i].rate+"_"+vars[j]+"_"+trk+"_pt"+Form("%3.1f",ptcuts[p])+extra+".png");
	
	delete leg;
	for (auto & graph : graphs) delete graph;
	delete canv;
      }
    }
  }
}

void StackValidation::MakeCMSSWKinematicDiffStacks(const TString & trk)
{
  // variables to plot
  std::vector<TString> diffs = {"nHits","invpt","phi","eta"};

  // diffferent reco collections
  std::vector<TString> coll = {"allmatch","bestmatch"};

  for (UInt_t c = 0; c < coll.size(); c++)
  {
    for (UInt_t d = 0; d < diffs.size(); d++)
    {
      for (UInt_t p = 0; p < ptcuts.size(); p++)
      {
	const Bool_t isLogy = true;
	TCanvas * canv = new TCanvas();
	canv->cd();
	canv->SetLogy(isLogy);

	TLegend * leg = new TLegend(0.85,0.80,1.0,1.0);
	
	// tmp min/max
	Double_t min =  1e9;
	Double_t max = -1e9;

	std::vector<TH1F*> hists(builds.size());
	for (UInt_t b = 0; b < builds.size(); b++)
        {
	  hists[b] = (TH1F*)files[b]->Get("kindiffs_cmssw/h_d"+diffs[d]+"_"+coll[c]+"_"+trk+"_pt"+Form("%3.1f",ptcuts[p]));
	  hists[b]->SetLineColor(builds[b].color);
	  hists[b]->SetMarkerColor(builds[b].color);

	  hists[b]->Scale(1.f/hists[b]->Integral());
	  hists[b]->GetYaxis()->SetTitle("Fraction of Tracks");
	  
	  GetMinMaxHist(hists[b],min,max);
	}

	for (UInt_t b = 0; b < builds.size(); b++)
	{
	  SetMinMaxHist(hists[b],min,max,isLogy);
	  
	  hists[b]->Draw(b>0?"EP SAME":"EP");
	  leg->AddEntry(hists[b],builds[b].label.Data(),"LEP");
	}
	
	leg->Draw("SAME");
	canv->SaveAs(label+"_"+coll[c]+"_d"+diffs[d]+"_"+trk+"_pt"+Form("%3.1f",ptcuts[p])+extra+".png");
	
	delete leg;
	for (auto & hist : hists) delete hist;
	delete canv;
      } // end pt cut loop
    } // end var loop
  } // end coll loop
}
