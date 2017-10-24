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
}

StackValidation::~StackValidation()
{
  for (auto & file : files) delete file;
}

void StackValidation::MakeValidationStacks()
{
  std::vector<TString> vars;
  vars.push_back("pt");
  vars.push_back("phi");
  vars.push_back("eta");
  std::vector<Float_t> ptcuts;
  ptcuts.push_back(0.f);
  ptcuts.push_back(0.9f);
  ptcuts.push_back(2.f);

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
	  graphs[b] = ((TEfficiency*)files[b]->Get(rates[i].dir+"/"+rates[i].rate+"_"+rates[i].sORr+"_"+vars[j]+"_build_pt"+Form("%3.1f",ptcuts[p])))->CreateGraph();
	  graphs[b]->SetLineColor(builds[b].color);
	  graphs[b]->SetMarkerColor(builds[b].color);

	  graphs[b]->Draw(b>0?"PZ SAME":"ALZ");

	  if (!rates[i].rate.Contains("ineff",TString::kExact)) graphs[b]->GetYaxis()->SetRangeUser(0.0,1.05);
	  else graphs[b]->GetYaxis()->SetRangeUser(0.0,0.25);
	  
	  leg->AddEntry(graphs[b],builds[b].label.Data(),"LEP");
	}
	
	leg->Draw("SAME");
	canv->SaveAs(label+"_"+rates[i].rate+"_"+vars[j]+"_pt"+Form("%3.1f",ptcuts[p])+extra+".png");
	
	delete leg;
	for (auto & graph : graphs) delete graph;
	delete canv;
      }
    }
  }
}
