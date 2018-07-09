#include "StackValidation.hh"

StackValidation::StackValidation(const TString & label, const TString & extra, const Bool_t cmsswComp) : label(label), extra(extra), cmsswComp(cmsswComp)
{
  gStyle->SetOptStat(0);

  // setup builds --> if simvalidation, add cmssw tracks to validation
  setupBuilds(!cmsswComp);

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
  // kinematic variables to plot
  std::vector<TString> vars = {"pt","eta","phi"};

  // indices for loops match PlotValidation.cpp
  for (UInt_t l = 0; l < rates.size(); l++)
  {
    for (UInt_t k = 0; k < ptcuts.size(); k++)
    {
      for (UInt_t i = 0; i < vars.size(); i++)
      {
	TCanvas * canv = new TCanvas();
	canv->cd();

	TLegend * leg = new TLegend(0.85,0.80,1.0,1.0);

	// tmp axis titles, not sure why ROOT is deleting them
	TString xtitle = "";
	TString ytitle = "";
	
	std::vector<TGraphAsymmErrors*> graphs(builds.size());
	for (UInt_t b = 0; b < builds.size(); b++)
	{
	  auto & graph = graphs[b];

	  graph = ((TEfficiency*)files[b]->Get(rates[l].dir+"/"+rates[l].rate+"_"+rates[l].sORr+"_"+vars[i]+"_"+trk+"_pt"+Form("%3.1f",ptcuts[k])))->CreateGraph();
	  graph->SetLineColor(builds[b].color);
	  graph->SetMarkerColor(builds[b].color);

	  // store tmp titles
	  if (b == 0)
	  {
	    xtitle = graph->GetXaxis()->GetTitle();
	    ytitle = graph->GetYaxis()->GetTitle();
	  }
	  
	  graph->Draw(b>0?"PZ SAME":"APZ");

	  if (!rates[l].rate.Contains("ineff",TString::kExact)) graph->GetYaxis()->SetRangeUser(0.0,1.05);
	  else graph->GetYaxis()->SetRangeUser(0.0,0.25);
	  
	  leg->AddEntry(graph,builds[b].label.Data(),"LEP");
	}

	// print standard plot for every rate/variable
	leg->Draw("SAME");
	canv->SaveAs(label+"_"+rates[l].rate+"_"+vars[i]+"_"+trk+"_pt"+Form("%3.1f",ptcuts[k])+extra+".png");
	
	// zoom in on pt range
	if (i == 0)
	{
	  std::vector<TGraphAsymmErrors*> zoomgraphs(builds.size());
	  for (UInt_t b = 0; b < builds.size(); b++)
	  {
	    zoomgraphs[b] = (TGraphAsymmErrors*)graphs[b]->Clone(Form("%s_zoom",graphs[b]->GetName()));
	    zoomgraphs[b]->GetXaxis()->SetRangeUser(0,10);
	    zoomgraphs[b]->Draw(b>0?"PZ SAME":"APZ");
	  }

	  leg->Draw("SAME");
	  canv->SaveAs(label+"_"+rates[l].rate+"_"+vars[i]+"_zoom_"+trk+"_pt"+Form("%3.1f",ptcuts[k])+extra+".png");

	  for (auto & zoomgraph : zoomgraphs) delete zoomgraph; 
	}

	// make logx plots for pt: causes LOTS of weird effects... workarounds for now
	if (i == 0)
	{
	  canv->SetLogx(1);

	  // apparently logx removes titles and ranges???
	  for (UInt_t b = 0; b < builds.size(); b++)
	  {
	    auto & graph = graphs[b];
	    graph->GetXaxis()->SetRangeUser(0.01,graph->GetXaxis()->GetBinUpEdge(graph->GetXaxis()->GetNbins()));

	    if (!rates[l].rate.Contains("ineff",TString::kExact)) graph->GetYaxis()->SetRangeUser(0.0,1.05);
	    else graph->GetYaxis()->SetRangeUser(0.0,0.25);

	    graph->GetXaxis()->SetTitle(xtitle);
	    graph->GetYaxis()->SetTitle(ytitle);

	    graph->Draw(b>0?"PZ SAME":"APZ");
	  }

	  leg->Draw("SAME");
	  canv->SaveAs(label+"_"+rates[l].rate+"_"+vars[i]+"_logx_"+trk+"_pt"+Form("%3.1f",ptcuts[k])+extra+".png");
	}

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
  std::vector<TString> diffs = {"nHits","invpt","eta","phi"};

  // diffferent reco collections
  std::vector<TString> coll = {"allmatch","bestmatch"};

  // indices for loops match PlotValidation.cpp
  for (UInt_t o = 0; o < coll.size(); o++)
  {
    for (UInt_t p = 0; p < diffs.size(); p++)
    {
      for (UInt_t k = 0; k < ptcuts.size(); k++)
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
	  hists[b] = (TH1F*)files[b]->Get("kindiffs_cmssw/h_d"+diffs[p]+"_"+coll[o]+"_"+trk+"_pt"+Form("%3.1f",ptcuts[k]));
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
	  
	  const TString mean = Form("%4.1f",hists[b]->GetMean());
	  leg->AddEntry(hists[b],builds[b].label+" "+" [#mu = "+mean+"]","LEP");
	}
	
	leg->Draw("SAME");
	canv->SaveAs(label+"_"+coll[o]+"_d"+diffs[p]+"_"+trk+"_pt"+Form("%3.1f",ptcuts[k])+extra+".png");
	
	delete leg;
	for (auto & hist : hists) delete hist;
	delete canv;
      } // end pt cut loop
    } // end var loop
  } // end coll loop
}
