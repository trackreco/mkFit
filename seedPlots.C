void makePlot(const TTree*& tree1, const TTree*& tree2, const TString label1, const TString label2, const TString outdir,
	      const TString plotname, const TString plottitle, const unsigned int nBins, const int bmin, const int bmax);
  
void seedPlots(){
  gStyle->SetOptStat("omenu");

  TString label1 = "ht_curve";
  TString label2 = "ht_approx";

  TString outdir = "comp/v_approx";
  // make output directory
  FileStat_t dummyFileStat;
  if (gSystem->GetPathInfo(outdir.Data(), dummyFileStat) == 1){
    TString mkDir = "mkdir -p ";
    mkDir += outdir.Data();
    gSystem->Exec(mkDir.Data());
  }

  TFile * file1 = TFile::Open(Form("%s.root",label1.Data()));
  TFile * file2 = TFile::Open(Form("%s.root",label2.Data()));

  TTree * tree1 = (TTree*)file1->Get("seedtree");
  TTree * tree2 = (TTree*)file2->Get("seedtree");

  makePlot(tree1,tree2,label1,label2,outdir,"nTkAll","All nTracks/event before filtering",100,6000,12000);
  makePlot(tree1,tree2,label1,label2,outdir,"nTkAllMC","MC nTracks/event before filtering",110,401,511);
  makePlot(tree1,tree2,label1,label2,outdir,"nTkCut","All nTracks/event after filtering",100,1000,5000);
  makePlot(tree1,tree2,label1,label2,outdir,"nTkCutMC","MC nTracks/event after filtering",110,401,511);
}

void makePlot(const TTree*& tree1, const TTree*& tree2, const TString label1, const TString label2, const TString outdir,
	      const TString plotname, const TString plottitle, const unsigned int nBins, const int bmin, const int bmax) {

  TCanvas * canvas = new TCanvas();
  canvas->cd();

  TPad * pad = new TPad("pad","pad",0.,0.,1.,1.);
  pad->Draw();
  pad->cd();

  TH1F * h1 = new TH1F("h1",plottitle.Data(),nBins,bmin,bmax);
  tree1->Draw(Form("%s>>h1",plotname.Data()),"","goff");
  h1->SetName(label1.Data());

  TH1F * h2 = new TH1F("h2",plottitle.Data(),nBins,bmin,bmax);
  tree2->Draw(Form("%s>>h2",plotname.Data()),"","goff");
  h2->SetName(label2.Data());

  //h1->SetMaximum((h1->GetMaximum()>h2->GetMaximum()?h1->GetMaximum*1.1:h2->GetMaximum()*1.1));
  h1->SetLineColor(kRed);
  h1->Draw();
  h2->SetLineColor(kBlue);
  h2->Draw("sames");

  pad->Update();

  TPaveStats * st1 = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
  st1->SetX1NDC(0.32);
  st1->SetX2NDC(0.52);
  st1->SetY1NDC(0.74);
  st1->SetY2NDC(0.94);

  TPaveStats * st2 = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
  st2->SetX1NDC(st1->GetX2NDC() + 0.02);
  st2->SetX2NDC(2.*st1->GetX2NDC() - st1->GetX1NDC() + 0.02);
  st2->SetY1NDC(st1->GetY1NDC());
  st2->SetY2NDC(st1->GetY2NDC());
  
  TLegend * leg = new TLegend(0.12,0.85,0.3,0.94);
  leg->AddEntry(h1, label1.Data(), "L" );
  leg->AddEntry(h2, label2.Data(), "L" );
  leg->Draw("same");

  pad->cd();
  if (plotname.Contains("nTkAll",TString::kExact) && !plotname.Contains("MC")) {
    pad->SetLogx(1);
  }
  else {
    pad->SetLogx(0);
  }
  
  if ((plotname.Contains("nTkAll",TString::kExact) && !plotname.Contains("MC")) || (plotname.Contains("nTkCutMC",TString::kExact))) {
    pad->SetLogy(1);
  }
  else {
    pad->SetLogy(0);
  }
  canvas->SaveAs(Form("%s/%s_%s_v_%s.png",outdir.Data(),plotname.Data(),label1.Data(),label2.Data()));

  delete leg;
  delete st2;
  delete st1;
  delete h2;
  delete h1;
  delete pad;
  delete canvas;
}
