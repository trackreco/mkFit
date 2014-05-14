{

  gROOT->Reset();

  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1011);

  bool doFit = 1;
  bool doBuild = 0;

  if (doFit) {

    TFile *_file0 = TFile::Open("validationtree.root");
    
    TCanvas c1;
    TH1F* pt_res = new TH1F("pt_res","pt resolution",100,-0.5,0.5);
    pt_res->GetXaxis()->SetTitle("(p_{T}^{MC} - p_{T}^{fit})/p_{T}^{MC}");
    tree->Draw("(pt_mc-pt_fit)/pt_mc>>pt_res");
    pt_res->Fit("gaus");
    c1.SaveAs("pt_res.png");  
    
    TCanvas c2;
    TH1F* pt_pul = new TH1F("pt_pul","pt pull",100,-10,10);
    pt_pul->GetXaxis()->SetTitle("(p_{T}^{MC} - p_{T}^{fit})/#sigma(p_{T}^{fit})");
    tree->Draw("(pt_mc-pt_fit)/pt_err>>pt_pul");
    pt_pul->Fit("gaus");
    c2.SaveAs("pt_pull.png");  

  }

  if (doBuild) {

    TFile *_file0 = TFile::Open("build_validationtree.root");

    TCanvas c1;
    TH1F* nhits = new TH1F("nhits","nhits",15,0,15);
    nhits->GetXaxis()->SetTitle("number of hits");
    tree->Draw("nhits>>nhits");
    c1.SaveAs("nhits.png");  

    TCanvas c2;
    TH1F* chi2 = new TH1F("chi2","normalized chi2",100,0,10);
    chi2->GetXaxis()->SetTitle("track #chi^{2}/ndof");
    tree->Draw("chi2/(3*nhits-5)>>chi2");//fixme
    c2.SaveAs("chi2.png");  

    
  }

}
