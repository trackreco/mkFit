{

  gROOT->Reset();

  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1011);

  bool doFit = 1;
  bool doBuild = 1;
  bool doSim = 1;

  if (doFit) {

    TFile *_file0 = TFile::Open("validationtree.root");
    
    TCanvas cf1;
    TH1F* pt_res = new TH1F("pt_res","pt resolution",100,-0.5,0.5);
    pt_res->GetXaxis()->SetTitle("(p_{T}^{MC} - p_{T}^{fit})/p_{T}^{MC}");
    tree->Draw("(pt_mc-pt_fit)/pt_mc>>pt_res");
    pt_res->Fit("gaus","","",-0.07,0.07);
    cf1.SaveAs("pt_res.png");  
    
    TCanvas cf2;
    TH1F* pt_pul = new TH1F("pt_pul","pt pull",100,-10,10);
    pt_pul->GetXaxis()->SetTitle("(p_{T}^{MC} - p_{T}^{fit})/#sigma(p_{T}^{fit})");
    tree->Draw("(pt_mc-pt_fit)/pt_err>>pt_pul");
    pt_pul->Fit("gaus","","",-3,3);
    cf2.SaveAs("pt_pull.png");  
    
    //TCanvas cf3;
    //TH1F* ipt_res = new TH1F("ipt_res","inverse pt resolution",100,-0.5,0.5);
    //ipt_res->GetXaxis()->SetTitle("(1./p_{T}^{MC} - 1./p_{T}^{fit})/(1./p_{T}^{MC})");
    //tree->Draw("(1./pt_mc-1./pt_fit)/(1./pt_mc)>>ipt_res");
    //ipt_res->Fit("gaus");
    //cf1.SaveAs("ipt_res.png");  

  }

  if (doBuild) {

    TFile *_file0 = TFile::Open("build_validationtree.root");

    TCanvas cb1;
    TH1F* nhits = new TH1F("nhits","nhits",15,0,15);
    nhits->GetXaxis()->SetTitle("number of hits");
    tree->Draw("nhits>>nhits");
    cb1.SaveAs("nhits.png");  

    TCanvas cb2;
    TH1F* chi2 = new TH1F("chi2","normalized chi2",100,0,10);
    chi2->GetXaxis()->SetTitle("track #chi^{2}/ndof");
    tree->Draw("chi2/(3*(nhits-3)-5)>>chi2");//fixme //nhits-3 not to count seed which is fake at this stage
    cb2.SaveAs("chi2.png");  
    
  }

  if (doSim) {

    TFile *_file0 = TFile::Open("build_validationtree.root");

    TH1F* h_gen_trk_Pt = (TH1F*) _file0->Get("h_gen_trk_Pt");	  
    TH1F* h_gen_trk_Px = (TH1F*) _file0->Get("h_gen_trk_Px");	  
    TH1F* h_gen_trk_Py = (TH1F*) _file0->Get("h_gen_trk_Py");	  
    TH1F* h_gen_trk_Pz = (TH1F*) _file0->Get("h_gen_trk_Pz");	  
    TH1F* h_gen_trk_dPhi = (TH1F*) _file0->Get("h_gen_trk_dPhi");	  
    TH1F* h_gen_trk_dR = (TH1F*) _file0->Get("h_gen_trk_dR");	  
    TH1F* h_gen_trk_eta = (TH1F*) _file0->Get("h_gen_trk_eta");	  
    TH1F* h_gen_trk_mindPhi = (TH1F*) _file0->Get("h_gen_trk_mindPhi");
    TH1F* h_gen_trk_mindR = (TH1F*) _file0->Get("h_gen_trk_mindR");  
    TH1F* h_gen_trk_phi = (TH1F*) _file0->Get("h_gen_trk_phi");    

    TH1F* h_gen_trk_Pt = (TH1F*) _file0->Get("h_gen_trk_Pt");	  
    TCanvas c1;
    //h_gen_trk_Pt->GetYaxis()->SetRangeUser(0.,h_gen_trk_Pt->GetMaximum()*1.1);
    h_gen_trk_Pt->Draw();
    c1.SaveAs("gen_trk_Pt.png");
    TH1F* h_gen_trk_Px = (TH1F*) _file0->Get("h_gen_trk_Px");	  
    TCanvas c2;
    h_gen_trk_Px->Draw();
    c2.SaveAs("gen_trk_Px.png");
    TH1F* h_gen_trk_Py = (TH1F*) _file0->Get("h_gen_trk_Py");	  
    TCanvas c3;
    h_gen_trk_Py->Draw();
    c3.SaveAs("gen_trk_Py.png");
    TH1F* h_gen_trk_Pz = (TH1F*) _file0->Get("h_gen_trk_Pz");	  
    TCanvas c4;
    h_gen_trk_Pz->Draw();
    c4.SaveAs("gen_trk_Pz.png");

    TH1F* h_gen_trk_eta = (TH1F*) _file0->Get("h_gen_trk_eta");	  
    TCanvas c5;
    h_gen_trk_eta->Draw();
    c5.SaveAs("gen_trk_eta.png");
    TH1F* h_gen_trk_phi = (TH1F*) _file0->Get("h_gen_trk_phi");	  
    TCanvas c6;
    h_gen_trk_phi->Draw();
    c6.SaveAs("gen_trk_phi.png");

    TH1F* h_gen_trk_dPhi = (TH1F*) _file0->Get("h_gen_trk_dPhi");	  
    TCanvas c7;
    h_gen_trk_dPhi->Draw();
    c7.SaveAs("gen_trk_dPhi.png");
    TH1F* h_gen_trk_dR = (TH1F*) _file0->Get("h_gen_trk_dR");	  
    TCanvas c8;
    h_gen_trk_dR->Draw();
    c8.SaveAs("gen_trk_dR.png");

    TH1F* h_gen_trk_mindPhi = (TH1F*) _file0->Get("h_gen_trk_mindPhi");	  
    TCanvas c9;
    h_gen_trk_mindPhi->Draw();
    c9.SaveAs("gen_trk_mindPhi.png");
    TH1F* h_gen_trk_mindR = (TH1F*) _file0->Get("h_gen_trk_mindR");	  
    TCanvas c10;
    h_gen_trk_mindR->Draw();
    c10.SaveAs("gen_trk_mindR.png");

  }

}
