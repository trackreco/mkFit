{

  gROOT->Reset();

  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1011);

  bool doFit = 1;
  bool doConformal = 0;
  bool doBuild = 1;
  bool doSim = 0;

  if (doFit) {

    TFile *_file0 = TFile::Open("validationtree.root");
    
    TCanvas cf1;
    TH1F* pt_res = new TH1F("pt_res","pt resolution",100,-0.5,0.5);
    pt_res->GetXaxis()->SetTitle("(p_{T}^{MC} - p_{T}^{fit})/p_{T}^{MC}");
    //tree->Draw("(1./pt_mc-1./pt_fit)/(1./pt_mc)>>pt_res");
    TTree * tree   = (TTree*) _file0->Get("ptTree");
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

    if (doConformal) {
      TCanvas ccf1;
      TH1F* x_res = new TH1F("x_res","x residual",100,-0.05,0.05);
      x_res->GetXaxis()->SetTitle("x^{cfit} - x^{MC}");
      tree->Draw("cfitHit0_x-simHit0_x>>x_res");
      x_res->Fit("gaus","","",-0.02,0.02);
      ccf1.SaveAs("cfit_x_res.png");  
      TCanvas ccf2;
      TH1F* y_res = new TH1F("y_res","y residual",100,-0.05,0.05);
      y_res->GetXaxis()->SetTitle("y^{cfit} - y^{MC}");
      tree->Draw("cfitHit0_y-simHit0_y>>y_res");
      y_res->Fit("gaus","","",-0.02,0.02);
      ccf2.SaveAs("cfit_y_res.png");  
      TCanvas ccf3;
      TH1F* z_res = new TH1F("z_res","z residual",100,-0.5,0.5);
      z_res->GetXaxis()->SetTitle("z^{cfit} - z^{MC}");
      tree->Draw("cfitHit0_z-simHit0_z>>z_res");
      z_res->Fit("gaus","","",-0.3,0.3);
      ccf3.SaveAs("cfit_z_res.png");  

      TCanvas ccf4;
      TH1F* px_res = new TH1F("px_res","px resolution",100,-5,5);
      px_res->GetXaxis()->SetTitle("px^{cfit} - px^{MC} / px^{MC}");
      tree->Draw("(cfitHit0_px-simHit0_px)/simHit0_px>>px_res");
      px_res->Fit("gaus","","",-1,1);
      ccf4.SaveAs("cfit_px_res.png");  
      TCanvas ccf5;
      TH1F* py_res = new TH1F("py_res","py resolution",100,-5,5);
      py_res->GetXaxis()->SetTitle("py^{cfit} - py^{MC} / py^{MC}");
      tree->Draw("(cfitHit0_py-simHit0_py)/simHit0_py>>py_res");
      py_res->Fit("gaus","","",-1,1);
      ccf5.SaveAs("cfit_py_res.png");  
      TCanvas ccf6;
      TH1F* pz_res = new TH1F("pz_res","pz resolution",100,-5,5);
      pz_res->GetXaxis()->SetTitle("pz^{cfit} - pz^{MC} / pz^{MC}");
      tree->Draw("(cfitHit0_pz-simHit0_pz)/simHit0_pz>>pz_res");
      pz_res->Fit("gaus","","",-1,1);
      ccf6.SaveAs("cfit_pz_res.png");  

      TCanvas ccf7;
      TH1F* cfit_pt_res = new TH1F("cfit_pt_res","1/pt residue",100,-0.05,0.05);
      cfit_pt_res->GetXaxis()->SetTitle("1./pt^{cfit} - 1./pt^{MC}");
      //tree->Draw("(sqrt(cfitHit0_px*cfitHit0_px+cfitHit0_py*cfitHit0_py)-sqrt(simHit0_px*simHit0_px+simHit0_py*simHit0_py))/(sqrt(simHit0_px*simHit0_px+simHit0_py*simHit0_py))>>cfit_pt_res");//,"sqrt(simHit0_px*simHit0_px+simHit0_py*simHit0_py)<5"
      tree->Draw("(1./sqrt(cfitHit0_px*cfitHit0_px+cfitHit0_py*cfitHit0_py)-1./sqrt(simHit0_px*simHit0_px+simHit0_py*simHit0_py))>>cfit_pt_res");// // /(1./sqrt(simHit0_px*simHit0_px+simHit0_py*simHit0_py))
      cfit_pt_res->Fit("gaus","","",-1,1);
      ccf7.SaveAs("cfit_pt_res.png");  
      
      TCanvas ccf8;
      TH1F* cfit_phi_res = new TH1F("cfit_phi_res","phi residue",100,-0.02,0.02);
      cfit_phi_res->GetXaxis()->SetTitle("phi^{cfit} - phi^{MC}");
      tree->Draw("atan2(cfitHit0_py,cfitHit0_px)-atan2(simHit0_py,simHit0_px)>>cfit_phi_res");
      cfit_phi_res->Fit("gaus","","",-0.1,0.1);
      ccf8.SaveAs("cfit_phi_res.png");  
      
      TCanvas ccf9;
      TH1F* cfit_theta_res = new TH1F("cfit_theta_res","theta residue",100,-0.03,0.03);
      cfit_theta_res->GetXaxis()->SetTitle("theta^{cfit} - theta^{MC}");
      tree->Draw("atan2(sqrt(cfitHit0_px*cfitHit0_px+cfitHit0_py*cfitHit0_py),cfitHit0_pz)-atan2(sqrt(simHit0_px*simHit0_px+simHit0_py*simHit0_py),simHit0_pz)>>cfit_theta_res");
      cfit_theta_res->Fit("gaus","","",-0.1,0.1);
      ccf9.SaveAs("cfit_theta_res.png");  



      TCanvas ccfp1;
      TH1F* x_pull = new TH1F("x_pull","x pull",100,-5,5);
      x_pull->GetXaxis()->SetTitle("(x^{cfit} - x^{MC})/#sigma(x^{cfit})");
      tree->Draw("(cfitHit0_x-simHit0_x)/cfitHit0_xe>>x_pull");
      x_pull->Fit("gaus","","",-2,2);
      ccfp1.SaveAs("cfit_x_pull.png");  
      TCanvas ccfp2;
      TH1F* y_pull = new TH1F("y_pull","y pull",100,-5,5);
      y_pull->GetXaxis()->SetTitle("(y^{cfit} - y^{MC})/#sigma(y^{cfit})");
      tree->Draw("(cfitHit0_y-simHit0_y)/cfitHit0_ye>>y_pull");
      y_pull->Fit("gaus","","",-2,2);
      ccfp2.SaveAs("cfit_y_pull.png");  
      TCanvas ccfp3;
      TH1F* z_pull = new TH1F("z_pull","z pull",100,-5,5);
      z_pull->GetXaxis()->SetTitle("(y^{cfit} - y^{MC})/#sigma(y^{cfit})");
      tree->Draw("(cfitHit0_z-simHit0_z)/cfitHit0_ze>>z_pull");
      z_pull->Fit("gaus","","",-2,2);
      ccfp3.SaveAs("cfit_z_pull.png");  

      TCanvas ccfp4;
      TH1F* px_pull = new TH1F("px_pull","px pull",100,-5,5);
      px_pull->GetXaxis()->SetTitle("(px^{cfit} - px^{MC})/#sigma(px^{cfit})");
      tree->Draw("(cfitHit0_px-simHit0_px)/cfitHit0_pxe>>px_pull");//,"sqrt(simHit0_px*simHit0_px+simHit0_py*simHit0_py)>2"
      px_pull->Fit("gaus","","",-2,2);
      ccfp4.SaveAs("cfit_px_pull.png");  
      TCanvas ccfp5;
      TH1F* py_pull = new TH1F("py_pull","py pull",100,-5,5);
      py_pull->GetXaxis()->SetTitle("(py^{cfit} - py^{MC})/#sigma(py^{cfit})");
      tree->Draw("(cfitHit0_py-simHit0_py)/cfitHit0_pye>>py_pull");
      py_pull->Fit("gaus","","",-2,2);
      ccfp5.SaveAs("cfit_py_pull.png");  
      TCanvas ccfp6;
      TH1F* pz_pull = new TH1F("pz_pull","pz pull",100,-5,5);
      pz_pull->GetXaxis()->SetTitle("(pz^{cfit} - pz^{MC})/#sigma(pz^{cfit})");
      tree->Draw("(cfitHit0_pz-simHit0_pz)/cfitHit0_pze>>pz_pull");
      pz_pull->Fit("gaus","","",-2,2);
      ccfp6.SaveAs("cfit_pz_pull.png");  

      TCanvas ccfr4;
      TH1F* px_residue = new TH1F("px_residue","px residue",100,-20,20);
      px_residue->GetXaxis()->SetTitle("px^{cfit} - px^{MC}");
      tree->Draw("cfitHit0_px-simHit0_px>>px_residue");//,"sqrt(simHit0_px*simHit0_px+simHit0_py*simHit0_py)>2"
      px_residue->Fit("gaus","","",-2,2);
      ccfr4.SaveAs("cfit_px_residue.png");  
      TCanvas ccfr5;
      TH1F* py_residue = new TH1F("py_residue","py residue",100,-20,20);
      py_residue->GetXaxis()->SetTitle("py^{cfit} - py^{MC}");
      tree->Draw("cfitHit0_py-simHit0_py>>py_residue");
      py_residue->Fit("gaus","","",-2,2);
      ccfr5.SaveAs("cfit_py_residue.png");  
      TCanvas ccfr6;
      TH1F* pz_residue = new TH1F("pz_residue","pz residue",100,-20,20);
      pz_residue->GetXaxis()->SetTitle("pz^{cfit} - pz^{MC}");
      tree->Draw("cfitHit0_pz-simHit0_pz>>pz_residue");
      pz_residue->Fit("gaus","","",-2,2);
      ccfr6.SaveAs("cfit_pz_residue.png");  
      

    }

  }

  if (doBuild) {

    TFile *_file0 = TFile::Open("build_validationtree.root");
    TTree * tree   = (TTree*)_file0->Get("tree");

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

    TCanvas c1;
    //h_gen_trk_Pt->GetYaxis()->SetRangeUser(0.,h_gen_trk_Pt->GetMaximum()*1.1);
    h_gen_trk_Pt->Draw();
    c1.SaveAs("gen_trk_Pt.png");

    TCanvas c2;
    h_gen_trk_Px->Draw();
    c2.SaveAs("gen_trk_Px.png");

    TCanvas c3;
    h_gen_trk_Py->Draw();
    c3.SaveAs("gen_trk_Py.png");

    TCanvas c4;
    h_gen_trk_Pz->Draw();
    c4.SaveAs("gen_trk_Pz.png");

    TCanvas c5;
    h_gen_trk_eta->Draw();
    c5.SaveAs("gen_trk_eta.png");

    TCanvas c6;
    h_gen_trk_phi->Draw();
    c6.SaveAs("gen_trk_phi.png");

    TCanvas c7;
    h_gen_trk_dPhi->Draw();
    c7.SaveAs("gen_trk_dPhi.png");

    TCanvas c8;
    h_gen_trk_dR->Draw();
    c8.SaveAs("gen_trk_dR.png");

    TCanvas c9;
    h_gen_trk_mindPhi->Draw();
    c9.SaveAs("gen_trk_mindPhi.png");

    TCanvas c10;
    h_gen_trk_mindR->Draw();
    c10.SaveAs("gen_trk_mindR.png");

  }

}
