void plotCompare_phiSmear(){

  gROOT->Reset();
  gStyle->SetOptStat(111111);

  TFile *_file0 = TFile::Open("build_validationtree_old.root");
  TFile *_file1 = TFile::Open("build_validationtree_new.root");

  TH1F* h_gen_trk_Pt_old = (TH1F*) _file0->Get("h_gen_trk_Pt");	  
  TH1F* h_gen_trk_Px_old = (TH1F*) _file0->Get("h_gen_trk_Px");	  
  TH1F* h_gen_trk_Py_old = (TH1F*) _file0->Get("h_gen_trk_Py");	  
  TH1F* h_gen_trk_Pz_old = (TH1F*) _file0->Get("h_gen_trk_Pz");	  
  TH1F* h_gen_trk_dPhi_old = (TH1F*) _file0->Get("h_gen_trk_dPhi");	  
  TH1F* h_gen_trk_dR_old = (TH1F*) _file0->Get("h_gen_trk_dR");	  
  TH1F* h_gen_trk_eta_old = (TH1F*) _file0->Get("h_gen_trk_eta");	  
  TH1F* h_gen_trk_mindPhi_old = (TH1F*) _file0->Get("h_gen_trk_mindPhi");
  TH1F* h_gen_trk_mindR_old = (TH1F*) _file0->Get("h_gen_trk_mindR");  
  TH1F* h_gen_trk_phi_old = (TH1F*) _file0->Get("h_gen_trk_phi");    
  TH1F* h_rec_trk_nHits_old = (TH1F*) _file0->Get("h_rec_trk_nHits");    
  TH1F* h_gen_hits_rad_old = (TH1F*) _file0->Get("h_gen_hits_rad");    
  TH1F* h_gen_hits_rad_lay3_old = (TH1F*) _file0->Get("h_gen_hits_rad_lay3");    
  TH1F* h_gen_hits_cov00_old = (TH1F*) _file0->Get("h_gen_hits_cov00");    
  TH1F* h_gen_hits_cov11_old = (TH1F*) _file0->Get("h_gen_hits_cov11");    

  TH1F* h_gen_trk_Pt = (TH1F*) _file1->Get("h_gen_trk_Pt");	  
  TH1F* h_gen_trk_Px = (TH1F*) _file1->Get("h_gen_trk_Px");	  
  TH1F* h_gen_trk_Py = (TH1F*) _file1->Get("h_gen_trk_Py");	  
  TH1F* h_gen_trk_Pz = (TH1F*) _file1->Get("h_gen_trk_Pz");	  
  TH1F* h_gen_trk_dPhi = (TH1F*) _file1->Get("h_gen_trk_dPhi");	  
  TH1F* h_gen_trk_dR = (TH1F*) _file1->Get("h_gen_trk_dR");	  
  TH1F* h_gen_trk_eta = (TH1F*) _file1->Get("h_gen_trk_eta");	  
  TH1F* h_gen_trk_mindPhi = (TH1F*) _file1->Get("h_gen_trk_mindPhi");
  TH1F* h_gen_trk_mindR = (TH1F*) _file1->Get("h_gen_trk_mindR");  
  TH1F* h_gen_trk_phi = (TH1F*) _file1->Get("h_gen_trk_phi");    
  TH1F* h_rec_trk_nHits = (TH1F*) _file1->Get("h_rec_trk_nHits");    
  TH1F* h_gen_hits_rad = (TH1F*) _file1->Get("h_gen_hits_rad");    
  TH1F* h_gen_hits_rad_lay3 = (TH1F*) _file1->Get("h_gen_hits_rad_lay3");    
  TH1F* h_gen_hits_cov00 = (TH1F*) _file0->Get("h_gen_hits_cov00");    
  TH1F* h_gen_hits_cov11 = (TH1F*) _file0->Get("h_gen_hits_cov11");    

  bool savePNG = true;
  TFile * rootout = new TFile("compare_plots/compareBuild.root","RECREATE");
  TCanvas * canvas = new TCanvas();
    
  // Do for all plots
  createPlot(canvas,h_gen_trk_Pt_old,h_gen_trk_Pt,rootout,savePNG);
  createPlot(canvas,h_gen_trk_Px_old,h_gen_trk_Px,rootout,savePNG);
  createPlot(canvas,h_gen_trk_Py_old,h_gen_trk_Py,rootout,savePNG);
  createPlot(canvas,h_gen_trk_Pz_old,h_gen_trk_Pz,rootout,savePNG);
  createPlot(canvas,h_gen_trk_dPhi_old,h_gen_trk_dPhi,rootout,savePNG);
  createPlot(canvas,h_gen_trk_dR_old,h_gen_trk_dR,rootout,savePNG);
  createPlot(canvas,h_gen_trk_eta_old,h_gen_trk_eta,rootout,savePNG);
  createPlot(canvas,h_gen_trk_mindPhi_old,h_gen_trk_mindPhi,rootout,savePNG);
  createPlot(canvas,h_gen_trk_mindR_old,h_gen_trk_mindR,rootout,savePNG);
  createPlot(canvas,h_gen_trk_phi_old,h_gen_trk_phi,rootout,savePNG);
  createPlot(canvas,h_rec_trk_nHits_old,h_rec_trk_nHits,rootout,savePNG);
  createPlot(canvas,h_gen_hits_rad_old,h_gen_hits_rad,rootout,savePNG);
  createPlot(canvas,h_gen_hits_rad_lay3_old,h_gen_hits_rad_lay3,rootout,savePNG);
  createPlot(canvas,h_gen_hits_cov00_old,h_gen_hits_cov00,rootout,savePNG);
  createPlot(canvas,h_gen_hits_cov11_old,h_gen_hits_cov11,rootout,savePNG);
  
  canvas->Close();
  rootout->Close();
}

void createPlot(TCanvas *canvas, TH1F * hist_old, TH1F * hist_new, TFile * rootout, bool savePNG) {  

  std::string hist_name_old = hist_old->GetName();
  std::string hist_name_new = hist_new->GetName();

  hist_name_old.append(" - XY smear");
  hist_name_new.append(" - Phi smear");

  hist_old->SetName(hist_name_old.c_str());
  hist_new->SetName(hist_name_new.c_str());

  hist_old->SetTitle("");
  hist_new->SetTitle("");

  bool DrawRatio = true;                                                
  canvas->cd();                                                                                                                                                                                                                         
  TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,0.8);                                                                                                                                                                        
  mainpad->Draw();  
  mainpad->cd();  

  double max = 0;                                                                                                                                                                                                       
  double V1max = hist_new->GetBinContent(hist_old->GetMaximumBin());  
  double V2max = hist_old->GetBinContent(hist_new->GetMaximumBin());

  max = (V1max>V2max) ? V1max : V2max;  

  hist_old->SetLineWidth(2);
  hist_old->SetLineStyle(1); 
  hist_old->SetLineColor(kRed);  
  hist_old->SetMaximum(max*(1.1));   
  hist_old->Draw(); 

  hist_new->SetLineWidth(2); 
  hist_new->SetLineStyle(1);
  hist_new->SetLineColor(kBlue);         
  hist_new->Draw("sames"); 

  mainpad->Update();


  if (hist_name_old != "h_rec_trk_nHits - XY smear"){
    TPaveStats *st1 = (TPaveStats*)(hist_old->GetListOfFunctions()->FindObject("stats"));            
    st1->SetX1NDC(0.77);                                                                                                                                                                                            
    st1->SetY1NDC(0.80);                                                                                                                                                                                               
    st1->SetX2NDC(0.98);                                                                                                                                               
    st1->SetY2NDC(0.97);

    Double_t defaulth = st1->GetY2NDC() - st1->GetY1NDC();
    Double_t gaph = 0.02;                                                                                                                                                                                           

    TPaveStats *st2 = (TPaveStats*)(hist_new->GetListOfFunctions()->FindObject("stats"));                                                                                               
    st2->SetX1NDC(0.77);                                                                                                                                                                         
    st2->SetY1NDC(st1->GetY1NDC() - 1.0*defaulth - gaph);                                                                                                                                        
    st2->SetX2NDC(0.98);                                                                                                                                                 
    st2->SetY2NDC(st1->GetY1NDC() - gaph);                                                                                                                                                                               
  }
  else {
    TPaveStats *st1 = (TPaveStats*)(hist_old->GetListOfFunctions()->FindObject("stats"));            
    st1->SetX1NDC(0.57);                                                                                                                                                                                            
    st1->SetY1NDC(0.60);                                                                                                                                                                                               
    st1->SetX2NDC(0.78);                                                                                                                                               
    st1->SetY2NDC(0.77);

    Double_t defaulth = st1->GetY2NDC() - st1->GetY1NDC();
    Double_t gaph = 0.02;                                                                                                                                                                                           

    TPaveStats *st2 = (TPaveStats*)(hist_new->GetListOfFunctions()->FindObject("stats"));                                                                                               
    st2->SetX1NDC(0.57);                                                                                                                                                                         
    st2->SetY1NDC(st1->GetY1NDC() - 1.0*defaulth - gaph);                                                                                                                                        
    st2->SetX2NDC(0.78);                                                                                                                                                 
    st2->SetY2NDC(st1->GetY1NDC() - gaph);                                                                                                                                                                               
  }

  TLegend *leg = new TLegend(0.32,0.86,0.76,0.97);                                                                                                                                                                    
  leg->SetTextSize(0.042);
  leg->SetTextFont(42);
  leg->SetFillColor(10); 
  leg->SetBorderSize(1);
  leg->AddEntry(hist_old, "X Y Smear", "L" );
  leg->AddEntry(hist_new, "Phi Smear", "L" );                                                                                                                                                            
  leg->Draw("SAME");                             

  if (DrawRatio){
    canvas->cd();
    TPad* respad = new TPad("respad","respad",0.0,0.78,1.0,0.95); 
    respad->SetTopMargin(1.05); 
    respad->Draw();
    respad->cd(); 

    TH1F* hratio = (TH1F*) hist_new->Clone("hratio");
    hratio->Divide(hist_old);
    hratio->SetMaximum(hratio->GetMaximum()*1.1);
    hratio->SetMinimum(hratio->GetMinimum()*1.1);
    //    if (hratio->GetMinimum()==0.0) hratio->SetMinimum(1.0/hratio->GetMaximum());
    //    hratio->SetMinimum(1.0/hratio->GetMaximum()); 
    hratio->GetYaxis()->SetLabelSize(0.1);                                                                                                                                                            
    //    hratio->GetYaxis()->SetRangeUser(0,2);                                                                                                                                                            
    hratio->GetXaxis()->SetLabelSize(0);                                                                                                                                                   
    hratio->GetXaxis()->SetTitleSize(0);                                                                                                                                                                         
    hratio->GetYaxis()->SetTitleSize(0.22);                                                                                                                                                                     
    hratio->GetYaxis()->SetTitleOffset(0.15);                                                                                                                                                                 
    hratio->GetYaxis()->SetLabelSize(0.2);                                                                                                                                                                    
    hratio->GetYaxis()->SetNdivisions(5);                                                                                                                                                
    hratio->GetYaxis()->SetTitle("Phi/XY");  
    hratio->SetStats(kFALSE);
    hratio->Draw(); 
  }                          
  
  rootout->cd();
  canvas->Write(Form("%s_compare",hist_name_new.c_str()));

  if (savePNG){
    canvas->Print(Form("compare_plots/%s_compare.png",hist_name_new.c_str()),"png");
  }

}
