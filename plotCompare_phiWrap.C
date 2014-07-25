void plotCompare_phiWrap(){

  // General Style
  gROOT->Reset();
  gStyle->SetOptStat(111111);
  gStyle->SetHatchesSpacing(2.0);  
  gStyle->SetHatchesLineWidth(1.0);

  // Open Input
  TFile *_file0 = TFile::Open("build_validationtree_old.root");
  TFile *_file1 = TFile::Open("build_validationtree_new.root");

  // Define Output
  Bool_t savePNG = true;

  FileStat_t dummyFileStat;
  TString outDir = "compare_plots_piPM";

  if (gSystem->GetPathInfo(outDir.Data(), dummyFileStat) == 1){
    TString mkDir = "mkdir -p ";
    mkDir += outDir.Data();
    gSystem->Exec(mkDir.Data());
  }
  
  TFile * rootout = new TFile(Form("%s/compareBuild.root",outDir.Data()),"RECREATE");
  TCanvas * canvas = new TCanvas();

  // Get all histos
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
  TH1F* h_rec_seed_ndPhiOB_lay3_old = (TH1F*) _file0->Get("h_rec_seed_ndPhiOB_lay3");
  TH1F* h_rec_seed_ndPhiOB_lay5_old = (TH1F*) _file0->Get("h_rec_seed_ndPhiOB_lay5");
  TH1F* h_rec_seed_ndPhiOB_lay7_old = (TH1F*) _file0->Get("h_rec_seed_ndPhiOB_lay7");
  TH1F* h_rec_seed_ndPhiOB_lay10_old = (TH1F*) _file0->Get("h_rec_seed_ndPhiOB_lay10");
  TH1F* h_rec_seed_ndPhiOB_lay3_fracCands_old = (TH1F*) _file0->Get("h_rec_seed_ndPhiOB_lay3_fracCands");
  TH1F* h_rec_seed_ndPhiOB_lay5_fracCands_old = (TH1F*) _file0->Get("h_rec_seed_ndPhiOB_lay5_fracCands");
  TH1F* h_rec_seed_ndPhiOB_lay7_fracCands_old = (TH1F*) _file0->Get("h_rec_seed_ndPhiOB_lay7_fracCands");
  TH1F* h_rec_seed_ndPhiOB_lay10_fracCands_old = (TH1F*) _file0->Get("h_rec_seed_ndPhiOB_lay10_fracCands");

  TH1F* h_gen_trk_Pt_new = (TH1F*) _file1->Get("h_gen_trk_Pt");	  
  TH1F* h_gen_trk_Px_new = (TH1F*) _file1->Get("h_gen_trk_Px");	  
  TH1F* h_gen_trk_Py_new = (TH1F*) _file1->Get("h_gen_trk_Py");	  
  TH1F* h_gen_trk_Pz_new = (TH1F*) _file1->Get("h_gen_trk_Pz");	  
  TH1F* h_gen_trk_dPhi_new = (TH1F*) _file1->Get("h_gen_trk_dPhi");	  
  TH1F* h_gen_trk_dR_new = (TH1F*) _file1->Get("h_gen_trk_dR");	  
  TH1F* h_gen_trk_eta_new = (TH1F*) _file1->Get("h_gen_trk_eta");	  
  TH1F* h_gen_trk_mindPhi_new = (TH1F*) _file1->Get("h_gen_trk_mindPhi");
  TH1F* h_gen_trk_mindR_new = (TH1F*) _file1->Get("h_gen_trk_mindR");  
  TH1F* h_gen_trk_phi_new = (TH1F*) _file1->Get("h_gen_trk_phi");    
  TH1F* h_rec_trk_nHits_new = (TH1F*) _file1->Get("h_rec_trk_nHits");    
  TH1F* h_gen_hits_rad_new = (TH1F*) _file1->Get("h_gen_hits_rad");    
  TH1F* h_gen_hits_rad_lay3_new = (TH1F*) _file1->Get("h_gen_hits_rad_lay3");    
  TH1F* h_rec_seed_ndPhiOB_lay3_new = (TH1F*) _file1->Get("h_rec_seed_ndPhiOB_lay3");
  TH1F* h_rec_seed_ndPhiOB_lay5_new = (TH1F*) _file1->Get("h_rec_seed_ndPhiOB_lay5");
  TH1F* h_rec_seed_ndPhiOB_lay7_new = (TH1F*) _file1->Get("h_rec_seed_ndPhiOB_lay7");
  TH1F* h_rec_seed_ndPhiOB_lay10_new = (TH1F*) _file1->Get("h_rec_seed_ndPhiOB_lay10");
  TH1F* h_rec_seed_ndPhiOB_lay3_fracCands_new = (TH1F*) _file1->Get("h_rec_seed_ndPhiOB_lay3_fracCands");
  TH1F* h_rec_seed_ndPhiOB_lay5_fracCands_new = (TH1F*) _file1->Get("h_rec_seed_ndPhiOB_lay5_fracCands");
  TH1F* h_rec_seed_ndPhiOB_lay7_fracCands_new = (TH1F*) _file1->Get("h_rec_seed_ndPhiOB_lay7_fracCands");
  TH1F* h_rec_seed_ndPhiOB_lay10_fracCands_new = (TH1F*) _file1->Get("h_rec_seed_ndPhiOB_lay10_fracCands");
    
  // Make Comparison plots
  createPlot(canvas,h_gen_trk_Pt_old,h_gen_trk_Pt_new,outDir,rootout,savePNG);
  createPlot(canvas,h_gen_trk_Px_old,h_gen_trk_Px_new,outDir,rootout,savePNG);
  createPlot(canvas,h_gen_trk_Py_old,h_gen_trk_Py_new,outDir,rootout,savePNG);
  createPlot(canvas,h_gen_trk_Pz_old,h_gen_trk_Pz_new,outDir,rootout,savePNG);
  createPlot(canvas,h_gen_trk_dPhi_old,h_gen_trk_dPhi_new,outDir,rootout,savePNG);
  createPlot(canvas,h_gen_trk_dR_old,h_gen_trk_dR_new,outDir,rootout,savePNG);
  createPlot(canvas,h_gen_trk_eta_old,h_gen_trk_eta_new,outDir,rootout,savePNG);
  createPlot(canvas,h_gen_trk_mindPhi_old,h_gen_trk_mindPhi_new,outDir,rootout,savePNG);
  createPlot(canvas,h_gen_trk_mindR_old,h_gen_trk_mindR_new,outDir,rootout,savePNG);
  createPlot(canvas,h_gen_trk_phi_old,h_gen_trk_phi_new,outDir,rootout,savePNG);
  createPlot(canvas,h_rec_trk_nHits_old,h_rec_trk_nHits_new,outDir,rootout,savePNG);
  createPlot(canvas,h_gen_hits_rad_old,h_gen_hits_rad_new,outDir,rootout,savePNG);
  createPlot(canvas,h_gen_hits_rad_lay3_old,h_gen_hits_rad_lay3_new,outDir,rootout,savePNG);
  createPlot(canvas,h_rec_seed_ndPhiOB_lay3_old,h_rec_seed_ndPhiOB_lay3_new,outDir,rootout,savePNG);
  createPlot(canvas,h_rec_seed_ndPhiOB_lay5_old,h_rec_seed_ndPhiOB_lay5_new,outDir,rootout,savePNG);
  createPlot(canvas,h_rec_seed_ndPhiOB_lay7_old,h_rec_seed_ndPhiOB_lay7_new,outDir,rootout,savePNG);
  createPlot(canvas,h_rec_seed_ndPhiOB_lay10_old,h_rec_seed_ndPhiOB_lay10_new,outDir,rootout,savePNG);
  createPlot(canvas,h_rec_seed_ndPhiOB_lay3_fracCands_old,h_rec_seed_ndPhiOB_lay3_fracCands_new,outDir,rootout,savePNG);
  createPlot(canvas,h_rec_seed_ndPhiOB_lay5_fracCands_old,h_rec_seed_ndPhiOB_lay5_fracCands_new,outDir,rootout,savePNG);
  createPlot(canvas,h_rec_seed_ndPhiOB_lay7_fracCands_old,h_rec_seed_ndPhiOB_lay7_fracCands_new,outDir,rootout,savePNG);
  createPlot(canvas,h_rec_seed_ndPhiOB_lay10_fracCands_old,h_rec_seed_ndPhiOB_lay10_fracCands_new,outDir,rootout,savePNG);
  
  // Close the output files
  canvas->Close();
  rootout->Close();
}

void createPlot(TCanvas *canvas, TH1F * hist_old, TH1F * hist_new, TString outDir, TFile * rootout, Bool_t savePNG) {  

  Bool_t DrawRatio = kTRUE;   
  canvas->cd();
  TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,0.8);                                                                                                                                                                        
  mainpad->Draw();  
  mainpad->cd();  

  Double_t max = 0;                                                                                                                                                                                                       
  Double_t V1max = hist_new->GetBinContent(hist_old->GetMaximumBin());  
  Double_t V2max = hist_old->GetBinContent(hist_new->GetMaximumBin());

  max = (V1max>V2max) ? V1max : V2max;  

  TString hist_name = hist_old->GetName();
  canvas->SetName(hist_name.Data());

  hist_old->SetName("No Periodicity");
  hist_new->SetName("Phi Wrap Included");

  //  hist_old->GetXaxis()->SetRangeUser(1,hist_old->GetXaxis()->GetLast());
  //  hist_new->GetXaxis()->SetRangeUser(1,hist_new->GetXaxis()->GetLast());

  hist_old->SetTitle("");
  hist_new->SetTitle("");

  //  hist_old->SetFillStyle(3215);
  //  hist_old->SetFillColor(kRed);

  hist_old->GetYaxis()->SetTitleSize(0.05);
  hist_old->GetYaxis()->SetTitleOffset(0.9);                                                                                                                                                                 
  
  hist_old->SetLineWidth(2);
  hist_old->SetLineStyle(1); 
  hist_old->SetLineColor(kRed);  
  hist_old->SetMaximum(max*(1.1));   
  hist_old->Draw("E"); 

  //  hist_new->SetFillStyle(3251);
  //  hist_new->SetFillColor(kBlue);
  hist_new->GetYaxis()->SetTitleOffset(1.15);                                                                                                                                                                 
  hist_new->SetLineWidth(2); 
  hist_new->SetLineStyle(1);
  hist_new->SetLineColor(kBlue);         
  hist_new->Draw("E sames"); 

  mainpad->Update();

  /*  
  if (hist_name.EqualTo("h_rec_seed_ndPhiOB_lay5",TString::kExact) == 1){
    cout << "Old 5: ," << hist_old->Integral(2,20) << " +/- " << TMath::Sqrt(hist_old->Integral(2,20)) << endl;
    cout << "New 5: ," << hist_new->Integral(2,20) << " +/- " << TMath::Sqrt(hist_new->Integral(2,20)) << endl;
  }
  else if (hist_name.EqualTo("h_rec_seed_ndPhiOB_lay7",TString::kExact) == 1){
    cout << "Old 7: ," << hist_old->Integral(2,20) << " +/- " << TMath::Sqrt(hist_old->Integral(2,20)) << endl;
    cout << "New 7: ," << hist_new->Integral(2,20) << " +/- " << TMath::Sqrt(hist_new->Integral(2,20)) << endl;
  }
  else if (hist_name.EqualTo("h_rec_seed_ndPhiOB_lay10",TString::kExact) == 1){
    cout << "Old 10: ," << hist_old->Integral(2,20) << " +/- " << TMath::Sqrt(hist_old->Integral(2,20)) << endl;
    cout << "New 10: ," << hist_new->Integral(2,20) << " +/- " << TMath::Sqrt(hist_new->Integral(2,20)) << endl;
  }
  */ 

  if (hist_name.EqualTo("h_rec_trk_nHits",TString::kExact) != 1){
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
  leg->AddEntry(hist_old, "No Periodicity", "L" );
  leg->AddEntry(hist_new, "Phi Wrap Included", "L" );                                                                                                                                                            
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
    //    hratio->GetYaxis()->SetRangeUser(0,2);                                                                                                                                                            
    hratio->GetXaxis()->SetLabelSize(0);                                                                                                                                                   
    hratio->GetXaxis()->SetTitleSize(0);                                                                                                                                                                         
    hratio->GetYaxis()->SetTitleSize(0.22);                                                                                                                                                                     
    hratio->GetYaxis()->SetTitleOffset(0.2);                                                                                                                                                                 
    hratio->GetYaxis()->SetLabelSize(0.2);                                                                                                                                                                    
    hratio->GetYaxis()->SetNdivisions(5);                                                                                                                                                
    hratio->GetYaxis()->SetTitle("New/Old");  
    hratio->SetStats(kFALSE);
    hratio->Draw(); 
  }                          
  
  rootout->cd();
  canvas->Write(Form("%s_compare",hist_name.Data()));

  if (savePNG){
    canvas->Print(Form("%s/%s_compare.png",outDir.Data(),hist_name.Data()),"png");
  }

}
