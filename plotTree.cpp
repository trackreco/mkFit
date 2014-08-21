#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"

#include <iostream>
#include <cmath>

void plotTree(){

  gROOT->Reset();

  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1011);

  bool doFit = 1;
  bool doBuild = 0;
  bool doSim = 0;

  if (doFit) {
    
    // ++++ Make Output Directory ++++ //
    
    FileStat_t dummyFileStat;
    TString outDir = "png_fittest_Solids";

    if (gSystem->GetPathInfo(outDir.Data(), dummyFileStat) == 1){
      TString mkDir = "mkdir -p ";
      mkDir += outDir.Data();
      gSystem->Exec(mkDir.Data());
    }

    // ++++ Get Root File, Trees ++++ //

    TString rootfile = "validationtree.root";
    
    TFile * file     = TFile::Open(Form("%s",rootfile.Data()));
    TTree * ptTree   = (TTree*)file->Get("ptTree");
    TTree * posTree  = (TTree*)file->Get("posTree");
    
    // ++++ Initialize Histos ++++ //
    
    /////////// pt

    TCanvas c_pt;

    TH1F* h_pt_res_fit = new TH1F("h_pt_res_fit","p_{T} Resolution MC-Fit",200,-1.0,1.0);
    h_pt_res_fit->GetXaxis()->SetTitle("(p_{T}^{MC} - p_{T}^{fit})/p_{T}^{MC}");
    h_pt_res_fit->GetYaxis()->SetTitle("Tracks");

    TH1F* h_pt_pull_fit = new TH1F("h_pt_pull_fit","p_{T} Pull MC-Fit",200,-10,10);
    h_pt_pull_fit->GetXaxis()->SetTitle("(p_{T}^{MC} - p_{T}^{fit})/#sigma(p_{T}^{fit})");
    h_pt_pull_fit->GetYaxis()->SetTitle("Tracks");

    /////////// xpos

    TCanvas c_xpos;

    TH1F* h_x_res_mc = new TH1F("h_x_res_mc","x Resolution MC Hit (#phi)-Initial Hit (XY)",200,-0.01,0.01);
    h_x_res_mc->GetXaxis()->SetTitle("(x^{MC} - x^{init})/x^{init}");    
    h_x_res_mc->GetYaxis()->SetTitle("Hits");

    TH1F* h_x_pull_mc = new TH1F("h_x_pull_mc","x Pull MC Hit (#phi)-Initial Hit (XY)",200,-10.0,10.0);
    h_x_pull_mc->GetXaxis()->SetTitle("(x^{MC} - x^{init})/#sigma(x^{MC})");    
    h_x_pull_mc->GetYaxis()->SetTitle("Hits");    

    TH1F* h_x_res_prop = new TH1F("h_x_res_prop","x Resolution Initial Hit-Prop Hit",200,-0.01,0.01);
    h_x_res_prop->GetXaxis()->SetTitle("(x^{init} - x^{prop})/x^{init}");    
    h_x_res_prop->GetYaxis()->SetTitle("Hits");    

    TH1F* h_x_pull_prop = new TH1F("h_x_pull_prop","x Pull Initial Hit-Prop Hit",200,-10.0,10.0);
    h_x_pull_prop->GetXaxis()->SetTitle("(x^{init} - x^{prop})/#sigma(x^{prop})");    
    h_x_pull_prop->GetYaxis()->SetTitle("Hits");

    TH1F* h_x_res_update = new TH1F("h_x_res_update","x Resolution Initial Hit-Update Hit",200,-0.01,0.01);
    h_x_res_update->GetXaxis()->SetTitle("(x^{init} - x^{update})/x^{init}");    
    h_x_res_update->GetYaxis()->SetTitle("Hits");    

    TH1F* h_x_pull_update = new TH1F("h_x_pull_update","x Pull Init Hit-Update Hit",100,-10.0,10.0);
    h_x_pull_update->GetXaxis()->SetTitle("(x^{init} - x^{update})/#sigma(x^{update})");    
    h_x_pull_update->GetYaxis()->SetTitle("Hits");

    /////////// phi

    TCanvas c_phi;
  
    TH1F* h_phi_res_update = new TH1F("h_phi_res_update","#phi Resolution Initial-Updated",200,-0.01,0.01);
    h_phi_res_update->GetXaxis()->SetTitle("(#phi^{init} - #phi^{updated})/#phi^{init}");
    h_phi_res_update->GetYaxis()->SetTitle("Hits");

    TH1F* h_phi_pull_update = new TH1F("h_phi_pull_update","#phi Pull",200,-10,10);
    h_phi_pull_update->GetXaxis()->SetTitle("(#phi^{init} - #phi^{updated})/#sigma(#phi^{udpated})");
    h_phi_pull_update->GetYaxis()->SetTitle("Hits");

    TH1F* h_phi_overflow = new TH1F("h_phi_overflow","#phi of x_update pull |#sigma| > 10",67,-3.35,3.35);
    h_phi_overflow->GetXaxis()->SetTitle("#phi");    
    h_phi_overflow->GetYaxis()->SetTitle("Hits");

    // ++++ Calculate Resolution/Sigma Fill Histos ++++ //

    /////////// pt

    Float_t pt_mc, pt_fit, pt_err;
    ptTree->SetBranchAddress("pt_mc",&pt_mc);
    ptTree->SetBranchAddress("pt_fit",&pt_fit);
    ptTree->SetBranchAddress("pt_err",&pt_err);

    Float_t ptRes = 0., ptSig = 0.;

    UInt_t n_entriesPt = (UInt_t)ptTree->GetEntries();
    for(UInt_t i = 0; i < n_entriesPt; i++){
      ptTree->GetEntry(i);

      ptRes = (pt_mc - pt_fit)/pt_mc;
      ptSig = (pt_mc - pt_fit)/pt_err;

      h_pt_res_fit->Fill(ptRes);      
      if (!isnan(ptSig)){
	h_pt_pull_fit->Fill(ptSig);      
      }
    }
    /////////// xpos

    Float_t x_init, x_mc, x_prop, x_update, x_mcerr, x_perr, x_uerr;
    posTree->SetBranchAddress("x_init",&x_init);
    posTree->SetBranchAddress("x_mc",&x_mc);
    posTree->SetBranchAddress("x_prop",&x_prop);
    posTree->SetBranchAddress("x_mcerr",&x_mcerr);
    posTree->SetBranchAddress("x_perr",&x_perr);
    posTree->SetBranchAddress("x_update",&x_update);
    posTree->SetBranchAddress("x_uerr",&x_uerr);

    /////////// phi

    Float_t phi_init, phi_update, phi_uerr;
    posTree->SetBranchAddress("phi_init",&phi_init);
    posTree->SetBranchAddress("phi_update",&phi_update);
    posTree->SetBranchAddress("phi_uerr",&phi_uerr);

    Float_t mcRes = 0., mcSig = 0., propRes = 0., propSig = 0., updateRes = 0., updateSig = 0.;
    Float_t phiRes = 0., phiSig = 0.;

    UInt_t n_entriesPos = (UInt_t)posTree->GetEntries();
    for(UInt_t i = 0; i < n_entriesPos; i++){
      posTree->GetEntry(i);

      /////////// xpos

      mcRes = (x_init - x_mc)/x_init;
      mcSig = (x_init - x_mc)/sqrt(x_mcerr);
            
      propRes = (x_init - x_prop)/x_init;
      propSig = (x_init - x_prop)/sqrt(x_perr);

      updateRes = (x_init - x_update)/x_init;
      updateSig = (x_init - x_update)/sqrt(x_uerr);

      h_x_res_mc->Fill(mcRes);      
      if (!isnan(mcSig)){
	h_x_pull_mc->Fill(mcSig);      
      }
      
      h_x_res_prop->Fill(propRes);      
      if (!isnan(propSig)){
	h_x_pull_prop->Fill(propSig);      
      }
      
      h_x_res_update->Fill(updateRes);      
      if (!isnan(updateSig)){
	h_x_pull_update->Fill(updateSig);      
      }

      /////////// phi

      phiRes = (phi_init - phi_update)/phi_init;
      phiSig = (phi_init - phi_update)/phi_uerr;

      h_phi_res_update->Fill(phiRes);
      
      if (!isnan(phiSig)){
	h_phi_pull_update->Fill(phiSig);
      }
      
      if ( (std::abs(updateSig) > 10.) && (!isnan(updateSig)) ){
	h_phi_overflow->Fill(phi_init);
      }
    }

    // ++++ Draw Histos ++++ //

    c_pt.cd();
    h_pt_res_fit->Draw();
    h_pt_res_fit->Fit("gaus","","",-0.3,0.3);
    c_pt.SaveAs(Form("%s/pt_res_fit.png",outDir.Data()));  
    
    h_pt_pull_fit->Draw();
    h_pt_pull_fit->Fit("gaus","","",-3.0,3.0);
    c_pt.SaveAs(Form("%s/pt_pull_fit.png",outDir.Data()));  

    c_phi.cd();
    h_phi_res_update->Draw();
    h_phi_res_update->Fit("gaus","","",-0.005,0.005);
    c_phi.SaveAs(Form("%s/phi_res_update.png",outDir.Data()));  

    h_phi_pull_update->Draw();
    h_phi_pull_update->Fit("gaus","","",-3.0,3.0);
    c_phi.SaveAs(Form("%s/phi_pull_update.png",outDir.Data()));  
 
    h_phi_overflow->Draw();
    c_phi.SaveAs(Form("%s/phi_overflow.png",outDir.Data()));  
    
    c_xpos.cd();
    h_x_res_mc->Draw();
    h_x_res_mc->Fit("gaus","","",-0.005,0.005);
    c_xpos.SaveAs(Form("%s/x_res_MC.png",outDir.Data()));  

    h_x_pull_mc->Draw();
    h_x_pull_mc->Fit("gaus","","",-3.0,3.0);
    c_xpos.SaveAs(Form("%s/x_pull_MC.png",outDir.Data()));  

    h_x_res_prop->Draw();
    h_x_res_prop->Fit("gaus","","",-0.005,0.005);
    c_xpos.SaveAs(Form("%s/x_res_prop.png",outDir.Data()));  

    h_x_pull_prop->Draw();
    h_x_pull_prop->Fit("gaus","","",-3.0,3.0);
    c_xpos.SaveAs(Form("%s/x_pull_prop.png",outDir.Data()));  

    h_x_res_update->Draw();
    h_x_res_update->Fit("gaus","","",-0.005,0.005);
    c_xpos.SaveAs(Form("%s/x_res_update.png",outDir.Data()));  
    
    h_x_pull_update->Draw();
    h_x_pull_update->Fit("gaus","","",-3.0,3.0);
    c_xpos.SaveAs(Form("%s/x_pull_update.png",outDir.Data()));  

    // ++++ Move Input Root File ++++ //

    TString mvFile = "mv ";
    mvFile += rootfile.Data();
    mvFile += " ";
    mvFile += outDir.Data();
    gSystem->Exec(mvFile.Data());
    
  }
  /*
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
  */
}
