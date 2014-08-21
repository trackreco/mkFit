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

  /*
  for (int i = 0; i < 8; i++){
    float pi4_i = TMath::TwoPi()*(float(i)/8.0);
    std::cout << pi4_i << std::endl;
  }
  */

  if (doFit) {
    
    FileStat_t dummyFileStat;
    TString outDir = "png_fittest";

    if (gSystem->GetPathInfo(outDir.Data(), dummyFileStat) == 1){
      TString mkDir = "mkdir -p ";
      mkDir += outDir.Data();
      gSystem->Exec(mkDir.Data());
    }

    TString rootfile = "validationtree.root";
    
    TFile * _file0   = TFile::Open(Form("%s",rootfile.Data()));
    TTree * ptTree  = (TTree*)_file0->Get("ptTree");
    TTree * posTree = (TTree*)_file0->Get("posTree");

    TCanvas c_pt;
    c_pt.cd();
    TH1F* h_pt_res_fit = new TH1F("h_pt_res_fit","p_{T} Resolution MC-Fit",200,-1.0,1.0);
    h_pt_res_fit->GetXaxis()->SetTitle("(p_{T}^{MC} - p_{T}^{fit})/p_{T}^{MC}");
    h_pt_res_fit->GetYaxis()->SetTitle("Tracks");
    ptTree->Draw("(pt_mc-pt_fit)/pt_mc>>h_pt_res_fit");
    h_pt_res_fit->Fit("gaus","","",-0.3,0.3);
    c_pt.SaveAs(Form("%s/pt_res_fit.png",outDir.Data()));  
    
    c_pt.cd();
    TH1F* h_pt_pull_fit = new TH1F("h_pt_pull_fit","p_{T} Pull MC-Fit",200,-10,10);
    h_pt_pull_fit->GetXaxis()->SetTitle("(p_{T}^{MC} - p_{T}^{fit})/#sigma(p_{T}^{fit})");
    h_pt_pull_fit->GetYaxis()->SetTitle("Tracks");
    ptTree->Draw("(pt_mc-pt_fit)/pt_err>>h_pt_pull_fit");
    h_pt_pull_fit->Fit("gaus","","",-3.0,3.0);
    c_pt.SaveAs(Form("%s/pt_pull_fit.png",outDir.Data()));  
  
    TCanvas c_phi;
    c_phi.cd();
    TH1F* h_phi_res_update = new TH1F("h_phi_res_update","#phi Resolution Initial-Updated",200,-0.01,0.01);
    h_phi_res_update->GetXaxis()->SetTitle("(#phi^{init} - #phi^{updated})/#phi^{init}");
    h_phi_res_update->GetYaxis()->SetTitle("Hits");
    posTree->Draw("(phi_init-phi_update)/phi_init>>h_phi_res_update");
    h_phi_res_update->Fit("gaus","","",-0.005,0.005);
    c_phi.SaveAs(Form("%s/phi_res_update.png",outDir.Data()));  

    c_phi.cd();
    TH1F* h_phi_pull_update = new TH1F("h_phi_pull_update","#phi Pull",200,-10,10);
    h_phi_pull_update->GetXaxis()->SetTitle("(#phi^{init} - #phi^{updated})/#sigma(#phi^{udpated})");
     h_phi_pull_update->GetYaxis()->SetTitle("Hits");
    posTree->Draw("(phi_mc-phi_update)/phi_uerr>>h_phi_pull_update");
    h_phi_pull_update->Fit("gaus","","",-3.0,3.0);
    c_phi.SaveAs(Form("%s/phi_pull_update.png",outDir.Data()));  

    TCanvas c_xpos;
    c_xpos.cd();
    TH1F* h_x_res_MC = new TH1F("h_x_res_MC","x Resolution MC Hit (#phi)-Initial Hit (XY)",200,-0.01,0.01);
    h_x_res_MC->GetXaxis()->SetTitle("(x^{MC} - x^{init})/x^{init}");    
    h_x_res_MC->GetYaxis()->SetTitle("Hits");
    posTree->Draw("(x_mc-x_init)/x_init>>h_x_res_MC");
    h_x_res_MC->Fit("gaus","","",-0.005,0.005);
    c_xpos.SaveAs(Form("%s/x_res_MC.png",outDir.Data()));  

    c_xpos.cd();
    TH1F* h_x_pull_MC = new TH1F("h_x_pull_MC","x Pull MC Hit (#phi)-Initial Hit (XY)",200,-10.0,10.0);
    h_x_pull_MC->GetXaxis()->SetTitle("(x^{MC} - x^{init})/#sigma(x^{MC})");    
    h_x_pull_MC->GetYaxis()->SetTitle("Hits");    
    posTree->Draw("(x_mc-x_init)/x_mcerr>>h_x_pull_MC");
    h_x_pull_MC->Fit("gaus","","",-3.0,3.0);
    c_xpos.SaveAs(Form("%s/x_pull_MC.png",outDir.Data()));  

    c_xpos.cd();
    TH1F* h_x_res_prop = new TH1F("h_x_res_prop","x Resolution Initial Hit-Prop Hit",200,-0.01,0.01);
    h_x_res_prop->GetXaxis()->SetTitle("(x^{init} - x^{prop})/x^{init}");    
    h_x_res_prop->GetYaxis()->SetTitle("Hits");    
    posTree->Draw("(x_init-x_prop)/x_init>>h_x_res_prop");
    h_x_res_prop->Fit("gaus","","",-0.005,0.005);
    c_xpos.SaveAs(Form("%s/x_res_prop.png",outDir.Data()));  

    c_xpos.cd();
    TH1F* h_x_pull_prop = new TH1F("h_x_pull_prop","x Pull Initial Hit-Prop Hit",200,-10.0,10.0);
    h_x_pull_prop->GetXaxis()->SetTitle("(x^{init} - x^{prop})/#sigma(x^{prop})");    
    h_x_pull_prop->GetYaxis()->SetTitle("Hits");
    posTree->Draw("(x_init-x_prop)/x_perr>>h_x_pull_prop");
    h_x_pull_prop->Fit("gaus","","",-3.0,3.0);
    c_xpos.SaveAs(Form("%s/x_pull_prop.png",outDir.Data()));  

    c_xpos.cd();
    TH1F* h_x_res_update = new TH1F("h_x_res_update","x Resolution Initial Hit-Update Hit",200,-0.01,0.01);
    h_x_res_update->GetXaxis()->SetTitle("(x^{init} - x^{update})/x^{init}");    
    h_x_res_update->GetYaxis()->SetTitle("Hits");    
    posTree->Draw("(x_init-x_update)/x_init>>h_x_res_update");
    h_x_res_update->Fit("gaus","","",-0.005,0.005);
    c_xpos.SaveAs(Form("%s/x_res_update.png",outDir.Data()));  

    TCanvas c_xpos;
    c_xpos.cd();
    TH1F* h_x_pull_update = new TH1F("h_x_pull_update","x Pull Init Hit-Update Hit",100,-10.0,10.0);
    h_x_pull_update->GetXaxis()->SetTitle("(x^{init} - x^{update})/#sigma(x^{update})");    
    h_x_pull_update->GetYaxis()->SetTitle("Hits");

    
    TH1F* h_phi_overflow = new TH1F("h_phi_overflow","#phi of x_update pull |#sigma| > 10",67,-3.35,3.35);
    h_phi_overflow->GetXaxis()->SetTitle("#phi");    
    h_phi_overflow->GetYaxis()->SetTitle("Hits");
    
    Float_t x_init, x_update, x_uerr, phi_init;
    posTree->SetBranchAddress("x_init",&x_init);
    posTree->SetBranchAddress("x_update",&x_update);
    posTree->SetBranchAddress("x_uerr",&x_uerr);
    //    posTree->SetBranchAddress("phi_init",&phi_init);

    Int_t n_entries = (Int_t)posTree->GetEntries();
    int bin =0;
    int bin_counter[103];
    for (int j = -1; j <102; j++){
      bin_counter[j] = 0;
    }

    Double_t updateSigma = 0;
    int ninside = 0, noutside = 0;
    int nover = 0, nunder = 0,ihateroot=0;
    std::cout << "before: " << h_x_pull_update->GetBinContent(101) << std::endl;
    for(Int_t i = 0; i<n_entries; i++){
      updateSigma =0.0;
      posTree->GetEntry(i);
      updateSigma = (x_init - x_update)/sqrt(x_uerr);
      //      bin = h_x_pull_update->Fill(updateSigma);      


      if (isnan(updateSigma)) {
	std::cout << "problem: updateSigma:" << updateSigma 
		  <<", " << x_init
		  <<", " << x_update 
		  <<", " << x_uerr 
		  << std::endl;
      }
      else 
	bin = h_x_pull_update->Fill(updateSigma);      




      if (std::abs(updateSigma) > 10.){
// 	 	std::cout << "x_init : " << x_init << " x_update: " << x_update << std::endl
//  		  << "sqrt(x_uerr): " << x_uerr << " Entry: " << i << std::endl
//  	      	  << "updateSigma: " << updateSigma << std::endl;
	
      //h_phi_overflow->Fill(phi_init);
	h_x_pull_update->Fill(updateSigma);      
	ninside++;
      }
      else
	noutside++;
      if ( updateSigma >= 10.0 ) 
	nover++;
      if ( updateSigma < -10.0 ) 
	nunder++;
      if ( (updateSigma < h_x_pull_update->GetBinLowEdge(1)) || (updateSigma>(h_x_pull_update->GetBinLowEdge(100)+h_x_pull_update->GetBinWidth(100))))
	++ihateroot;

    }

    for (int j = -1; j <102; j++){
      std::cout << "bin: " << j << " counts: " << bin_counter[j] << std::endl;
    }




    std::cout <<" n_entries = " << n_entries << std::endl;
    std::cout <<" noutside, ninside = " << noutside << ", " << ninside << std::endl;
    std::cout <<" nover, nunder = " << nover << ", " << nunder << std::endl;
    std::cout << "after: " << h_x_pull_update->GetBinContent(101) << std::endl;
    std::cout << "ihr = " << ihateroot << std::endl;
    
    //    posTree->Draw("(x_init-x_update)/sqrt(x_uerr)>>h_x_pull_update");
    h_x_pull_update->Draw();
    h_x_pull_update->Fit("gaus","","",-3.0,3.0);
    c_xpos.SaveAs(Form("%s/x_pull_update.png",outDir.Data()));  
    /*
    c_xpos.cd();
    h_phi_overflow->Draw();
    c_xpos.SaveAs(Form("%s/phi_overflow.png",outDir.Data()));  
    */
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
