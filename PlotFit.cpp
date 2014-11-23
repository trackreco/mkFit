#include "PlotFit.hh"

PlotFit::PlotFit(TString inName, TString outName, TString outType){

  // ++++ Get Root File ++++ //

  fInRoot = TFile::Open(Form("%s",inName.Data()));

  // ++++ Define Output Parameters, Make Directory/File ++++ //

  fOutName = outName;
  fOutType = outType;
  
  FileStat_t dummyFileStat;
    
  if (gSystem->GetPathInfo(fOutName.Data(), dummyFileStat) == 1){
    TString mkDir = "mkdir -p ";
    mkDir += fOutName.Data();
    gSystem->Exec(mkDir.Data());
  }

  fOutRoot = new TFile(Form("%s/%s.root",fOutName.Data(),fOutName.Data()), "RECREATE");

  // General style

  gROOT->Reset();
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1011);

  fTH1Canv = new TCanvas();
  fTH2Canv = new TCanvas();
}

PlotFit::~PlotFit(){
  fInRoot->Delete();
  fOutRoot->Delete();
  fTH1Canv->Close();
  fTH2Canv->Close();
}

void PlotFit::PlotAllHistos(){
  PlotFit::PlotPosResPull();
  PlotFit::PlotOverFlow();
  PlotFit::PlotPtResPull();
  PlotFit::PlotGeo();
}

void PlotFit::PlotPosResPull(){

  // Loop over entries in tree first
  // Then loop over x,y,z,phi
  // Then MC,prop,update to make pull/res plots

  // Get trees
  TTree * posTree  = (TTree*)fInRoot->Get("posTree");

  //Declare strings for branches and plots
  TString pos[4]   = {"x","y","z","phi"};
  TString spos[4]  = {"x","y","z","#phi"};
  TString step[3]  = {"mc","prop","update"};
  TString sstep[3] = {"MC","Prop","Update"};
  TString estep[3] = {"mc","p","u"}; 

  // Floats to be filled for trees
  Float_t init_val[4];
  Float_t step_val[4][3];
  Float_t step_err[4][3];
  Float_t step_out[2] = {0.,0.}; // res/pull output
  
  // Create pos plots
  TH1F * posResPlot[4][3];
  TH1F * posPullPlot[4][3];
  for (UInt_t i = 0; i < 4; i++){ // loop over pos index
    for (UInt_t j = 0; j < 3; j++){ // loop over step index
      //Res
      posResPlot[i][j] = new TH1F(Form("h_%s_res_%s",pos[i].Data(),step[j].Data()),Form("%s Resolution (%s Hit vs. Initial Hit)",spos[i].Data(),sstep[j].Data()),100,-0.01,0.01);
      posResPlot[i][j]->GetXaxis()->SetTitle(Form("(%s^{%s} - %s^{init})/%s^{init}",spos[i].Data(),step[j].Data(),spos[i].Data(),spos[i].Data()));
      posResPlot[i][j]->GetYaxis()->SetTitle("Hits");
      
      //Pull
      posPullPlot[i][j] = new TH1F(Form("h_%s_pull_%s",pos[i].Data(),step[j].Data()),Form("%s Pull (%s Hit vs. Initial Hit)",spos[i].Data(),sstep[j].Data()),100,-10.0,10.0);
      posPullPlot[i][j]->GetXaxis()->SetTitle(Form("(%s^{%s} - %s^{init})/#sigma(%s^{%s})",spos[i].Data(),step[j].Data(),spos[i].Data(),spos[i].Data(),sstep[j].Data()));
      posPullPlot[i][j]->GetYaxis()->SetTitle("Hits");    
    }
  }

  //Initialize step_val arrays, SetBranchAddress
  for (UInt_t i = 0; i < 4; i++){ // loop over pos index
    //Initialize init val array
    init_val[i] = 0.;
    
    //Set init_val branch
    posTree->SetBranchAddress(Form("%s_init",pos[i].Data()),&(init_val[i]));
    for (UInt_t j = 0; j < 3; j++){ // loop over step index
      //Step value + errors arrays  
      step_val[i][j] = 0.;
      step_err[i][j] = 0.;
      
      //Set step_val/err branch
      posTree->SetBranchAddress(Form("%s_%s",pos[i].Data(),step[j].Data()),&(step_val[i][j]));
      posTree->SetBranchAddress(Form("%s_%serr",pos[i].Data(),estep[j].Data()),&(step_err[i][j]));
    }
  }

  // Fill histos, compute res/pull from tree branches 
  for (UInt_t k = 0; k < (UInt_t)posTree->GetEntries(); k++){
    posTree->GetEntry(k);
    for (UInt_t i = 0; i < 4; i++){  // loop over pos index
      for (UInt_t j = 0; j < 3; j++){ // loop over step index
	PlotFit::ComputeResPull(init_val[i],step_val[i][j],step_err[i][j],step_out);
	if (!isnan(step_out[0])){ // fill if not nan
	  posResPlot[i][j]->Fill(step_out[0]);
	}
	if (!isnan(step_out[1])){ // fill if not nan
	  posPullPlot[i][j]->Fill(step_out[1]);
	}
      } // end loop over steps
    } // end loop over pos
  } // end loop over entry in tree

  // Draw, fit, and save plots
  for (UInt_t i = 0; i < 4; i++){
    for (UInt_t j = 0; j < 3; j++){
      PlotFit::DrawFitSaveTH1Plot(posResPlot[i][j],0.005,Form("%s_res_%s",pos[i].Data(),step[j].Data()));
      PlotFit::DrawFitSaveTH1Plot(posPullPlot[i][j],3.0,Form("%s_pull_%s",pos[i].Data(),step[j].Data()));
    }
  }  
  
  posTree->Delete();
}

void PlotFit::PlotOverFlow(){

  // Difference here is after calculating x_pull at any step, save phi value
  // Can make overflow plots general...but just restricting them to phi only

  // Get trees
  TTree * posTree  = (TTree*)fInRoot->Get("posTree");

  //Strings for branches and plots
  TString step[3]   = {"mc","prop","update"};
  TString sstep[3]  = {"MC","Prop","Update"};
  TString estep[3]  = {"mc","p","u"}; 

  // Floats to be filled for trees
  Float_t pos_init_val = 0.;
  Float_t pos_step_val[3];
  Float_t pos_step_err[3];
  Float_t pos_step_out[2] = {0.,0.};

  Float_t phi_init_val = 0.;
  Float_t phi_step_val[3];

  // Create overflow plots
  TH1F * posOverFlow_stepPhi[3];
  TH1F * posOverFlow_initPhi[3];
  for (UInt_t j = 0; j < 3; j++){ // loop over step index
    // Plots to save computed step phi value
    posOverFlow_stepPhi[j] = new TH1F(Form("h_%s_phi_overflow_%s",step[j].Data(),step[j].Data()),Form("%s #phi for x %s Pull |#sigma| > 10",sstep[j].Data(),sstep[j].Data()),67,-3.35,3.35);
    posOverFlow_stepPhi[j]->GetXaxis()->SetTitle(Form("#phi^{%s}",sstep[j].Data())); 
    posOverFlow_stepPhi[j]->GetYaxis()->SetTitle("Hits");
    
    // Plots to save generated phi value, no smear
    posOverFlow_initPhi[j] = new TH1F(Form("h_init_phi_overflow_%s",step[j].Data()),Form("Initial #phi for x %s Pull |#sigma| > 10",sstep[j].Data()),67,-3.35,3.35);
    posOverFlow_initPhi[j]->GetXaxis()->SetTitle("#phi^{init}"); 
    posOverFlow_initPhi[j]->GetYaxis()->SetTitle("Hits");
  }

  //Initialize step_val arrays, SetBranchAddress
  //Set init_val branches
  posTree->SetBranchAddress("x_init",&pos_init_val);
  posTree->SetBranchAddress("phi_init",&phi_init_val);
    
  //Set step_val branches
  for (UInt_t j = 0; j < 3; j++){ // loop over step index
    //Step value + errors arrays (pos) initialize + set address
    pos_step_val[j] = 0.;
    pos_step_err[j] = 0.;
    
    posTree->SetBranchAddress(Form("x_%s",step[j].Data()),&(pos_step_val[j]));
    posTree->SetBranchAddress(Form("x_%serr",estep[j].Data()),&(pos_step_err[j]));

    //Step val for phi array initialize + set address
    phi_step_val[j] = 0.;
    posTree->SetBranchAddress(Form("phi_%s",step[j].Data()),&(phi_step_val[j]));
  }

  // Fill histos, compute res/pull from tree branches 
  for (UInt_t k = 0; k < (UInt_t)posTree->GetEntries(); k++){
    posTree->GetEntry(k);
    for (UInt_t j = 0; j < 3; j++){ // loop over step index
      PlotFit::ComputeResPull(pos_init_val,pos_step_val[j],pos_step_err[j],pos_step_out);
      if ( (std::abs(pos_step_out[1]) > 10.) && (!isnan(pos_step_out[1])) ){
	posOverFlow_stepPhi[j]->Fill(phi_step_val[j]);
	posOverFlow_initPhi[j]->Fill(phi_init_val);
      }
    } // end loop over steps
  } // end loop over entry in tree
 
  // Draw and save plots
  for (UInt_t j = 0; j < 3; j++){
    PlotFit::DrawSaveTH1Plot(posOverFlow_stepPhi[j],Form("%s_phi_%s_xPullOverFlow",step[j].Data(),step[j].Data()));
    PlotFit::DrawSaveTH1Plot(posOverFlow_initPhi[j],Form("init_phi_%s_xPullOverFlow",step[j].Data()));
  }  
  posTree->Delete();
}

void PlotFit::PlotPtResPull(){

  // Get the trees
  TTree *  ptTree = (TTree*)fInRoot->Get("fittree");

  // Declare/initialize branch variables
  Float_t init_val = 0.;
  Float_t step_val = 0.;
  Float_t step_err = 0.;
  Float_t step_out[2] = {0.,0.};
  
  // Create pt plots
  //Res
  TH1F * ptResPlot = new TH1F("h_pt_res_fit","p_{T} Resolution (MC-Fit Track)",100,-0.5,0.5);
  ptResPlot->GetXaxis()->SetTitle("(p_{T}^{MC} - p_{T}^{fit})/p_{T}^{MC}");
  ptResPlot->GetYaxis()->SetTitle("Tracks");
      
  //Pull
  TH1F * ptPullPlot = new TH1F("h_pt_pull_fit","p_{T} Pull (MC-Fit Track)",100,-10.0,10.0);
  ptPullPlot->GetXaxis()->SetTitle("(p_{T}^{MC} - p_{T}^{fit})/#sigma(p_{T}^{fit})");
  ptPullPlot->GetYaxis()->SetTitle("Tracks");

  //Set branch addresses for pt plots
  ptTree->SetBranchAddress("pt_mc",&init_val);
  ptTree->SetBranchAddress("pt_fit",&step_val);
  ptTree->SetBranchAddress("pt_err",&step_err);
  
  for (UInt_t k = 0; k < (UInt_t)ptTree->GetEntries(); k++){
    ptTree->GetEntry(k);
    PlotFit::ComputeResPull(init_val,step_val,(step_err*step_err),step_out); // step_err is already in the right units, but need it to match the function call
    if (!isnan(step_out[0])){ // fill if not nan
      ptResPlot->Fill(step_out[0]);
    }
    if (!isnan(step_out[1])){ // fill if not nan
      ptPullPlot->Fill(step_out[1]);
    }
  } // end loop over entry in tree

  // Draw, fit, and save plots
  PlotFit::DrawFitSaveTH1Plot(ptResPlot,0.3,"pt_res_fit");
  PlotFit::DrawFitSaveTH1Plot(ptPullPlot,3.0,"pt_pull_fit");

  ptTree->Delete();
}

void PlotFit::PlotGeo(){

  // Get trees
  TTree * posTree  = (TTree*)fInRoot->Get("posTree");

  // Declare strings for plots and branches
  TString coord[2][2]  = {{"x","y"},{"z","r"}}; 
  TString scoord[2][2] = {{"x","y"},{"z","Radius"}}; 
  TString step[3]      = {"init","mc","update"};
  TString sstep[3]     = {"Initial","MC","Update"};
  
  // Declare step_val for all sets, last two are each x,y
  Float_t step_val[2][3][2];   
    
  // Settings for geoPlots
  Int_t nBins[2][2]      = {{840,840}, {2200,420}};
  Float_t xRange[2][2]   = {{-42.0,42.0},{-110.0,110.0}};
  Float_t yRange[2][2]   = {{-42.0,42.0},{0.0,42.0}};
  Int_t pixelRange[2][2] = {{2000,2000},{2000,1500}};

  // Create geo plots
  TH2F * geoPlot[2][3];
  for (UInt_t i = 0; i < 2; i++){ // loop over coord index
    for (UInt_t j = 0; j < 3; j++){ // loop over step index
      //Geo
      geoPlot[i][j] = new TH2F(Form("h_%s_vs_%s_%s",coord[i][0].Data(),coord[i][1].Data(),step[j].Data()),Form("%s vs %s %s",scoord[i][0].Data(),scoord[i][1].Data(),sstep[j].Data()),nBins[i][0],xRange[i][0],xRange[i][1],nBins[i][1],yRange[i][0],yRange[i][1]);
      geoPlot[i][j]->GetXaxis()->SetTitle(Form("%s^{%s} (cm)",scoord[i][0].Data(),step[j].Data()));
      geoPlot[i][j]->GetYaxis()->SetTitle(Form("%s^{%s} (cm)",scoord[i][1].Data(),step[j].Data()));
    }
  }

  //Initialize step_vals, set branch addresses
  for (UInt_t i = 0; i < 2; i++){ // loop over coord index
    for (UInt_t j = 0; j < 3; j++){ // loop over step index
      step_val[i][j][0] = 0.;
      step_val[i][j][1] = 0.;

      //Set Branch Address
      posTree->SetBranchAddress(Form("%s_%s",coord[i][0].Data(),step[j].Data()),&(step_val[i][j][0]));
      posTree->SetBranchAddress(Form("%s_%s",coord[i][1].Data(),step[j].Data()),&(step_val[i][j][1]));
    }
  }
      
  // Get values, fill geo plot
  for (UInt_t k = 0; k < (UInt_t)posTree->GetEntries(); k++){
    posTree->GetEntry(k);
    for (UInt_t i = 0; i < 2; i++){  // loop over coord index
      for (UInt_t j = 0; j < 3; j++){ // loop over step index
	geoPlot[i][j]->Fill(step_val[i][j][0],step_val[i][j][1]);
      }
    }
  }
  
  //Save Geo Plots
  for (UInt_t i = 0; i < 2; i++){ // loop over coord index
    for (UInt_t j = 0; j < 3; j++){ // loop over step index
      //Geo
      PlotFit::DrawSaveTH2Plot(geoPlot[i][j],Form("%s_vs_%s_%s",coord[i][0].Data(),coord[i][1].Data(),step[j].Data()),pixelRange[i]);
    }
  }
  posTree->Delete();
}

void PlotFit::ComputeResPull(const Float_t& init_val, const Float_t& step_val, const Float_t& step_err, Float_t step_out[]){
  step_out[0] = (step_val - init_val)/init_val;
  step_out[1] = (step_val - init_val)/sqrt(step_err);
}

void PlotFit::DrawFitSaveTH1Plot(TH1F * hist, Float_t fitRange, TString plotName){
  fTH1Canv->cd();
  hist->Draw();
  fOutRoot->cd();
  hist->Write();
  fTH1Canv->cd();
  hist->Fit("gaus","","",-fitRange,fitRange);
  fTH1Canv->SaveAs(Form("%s/%s_%s.%s",fOutName.Data(),plotName.Data(),fOutName.Data(),fOutType.Data()));  
}

void PlotFit::DrawSaveTH1Plot(TH1F * hist, TString plotName){
  fTH1Canv->cd();
  hist->Draw();
  fOutRoot->cd();
  hist->Write();
  fTH1Canv->cd();
  fTH1Canv->SaveAs(Form("%s/%s_%s.%s",fOutName.Data(),plotName.Data(),fOutName.Data(),fOutType.Data()));  
}

void PlotFit::DrawSaveTH2Plot(TH2F * hist, TString plotName, Int_t pixelRange[]){
  fTH2Canv->cd();
  fTH2Canv->SetCanvasSize(pixelRange[0],pixelRange[1]);
  hist->SetStats(0);
  hist->Draw();
  fOutRoot->cd();
  hist->Write();
  fTH2Canv->cd();
  fTH2Canv->SaveAs(Form("%s/%s_%s.%s",fOutName.Data(),plotName.Data(),fOutName.Data(),fOutType.Data()));  
}


