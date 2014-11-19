void runOverplot(){

  gStyle->SetOptStat(0);
  gStyle->SetFillColor(0);

  // Input parameters
  TString inPlotName = "";
  TString name_type[2]  = {"h_rec_eff_vs_Sim","h_rec_fake_vs_Rec"};
  TString name_var[3]   = {"Pt","Phi","Eta"};
  TString name_assoc[2]  = {"_RD","_SD"};

  // open files
  TFile * file0[10];
  for (int nhits = 3; nhits <= 10; nhits++){
    file0[nhits] = TFile::Open(Form("eff_nhits%i/eff_nhits%i.root",nhits,nhits));
  }

  // Output parameters
  Int_t fColors[11]={1,600,632,418,880,401,433,807,420,884,793};
  TLegend * leg[2][3][2];
  TCanvas * c[2][3][2];

  TH1F * overplot[2][3][2];
  TString outPlotName  = "";
  TString outputName_type[2] = {"Eff_vs_Sim","Fake_vs_Rec"};  

  // make out dir
  TString outDir = "overplot_nHits";
  FileStat_t dummyFileStat;
    
  if (gSystem->GetPathInfo(outDir.Data(), dummyFileStat) == 1){
    TString mkDir = "mkdir -p ";
    mkDir += outDir.Data();
    gSystem->Exec(mkDir.Data());
  }

  // Do the filling

  for (UInt_t i = 0; i < 2; i++){ // loop over efficiency then fake rate
    std::cout << "i:" << i << std::endl;
    for (UInt_t j = 0; j < 3; j++){ // loop over track variable
      std::cout << "j:" << j << std::endl;
      for (UInt_t k = 0; k < 2; k++){ // loop over associator
	std::cout << "k:" << k << std::endl;
	c[i][j][k] = new TCanvas();
	leg[i][j][k] = new TLegend(.95, .7, 1.0, .9);
	overplot[i][j][k] = new TH1F();

	c[i][j][k]->cd();

	// Make one long string -- reset each time
	inPlotName = "";
	inPlotName+=name_type[i];
	inPlotName+=name_var[j];
	inPlotName+=name_assoc[k];

	outPlotName = "";	
	outPlotName+=outputName_type[i];
	outPlotName+=name_var[j];
	outPlotName+=name_assoc[k];

	std::cout << inPlotName.Data() << std::endl;

      	for (int nhits = 3; nhits <= 10; nhits++){ // loop over nhits 
	  std::cout << "nhits: " << nhits << std::endl;

	  overplot[i][j][k] = (TH1F*) file0[nhits]->Get(Form("%s",inPlotName.Data()));
	  overplot[i][j][k]->SetLineColor(fColors[nhits]);
	  overplot[i][j][k]->SetMarkerColor(fColors[nhits]);
	  overplot[i][j][k]->Draw((nhits>3)?"SAME":"");
	
	  //PUT IN COLORs
	  leg[i][j][k]->AddEntry(overplot[i][j][k],Form("%i",nhits),"l");

	} // END LOOP OVER NHITS
	c[i][j][k]->cd();
	leg[i][j][k]->Draw();
	c[i][j][k]->SaveAs(Form("%s/%s_overplot_nhits.pdf",outDir.Data(),outPlotName.Data()));  	
      }
    }
  }

}


