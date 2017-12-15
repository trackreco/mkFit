/*
  Investigating barrel, mostly pixel layers

Pixels, lay < 4
  ex_t, ey_t ~ 1 mum,    ez_t ~ 3 - 5 mum
  ex_h, ey_h ~ 0.01 mum, ez_h ~ 4 - 8 mum

Lay = 4
  ex_h ey_h ~ 5 nm - 6 mum; ez_h = 11.285
 */

#include <initializer_list>

TTree hps, hpf, hss, hsf;
TTree tps, tpf, tss, tsf;

TTree *a, *b;

TString base_cut = "chi2<500 && pt>0.5 && pt<20 && abs(phi)<7 && theta>=0";
TString inbrl_cut = "layer<7";
TString inecp_cut = "layer>17 && layer<30";

TString smerge(std::initializer_list<TString> l)
{
  auto i = l.begin();
  TString r(*i);
  while ((++i) != l.end())
  {
     r += " && ";
     r += *i;
  }
  return r;
}

// Interesting plots:
// "chi2:theta"
// "chi2:layer" per region
// "z_h:theta" per layer
// "pt:theta" per layer, more pronounced on inner pixels

void plot(TTree& t, const char *name)
{
  TCanvas *c = new TCanvas(name, name, 1200, 800);
  c->Divide(2,2);

  c->cd(1);
  t.Draw("log10(chi2)");

  c->cd(2);
  t.Draw("chi2", "chi2<500");

  c->cd(3);
  t.Draw("layer");

  c->cd(4);
  t.Draw("chi2:layer", "chi2<500", "col");
}

void plot_brl(TTree& t)
{
  TString sel = smerge({base_cut, inbrl_cut});

  TCanvas *c = new TCanvas("brl_all_lyrs", "barrel all layers", 1600, 800);
  c->Divide(3,2);
  int i = 0;

  c->cd(++i);  t.Draw("y_h:x_h", "chi2<500 && layer<10", "");
  c->cd(++i);  t.Draw("chi2:theta", sel, "col");
  c->cd(++i);  t.Draw("chi2:pt", sel, "col");
  c->cd(++i);  t.Draw("chi2:z_h", sel, "col");
  c->cd(++i);  t.Draw("pt:theta", sel, "col");
  c->cd(++i);  t.Draw("z_h:theta", sel, "col");

  c = new TCanvas("brl_per_lyr", "barrel separated by layer", 1600, 800);
  c->Divide(3,2);
  i = 0;
  
  c->cd(++i);  t.Draw("z_h:theta", sel + " && layer==0", "col");
  c->cd(++i);  t.Draw("z_h:theta", sel + " && layer==1", "col");
  c->cd(++i);  t.Draw("z_h:theta", sel + " && layer==2", "col");
  c->cd(++i);  t.Draw("z_h:theta", sel + " && layer==3", "col");
  c->cd(++i);  t.Draw("z_h:theta", sel + " && layer==4", "col");
  c->cd(++i);  t.Draw("z_h:theta", sel + " && layer==5", "col");
}

void plot_ecp(TTree& t)
{
  TString sel = smerge({base_cut, inecp_cut});

  TCanvas *c = new TCanvas("ecp_all_lyrs", "endcap pos all layers", 1600, 800);
  c->Divide(3,2);
  int i = 0;

  c->cd(++i);  t.Draw("hypot(x_h,y_h):z_h", "chi2<500 && layer>17 && layer < 24", "");
  c->cd(++i);  t.Draw("chi2:theta", sel, "col");
  c->cd(++i);  t.Draw("chi2:pt", sel, "col");
  c->cd(++i);  t.Draw("chi2:hypot(x_h,y_h)", sel, "col");
  c->cd(++i);  t.Draw("pt:theta", sel, "col");
  c->cd(++i);  t.Draw("hypot(x_h,y_h):theta", sel, "col");

  c = new TCanvas("ecp_per_lyr", "endcap pos separated by layer", 1600, 800);
  c->Divide(3,2);
  i = 0;
  
  c->cd(++i);  t.Draw("hypot(x_h,y_h):theta", sel + " && layer==18", "col");
  c->cd(++i);  t.Draw("z_h:theta", sel + " && layer==19", "col");
  c->cd(++i);  t.Draw("z_h:theta", sel + " && layer==20", "col");
  c->cd(++i);  t.Draw("z_h:theta", sel + " && layer==21", "col");
  c->cd(++i);  t.Draw("z_h:theta", sel + " && layer==22", "col");
  c->cd(++i);  t.Draw("z_h:theta", sel + " && layer==23", "col");
}

void plot_res(TTree& t)
{
  TCanvas *c = new TCanvas("hit_errors_brl", "hit errors barrel", 1200, 800);
  c->Divide(2,2);
  int i = 0;

  c->cd(++i);  t.Draw("1e4*sqrt(ex_h + ey_h)", "1e4*sqrt(ex_h + ey_h)<200 && layer < 4");
  c->cd(++i);  t.Draw("1e4*sqrt(ex_h + ey_h)", "1e4*sqrt(ex_h + ey_h)<200 && layer >= 4 && layer < 18");
  c->cd(++i);  t.Draw("1e4*sqrt(ez_h)", "1e4*sqrt(ez_h)<200 && layer < 4");

  c = new TCanvas("hit_errors_endcap", "hit errors endcap", 1200, 400);
  c->Divide(2,1);
  i = 0;

  c->cd(++i);  t.Draw("1e4*sqrt(ex_h + ey_h)", "1e4*sqrt(ex_h + ey_h)<200 && ((layer >= 18 && layer < 21) || (layer >= 45 && layer < 47))");
  c->cd(++i);  t.Draw("1e4*sqrt(ex_h + ey_h)", "1e4*sqrt(ex_h + ey_h)<200 && ((layer >= 21 && layer < 45) || (layer >= 47))");
}

void plot_res_theta(TTree& t)
{
  TCanvas *c = new TCanvas("hit_errors_brl_theta", "hit errors barrel vs theta", 1200, 800);
  c->Divide(2,2);
  int i = 0;

  c->cd(++i);  t.Draw("1e4*sqrt(ex_h + ey_h):theta", "1e4*sqrt(ex_h + ey_h)<200 && layer < 4", "col");
  c->cd(++i);  t.Draw("1e4*sqrt(ex_h + ey_h):theta", "1e4*sqrt(ex_h + ey_h)<200 && layer >= 4 && layer < 18", "col");
  c->cd(++i);  t.Draw("1e4*sqrt(ez_h):theta", "1e4*sqrt(ez_h)<200 && layer < 4", "col");

  c = new TCanvas("hit_errors_endcap_theta", "hit errors endcap vs theta", 1200, 400);
  c->Divide(2,1);
  i = 0;

  c->cd(++i);  t.Draw("1e4*sqrt(ex_h + ey_h):theta", "1e4*sqrt(ex_h + ey_h)<200 && ((layer >= 18 && layer < 21) || (layer >= 45 && layer < 47))", "col");
  c->cd(++i);  t.Draw("1e4*sqrt(ex_h + ey_h):theta", "1e4*sqrt(ex_h + ey_h)<200 && ((layer >= 21 && layer < 45) || (layer >= 47))", "col");
}

//------------------------------------------------------------------------------

void draw_w_ovrflw(TTree *t, TH1* h, const char* var, const char* sel, const char* opt, int col)
{
  TString vv;
  vv.Form("%s>>+%s", var, h->GetName());
  t->Draw(vv, sel, opt);

  int nb = h->GetNbinsX();
  h->SetBinContent(nb, h->GetBinContent(nb+1));
  h->SetLineColor(col);

  gPad->SetLogy(1);
}

//------------------------------------------------------------------------------


void chi()
{
  //hps.ReadFile("std-hit-pure.rtt"); hpf.ReadFile("wfix-hit-pure.rtt");
  //hss.ReadFile("std-hit-sim.rtt");  hsf.ReadFile("wfix-hit-sim.rtt");

  tps.ReadFile("std-trk-pure.rtt"); tpf.ReadFile("wfix-trk-pure.rtt");
  a = &tps; b = &tpf;

  // tss.ReadFile("std-trk-sim.rtt");  tsf.ReadFile("wfix-trk-sim.rtt");
  // a = &tss; b = &tsf;


  TCanvas *c = new TCanvas("chi2trk","chi2 of backward-fitted tracks",1600, 1000);
  c->Divide(3,2);
  int i = 0;

  c->cd(++i); draw_w_ovrflw(a, new TH1F("stdchi2",  "", 100, 0, 1000), "chi2", "chi2>0", "", kBlue);
  //c->cd(++i);
  draw_w_ovrflw(b, new TH1F("wfixchi2", "", 100, 0, 1000), "chi2", "chi2>0", "same", kRed);

  c->cd(++i);
  a->Draw("chi2pdof:pt", "chi2pdof > 0 && chi2pdof < 100 && pt < 10", "colz");

  c->cd(++i);
  a->Draw("chi2pdof:theta", "chi2pdof > 0 && chi2pdof < 100 && pt < 10", "colz");

  // --------------------------------
  
  c->cd(++i); draw_w_ovrflw(a, new TH1F("stdchi2pdof",  "", 100, 0, 100), "chi2pdof", "chi2pdof>0", "", kBlue);
  //c->cd(++i);
  draw_w_ovrflw(b, new TH1F("wfixchi2pdof", "", 100, 0, 100), "chi2pdof", "chi2pdof>0", "same", kRed);

  c->cd(++i);
  a->Draw("chi2pdof:pt", "chi2pdof > 0 && chi2pdof < 100 && pt < 10", "colz");

  c->cd(++i);
  a->Draw("chi2pdof:theta", "chi2pdof > 0 && chi2pdof < 100 && pt < 10", "colz");

  /*
  // plot(p, "pure");
  plot(s, "sim");

  plot_brl(s);
  plot_ecp(s);

  plot_res(s);
  plot_res_theta(s);
  */
}
