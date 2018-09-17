### 
#Prior running, make sure to have sourced the necessary environment:
#source xeon_scripts/common-variables.sh
#source xeon_scripts/init-env.sh
#source val_scripts/init-root.sh
###

import os,sys
import ROOT
import copy

def getCanvasMainPad( logY ):
  pad1 = ROOT.TPad("pad1", "pad1", 0, 0.2, 1, 1)
  pad1.SetBottomMargin(0.15)
  if( logY ):
    pad1.SetLogy()
  return pad1

def getCanvasRatioPad( logY ):
  pad2 = ROOT.TPad("pad2", "pad2", 0, 0, 1, 0.21)
  pad2.SetTopMargin(0.05)
  pad2.SetBottomMargin(0.1)
  return pad2

def getRatioAxes( xMin, xMax, yMin, yMax ):
  h2_axes_ratio = ROOT.TH2D("axes_ratio", "", 10, xMin, xMax, 10, yMin, yMax )
  h2_axes_ratio.SetStats(0)
  h2_axes_ratio.GetXaxis().SetLabelSize(0.00)
  h2_axes_ratio.GetXaxis().SetTickLength(0.09)
  h2_axes_ratio.GetYaxis().SetNdivisions(5,5,0)
  h2_axes_ratio.GetYaxis().SetTitleSize(0.13)
  h2_axes_ratio.GetYaxis().SetTitleOffset(0.37)
  h2_axes_ratio.GetYaxis().SetLabelSize(0.13)
  h2_axes_ratio.GetYaxis().SetTitle("Ratio")
  return h2_axes_ratio


if(len(sys.argv)!=8):
  print "python [-i] compareTotals.py <original-dir (containing original validation dir)> <updated-dir (containing updated validation dir)> <validation-dir> <output-dir> <what-to-compare [eff;fr;dr;]> <validation-type [sim;cmssw;]> <vsvar [eta;pt;]>"
  exit()

dirnames = []
dirnames.append(sys.argv[1]) #original (1)
dirnames.append(sys.argv[2]) #updated (2)

stddir = sys.argv[3]

outdir = sys.argv[4] 

if not os.path.exists(outdir):
    os.mkdir(outdir)

iseff = False
isfr  = False
isdr  = False
whattot = sys.argv[5]
if(whattot=="eff"):
  iseff=True
elif(whattot=="fr"):
  isfr =True
elif(whattot=="dr"):
  isdr =True
else:
  print "What do you want to compare?"
  print "python [-i] compareTotals.py <original-dir (containing original validation dir)> <updated-dir (containing updated validation dir)> <validation-dir> <output-dir> <what-to-compare [eff;fr;dr;]> <validation-type [sim;cmssw;]> <vsvar [eta;pt;]>"
  exit()  

issim   = False
iscmssw = False
valtype = sys.argv[6]
if(valtype=="sim"):
  issim   =True
elif(valtype=="cmssw"):
  iscmssw =True
else:
  print "SIM or CMSSW validaiton?"
  print "python [-i] compareTotals.py <original-dir (containing original validation dir)> <updated-dir (containing updated validation dir)> <validation-dir> <output-dir> <what-to-compare [eff;fr;dr;]> <validation-type [sim;cmssw;]> <vsvar [eta;pt;]>"
  exit()  

vseta = False
vspt  = False
vsvar = sys.argv[7]
if(vsvar=="eta"):
  vseta = True
elif(vsvar=="pt"):
  vspt  = True
else:
  print "vs eta or pT?"
  print "python [-i] compareTotals.py <original-dir (containing original validation dir)> <updated-dir (containing updated validation dir)> <validation-dir> <output-dir> <what-to-compare [eff;fr;dr;]> <validation-type [sim;cmssw;]> <vsvar [eta;pt;]>"
  exit()  
  
names = ["Original", "Updated"]
colors = [1,2]
titles = ["p_{T}(MC)>0 GeV", "p_{T}(MC)>0.9 GeV", "p_{T}(MC)>2.0 GeV"]
#titles = ["0.0 < p_{T}(MC) < 0.9 GeV", "0.9 < p_{T}(MC) < 2.0 GeV", "p_{T}(MC)>2.0 GeV"]

fnames = []
fs = []
for d in range(0, len(dirnames)):
    fnames.append(dirnames[d]+"/"+stddir+"/plots.root")
    fs.append(ROOT.TFile.Open(fnames[d]))

eff = []
eff_pass = []
eff_tot = []
eff_ratio = []

toget=[]
if(iseff):
  if(issim): 
    if(vseta):
      toget.append("efficiency/eff_sim_eta_build_pt0.0")
      toget.append("efficiency/eff_sim_eta_build_pt0.9")
      toget.append("efficiency/eff_sim_eta_build_pt2.0")
    else:
      toget.append("efficiency/eff_sim_pt_build_pt0.0")
      toget.append("efficiency/eff_sim_pt_build_pt0.9")
      toget.append("efficiency/eff_sim_pt_build_pt2.0")
  else:
    if(vseta):
      toget.append("efficiency_cmssw/eff_cmssw_eta_build_pt0.0")
      toget.append("efficiency_cmssw/eff_cmssw_eta_build_pt0.9")
      toget.append("efficiency_cmssw/eff_cmssw_eta_build_pt2.0")    
    else:
      toget.append("efficiency_cmssw/eff_cmssw_pt_build_pt0.0")
      toget.append("efficiency_cmssw/eff_cmssw_pt_build_pt0.9")
      toget.append("efficiency_cmssw/eff_cmssw_pt_build_pt2.0")    

elif(isdr):
  if(issim): 
    if(vseta):
      toget.append("duplicaterate/dr_sim_eta_build_pt0.0")
      toget.append("duplicaterate/dr_sim_eta_build_pt0.9")
      toget.append("duplicaterate/dr_sim_eta_build_pt2.0")
    else:
      toget.append("duplicaterate/dr_sim_pt_build_pt0.0")
      toget.append("duplicaterate/dr_sim_pt_build_pt0.9")
      toget.append("duplicaterate/dr_sim_pt_build_pt2.0")
  else:
    if(vseta):
      toget.append("duplicaterate_cmssw/dr_cmssw_eta_build_pt0.0")
      toget.append("duplicaterate_cmssw/dr_cmssw_eta_build_pt0.9")
      toget.append("duplicaterate_cmssw/dr_cmssw_eta_build_pt2.0")    
    else:
      toget.append("duplicaterate_cmssw/dr_cmssw_pt_build_pt0.0")
      toget.append("duplicaterate_cmssw/dr_cmssw_pt_build_pt0.9")
      toget.append("duplicaterate_cmssw/dr_cmssw_pt_build_pt2.0")    

elif(isfr):
  if(issim): 
    if(vseta):
      toget.append("fakerate/fr_reco_eta_build_pt0.0")
      toget.append("fakerate/fr_reco_eta_build_pt0.9")
      toget.append("fakerate/fr_reco_eta_build_pt2.0")
    else:
      toget.append("fakerate/fr_reco_pt_build_pt0.0")
      toget.append("fakerate/fr_reco_pt_build_pt0.9")
      toget.append("fakerate/fr_reco_pt_build_pt2.0")
  else:
    if(vseta):
      toget.append("fakerate_cmssw/fr_reco_eta_build_pt0.0")
      toget.append("fakerate_cmssw/fr_reco_eta_build_pt0.9")
      toget.append("fakerate_cmssw/fr_reco_eta_build_pt2.0")    
    else:
      toget.append("fakerate_cmssw/fr_reco_pt_build_pt0.0")
      toget.append("fakerate_cmssw/fr_reco_pt_build_pt0.9")
      toget.append("fakerate_cmssw/fr_reco_pt_build_pt2.0")    


for d in range(0, len(fnames)):
    thiseff = []
    thiseff_pass = []
    thiseff_tot = []
    aux_ratio = []
    thiseff_ratio = []

    fs[d].cd()

    thiseff.append(fs[d].Get(toget[0]))
    thiseff_pass.append(thiseff[0].GetCopyPassedHisto())
    thiseff_tot.append(thiseff[0].GetCopyTotalHisto())
    aux_ratio = copy.deepcopy(thiseff_pass)
    aux_ratio[0].Sumw2()
    aux_ratio[0].Divide(thiseff_tot[0])
    thiseff_ratio.append(aux_ratio[0])
    thiseff_ratio[0].SetLineColor(colors[d])
    thiseff_ratio[0].SetMarkerColor(colors[d])
    thiseff_ratio[0].SetMarkerSize(0.3)
    thiseff_ratio[0].SetMarkerStyle(20)
    thiseff_ratio[0].SetStats(0)
    thiseff_ratio[0].SetTitle(titles[0])

    thiseff.append(fs[d].Get(toget[1]))
    thiseff_pass.append(thiseff[1].GetCopyPassedHisto())
    thiseff_tot.append(thiseff[1].GetCopyTotalHisto())
    aux_ratio = copy.deepcopy(thiseff_pass)
    aux_ratio[1].Sumw2()
    aux_ratio[1].Divide(thiseff_tot[1])
    thiseff_ratio.append(aux_ratio[1])
    thiseff_ratio[1].SetLineColor(colors[d])
    thiseff_ratio[1].SetMarkerColor(colors[d])
    thiseff_ratio[1].SetMarkerStyle(20)
    thiseff_ratio[1].SetMarkerSize(0.3)
    thiseff_ratio[1].SetStats(0)
    thiseff_ratio[1].SetTitle(titles[1])

    thiseff.append(fs[d].Get(toget[2]))
    thiseff_pass.append(thiseff[2].GetCopyPassedHisto())
    thiseff_tot.append(thiseff[2].GetCopyTotalHisto())
    aux_ratio = copy.deepcopy(thiseff_pass)
    aux_ratio[2].Sumw2()
    aux_ratio[2].Divide(thiseff_tot[2])
    thiseff_ratio.append(aux_ratio[2])
    thiseff_ratio[2].SetLineColor(colors[d])
    thiseff_ratio[2].SetMarkerColor(colors[d])
    thiseff_ratio[2].SetMarkerSize(0)
    thiseff_ratio[2].SetStats(0)
    thiseff_ratio[2].SetTitle(titles[2])

    fs[d].Close()

    eff.append(thiseff)
    eff_pass.append(thiseff_pass)
    eff_tot.append(thiseff_tot)
    eff_ratio.append(thiseff_ratio)

eff_ratio_xratio = copy.deepcopy(eff_ratio)
eff_ratio_xratio_denom = copy.deepcopy(eff_ratio)

ratios_noclean = []
for d in range(0, len(fnames)):
    thisratios_noclean = []
    for p in range(0,3):
        for b in range(1,eff_ratio_xratio_denom[d][p].GetNbinsX()+1):
            eff_ratio_xratio_denom[d][p].SetBinError(b,0.0)
        thisratios_noclean.append(eff_ratio_xratio[d][p])
        thisratios_noclean[p].Divide(eff_ratio_xratio_denom[0][p])
        thisratios_noclean[p].GetYaxis().SetTitle("Updated / Original")
        thisratios_noclean[p].SetLineColor(colors[d])
        thisratios_noclean[p].SetMarkerColor(colors[d])
        thisratios_noclean[p].SetMarkerSize(0)
        thisratios_noclean[p].SetStats(0)
        ### Suppress y-error bars for ratio
        for b in range(1,thisratios_noclean[p].GetNbinsX()+1):
            thisratios_noclean[p].SetBinError(b,1e-7)
        
    ratios_noclean.append(thisratios_noclean)


### Drawing
legend = ROOT.TLegend(0.7,0.7, 0.87, 0.87);
legend.SetLineColor(0)
legend.SetFillColor(0)
for d in range(0,len(fnames)):
    legend.AddEntry(eff_ratio[d][0], names[d], "L")

ROOT.gStyle.SetOptStat(0)

outname=[outdir+"/"+whattot+"_"+valtype+"_vs"+vsvar+"_pT0.png",outdir+"/"+whattot+"_"+valtype+"_vs"+vsvar+"_pT0p9.png",outdir+"/"+whattot+"_"+valtype+"_vs"+vsvar+"_pT2p0.png"]
for p in range(0,3):
  can = ROOT.TCanvas(titles[p], "", 600, 600)
  can.cd()
  
  pad1 = getCanvasMainPad(0)
  pad1.SetTickx()
  pad1.SetTicky()
  
  pad2 = getCanvasRatioPad(0)
  pad2.SetTickx()
  pad2.SetTicky()
  
  can.cd()
  pad1.Draw()
  pad1.cd()
  
  eff_ratio[0][p].GetYaxis().SetRangeUser(0.0, 1.5)
  eff_ratio[0][p].GetYaxis().SetTitleOffset(1.4)
  eff_ratio[0][p].Draw("PE")
  eff_ratio[1][p].Draw("PE,same")
  
  legend.Draw("same")
  
  can.cd()
  pad2.Draw()
  pad2.cd()
  
  xmin = ratios_noclean[0][p].GetXaxis().GetBinLowEdge(1)
  xmax = ratios_noclean[0][p].GetXaxis().GetBinLowEdge(ratios_noclean[0][p].GetNbinsX())+ratios_noclean[0][p].GetXaxis().GetBinWidth(ratios_noclean[0][p].GetNbinsX())
  ymin = 0.90
  ymax = 1.10
  
  hraxes = getRatioAxes(xmin,xmax,ymin,ymax)
  
  hraxes.Draw("")
  ratios_noclean[0][p].Draw("PE,same")
  ratios_noclean[1][p].Draw("PE,same")
  
  can.cd()
  pad1.Draw()
  pad2.Draw()
  
  can.SaveAs(outname[p]);
  
  del hraxes
  del pad1
  del pad2
  del can

  tot = [eff_tot[0][p].Integral(), eff_tot[1][p].Integral()]
  passing = [eff_pass[0][p].Integral(), eff_pass[1][p].Integral()]
  efficiency = []
  relloss = []
  for d in range(0,len(tot)):
    efficiency.append(passing[d]/tot[d])
    relloss.append(efficiency[d]/efficiency[0])

  print titles[p]
  print "Totals:"
  for d in range(0,len(tot)):
    print int(tot[d]),
  print "\nPassing:"
  for d in range(0,len(tot)):
    print int(passing[d]),
  print "\nEfficiency:"
  for d in range(0,len(tot)):
    print "%0.4f" % efficiency[d],
  print "\nRelative loss (/NO):"
  for d in range(0,len(tot)):
    print "%0.4f" % relloss[d],
  print "\n"
