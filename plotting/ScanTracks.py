print "script instructions:  ",
print "python scanTracks.py path_and_file1 path_and_file2 output_text_file (build / fit) (outdir)"

from ROOT import *
import sys,os,copy

gStyle.SetOptStat(0000)

path_and_file1=sys.argv[1]
path_and_file2=sys.argv[2]
ofilename=sys.argv[3]

suffix="build"
try:
    suffix=sys.argv[4]
except:
    pass

outdir="."
try:
    outdir=sys.argv[5]
except:
    pass

if not os.path.exists(outdir):
    os.mkdir(outdir)

print "running: python scanTracks.py", sys.argv[1], sys.argv[2], sys.argv[3], suffix, outdir

a=TFile.Open(path_and_file1)
tree_ref=a.Get("efftree")

b=TFile.Open(path_and_file2)
tree_compare=b.Get("efftree")

tree_compare.BuildIndex("evtID","mcID")
tree_ref.AddFriend(tree_compare,"compare")

print "check mc_pt match", tree_ref.GetEntries("pt_mc_gen-compare.pt_mc_gen!=0")
print "check mc_phi match", tree_ref.GetEntries("phi_mc_gen-compare.phi_mc_gen!=0")

print "check evt match",  tree_ref.GetEntries("evtID-compare.evtID!=0")
print "check mcID match",  tree_ref.GetEntries("mcID-compare.mcID!=0")

#masks=["mcmask_seed","mcmask_build","mcmask_fit",
#"mcTSmask_seed","mcTSmask_build","mcTSmask_fit"]

#kinematic_vars=["xhit","yhit","zhit","pt","phi","eta",
#"nHits","nLayers","nHitsMatched","fracHitsMatched",
#"lastlyr","dphi","hitchi2","score","helixchi2"]

mask="mcmask"
cutstring_all="("+mask+"_"+suffix+">=0)&&(pt_mc_gen-compare.pt_mc_gen==0)"
cutstring1="("+mask+"_"+suffix+"==1)&&(compare."+mask+"_"+suffix+"==0)&&(pt_mc_gen-compare.pt_mc_gen==0)"
cutstring2="("+mask+"_"+suffix+"==0)&&(compare."+mask+"_"+suffix+"==1)&&(pt_mc_gen-compare.pt_mc_gen==0)"

print cutstring1, "--> eff in ref - not in compare tree"
print cutstring2, "--> eff in compare - not in ref tree"

print "tracks in ref tree- not in compare tree:", tree_ref.GetEntries(cutstring1)
print "tracks in compare tree- not in ref tree:", tree_ref.GetEntries(cutstring2)

scanstring="evtID:mcID:seedID_seed:seedID_"+suffix+":compare.seedID_seed:compare.seedID_"+suffix+":pt_mc_gen:eta_mc_gen:phi_mc_gen:compare.pt_mc_gen:compare.eta_mc_gen:compare.phi_mc_gen:pt_"+suffix+":eta_"+suffix+":phi_"+suffix+":compare.pt_"+suffix+":compare.eta_"+suffix+":compare.phi_"+suffix+":fracHitsMatched_"+suffix+":compare.fracHitsMatched_"+suffix

print "scan these variables: ", scanstring

tree_ref.GetPlayer().SetScanRedirect(kTRUE)

tree_ref.GetPlayer().SetScanFileName(ofilename+"ref_no_compare.txt")
tree_ref.Scan(scanstring, cutstring1)
tree_ref.GetPlayer().SetScanFileName(ofilename+"compare_no_ref.txt")
tree_ref.Scan(scanstring, cutstring2)

c=TCanvas("c","c")

def make_scanplots(histo1, histo2, varname, varfilename, cutstring1, cutstring2, logx=False, logy=False): # canvas and tree not passed!
    
    histo1.SetLineWidth(2)
    histo2.SetLineWidth(2)
    histo1.SetLineColor(kRed)
    histo2.SetLineColor(kGreen)
    
    tree_ref.Draw(varname+">>"+histo1.GetName(), cutstring1)
    tree_ref.Draw(varname+">>"+histo2.GetName(), cutstring2)
        
    if logy: 
        histo1.GetYaxis().SetRangeUser(0.001,1.5*max(histo1.GetMaximum(),histo2.GetMaximum()))
    else:
        histo1.GetYaxis().SetRangeUser(0,1.5*max(histo1.GetMaximum(),histo2.GetMaximum()))
    
    
    legend=TLegend(0.7,0.75,0.89,0.88)
    legend.AddEntry(histo1,"ref_no_compare")
    legend.AddEntry(histo2,"compare_no_ref")

    if logx: c.SetLogx()
    if logy: c.SetLogy()
    histo1.Draw()
    histo2.Draw("same")
    legend.Draw()
    c.SaveAs(ofilename+"_"+varfilename+"_compare.png")
    c.SetLogx(0)
    c.SetLogy(0)
 

def make_scanplots_eff(histo_denom, histo1, histo2, varname, varfilename, cutstring_all, logx=False, logy=False): # canvas and tree not passed!
    
    print cutstring_all
    tree_ref.Draw(varname+">>"+histo_denom.GetName(), cutstring_all)    
    print histo_denom.GetEntries(), histo1.GetEntries(), histo2.GetEntries()
    eff1=TEfficiency(histo1, histo_denom)
    eff2=TEfficiency(histo2, histo_denom)

    eff1.SetLineWidth(2)
    eff2.SetLineWidth(2)
    eff1.SetLineColor(kRed)
    eff2.SetLineColor(kGreen)
    eff1.SetTitle(histo1.GetTitle()+";"+histo1.GetXaxis().GetTitle()+";Eff.")

    legend=TLegend(0.7,0.75,0.89,0.88)
    legend.AddEntry(histo1,"ref_no_compare")
    legend.AddEntry(histo2,"compare_no_ref")

    if logx: c.SetLogx()
    if logy: c.SetLogy()
    eff1.Draw("")
    eff2.Draw("same")
    legend.Draw()
    c.SaveAs(ofilename+"_"+varfilename+"_effcompare.png")
    c.SetLogx(0)
    c.SetLogy(0)
 
    
 
# pt eta phi for the scanTracks
ptbins=[0.01, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 40, 50, 100, 200, 500, 1000]
ptbins2=[0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 1, 1.25, 1.5, 1.75, 2. ]
from array import array

pt_1=TH1F("pt_1","compare tracks by sim p_{T}; sim p_{T}; N. tracks",len(ptbins)-1,array('d', ptbins))
pt_2=TH1F("pt_2","compare tracks by sim p_{T}; sim p_{T}; N. tracks",len(ptbins)-1,array('d', ptbins))  
pt_D=TH1F("pt_D","compare tracks by sim p_{T}; sim p_{T}; N. tracks",len(ptbins)-1,array('d', ptbins))  
pt_11=TH1F("pt_11","compare tracks by sim p_{T}; sim p_{T}; N. tracks",len(ptbins2)-1,array('d', ptbins2))
pt_21=TH1F("pt_21","compare tracks by sim p_{T}; sim p_{T}; N. tracks",len(ptbins2)-1,array('d', ptbins2))  
pt_D1=TH1F("pt_D1","compare tracks by sim p_{T}; sim p_{T}; N. tracks",len(ptbins2)-1,array('d', ptbins2))  
eta_1=TH1F("eta_1","compare tracks by sim #eta; sim #eta; N. tracks",60,-3,3)
eta_2=TH1F("eta_2","compare tracks by sim #eta; sim #eta; N. tracks",60,-3,3)  
eta_D=TH1F("eta_D","compare tracks by sim #eta; sim #eta; N. tracks",60,-3,3)  
phi_1=TH1F("phi_1","compare tracks by sim #phi; sim #phi; N. tracks",70,-3.5,3.5)
phi_2=TH1F("phi_2","compare tracks by sim #phi; sim #phi; N. tracks",70,-3.5,3.5)  
phi_D=TH1F("phi_D","compare tracks by sim #phi; sim #phi; N. tracks",70,-3.5,3.5)  
nL_1=TH1F("nL_1","compare tracks by sim nLayers; sim nLayers; N. tracks",26,-0.5,25.5)
nL_2=TH1F("nL_2","compare tracks by sim nLayers; sim nLayers; N. tracks",26,-0.5,25.5)
nL_D=TH1F("nL_D","compare tracks by sim nLayers; sim nLayers; N. tracks",26,-0.5,25.5)
nH_1=TH1F("nH_1","compare tracks by sim nHits; sim nHits; N. tracks",30,-0.5,29.5)
nH_2=TH1F("nH_2","compare tracks by sim nHits; sim nHits; N. tracks",30,-0.5,29.5)      
nH_D=TH1F("nH_D","compare tracks by sim nHits; sim nHits; N. tracks",30,-0.5,29.5)    
 
make_scanplots(pt_1, pt_2, "pt_mc_gen", "pt", cutstring1, cutstring2, logx=True, logy=True)
make_scanplots(pt_11, pt_21, "pt_mc_gen", "ptlow", cutstring1, cutstring2, logx=False, logy=False)
make_scanplots(eta_1, eta_2, "eta_mc_gen", "eta", cutstring1, cutstring2, logx=False, logy=False)
make_scanplots(phi_1, phi_2, "phi_mc_gen", "phi", cutstring1, cutstring2, logx=False, logy=False)
make_scanplots(nL_1, nL_2, "nLayers_mc", "nLayers", cutstring1, cutstring2, logx=False, logy=False)    
make_scanplots(nH_1, nH_2, "nHits_mc", "nHits", cutstring1, cutstring2, logx=False, logy=False)   

make_scanplots_eff(pt_D, pt_1, pt_2, "pt_mc_gen", "pt", cutstring_all, logx=True, logy=True)
make_scanplots_eff(pt_D1, pt_11, pt_21, "pt_mc_gen", "ptlow", cutstring_all, logx=False, logy=False)
make_scanplots_eff(eta_D, eta_1, eta_2, "eta_mc_gen", "eta", cutstring_all, logx=False, logy=False)
make_scanplots_eff(phi_D, phi_1, phi_2, "phi_mc_gen", "phi", cutstring_all, logx=False, logy=False)
make_scanplots_eff(nL_D, nL_1, nL_2, "nLayers_mc", "nLayers", cutstring_all, logx=False, logy=False)    
make_scanplots_eff(nH_D, nH_1, nH_2, "nHits_mc", "nHits", cutstring_all, logx=False, logy=False)   
