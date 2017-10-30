import os.path, glob, sys
import ROOT

arch   = sys.argv[1]
sample = sys.argv[2]
build  = sys.argv[3]
suffix = sys.argv[4]

g = ROOT.TFile("test_"+arch+"_"+sample+"_"+build+"_"+suffix+".root","recreate")

# declare hists
h_MXNH  = ROOT.TH1F("h_MXNH", "h_MXNH", 35, 0, 35)
h_MXNH.GetXaxis().SetTitle("nHits_{found}")

h_MXC2  = ROOT.TH1F("h_MXC2", "h_MXC2", 40, 0, 20)
h_MXC2.GetXaxis().SetTitle("#chi^{2}")

h_MXPT  = ROOT.TH1F("h_MXPT", "h_MXPT", 20, 0, 20)
h_MXPT.GetXaxis().SetTitle("p_{T}^{rec}")

h_MXPTm = ROOT.TH1F("h_MXPTm", "h_MXPTm", 20, 0, 20)
h_MXPTm.GetXaxis().SetTitle("p_{T}^{sim}")

h_MXPTr = ROOT.TH1F("h_MXPTr", "h_MXPTr", 40, -0.4, 0.4)
h_MXPTr.GetXaxis().SetTitle("(p_{T}^{rec}-p_{T}^{sim})/p_{T}^{sim}")

h_NHDS  = ROOT.TH1F("h_NHDS", "h_NHDS", 15, -7.5, 7.5)
h_NHDS.GetXaxis().SetTitle("nHits_{found}-nHits_{sim}")

with open('log_'+arch+'_'+sample+'_'+build+'_'+suffix+'_DumpForPlots.txt') as f :
    for line in f :
        if "MX - found track with nFoundHits" in line :
            lsplit = line.split()

            NH = float(lsplit[6])
            h_MXNH.Fill(NH)

            C2 = float(lsplit[8])
            h_MXC2.Fill(C2/(3.*(NH-3.)-5.))

            PT = float(lsplit[10])
            h_MXPT.Fill(PT)

            PTMC = float(lsplit[12])
            NHS = float(lsplit[12])
            if PTMC > 0.01 and NHS >= 3 :
                h_MXPTm.Fill(PTMC)
                h_MXPTr.Fill((PT-PTMC)/PTMC)
                h_NHDS.Fill(NH-NHS)

g.Write()
g.Close()
