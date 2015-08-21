import os.path, glob, sys
import ROOT

suffix = "ttm2"

g = ROOT.TFile("test"+suffix+"_bh.root","recreate")

h_SM3 = ROOT.TH1F("h_SM3", "h_SM3", 300, 0, 300)
h_SM6 = ROOT.TH1F("h_SM6", "h_SM6", 300, 0, 300)
h_SM9 = ROOT.TH1F("h_SM9", "h_SM9", 300, 0, 300)

h_MX3 = ROOT.TH1F("h_MX3", "h_MX3", 300, 0, 300)
h_MX6 = ROOT.TH1F("h_MX6", "h_MX6", 300, 0, 300)
h_MX9 = ROOT.TH1F("h_MX9", "h_MX9", 300, 0, 300)

h_MXNH = ROOT.TH1F("h_MXNH", "h_MXNH", 11, 0, 11)
h_MXC2 = ROOT.TH1F("h_MXC2", "h_MXC2", 40, 0, 20)
h_MXPT = ROOT.TH1F("h_MXPT", "h_MXPT", 20, 0, 20)
h_MXPTm = ROOT.TH1F("h_MXPTm", "h_MXPTm", 20, 0, 20)
h_MXPTr = ROOT.TH1F("h_MXPTr", "h_MXPTr", 40, -2, 2)

# simtrack only validation
h_SIMNH = ROOT.TH1F("h_SIMNH", "h_SIMNH", 11, 0, 11)
h_SIMC2 = ROOT.TH1F("h_SIMC2", "h_SIMC2", 40, 0, 20)
h_SIMPT = ROOT.TH1F("h_SIMPT", "h_SIMPT", 20, 0, 20)
h_SIMPHI = ROOT.TH1F("h_SIMPHI", "h_SIMPHI", 32, -3.2, 3.2)
h_SIMETA = ROOT.TH1F("h_SIMETA", "h_SIMETA", 20, -2, 2)

with open('log'+suffix+'.txt') as f:
    for line in f:
        if "number of hits in window in layer 3" in line:
            val = line.split()[10]
            if "SM" in line:
                h_SM3.Fill(float(val))
            if "MX" in line:
                h_MX3.Fill(float(val))
        elif "number of hits in window in layer 6" in line:
            val = line.split()[10]
            if "SM" in line:
                h_SM6.Fill(float(val))
            if "MX" in line:
                h_MX6.Fill(float(val))
        elif "number of hits in window in layer 9" in line:
            val = line.split()[10]
            if "SM" in line:
                h_SM9.Fill(float(val))
            if "MX" in line:
                h_MX9.Fill(float(val))
        elif "MXBH - found track with nHitIdx" in line:
            NH = line.split()[5].split('=')[1]
            h_MXNH.Fill(float(NH))
            C2 = line.split()[6].split('=')[1]
            h_MXC2.Fill(float(C2)/(3.*(float(NH)-3.)-5.))
            if (float(NH)>9):
                PT = line.split()[7].split('=')[1]
                h_MXPT.Fill(float(PT))
                PTMC = line.split()[8].split('=')[1]
                h_MXPTm.Fill(float(PTMC))
                h_MXPTr.Fill((float(PT)-float(PTMC))/float(PTMC))
        elif "MX - simtrack with nHits" in line:
            NH = line.split()[4].split('=')[1]
            h_SIMNH.Fill(float(NH))
            C2 = line.split()[5].split('=')[1]
            h_SIMC2.Fill(float(C2))
            PT = line.split()[6].split('=')[1]
            h_SIMPT.Fill(float(PT))
            PHI = line.split()[7].split('=')[1]
            h_SIMPHI.Fill(float(PHI))
            ETA = line.split()[8].split('=')[1]
            h_SIMETA.Fill(float(ETA))
            

h_SM3.Draw()
h_MX3.Draw()
h_SM6.Draw()
h_MX6.Draw()
h_SM9.Draw()
h_MX9.Draw()
h_MXNH.Draw()
h_MXC2.Draw()
h_MXPT.Draw()
h_MXPTm.Draw()
h_MXPTr.Draw()
h_SIMNH.Draw()
h_SIMC2.Draw()
h_SIMPT.Draw()
h_SIMPHI.Draw()
h_SIMETA.Draw()
g.Write()
g.Close()
