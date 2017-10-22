import os.path, glob, sys
import ROOT

if len(sys.argv)!=6: exit

hORm   = sys.argv[1]
sample = sys.argv[2]
region = sys.argv[3]
test   = sys.argv[4]
suffix = sys.argv[5]

g = ROOT.TFile("test_"+hORm+"_"+sample+"_"+region+"_"+test+"_"+suffix+".root","recreate")

h_SM3 = ROOT.TH1F("h_SM3", "h_SM3", 300, 0, 300)
h_SM6 = ROOT.TH1F("h_SM6", "h_SM6", 300, 0, 300)
h_SM9 = ROOT.TH1F("h_SM9", "h_SM9", 300, 0, 300)

h_MX3 = ROOT.TH1F("h_MX3", "h_MX3", 300, 0, 300)
h_MX6 = ROOT.TH1F("h_MX6", "h_MX6", 300, 0, 300)
h_MX9 = ROOT.TH1F("h_MX9", "h_MX9", 300, 0, 300)

h_MXNH  = ROOT.TH1F("h_MXNH", "h_MXNH", 26, 0, 26)
h_MXNH.GetXaxis().SetTitle("nHits_{found}")
h_MXC2  = ROOT.TH1F("h_MXC2", "h_MXC2", 40, 0, 20)
h_MXPT  = ROOT.TH1F("h_MXPT", "h_MXPT", 20, 0, 20)
h_MXPT.GetXaxis().SetTitle("p_{T}^{rec}")
h_MXPTm = ROOT.TH1F("h_MXPTm", "h_MXPTm", 20, 0, 20)
h_MXPTm.GetXaxis().SetTitle("p_{T}^{sim}")
h_MXPTr = ROOT.TH1F("h_MXPTr", "h_MXPTr", 40, -0.4, 0.4)
h_MXPTr.GetXaxis().SetTitle("(p_{T}^{rec}-p_{T}^{sim})/p_{T}^{sim}")

h_NHDS  = ROOT.TH1F("h_NHDS", "h_NHDS", 15, -7.5, 7.5)
h_NHDS.GetXaxis().SetTitle("nHits_{found}-nHits_{sim}")
h_NHDR  = ROOT.TH1F("h_NHDR", "h_NHDR", 15, -7.5, 7.5)
h_NHDR.GetXaxis().SetTitle("nHits_{found}-nHits_{cmssw}")

# simtrack only validation
h_SIMNH  = ROOT.TH1F("h_SIMNH", "h_SIMNH", 26, 0, 26)
h_SIMC2  = ROOT.TH1F("h_SIMC2", "h_SIMC2", 40, 0, 20)
h_SIMPT  = ROOT.TH1F("h_SIMPT", "h_SIMPT", 20, 0, 20)
h_SIMPHI = ROOT.TH1F("h_SIMPHI", "h_SIMPHI", 32, -3.2, 3.2)
h_SIMETA = ROOT.TH1F("h_SIMETA", "h_SIMETA", 20, -2, 2)

# recseed only validation
h_SEEDNH     = ROOT.TH1F("h_SEEDNH", "h_SEEDNH", 26, 0, 26)
h_SEEDC2     = ROOT.TH1F("h_SEEDC2", "h_SEEDC2", 40, 0, 20)
h_SEEDPOSETA = ROOT.TH1F("h_SEEDPOSETA", "h_SEEDPOSETA", 20, -2, 2)
h_SEEDPOSPHI = ROOT.TH1F("h_SEEDPOSPHI", "h_SEEDPOSPHI", 32, -3.2, 3.2)
h_SEEDPOSR   = ROOT.TH1F("h_SEEDPOSR", "h_SEEDPOSR", 20, 10, 14)
h_SEEDPT     = ROOT.TH1F("h_SEEDPT", "h_SEEDPT", 20, 0, 20)

with open('log_'+hORm+'_'+sample+'_'+region+'_'+test+'_'+suffix+'.txt') as f:
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
        elif "MX - found track with nFoundHits" in line:
            NH = line.split()[5].split('=')[1]
            h_MXNH.Fill(float(NH))
            C2 = line.split()[6].split('=')[1]
            h_MXC2.Fill(float(C2)/(3.*(float(NH)-3.)-5.))
            PT = line.split()[7].split('=')[1]
            h_MXPT.Fill(float(PT))
            PTMC = line.split()[8].split('=')[1]
            NHS = line.split()[9].split('=')[1]
            NHR = line.split()[10].split('=')[1]
            if float(PTMC)>0.01 and float(NHS)>=3 and float(NHR)>=3:
                h_MXPTm.Fill(float(PTMC))
                h_MXPTr.Fill((float(PT)-float(PTMC))/float(PTMC))
                h_NHDS.Fill(float(NH)-float(NHS))
                h_NHDR.Fill(float(NH)-float(NHR))
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
        elif "MX - found seed with nHits" in line:
            NH = line.split()[5].split('=')[1]
            h_SEEDNH.Fill(float(NH))
            C2 = line.split()[6].split('=')[1]
            h_SEEDC2.Fill(float(C2))
            POSETA = line.split()[7].split('=')[1]
            h_SEEDPOSETA.Fill(float(POSETA))
            POSPHI = line.split()[8].split('=')[1]
            h_SEEDPOSPHI.Fill(float(POSPHI))
            POSR = line.split()[9].split('=')[1]
            h_SEEDPOSR.Fill(float(POSR))
            PT = line.split()[10].split('=')[1]
            h_SEEDPT.Fill(float(PT))

#h_SM3.Draw()
#h_MX3.Draw()
#h_SM6.Draw()
#h_MX6.Draw()
#h_SM9.Draw()
#h_MX9.Draw()
#h_MXNH.Draw()
#h_NHDS.Draw()
#h_NHDR.Draw()
#h_MXC2.Draw()
#h_MXPT.Draw()
#h_MXPTm.Draw()
#h_MXPTr.Draw()
#h_SIMNH.Draw()
#h_SIMC2.Draw()
#h_SIMPT.Draw()
#h_SIMPHI.Draw()
#h_SIMETA.Draw()
#h_SEEDNH.Draw()
#h_SEEDC2.Draw()
#h_SEEDPOSETA.Draw()
#h_SEEDPOSPHI.Draw()
#h_SEEDPOSR.Draw()
#h_SEEDPT.Draw()
g.Write()
g.Close()
