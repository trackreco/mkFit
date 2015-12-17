import os.path, glob, sys
import ROOT
import array

g = ROOT.TFile("benchmark.root","recreate")

for test in {'BH','CE','CEST'}:
    pos = 15
    if 'BH' in test: pos = 9
    g_HOST_VU = ROOT.TGraph(4)
    g_HOST_VU_speedup = ROOT.TGraph(4)
    #os.system('grep Total log_10x20k_'+test+'_NVU*_NTH1.txt >& log_10x20k_'+test+'_VU.txt')
    with open('log_10x20k_'+test+'_VU.txt') as f:
        point = 0
        for line in f:
            lsplit = line.split()
            yval = float(lsplit[pos])
            for sub in lsplit: 
            #print val
                if "NVU" in sub:
                    for subsub in sub.split('_'):
                        if "NVU" in subsub:
                            xval = float(subsub.replace("NVU",""))
                            print xval, yval
                            g_HOST_VU.SetPoint(point,xval,yval)
                            point = point+1
    x0 = array.array('d',[0])
    y0 = array.array('d',[0])
    g_HOST_VU.GetPoint(0,x0,y0)
    point = 0
    for vu in {1,2,4,8}:
        xval = array.array('d',[0])
        yval = array.array('d',[0])
        g_HOST_VU.GetPoint(point,xval,yval)
        g_HOST_VU_speedup.SetPoint(point,xval[0],y0[0]/yval[0])
        point = point+1
    g_HOST_VU.Write("g_HOST_"+test+"_VU")
    g_HOST_VU_speedup.Write("g_HOST_"+test+"_VU_speedup")

    g_HOSTTH = ROOT.TGraph(4)
    g_HOSTTH_speedup = ROOT.TGraph(4)
    #os.system('grep Total log_10x20k_'+test+'_NVU8_NTH*.txt >& log_10x20k_'+test+'_TH.txt')
    with open('log_10x20k_'+test+'_TH.txt') as f:
        point = 0
        for line in f:
            lsplit = line.split()
            yval = float(lsplit[pos])
            for sub in lsplit: 
            #print val
                if "NTH" in sub:
                    for subsub in sub.split('_'):
                        if "NTH" in subsub:
                            print subsub
                            xval = float(subsub.replace("NTH","").replace(".txt:Total",""))
                            print xval, yval
                            g_HOSTTH.SetPoint(point,xval,yval)
                            point = point+1
    x0 = array.array('d',[0])
    y0 = array.array('d',[0])
    g_HOSTTH.GetPoint(0,x0,y0)
    point = 0
    for vu in {1,3,7,21}:
        xval = array.array('d',[0])
        yval = array.array('d',[0])
        g_HOSTTH.GetPoint(point,xval,yval)
        g_HOSTTH_speedup.SetPoint(point,xval[0],y0[0]/yval[0])
        point = point+1
    g_HOSTTH.Write("g_HOST_"+test+"_TH")
    g_HOSTTH_speedup.Write("g_HOST_"+test+"_TH_speedup")



g.Write()
g.Close()
