import os.path, glob, sys
import ROOT
import array

if len(sys.argv)<2 or len(sys.argv)>3: exit

hORm = sys.argv[1]

isCMSSW = False
if len(sys.argv)>2:
    if sys.argv[2]=="cmssw": isCMSSW=True
    else: exit

if hORm!='snb' and hORm!='snb_endcap' and hORm!='knc' and hORm!='knc_endcap': exit

g = ROOT.TFile('benchmark_'+hORm+'.root',"recreate")

for test in ['BH','COMB','FIT']:
    if isCMSSW and test=='FIT': continue
    if 'endcap' in hORm and not isCMSSW and 'FIT' not in test: continue
    print test

    if 'FIT'  in test: pos = 3  
    if 'BH'   in test: pos = 8  
    if 'COMB' in test: pos = 11 

    if isCMSSW:
        ntks = '100xTTbarPU35'
        nevt = 100.
    else :
        if test=='FIT' :
            ntks = '1kx10k'
            nevt = 1000.
        else :
            ntks = '20x10k' 
            nevt = 20.

    # Vectorization data points
    vuvals = ['1','2','4','8']
    if 'knc' in hORm: 
        nth = '1'
        vuvals.append('16')
        vuvals.append('16int')
    else:
        nth = '1'
        vuvals.append('8int')

    #Vectorization time    
    print "Vectorization"
    g_VU = ROOT.TGraph(4)
    g_VU_speedup = ROOT.TGraph(4)
    point = 0
    for vu in vuvals:
        if    vu == '16int': xval = 16.0
        elif  vu == '8int' : xval = 8.0
        else               : xval = float(vu)
        yval = 0.
        firstFound = False
        
        os.system('grep Matriplex log_'+hORm+'_'+ntks+'_'+test+'_NVU'+vu+'_NTH'+nth+'.txt >& log_'+hORm+'_'+ntks+'_'+test+'_VU.txt')
        with open('log_'+hORm+'_'+ntks+'_'+test+'_VU.txt') as f:
            for line in f:
                if 'Matriplex' not in line: continue
                if 'Total' in line: continue
                lsplit = line.split()
                if not firstFound:
                    firstFound = True
                    continue
                yval = yval+float(lsplit[pos])
        
        yval = nevt*yval/(nevt-1.)
        print xval, yval
        g_VU.SetPoint(point,xval,yval)
        point = point+1
    g_VU.Write("g_"+test+"_VU")

    # Vectorization speedup
    x0 = array.array('d',[0])
    y0 = array.array('d',[0])
    g_VU.GetPoint(0,x0,y0)
    point = 0
    for vu in vuvals:
        xval = array.array('d',[0])
        yval = array.array('d',[0])
        g_VU.GetPoint(point,xval,yval)
        speedup = 0.
        if yval[0]>0.: speedup = y0[0]/yval[0]
        g_VU_speedup.SetPoint(point,xval[0],speedup)
        point = point+1
    g_VU_speedup.Write("g_"+test+"_VU_speedup")

    # Parallelization datapoints
    if 'knc' in hORm: 
        nvu = '16int'
        thvals = [1,2,4,8,15,30,60,90,120,150,180,210,240]
    else :
        nvu = '8int'
        thvals = [1,2,4,6,8,12,16,20,24]
    
    # Parallelization time
    print "Parallelization"
    g_TH = ROOT.TGraph(len(thvals))
    g_TH_speedup = ROOT.TGraph(len(thvals))
    point = 0
    for th in thvals:
        xval = float(th)
        yval = 0.
        firstFound = False

        os.system('grep Matriplex log_'+hORm+'_'+ntks+'_'+test+'_NVU'+nvu+'_NTH'+str(th)+'.txt >& log_'+hORm+'_'+ntks+'_'+test+'_TH.txt')
        with open('log_'+hORm+'_'+ntks+'_'+test+'_TH.txt') as f:
            for line in f:
                if 'Matriplex' not in line: continue
                if 'Total' in line: continue
                lsplit = line.split()
                if not firstFound:
                    firstFound = True
                    continue
                yval = yval+float(lsplit[pos])

        yval = nevt*yval/(nevt-1.)
        print xval, yval
        g_TH.SetPoint(point,xval,yval)
        point = point+1
    g_TH.Write("g_"+test+"_TH")

    # Parallelization speedup
    x0 = array.array('d',[0])
    y0 = array.array('d',[0])
    g_TH.GetPoint(0,x0,y0)
    point = 0
    for th in thvals:
        xval = array.array('d',[0])
        yval = array.array('d',[0])
        g_TH.GetPoint(point,xval,yval)
        speedup = 0.
        if yval[0]>0.: speedup = y0[0]/yval[0]
        g_TH_speedup.SetPoint(point,xval[0],speedup)
        point = point+1
    g_TH_speedup.Write("g_"+test+"_TH_speedup")

g.Write()
g.Close()
