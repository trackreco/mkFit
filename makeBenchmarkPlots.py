import os.path, glob, sys
import ROOT
import array

if len(sys.argv)!=2: exit

hORm = sys.argv[1]

if hORm!='host' and hORm!='mic': exit

g = ROOT.TFile('benchmark_'+hORm+'.root',"recreate")

for test in ['BH','CE','CEST','ST','TBBST','FIT']:
    print test
    pos = 14
    ntks = '20k'
    if 'BH' in test: pos = 8
    if 'TBB' in test: pos = 17
    if 'ST' == test: pos = 11
    if 'FIT' in test: 
        pos = 3
        ntks = '1M'
    g_VU = ROOT.TGraph(4)
    g_VU_speedup = ROOT.TGraph(4)
    point = 0
    vuvals = ['1','2','4','8']
    if hORm == 'mic': 
        vuvals.append('16')
        vuvals.append('16int')
    else:
        vuvals.append('8int')
    for vu in vuvals:
        os.system('grep Matriplex log_'+hORm+'_10x'+ntks+'_'+test+'_NVU'+vu+'_NTH1.txt >& log_'+hORm+'_10x'+ntks+'_'+test+'_VU.txt')
        if vu == '16int':
            xval = 16.0
        elif vu == '8int':
            xval = 8.0
        else:
            xval = float(vu)
        yval = 0.
        firstFound = False
        with open('log_'+hORm+'_10x'+ntks+'_'+test+'_VU.txt') as f:
            for line in f:
                if 'Matriplex' not in line: continue
                if 'Total' in line: continue
                lsplit = line.split()
                if 'FIT' in test: 
                    tmp = float(lsplit[pos])
                    if firstFound is False or tmp<yval: yval=tmp
                    if not firstFound:
                        firstFound = True
                else:
                    if not firstFound:
                        firstFound = True
                        continue
                    yval = yval+float(lsplit[pos])
        if 'FIT' not in test: yval = 10.*yval/9.
        print xval, yval
        g_VU.SetPoint(point,xval,yval)
        point = point+1
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
    g_VU.Write("g_"+test+"_VU")
    g_VU_speedup.Write("g_"+test+"_VU_speedup")

    point = 0
    nvu = '8int'
    if hORm == 'mic': nvu = '16int'
    thvals = [1,3,7,21]
    if 'TBB' in test : thvals = [1,3,7,10,12,14,16,21]
    if hORm == 'mic': thvals = [1,3,7,21,42,63,84,105,126,147,168,189,210]
    g_TH = ROOT.TGraph(len(thvals))
    g_TH_speedup = ROOT.TGraph(len(thvals))
    for th in thvals:
        os.system('grep Matriplex log_'+hORm+'_10x'+ntks+'_'+test+'_NVU'+nvu+'_NTH'+str(th)+'.txt >& log_'+hORm+'_10x'+ntks+'_'+test+'_TH.txt')
        xval = float(th)
        yval = 0.
        firstFound = False
        with open('log_'+hORm+'_10x'+ntks+'_'+test+'_TH.txt') as f:
            for line in f:
                if 'Matriplex' not in line: continue
                if 'Total' in line: continue
                lsplit = line.split()
                if 'FIT' in test: 
                    tmp = float(lsplit[pos])
                    if firstFound is False or tmp<yval: yval=tmp
                    if not firstFound:
                        firstFound = True
                else:
                    if not firstFound:
                        firstFound = True
                        continue
                    yval = yval+float(lsplit[pos])
        if 'FIT' not in test: yval = 10.*yval/9.
        print xval, yval
        g_TH.SetPoint(point,xval,yval)
        point = point+1
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
    g_TH.Write("g_"+test+"_TH")
    g_TH_speedup.Write("g_"+test+"_TH_speedup")

g.Write()
g.Close()
