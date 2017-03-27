import os.path, glob, sys
import ROOT
import array
import math

if len(sys.argv)!=4: exit

hORm   = sys.argv[1] # host == Xeon SNB, KNC, KNL
sample = sys.argv[2] # toymc or cmssw
region = sys.argv[3] # barrel or endcap

g = ROOT.TFile('benchmark_'+hORm+'_'+sample+'_'+region+'.root',"recreate")

for test in ['BH','STD','CE','FIT']:
    if sample == 'CMSSW' and test == 'FIT' : continue
    if region == 'Endcap' and not sample == 'CMSSW' and not test == 'FIT' : continue
    print hORm,sample,region,test

    if test == 'FIT' : pos = 3  
    if test == 'BH'  : pos = 8  
    if test == 'STD' : pos = 11  
    if test == 'CE'  : pos = 14 

    if sample == 'CMSSW' :
        nevt = 100.
    else :
        if test == 'FIT' :
            nevt = 1000.
        else :
            nevt = 20.

    # Vectorization data points
    vuvals = ['1','2','4','8']
    if hORm == 'KNC' or hORm == 'KNL':
        nth = '1'
        vuvals.append('16')
        vuvals.append('16int')
    else:
        nth = '1'
        vuvals.append('8int')

    #Vectorization time    
    print "Vectorization"
    g_VU = ROOT.TGraphErrors(len(vuvals))
    g_VU_speedup = ROOT.TGraphErrors(len(vuvals))
    point = 0
    for vu in vuvals:
        if    vu == '16int': xval = 16.0
        elif  vu == '8int' : xval = 8.0
        else               : xval = float(vu)
        yvals = array.array('d');
        firstFound = False

        # Read in times from log files, store into yvals
        os.system('grep Matriplex log_'+hORm+'_'+sample+'_'+region+'_'+test+'_NVU'+vu+'_NTH'+nth+'.txt >& log_'+hORm+'_'+sample+'_'+region+'_'+test+'_VU.txt')
        with open('log_'+hORm+'_'+sample+'_'+region+'_'+test+'_VU.txt') as f:
            for line in f:
                if 'Matriplex' not in line: continue
                if 'Total' in line: continue
                lsplit = line.split()
                if not firstFound:
                    firstFound = True
                    continue
                yvals.append(float(lsplit[pos]))

        # Compute mean and uncertainty on mean
        sum = 0.;
        for yval in range(0,len(yvals)):
            sum = sum + yvals[yval]
        mean = sum/len(yvals)
        emean = 0.;
        for yval in range(0,len(yvals)):
            emean = emean + ((yvals[yval] - mean) * (yvals[yval] - mean))
        emean = math.sqrt(emean / (len(yvals) - 1))
        emean = emean/math.sqrt(len(yvals))

        print xval,mean,'+/-',emean
        g_VU.SetPoint(point,xval,mean)
        g_VU.SetPointError(point,0,emean)
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
    if hORm == 'KNC' or hORm == 'KNL':
        nvu = '16int'
        thvals = [1,2,4,8,15,30,60,90,120,150,180,210,240]
    else :
        nvu = '8int'
        thvals = [1,2,4,6,8,12,16,20,24]
    
    # Parallelization time
    print "Parallelization"
    g_TH = ROOT.TGraphErrors(len(thvals))
    g_TH_speedup = ROOT.TGraphErrors(len(thvals))
    point = 0
    for th in thvals:
        xval = float(th)
        yvals = array.array('d');
        firstFound = False

        os.system('grep Matriplex log_'+hORm+'_'+sample+'_'+region+'_'+test+'_NVU'+nvu+'_NTH'+str(th)+'.txt >& log_'+hORm+'_'+sample+'_'+region+'_'+test+'_TH.txt')
        with open('log_'+hORm+'_'+sample+'_'+region+'_'+test+'_TH.txt') as f:
            for line in f:
                if 'Matriplex' not in line: continue
                if 'Total' in line: continue
                lsplit = line.split()
                if not firstFound:
                    firstFound = True
                    continue
                yvals.append(float(lsplit[pos]))
        
        # Compute mean and uncertainty on mean
        sum = 0.;
        for yval in range(0,len(yvals)):
            sum = sum + yvals[yval]
        mean = sum/len(yvals)
        emean = 0.;
        for yval in range(0,len(yvals)):
            emean = emean + ((yvals[yval] - mean) * (yvals[yval] - mean))
        emean = math.sqrt(emean / (len(yvals) - 1))
        emean = emean/math.sqrt(len(yvals))

        print xval,mean,'+/-',emean
        g_TH.SetPoint(point,xval,mean)
        g_TH.SetPointError(point,0,emean)
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
