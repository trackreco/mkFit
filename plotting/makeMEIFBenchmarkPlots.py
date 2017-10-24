import os.path, glob, sys
import ROOT
import array
import math

arch   = sys.argv[1] # SNB, KNC, KNL
sample = sys.argv[2] 

g = ROOT.TFile('benchmarkMEIF_'+arch+'_'+sample+'.root','recreate')

# Parallelization datapoints
if arch is 'KNC' or arch is 'KNL' :
    nvu = '16int'
    evvals = ['1','2','4','8','16','32','64','128']
    if arch is 'KNC' :
        thvals = ['1','2','4','8','15','30','60','90','120','150','180','210','240']
    else : # KNL
        thvals = ['1','2','4','8','16','32','64','96','128','160','192','224','256']
else : # SNB
    nvu = '8int'
    evvals = ['1','2','4','8','12']
    thvals = ['1','2','4','6','8','12','16','20','24']

# extra text label
text = 'MEIF'

# text for grepping
grepnEV  = '=== TOTAL for'
grepTime = 'Total event loop time'

# needed for speedups
xval0 = array.array('d',[0])
yval0 = array.array('d',[0])

# time    
for evval in evvals :
    print arch,sample,"nEV: ",evval
    
    # define event float
    ev = float(evval)
        
    # define tgraphs vs absolute time and speedup
    g_time    = ROOT.TGraph()
    g_speedup = ROOT.TGraph()

    point = 0
    for thval in thvals :
        xval = float(thval)
        if ev > xval: continue;
            
        # extracted time
        yval = float(0)

        # open log file, grep for relevant lines
        with open('log_'+arch+'_'+sample+'_CE_NVU'+nvu+'_NTH'+thval+'_NEV'+evval+'.txt') as f :
            for line in f :
                if grepnEV in line :
                    lsplit = line.split()                
                    nev  = float(lsplit[3])
                elif grepTime in line :
                    lsplit = line.split()                
                    yval = float(lsplit[4])

        yval /= nev

        # Printout value for good measure
        print xval,yval

        # store val
        g_time.SetPoint(point,xval,yval)
        point = point+1

    # write out the plot
    g_time.Write('g_CE_'+text+'_nEV'+evval+'_time')

    # needed for speedup calculation
    if evval is '1' :
        g_time.GetPoint(0,xval0,yval0)        
        
    # speedup plots
    point = 0
    for thval in thvals :
        xval = float(thval)
        if ev > xval: continue;

        # set up inputs
        xval = array.array('d',[0])
        yval = array.array('d',[0])
        
        # get point from time
        g_time.GetPoint(point,xval,yval)

        speedup  = 0.
        if yval[0] > 0. : 
            speedup  = yval0[0]/yval[0]
                
        # store in speedup plot
        g_speedup.SetPoint(point,xval[0],speedup)
        point = point+1

    # always write out speedup
    g_speedup.Write('g_CE_'+text+'_nEV'+evval+'_speedup')

g.Write()
g.Close()
