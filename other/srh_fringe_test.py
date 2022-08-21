# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 01:47:01 2016

@author: Sergey
"""
import os, fnmatch;
from srhFitsFile import SrhFitsFile
import numpy as NP;
import datetime as DT;
import pylab as PL;
from PyQt5 import QtGui, QtCore, QtWidgets
from PyQt5.QtWidgets import QFileDialog
import sys;
from optparse import OptionParser;
from astropy.io import fits


dt_major = 3600.;
dt_minor = 900.;
freq_first = 2;
freq_intervals = 5;

def hhmm_format(t, pos):
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  return '%02d:%02d' % (hh,mm);

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;

def smooth(y, box_pts):
    box = NP.ones(box_pts)/box_pts
    y_smooth = NP.convolve(y, box, mode='same')
    return y_smooth

application = QtWidgets.QApplication.instance();
if not application:
    application = QtWidgets.QApplication(sys.argv);

dateName = DT.datetime.now().strftime("%Y%m%d");
fitNames = findFits('/home/svlesovoi/SRH/2019/04/14/','mf_' + '*.fit')

fitNames.sort();
fitNames = fitNames[1:];
srhFits = SrhFitsFile(fitNames[0], 256)

for fName in fitNames[1:]:
    srhFits.append(fName)

ant77Index = 44
ant78Index = 45
ant79Index = 46
ant80Index = 47
cor77_78Index = 512 + 15 + 28
cor78_79Index = 512 + 15 + 29

#ant77Index = 0
#ant78Index = 1
#ant79Index = 2
#ant80Index = 3
#cor77_78Index = 512
#cor78_79Index = 512 + 1
freqIndex = 1

time = srhFits.freqTime[freqIndex,:]
dataLength = time.shape[0] * srhFits.freqListLength

times = NP.zeros(dataLength)
LCP77 = NP.zeros(dataLength)
RCP77 = NP.zeros(dataLength)
LCP78 = NP.zeros(dataLength)
RCP78 = NP.zeros(dataLength)
LCP79 = NP.zeros(dataLength)
RCP79 = NP.zeros(dataLength)
LCP80 = NP.zeros(dataLength)
RCP80 = NP.zeros(dataLength)
LCP77_78 = NP.zeros(dataLength)
RCP77_78 = NP.zeros(dataLength)
LCP78_79 = NP.zeros(dataLength)
RCP78_79 = NP.zeros(dataLength)
for l in range(time.shape[0]):
    for f in range(srhFits.freqListLength):
        times[l*srhFits.freqListLength + f] = srhFits.freqTime[f,l]
        LCP77[l*srhFits.freqListLength + f] = srhFits.ampLcp[f,l,ant77Index]
        RCP77[l*srhFits.freqListLength + f] = srhFits.ampRcp[f,l,ant77Index]
        LCP78[l*srhFits.freqListLength + f] = srhFits.ampLcp[f,l,ant78Index]
        RCP78[l*srhFits.freqListLength + f] = srhFits.ampRcp[f,l,ant78Index]
        LCP79[l*srhFits.freqListLength + f] = srhFits.ampLcp[f,l,ant79Index]
        RCP79[l*srhFits.freqListLength + f] = srhFits.ampRcp[f,l,ant79Index]
        LCP80[l*srhFits.freqListLength + f] = srhFits.ampLcp[f,l,ant80Index]
        RCP80[l*srhFits.freqListLength + f] = srhFits.ampRcp[f,l,ant80Index]
        LCP77_78[l*srhFits.freqListLength + f] = srhFits.visLcp[f,l,cor77_78Index].real
        RCP77_78[l*srhFits.freqListLength + f] = srhFits.visRcp[f,l,cor77_78Index].real
        LCP78_79[l*srhFits.freqListLength + f] = srhFits.visLcp[f,l,cor78_79Index].real
        RCP78_79[l*srhFits.freqListLength + f] = srhFits.visRcp[f,l,cor78_79Index].real
        
I77 = (LCP77 + RCP77)*0.5
I78 = (LCP78 + RCP78)*0.5
I79 = (LCP79 + RCP79)*0.5
I80 = (LCP80 + RCP80)*0.5

t0 = 500*srhFits.freqListLength
t1 = time.shape[0]*srhFits.freqListLength
fig = PL.figure()
pl0 = PL.subplot(211)
pl1 = PL.subplot(212)
fig.suptitle('%s,  %d MHz'%(srhFits.dateObs.split('T')[0], 6000))
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
pl0.plot(times[t0:t1], I77[t0:t1], label='A77(81)')
pl0.plot(times[t0:t1], I78[t0:t1], label='A78(83)')
pl0.plot(times[t0:t1], I79[t0:t1], label='A79')
pl0.plot(times[t0:t1], I80[t0:t1], label='A80')

pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
pl0.set_xlabel('UT')
pl1.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
pl0.set_yscale('log')
pl1.set_xlabel('UT')
pl0.grid()
pl0.legend()

pl1.plot(times[t0:t1], LCP77_78[t0:t1], label='LCP77_78')
pl1.plot(times[t0:t1], RCP77_78[t0:t1], label='RCP77_78')
pl1.plot(times[t0:t1], LCP78_79[t0:t1], label='LCP78_79')
pl1.plot(times[t0:t1], RCP78_79[t0:t1], label='RCP78_79')
pl1.set_ylim(-0.25,0.25)

t_sen0 = t0 + 100
t_sen1 = t0 + 1100
pl1.plot(times[t_sen0:t_sen1], LCP77_78[t_sen0:t_sen1], 'o')

smooth_LCP = smooth(LCP77_78[t_sen0:t_sen1], 9)
dev_LCP = LCP77_78[t_sen0:t_sen1] - smooth_LCP
smooth_RCP = smooth(RCP77_78[t_sen0:t_sen1], 9)
dev_RCP = RCP77_78[t_sen0:t_sen1] - smooth_RCP

NP.max(smooth_LCP[10:-10]) / NP.std(dev_LCP[10:-10])
NP.max(smooth_RCP[10:-10]) / NP.std(dev_RCP[10:-10])

pl1.grid()
pl1.legend()
