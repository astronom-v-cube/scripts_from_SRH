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

application = QtWidgets.QApplication.instance();
if not application:
    application = QtWidgets.QApplication(sys.argv);

dateName = DT.datetime.now().strftime("%Y%m%d");
fitNames = findFits('/home/svlesovoi/SRH/2019/03/PowerTest/26/four/','mf_' + '*.fit')
#fitNames = findFits('/home/svlesovoi/SRH/2019/03/PowerTest/27/','mf_' + '*.fit')
#fitNames = findFits('/home/svlesovoi/SRH/2019/03/PowerTest/27/sun/','mf_' + '*.fit')
#fitNames = findFits('/home/svlesovoi/SRH/2019/03/PowerTest/28/','mf_' + '*.fit')
#fitNames = findFits('/home/svlesovoi/SRH/2019/03/PowerTest/29/4','mf_' + '*.fit')

fitNames.sort();
fitNames = fitNames[1:];
srhFits = SrhFitsFile(fitNames[0], 256)

for fName in fitNames[1:]:
    srhFits.append(fName)

ant77Index = 44
ant78Index = 45
ant79Index = 46
ant80Index = 47

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
        
I77 = (LCP77 + RCP77)*0.5
I78 = (LCP78 + RCP78)*0.5
I79 = (LCP79 + RCP79)*0.5
I80 = (LCP80 + RCP80)*0.5

t0 = 0*srhFits.freqListLength
t1 = t0 + 2500
PL.figure()
pl0 = PL.subplot(111)
PL.title('%s ,  %d MHz'%(srhFits.dateObs.split('T')[0], srhFits.freqList[freqIndex]))
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
pl0.plot(times[t0:t1], I77[t0:t1])
pl0.plot(times[t0:t1], I78[t0:t1])
pl0.plot(times[t0:t1], I79[t0:t1])
pl0.plot(times[t0:t1], I80[t0:t1])

pl0.set_yscale('log')
pl0.grid()
