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
fitNames = findFits('C:/SRH/2019/03/01/beginDay/','mf_' + '*.fit')

fitNames.sort();
fitNames = fitNames[1:];
srhFits = SrhFitsFile(fitNames[0], 256)

for fName in fitNames:
    srhFits.append(fName)

freqIndex = 25
ant72Index = 39
ant73Index = 40
ant74Index = 41
time = srhFits.freqTime[freqIndex,:]
I72 = srhFits.ampLcp[freqIndex,:,ant72Index] + srhFits.ampRcp[freqIndex,:,ant72Index]
I73 = srhFits.ampLcp[freqIndex,:,ant73Index] + srhFits.ampRcp[freqIndex,:,ant73Index]
I74 = srhFits.ampLcp[freqIndex,:,ant74Index] + srhFits.ampRcp[freqIndex,:,ant74Index]

t0 = 75
t1 = 81
t2 = 95
t3 = 115
t4 = 555
t5 = 570
t6 = 590
t7 = 605
t8 = 666
t9 = 684

stdSky72 = NP.std(I72[t0:t1])
stdSky73 = NP.std(I73[t0:t1])
stdSun72 = NP.std(I72[t2:t3])
stdSun73 = NP.std(I73[t2:t3])

stdSun72_1 = NP.std(I72[t4:t5])
stdSky72_1 = NP.std(I72[t6:t7])
stdSun73_1 = NP.std(I73[t4:t5])
stdSky73_1 = NP.std(I73[t6:t7])

stdLoad = NP.std(I72[t8:t9])

sun72   = NP.mean(I72[t2:t3]) - NP.mean(I72[t0:t1])
sun73   = NP.mean(I73[t2:t3]) - NP.mean(I73[t0:t1])
sun72_1 = NP.mean(I72[t4:t5]) - NP.mean(I72[t6:t7])
sun73_1 = NP.mean(I73[t4:t5]) - NP.mean(I73[t6:t7])

PL.figure()
pl0 = PL.subplot(111)
PL.title('2019 03 01, intensity at %d MHz'%(srhFits.freqList[freqIndex]))
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
pl0.plot(time[t0:t1], I72[t0:t1], 'o')
pl0.plot(time[t0:t1], I73[t0:t1], 'o')
pl0.plot(time[t2:t3], I72[t2:t3], 'o')
pl0.plot(time[t2:t3], I73[t2:t3], 'o')
pl0.plot(time[t4:t5], I72[t4:t5], 'o')
pl0.plot(time[t6:t7], I72[t6:t7], 'o')
pl0.plot(time[t4:t5], I73[t4:t5], 'o')
pl0.plot(time[t6:t7], I73[t6:t7], 'o')

pl0.plot(time, I72, label='ant72, RMS_sky %2.2f, RMS_sun %2.2f, sun %2.2f' % (stdSky72, stdSun72, sun72))
pl0.plot(time, I73, label='ant73, RMS_sky %2.2f, RMS_sun %2.2f, sun %2.2f' % (stdSky73, stdSun73, sun73))
pl0.plot(time[t4:], I72[t4:],label='ant72_1, RMS_sky %2.2f, RMS_sun %2.2f, sun %2.2f' % (stdSky72_1, stdSun72_1, sun72_1))
pl0.plot(time[t4:], I73[t4:],label='ant73_1, RMS_sky %2.2f, RMS_sun %2.2f, sun %2.2f' % (stdSky73_1, stdSun73_1, sun73_1))
pl0.plot(time[t8:t9], I72[t8:t9],'o',label='RMS_load %2.2f' % (stdLoad))

pl0.legend()
