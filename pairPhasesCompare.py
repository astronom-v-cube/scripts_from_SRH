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
fitNames = findFits('/home/svlesovoi/SRH/2019/07/26','mf_' + '*.fit')

fitNames.sort()
srhFits = SrhFitsFile(fitNames[0], 256)
for fName in fitNames[1:]:
    srhFits.append(fName)

freqIndex = 15
pairIndex = 512 + 13
time = srhFits.freqTime[freqIndex,:]

PL.figure()
if pairIndex < 512 + 15:
    PL.title('Pair %d, %d '%(192 - (pairIndex - 512), 192 - (pairIndex - 512 + 1)))
else:
    PL.title('Pair %d, %d '%(49 + (pairIndex - 15 - 512), 49 + (pairIndex - 15 - 512 + 1)))
pl0 = PL.subplot(111)
pl0.set_ylim(-4,4)
#pl0.set_ylim(0,.1)
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))

#pl0.plot(time, 10*NP.abs(srhFits.visLcp[freqIndex, :, pairIndex + 1]))
#pl0.plot(time, NP.angle(srhFits.visLcp[freqIndex, :, pairIndex + 1]))

pl0.plot(time, 10*srhFits.visLcp[freqIndex, :, pairIndex].real)
pl0.plot(time, 10*srhFits.visLcp[freqIndex, :, pairIndex].imag)
pl0.plot(time, 10*srhFits.visRcp[freqIndex, :, pairIndex].real)
pl0.plot(time, 10*srhFits.visRcp[freqIndex, :, pairIndex].imag)
pl0.grid()