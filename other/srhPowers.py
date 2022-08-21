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
fitNames = findFits('/home/svlesovoi/SRH/2018/10/26/','mf_' + '*.fit')

fitNames.sort();
fitNames = fitNames[1:];
srhFits = SrhFitsFile(fitNames[0], 256)

for fName in fitNames[1:]:
    srhFits.append(fName)

freqIndex = 1

time = srhFits.freqTime[freqIndex,:]
dataLength = time.shape[0]
times = NP.zeros(dataLength)
powers = NP.zeros((srhFits.freqListLength, dataLength))
for t in range(dataLength):
    for f in range(srhFits.freqListLength):
        for ant in range(srhFits.antennaNumbers.shape[0]):
            powers[f, t] += srhFits.ampLcp[f, t, ant]
        

PL.figure()
pl0 = PL.subplot(111)
PL.title('%s ,  %d MHz'%(srhFits.dateObs.split('T')[0], srhFits.freqList[freqIndex]))
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
pl0.imshow(powers, aspect=40, origin='lower')
