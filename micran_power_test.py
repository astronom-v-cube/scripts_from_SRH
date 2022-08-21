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
fitNames = findFits('C:/SRH/2019/03/powerTest/mod10','mf_' + '*.fit')

fitNames.sort();
fitNames = fitNames[1:];
srhFits = SrhFitsFile(fitNames[0], 256)

for fName in fitNames:
    srhFits.append(fName)

freqIndex = 25
ant72Index = 40
ant73Index = 41
ant74Index = 42
ant75Index = 43

time = srhFits.freqTime[freqIndex,:]
I72 = srhFits.ampLcp[freqIndex,:,ant72Index] + srhFits.ampRcp[freqIndex,:,ant72Index]
I73 = srhFits.ampLcp[freqIndex,:,ant73Index] + srhFits.ampRcp[freqIndex,:,ant73Index]
I74 = srhFits.ampLcp[freqIndex,:,ant74Index] + srhFits.ampRcp[freqIndex,:,ant74Index]
I75 = srhFits.ampLcp[freqIndex,:,ant75Index] + srhFits.ampRcp[freqIndex,:,ant75Index]

PL.figure()
pl0 = PL.subplot(111)
PL.title('2019 03 17, intensity at %d MHz'%(srhFits.freqList[freqIndex]))
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
pl0.plot(time, I72, '.', label='A72')
pl0.plot(time, I73, '.', label='A73')
pl0.plot(time, I74, '.', label='A74')
pl0.plot(time, I75, '.', label='A75')

t_index0 = 150

val72 = []
val73 = []
val74 = []
val75 = []
for step in range(26):
    t_index = int(t_index0 + step*6.7 + .5)
    t_step = time[t_index]

    val_step = NP.mean(I72[t_index - 2:t_index + 2])
    val72.append(val_step)
    pl0.plot([t_step, t_step],[val_step, val_step], 'o' , color='black')

    val_step = NP.mean(I73[t_index - 2:t_index + 2])
    val73.append(val_step)
    pl0.plot([t_step, t_step],[val_step, val_step], 'o' , color='black')

    val_step = NP.mean(I74[t_index - 2:t_index + 2])
    val74.append(val_step)
    pl0.plot([t_step, t_step],[val_step, val_step], 'o' , color='black')

    val_step = NP.mean(I75[t_index - 2:t_index + 2])
    val75.append(val_step)
    pl0.plot([t_step, t_step],[val_step, val_step], 'o' , color='black')

pl0.legend()

inPower = NP.linspace(-3,22,26)
inPower = 10**(0.1*inPower)/1e3
PL.figure()
PL.title('2019 03 17, digital receiver 10')
PL.xlabel('input power [mW]')
PL.ylabel('output power [arbitrary]')
PL.plot(inPower, val72, label='A72')
PL.plot(inPower, val73, label='A73')
PL.plot(inPower, val74, label='A74')
PL.plot(inPower, val75, label='A75')
PL.legend()
PL.grid()

