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
  hh = (int)(t / 3600.)
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
#fitNames = findFits('/home/svlesovoi/SRH/2019/03/PowerTest/29/2','mf_' + '*.fit')

fitNames.sort();
fitNames = fitNames[1:];
srhFits = SrhFitsFile(fitNames[0], 256)

for fName in fitNames[1:]:
    srhFits.append(fName)

freqIndex = 2
ant77Index = 44
ant78Index = 45
ant79Index = 46
ant80Index = 47

time = srhFits.freqTime[freqIndex,:]
I77 = (srhFits.ampLcp[freqIndex,:,ant77Index] + srhFits.ampRcp[freqIndex,:,ant77Index])*.5
I78 = (srhFits.ampLcp[freqIndex,:,ant78Index] + srhFits.ampRcp[freqIndex,:,ant78Index])*.5
#I79 = srhFits.ampLcp[freqIndex,:,ant74Index] + srhFits.ampRcp[freqIndex,:,ant79Index]
#I80 = srhFits.ampLcp[freqIndex,:,ant75Index] + srhFits.ampRcp[freqIndex,:,ant80Index]
LCP78 = srhFits.ampLcp[freqIndex,:,ant73Index]
RCP78 = srhFits.ampRcp[freqIndex,:,ant73Index]

t0 = 1350
t_sun = t0 + 500
t_sky_0 = t0 + 150
t_sky_1 = t0 + 195
PL.figure()
pl0 = PL.subplot(111)
PL.title('%s , intensity at %d MHz'%(srhFits.dateObs.split('T')[0], srhFits.freqList[freqIndex]))
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
#pl0.plot(time, NP.log10(I77), '.', label='A77')
#pl0.plot(time, NP.log10(I78), '.', label='A78')
#pl0.plot(time, NP.log10(I79), '.', label='A79')
#pl0.plot(time, NP.log10(I80), '.', label='A80')
sky_mean = NP.mean(I78[t_sky_0:t_sky_1])
sun_mean = NP.mean(I78[t_sun:]) - sky_mean
sun_std = NP.std(I78[t_sun:])
pl0.plot(time[t0:], I77[t0:], label='A81')
pl0.plot(time[t0:], I78[t0:], label='A83')
pl0.plot(time[t_sky_0:t_sky_1], I78[t_sky_0:t_sky_1], 'o', label=('sky %0.2f'%sky_mean))
pl0.plot(time[t_sun:], I78[t_sun:], 'o', label=('Sun %0.2f, N/S %0.3f'%(sun_mean, sun_std/sun_mean)))

pl0.plot(time[t0:], LCP78[t0:])
pl0.plot(time[t0:], RCP78[t0:])

#pl0.plot(time, I79, label='A79')
#pl0.plot(time, I80, label='A80')

#t_index0 = 150

#val72 = []
#val73 = []
#val74 = []
#val75 = []
#for step in range(26):
#    t_index = int(t_index0 + step*6.7 + .5)
#    t_step = time[t_index]
#
#    val_step = NP.mean(I72[t_index - 2:t_index + 2])
#    val72.append(val_step)
#    pl0.plot([t_step, t_step],[val_step, val_step], 'o' , color='black')
#
#    val_step = NP.mean(I73[t_index - 2:t_index + 2])
#    val73.append(val_step)
#    pl0.plot([t_step, t_step],[val_step, val_step], 'o' , color='black')
#
#    val_step = NP.mean(I74[t_index - 2:t_index + 2])
#    val74.append(val_step)
#    pl0.plot([t_step, t_step],[val_step, val_step], 'o' , color='black')
#
#    val_step = NP.mean(I75[t_index - 2:t_index + 2])
#    val75.append(val_step)
#    pl0.plot([t_step, t_step],[val_step, val_step], 'o' , color='black')
#
pl0.legend()
pl0.grid()

#inPower = NP.linspace(-3,22,26)
#inPower = 10**(0.1*inPower)/1e3
#PL.figure()
#PL.title('2019 03 26, digital receiver 11')
#PL.xlabel('input power [mW]')
#PL.ylabel('output power [arbitrary]')
#PL.plot(inPower, val72, label='A77')
#PL.plot(inPower, val73, label='A78')
#PL.plot(inPower, val74, label='A79')
#PL.plot(inPower, val75, label='A80')
#PL.legend()
#PL.grid()

