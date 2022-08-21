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
  t -= mm*60.;
  ss = (int)(t);
  return '%02d:%02d:%02d' % (hh,mm,ss);

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
#fitNames = findFits('/home/svlesovoi/SRH/2019/03/PowerTest/26/four/','mf_' + '*.fit')
#fitNames = findFits('/home/svlesovoi/SRH/2019/03/PowerTest/27/','mf_' + '*.fit')
#fitNames = findFits('/home/svlesovoi/SRH/2019/03/PowerTest/27/sun/','mf_' + '*.fit')
#fitNames = findFits('/home/svlesovoi/SRH/2019/03/PowerTest/28/','mf_' + '*.fit')
fitNames = findFits('/home/svlesovoi/SRH/2019/03/PowerTest/29/4','mf_' + '*.fit')

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
LCP = NP.zeros(dataLength)
RCP = NP.zeros(dataLength)
for l in range(time.shape[0]):
    for f in range(srhFits.freqListLength):
        times[l*srhFits.freqListLength + f] = srhFits.freqTime[f,l]
        LCP[l*srhFits.freqListLength + f] = srhFits.ampLcp[f,l,ant78Index]
        RCP[l*srhFits.freqListLength + f] = srhFits.ampRcp[f,l,ant78Index]
        
I83 = (LCP + RCP)*.5
V83 = (LCP - RCP)*.5

t0 = 300*srhFits.freqListLength
t1 = t0 + 2000
t_sun_0 = t0 + 0*srhFits.freqListLength
t_sun_1 = t0 + 30*srhFits.freqListLength
t_sky_0 = t0 + 110*srhFits.freqListLength
t_sky_1 = t0 + 140*srhFits.freqListLength
t_zero_0 = t0 + 380*srhFits.freqListLength
t_zero_1 = t0 + 390*srhFits.freqListLength
PL.figure()
pl0 = PL.subplot(111)
PL.title('%s , Stokes I, V at %d MHz'%(srhFits.dateObs.split('T')[0], srhFits.freqList[freqIndex]))
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
sky_mean = NP.mean(I83[t_sky_0:t_sky_1])
zero_mean = NP.mean(I83[t_zero_0:t_zero_1])
sun_mean = NP.mean(I83[t_sun_0:t_sun_1]) - sky_mean
sun_std = NP.std(I83[t_sun_0:t_sun_1])
pl0.plot(times[t0:t1], I83[t0:t1], label='I83')
pl0.plot(times[t0:t1], V83[t0:t1], label='V83')
pl0.plot(times[t_sky_0:t_sky_1], I83[t_sky_0:t_sky_1], 'o', label=('sky %0.2f'%sky_mean))
pl0.plot(times[t_sun_0:t_sun_1], I83[t_sun_0:t_sun_1], 'o', label=('Sun %0.2f, N/S %0.3f'%(sun_mean, sun_std/sun_mean)))
pl0.plot(times[t_zero_0:t_zero_1], I83[t_zero_0:t_zero_1], 'o', label=('zero %0.2f'%zero_mean))

#pl0.plot(times[t0:], LCP[t0:])
#pl0.plot(times[t0:], RCP[t0:])

pl0.legend()
pl0.grid()
