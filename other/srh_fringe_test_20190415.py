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
fitNames = findFits('/home/svlesovoi/SRH/2019/04/15/','mf_' + '*.fit')

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
fig.suptitle('%s, 4, 5, 6 GHz, RCP, LCP'%(srhFits.dateObs.split('T')[0]))
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))

pl0.plot(times[t0:t1], LCP77_78[t0:t1], label='LCP77_78')
pl0.plot(times[t0:t1], RCP77_78[t0:t1], label='RCP77_78', color='gray', linewidth=0.5)

pl0.set_ylim(-0.3,0.3)
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
pl0.set_xlabel('UT')
pl0.set_ylabel('correlation')
pl1.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
pl1.set_xlabel('UT')
pl1.set_ylabel('smoothed correlation')
pl0.grid()
pl0.legend()


smoothInterval = 900
t_f4_sen0 = t0 + 10
t_f4_sen1 = t_f4_sen0 + smoothInterval
pl0.plot(times[t_f4_sen0:t_f4_sen1], LCP77_78[t_f4_sen0:t_f4_sen1], color='red')

t_f5_sen0 = t0 + 1100
t_f5_sen1 = t_f5_sen0 + smoothInterval
pl0.plot(times[t_f5_sen0:t_f5_sen1], LCP77_78[t_f5_sen0:t_f5_sen1], color='green')

t_f6_sen0 = t0 + 2200
t_f6_sen1 = t_f6_sen0 + smoothInterval
pl0.plot(times[t_f6_sen0:t_f6_sen1], LCP77_78[t_f6_sen0:t_f6_sen1], color='magenta')

smooth_LCP_f4 = smooth(LCP77_78[t_f4_sen0:t_f4_sen1], 9)
dev_LCP_f4 = LCP77_78[t_f4_sen0:t_f4_sen1] - smooth_LCP_f4
smooth_RCP_f4 = smooth(RCP77_78[t_f4_sen0:t_f4_sen1], 9)
dev_RCP_f4 = RCP77_78[t_f4_sen0:t_f4_sen1] - smooth_RCP_f4

smooth_LCP_f5 = smooth(LCP77_78[t_f5_sen0:t_f5_sen1], 9)
dev_LCP_f5 = LCP77_78[t_f5_sen0:t_f5_sen1] - smooth_LCP_f5
smooth_RCP_f5 = smooth(RCP77_78[t_f5_sen0:t_f5_sen1], 9)
dev_RCP_f5 = RCP77_78[t_f5_sen0:t_f5_sen1] - smooth_RCP_f5

smooth_LCP_f6 = smooth(LCP77_78[t_f6_sen0:t_f6_sen1], 9)
dev_LCP_f6 = LCP77_78[t_f6_sen0:t_f6_sen1] - smooth_LCP_f6
smooth_RCP_f6 = smooth(RCP77_78[t_f6_sen0:t_f6_sen1], 9)
dev_RCP_f6 = RCP77_78[t_f6_sen0:t_f6_sen1] - smooth_RCP_f6

badValuesIndex = 10

s_n_LCP_f4 = NP.max(smooth_LCP_f4[badValuesIndex:-badValuesIndex]) / NP.std(dev_LCP_f4[badValuesIndex:-badValuesIndex])
s_n_RCP_f4 = NP.max(smooth_RCP_f4[badValuesIndex:-badValuesIndex]) / NP.std(dev_RCP_f4[badValuesIndex:-badValuesIndex])
s_n_LCP_f5 = NP.max(smooth_LCP_f5[badValuesIndex:-badValuesIndex]) / NP.std(dev_LCP_f5[badValuesIndex:-badValuesIndex])
s_n_RCP_f5 = NP.max(smooth_RCP_f5[badValuesIndex:-badValuesIndex]) / NP.std(dev_RCP_f5[badValuesIndex:-badValuesIndex])
s_n_LCP_f6 = NP.max(smooth_LCP_f6[badValuesIndex:-badValuesIndex]) / NP.std(dev_LCP_f6[badValuesIndex:-badValuesIndex])
s_n_RCP_f6 = NP.max(smooth_RCP_f6[badValuesIndex:-badValuesIndex]) / NP.std(dev_RCP_f6[badValuesIndex:-badValuesIndex])

pl1.plot(times[t0:t1], LCP77_78[t0:t1]*0)
pl1.plot(times[t_f4_sen0+badValuesIndex:t_f4_sen1-badValuesIndex], smooth_LCP_f4[badValuesIndex:-badValuesIndex], label=('4 GHz, S/N = %0.2f'%s_n_LCP_f4), color='red')
pl1.plot(times[t_f4_sen0+badValuesIndex:t_f4_sen1-badValuesIndex], dev_LCP_f4[badValuesIndex:-badValuesIndex], color='red')
pl1.plot(times[t_f5_sen0+badValuesIndex:t_f5_sen1-badValuesIndex], smooth_LCP_f5[badValuesIndex:-badValuesIndex], label=('5 GHz, S/N = %0.2f'%s_n_LCP_f5), color='green')
pl1.plot(times[t_f5_sen0+badValuesIndex:t_f5_sen1-badValuesIndex], dev_LCP_f5[badValuesIndex:-badValuesIndex], color='green')
pl1.plot(times[t_f6_sen0+badValuesIndex:t_f6_sen1-badValuesIndex], smooth_LCP_f6[badValuesIndex:-badValuesIndex], label=('6 GHz, S/N = %0.2f'%s_n_LCP_f6), color='magenta')
pl1.plot(times[t_f6_sen0+badValuesIndex:t_f6_sen1-badValuesIndex], dev_LCP_f6[badValuesIndex:-badValuesIndex], color='magenta')

pl1.grid()
pl1.legend()
