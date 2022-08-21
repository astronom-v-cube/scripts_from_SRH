# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 01:47:01 2016

@author: Sergey
"""
import os, fnmatch;
import numpy as NP;
import datetime as DT;
import pylab as PL;
import sys;
from astropy.io import fits
from PyQt5 import QtGui, QtCore, QtWidgets, QtNetwork
from PyQt5.QtWidgets import QFileDialog
from matplotlib.ticker import MultipleLocator, FormatStrFormatter;



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
    
if sys.platform == 'linux':
    font = QtGui.QFont()
    application.setFont(QtGui.QFont(font.defaultFamily(),8));

freqRange = '12-24'
fitNames = findFits('/home/svlesovoi/mSSRT/20190530/' + freqRange + '/','srh_20190530T*' + '*.fit')
#fitNames = findFits('/home/svlesovoi/mSSRT/20190529/skySun','srh_20190529T*' + '*.fit')
fitNames.sort();
fitNames = fitNames[1:];
freqAmount = 0
scanAmount = 0
visAmount = 0
for fName in fitNames:
    testFits = fits.open(fName)
    if (fName == fitNames[0]):
        freqAmount  = testFits[1].data['frequency'].shape[0]
        times = testFits[1].data['time']
        scanAmount = times.shape[1]
        vis = testFits[1].data['vis']
        visAmount = vis.shape[1] // scanAmount
        lcpVis = NP.reshape(vis,(freqAmount,scanAmount,visAmount))
        vis = 0
    else:
        times = NP.concatenate((times,testFits[1].data['time']),axis=1);
        vis = testFits[1].data['vis']
        lcpVis = NP.concatenate((lcpVis,NP.reshape(vis,(freqAmount,scanAmount,visAmount))),axis=1)
    testFits.close();

ind = NP.where(lcpVis[0,:,5].real > 100)
time81 = times[0,ind[0]]
vis81 = lcpVis[0,ind[0],4]
ind = NP.where(lcpVis[0,:,4].real > 100)
time83 = times[0,ind[0]]
vis83 = lcpVis[0,ind[0],5]

dt_major = 120
dt_minor = 30

#PL.plot(time81,vis81)
#PL.plot(time83,vis83)
#
sunInd0 = 10300
sunInd1 = 10550
skyInd0 = 11300
skyInd1 = 11550
#PL.plot(time81[sunInd0:sunInd1],vis81[sunInd0:sunInd1],'o')
#PL.plot(time83[sunInd0:sunInd1],vis83[sunInd0:sunInd1],'o')

fig1 = PL.figure(1,figsize=(20,12));
fig1.suptitle('MICRAN equipment test',fontsize='large');
sp1 = fig1.add_subplot(1,1,1);
sp1.set_ylabel('arbitrary');
sp1.set_xlabel('UT');
sp1.xaxis.set_major_locator(MultipleLocator(dt_major));
sp1.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
sp1.xaxis.set_minor_locator(MultipleLocator(dt_minor));
sp1.set_ylim(0,100)
sp1.plot(times[0,:],lcpVis[0,:,4].real)
sp1.plot(times[0,:],lcpVis[0,:,5].real)
sp1.plot(times[0,sunInd0:sunInd1],lcpVis[0,sunInd0:sunInd1,5].real,'.')
sp1.plot(times[0,skyInd0:skyInd1],lcpVis[0,skyInd0:skyInd1,5].real,'.')
#sp1.plot(times[0,:],lcpVis[0,:,0].real)
#sp1.plot(times[0,:],lcpVis[0,:,0].imag,'.')
skySun = lcpVis[0,sunInd0:sunInd1,5].real.mean() - lcpVis[0,skyInd0:skyInd1,5].real.mean()
skyNoise = lcpVis[0,skyInd0:skyInd1,5].real.std()
sunNoise = lcpVis[0,sunInd0:sunInd1,5].real.std()
