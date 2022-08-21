#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 10:53:41 2021

@author: sergeyvlesovoi
"""

import os, fnmatch;
from srhFitsFile36 import SrhFitsFile
import numpy as NP;
import datetime as DT;
import pylab as PL;
import sys;
import matplotlib;
import phaMatrixGen

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;

dateName = DT.datetime.now().strftime("%Y%m%d");
fitNames = findFits('/home/sergeyvlesovoi/SRH36/20210329','*.fit')
fitNames.sort();
fitNames = fitNames[1:];

visScale = 1/(2e6*49)
ampScale = visScale / 128

inputFits = SrhFitsFile(fitNames[-3],2049)
inputFits.setCalibIndex(5)
#inputFits.useNonlinearApproach = False
inputFits.flagAntennasByPhaseNoise(0.1)
inputFits.calculatePhaseCalibration(3)

inputFits.setFrequencyChannel(0)
inputFits.vis2uv(0,True,False)
inputFits.uv2lmImage()
PL.figure()
PL.imshow(inputFits.lcp.real,cmap='hot')
PL.title('lm 0')

inputFits.setFrequencyChannel(5)
inputFits.vis2uv(0,True,False)
inputFits.uv2lmImage()
PL.figure()
PL.imshow(inputFits.lcp.real,cmap='hot')
PL.title('lm 5')

inputFits.vis2uvDoubleFrequency(0,True,False,zeroChannelScale = 0.6)
PL.figure()
PL.imshow(NP.sqrt(NP.abs(inputFits.uvLcp)),cmap='hot')
PL.title('uv 0+5')

inputFits.uv2lmImage()
PL.figure()
PL.imshow(inputFits.lcp.real,cmap='hot')
PL.title('lm 0+5')

