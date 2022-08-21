#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 02:49:34 2021

@author: svlesovoi
"""
import srh36Utils
from srhFitsFile36 import SrhFitsFile
import pylab as PL
import numpy as NP

fitPath = 'SRH36_temp_20210914_1/'
fitNames = srh36Utils.findFits(fitPath,'*.fit')

for fileName in fitNames:
    if (fileName == fitNames[0]):
        fitFiles = SrhFitsFile(fileName,1025)
    else:
        fitFiles.append(fileName)

nRedLcpVis = fitFiles.visLcp[:,:,3007:3007+31]
nRedRcpVis = fitFiles.visRcp[:,:,3007:3007+31]
nLrRatio = nRedLcpVis / nRedRcpVis

flag = NP.zeros(31)
for pair in range(30):
    if (NP.abs(nRedLcpVis[0,:,pair]).mean() < 0.05 or NP.abs(nRedRcpVis[0,:,pair]).mean() < 0.05):
        flag[pair] = 1

PL.figure()
for pair in range(30):
    PL.plot(NP.abs(nRedLcpVis[0,:,pair]),NP.abs(nRedRcpVis[0,:,pair]),'.')

PL.figure()
for pair in range(30):
    if (not flag[pair]):
        PL.plot(NP.angle(nRedLcpVis[0,:,pair]),NP.angle(nRedRcpVis[0,:,pair]),'.')


