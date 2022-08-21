#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 28 01:06:33 2021

@author: svlesovoi
"""

import os, fnmatch;
from srhFitsFile36_amp import SrhFitsFile
import numpy as NP;
import pylab as PL;
from scipy.stats import linregress
import scipy.constants

Cf = scipy.constants.c / 1.467
L0 = 800.6
S0 = 795

eastFiberLength = NP.zeros(1+64)
eastFiberLength[2] = 801.29
eastFiberLength[3] = 802.57
eastFiberLength[4] = 808.1
eastFiberLength[5] = 807.88
eastFiberLength[6] = 808.0
eastFiberLength[7] = 808.3

eastFiberLength[9] = 801.5
eastFiberLength[10] = 800.4
eastFiberLength[11] = 802.57
eastFiberLength[12] = 807.9
eastFiberLength[13] = 801.94
eastFiberLength[14] = 801.94
eastFiberLength[15] = 807.25
eastFiberLength[16] = 807.67
eastFiberLength[17] = 803.56

eastFiberLength[35] = 807.1
eastFiberLength[36] = 802.14
eastFiberLength[37] = 807.78

eastFiberLength[50] = 807.46
eastFiberLength[51] = 807.25
eastFiberLength[52] = 808.1


northFiberLength = NP.zeros(1+31)
northFiberLength[30] = 800.78
northFiberLength[31] = 808.31

westFiberLength = NP.zeros(1+32)
westFiberLength[1] = 806.23
westFiberLength[2] = 806.39

westFiberLength[7] = 806.61
westFiberLength[9] = 807.35


def northAntennaFormat(x, pos):
    if x == 0:
        return 'C0'
    else:
        return 'N%d'%(x)

def eastAntennaFormat(x, pos):
    if x == 0:
        return 'C0'
    else:
        return 'E%d'%(x)

def westAntennaFormat(x, pos):
    if x == 0:
        return 'C0'
    else:
        return 'W%d'%(x)

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;

#fitNames = findFits('SRH36_temp_20210527_60/','*.fit')  
fitNames = findFits('SRH36_temp_20210529_60/','*.fit')  
fitNames.sort();
fitNames = fitNames[1:];

for fName in fitNames:
    if (fName == fitNames[0]):
        sF = SrhFitsFile(fName, 256)
    else:
        sF.append(fName)

northSlopes = NP.zeros(31)
eastSlopes = NP.zeros(64)# C0E1..E64
westSlopes = NP.zeros(32)#C0W1..W32

northAntDelay = NP.zeros(32)
eastAntDelay = NP.zeros(65)
westAntDelay = NP.zeros(33)

northAntLength = NP.zeros(32)
eastAntLength = NP.zeros(65)
westAntLength = NP.zeros(33)

samples = sF.visRcp.shape[1]

for scan in range(samples):
    phaSlope, intercept, r_value, p_value, std_err = linregress(sF.freqList*1e3, NP.unwrap(NP.angle(sF.visLcp[:,scan,32])))
    northSlopes[0] += phaSlope/(2*NP.pi)
    phaSlope, intercept, r_value, p_value, std_err = linregress(sF.freqList*1e3, NP.unwrap(NP.angle(sF.visRcp[:,scan,32])))
    northSlopes[0] += phaSlope/(2*NP.pi)
    for i in range(30):
        phaSlope, intercept, r_value, p_value, std_err = linregress(sF.freqList*1e3, NP.unwrap(NP.angle(sF.visLcp[:,scan,3007+i])))
        northSlopes[i+1] += phaSlope/(2*NP.pi)
        phaSlope, intercept, r_value, p_value, std_err = linregress(sF.freqList*1e3, NP.unwrap(NP.angle(sF.visRcp[:,scan,3007+i])))
        northSlopes[i+1] += phaSlope/(2*NP.pi)
    for i in range(64):
        phaSlope, intercept, r_value, p_value, std_err = linregress(sF.freqList*1e3, NP.unwrap(NP.angle(sF.visLcp[:,scan,3504+i])))
        eastSlopes[i] += phaSlope/(2*NP.pi)
        phaSlope, intercept, r_value, p_value, std_err = linregress(sF.freqList*1e3, NP.unwrap(NP.angle(sF.visRcp[:,scan,3504+i])))
        eastSlopes[i] += phaSlope/(2*NP.pi)
    for i in range(32):
        phaSlope, intercept, r_value, p_value, std_err = linregress(sF.freqList*1e3, NP.unwrap(NP.angle(sF.visLcp[:,scan,3503-i])))
        westSlopes[i] += phaSlope/(2*NP.pi)
        phaSlope, intercept, r_value, p_value, std_err = linregress(sF.freqList*1e3, NP.unwrap(NP.angle(sF.visRcp[:,scan,3503-i])))
        westSlopes[i] += phaSlope/(2*NP.pi)
            
northSlopes /= (2*samples)
eastSlopes /= (2*samples)
westSlopes /= (2*samples)

northAntLength[0] = L0
northAntDelay[0] = northAntLength[0]/Cf

for ant in range(31):
    northAntDelay[ant + 1] = northAntDelay[ant] + northSlopes[ant]
for ant in range(32):
    northAntLength[ant] = Cf*northAntDelay[ant]

fig = PL.figure(figsize = (6,4));
fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
fig.tight_layout()
fig.suptitle('C0N1..31');
pl = fig.add_subplot(1,1,1);
pl.set_ylabel('meter');
pl.set_xlabel('antenna');
pl.xaxis.set_major_locator(PL.MultipleLocator(4));
pl.xaxis.set_major_formatter(PL.FuncFormatter(northAntennaFormat));
pl.xaxis.set_minor_locator(PL.MultipleLocator(1));
pl.plot(northAntLength)
pl.plot(northAntLength,'o')
pl.plot(northAntDelay,'o')
pl.plot(northFiberLength,'.')
pl.set_ylim(790,810)
pl.grid()

eastAntDelay[0] = L0/Cf
eastAntLength[0] = L0

for ant in range(64):
    eastAntDelay[ant + 1] = eastAntDelay[ant] + eastSlopes[ant]
for ant in range(65):
    eastAntLength[ant] = Cf*eastAntDelay[ant]
#east correction
#eastAntLength[18:] += 10.5
#eastAntLength[51] += 5
#eastAntLength[52:] += 8.5
fig = PL.figure(figsize = (8,4));
fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
fig.tight_layout()
fig.suptitle('C0E1..64');
pl = fig.add_subplot(1,1,1);
pl.set_ylabel('meter');
pl.set_xlabel('antenna');
pl.xaxis.set_major_locator(PL.MultipleLocator(4));
pl.xaxis.set_major_formatter(PL.FuncFormatter(eastAntennaFormat));
pl.xaxis.set_minor_locator(PL.MultipleLocator(1));
pl.plot(eastAntLength)
pl.plot(eastAntLength,'o')
pl.plot(eastFiberLength,'.')
pl.set_ylim(790,810)
pl.grid()

westAntDelay[0] = L0/Cf
westAntLength[0] = L0

for ant in range(32):
    westAntDelay[ant + 1] = westAntDelay[ant] - westSlopes[ant]
for ant in range(33):
    westAntLength[ant] = Cf*westAntDelay[ant]

fig = PL.figure(figsize = (6,4));
fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
fig.tight_layout()
fig.suptitle('C0W1..32');
pl = fig.add_subplot(1,1,1);
pl.set_ylabel('meter');
pl.set_xlabel('antenna');
pl.xaxis.set_major_locator(PL.MultipleLocator(4));
pl.xaxis.set_major_formatter(PL.FuncFormatter(westAntennaFormat));
pl.xaxis.set_minor_locator(PL.MultipleLocator(1));
pl.plot(westAntLength)
pl.plot(westAntLength,'o')
pl.plot(westFiberLength,'.')
pl.set_ylim(790,810)
pl.grid()

wholeLengths = NP.concatenate((northAntLength, westAntLength, eastAntLength))
commonMin = wholeLengths.min()
northAntLength -= commonMin
westAntLength -= commonMin
eastAntLength -= commonMin
wholeLengths = NP.concatenate((northAntLength, westAntLength, eastAntLength))
commonMax = wholeLengths.max()
#northAntLength = commonMax - northAntLength
#westAntLength = commonMax - westAntLength
#eastAntLength = commonMax - eastAntLength

northAntDelay = NP.ceil((northAntLength/Cf)*1e12)
westAntDelay = NP.ceil((westAntLength/Cf)*1e12)
eastAntDelay = NP.ceil((eastAntLength/Cf)*1e12)

#T0 = 4.042e-6
#northAntDelay = NP.ceil((T0 - northAntDelay)*1e9)
#westAntDelay = NP.ceil((T0 - westAntDelay)*1e9)
#eastAntDelay = NP.ceil((T0 - eastAntDelay)*1e9)

#for ant in range(31):
#    print('N%d = %d'%(ant+1,northAntDelay[ant+1]))
#
#for ant in range(32):
#    print('W%d = %d'%(ant+1,westAntDelay[ant+1]))
#
#for ant in range(64):
#    print('E%d = %d'%(ant+1,eastAntDelay[ant+1]))
