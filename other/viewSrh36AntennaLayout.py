#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 01:16:47 2020

@author: svlesovoi
"""

import numpy as NP
import pylab as PL
from matplotlib.ticker import (MultipleLocator)

antennaLayoutFile = open('srh36antennaLayout.txt')
antennaLayoutText = antennaLayoutFile.readlines()
eastWestBases = []
eastWestOffsets = []
northBases = []
northOffsets = []
for i in range(32):
    eastWestBases.append(float(antennaLayoutText[i + 1].split(' ')[1]))
    eastWestOffsets.append(float(antennaLayoutText[i + 1].split(' ')[2]))
    
for i in range(64):
    eastWestBases.append(float(antennaLayoutText[i + 34].split(' ')[1]))
    eastWestOffsets.append(float(antennaLayoutText[i + 34].split(' ')[2]))
    
for i in range(33):
    northOffsets.append(-float(antennaLayoutText[i + 99].split(' ')[1]))
    northBases.append(float(antennaLayoutText[i + 99].split(' ')[2]))
    
northOffsetConst = northOffsets[0]
for i in range(33):
    northOffsets[i] -= northOffsetConst

westArm = []
eastArm = []
northArm = []
westArmIdeal = []
eastArmIdeal = []
northArmIdeal = []

westArm.append(-eastWestBases[31])
westArmIdeal.append(-9800.)
for i in NP.linspace(1,31,31, dtype='int'):
    westArm.append(westArm[i - 1] - eastWestBases[31 - i])
    westArmIdeal.append(westArmIdeal[i - 1] - 9800.)

eastArm.append(eastWestBases[32])
eastArmIdeal.append(9800.)
for i in NP.linspace(33,95,63, dtype='int'):
    eastArm.append(eastArm[i - 33] + eastWestBases[i])
    eastArmIdeal.append(eastArmIdeal[i - 33] + 9800.)

northArm.append(northBases[0])
northArmIdeal.append(9800.)
for i in NP.linspace(1,32,32, dtype='int'):
    northArm.append(northArm[i - 1] + northBases[i])
    northArmIdeal.append(northArmIdeal[i - 1] + 9800.)

westOffsets = list(reversed(eastWestOffsets[0:32]))
eastOffsets = eastWestOffsets[32:]

#northArmOffset = -30
#westArmOffset = -60
#eastArmOffset = -80
northArmOffset = 0
westArmOffset = 0
eastArmOffset = 0

fig, pl = PL.subplots()
pl.xaxis.set_major_locator(MultipleLocator(9800.))
pl.yaxis.set_major_locator(MultipleLocator(10.))
pl.grid()
#pl.plot(westArmIdeal, westArmIdeal, '+', color = 'black')
#pl.plot(eastArmIdeal, eastArmIdeal, '+', color = 'black')
#pl.plot(westArmIdeal, westArm, '.')
#pl.plot(eastArmIdeal, eastArm, '.')
pl.plot(westArmIdeal, NP.array(westArmIdeal) - NP.array(westArm) - northArmOffset - westArmOffset, '.')
pl.plot(eastArmIdeal, NP.array(eastArmIdeal) - NP.array(eastArm) - northArmOffset - eastArmOffset, '.')
pl.plot([eastArmIdeal[-1], westArmIdeal[-1]],[100,100], color='gray')
pl.plot([eastArmIdeal[-1], westArmIdeal[-1]],[-100,-100], color='gray')

fig, pl1 = PL.subplots()
pl1.xaxis.set_major_locator(MultipleLocator(9800.))
pl1.yaxis.set_major_locator(MultipleLocator(9800.))
pl1.grid()
pl1.plot(westArm, westOffsets, '.')
pl1.plot(eastArm, eastOffsets, '.')
pl1.plot(northOffsets, northArm,  '.')
pl1.plot([-northOffsetConst, -northOffsetConst],[0,northArm[-1]])
pl1.plot(NP.array(westArm) + westArmOffset, westOffsets, '+')
pl1.plot(NP.array(eastArm) + eastArmOffset, eastOffsets, '+')
pl1.plot(NP.array(northOffsets) - northArmOffset, northArm, '+')
