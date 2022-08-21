#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 01:20:42 2022

@author: sergey_lesovoi
"""

import json
import pylab as PL
import numpy as NP

EW_base_names = []

for ant in range(24):
    EW_base_names.append('%02d-%02d'%(24-ant,24-ant-1))

#EW_names.append('C1')
    
for ant in range(24):
    EW_base_names.append('%02d-%02d'%(ant,ant + 1))

antPosFiles = ['2022-04-03', '2022-04-04', '2022-04-24', '2022-04-26']

antPosDX = []
antPosDY = []
for filNam in antPosFiles:
    fil = open(filNam + '.json')
    jAntPos = json.load(fil)
    antPosDX.append(jAntPos.get('dx'))
    antPosDY.append(jAntPos.get('dy'))
    fil.close()

fig = PL.figure()
pl = fig.subplots(nrows=2,ncols=1)
for day in range(len(antPosDX)):
    pl[0].plot(EW_base_names, antPosDX[day],label=antPosFiles[day])
    pl[1].plot(EW_base_names, antPosDY[day],label=antPosFiles[day])
    pl[0].plot(EW_base_names, antPosDX[day],'.')
    pl[1].plot(EW_base_names, antPosDY[day],'.')
pl[0].set_ylabel('dx [mm]')
pl[0].set_ylim(-50,50)
pl[0].legend()
pl[0].grid()
pl[0].tick_params(axis='x',rotation=45)

pl[1].set_ylabel('dy [mm]')
pl[1].set_ylim(-50,50)
pl[1].legend()
pl[1].grid()
pl[1].tick_params(axis='x',rotation=45)
    