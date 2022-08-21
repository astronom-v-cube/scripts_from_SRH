#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 00:32:46 2022

@author: sergeyvlesovoi
"""

from srhAntennaPositions import SrhAntennaPositions
import pylab as PL

EW_names = []

for ant in range(24):
    EW_names.append('%02d-%02d'%(24-ant,24-ant-1))

for ant in range(24):
    EW_names.append('%02d-%02d'%(ant,ant + 1))

ap0424 = SrhAntennaPositions('2022-04-24')
bl0424_minus = ap0424.extractBaselines('1224','minus_pi_4')
bl0224_plus = ap0424.extractBaselines('1224','plus_pi_4')

ap0425 = SrhAntennaPositions('2022-04-25')
bl0425_minus = ap0425.extractBaselines('1224','minus_pi_4')
bl0425_plus = ap0425.extractBaselines('1224','plus_pi_4')

ap0426 = SrhAntennaPositions('2022-04-26')
bl0426_minus = ap0426.extractBaselines('1224','minus_pi_4')
bl0426_plus = ap0426.extractBaselines('1224','plus_pi_4')

slopeScale = 1
#sl0424_minus = ap0424.calcSlopes(bl0424_minus)*slopeScale
sl0424_plus = ap0424.calcSlopes(bl0424_plus)*slopeScale

sl0425_minus = ap0425.calcSlopes(bl0425_minus)*slopeScale
sl0425_plus = ap0425.calcSlopes(bl0425_plus)*slopeScale

sl0426_minus = ap0426.calcSlopes(bl0426_minus)*slopeScale
sl0426_plus = ap0426.calcSlopes(bl0426_plus)*slopeScale

dX0425 = (sl0425_minus - sl0425_plus)*0.707
dX0426 = (sl0426_minus - sl0426_plus)*0.707
dXmean0 = (dX0425 + dX0426) / 2

dY0425 = (sl0425_minus + sl0425_plus)*0.707
dY0426 = (sl0426_minus + sl0426_plus)*0.707
dYmean0 = (dY0425 + dY0426) / 2

fig = PL.figure()
pl = fig.subplots(nrows=2,ncols=1)

pl[0].plot(EW_names,dX0425,'.',label='0425')
pl[0].plot(EW_names,dX0426,'.',label='0426')
pl[0].plot(EW_names,dXmean0)
pl[0].plot([EW_names[0],EW_names[-1]],[0.005,0.005], '--', color='black',linewidth=0.6)
pl[0].plot([EW_names[0],EW_names[-1]],[-0.005,-0.005], '--', color='black',linewidth=0.6)
pl[0].set_xticklabels(EW_names, rotation=45)
pl[0].set_ylabel('dX [meter]')
pl[0].set_xlabel('baseline')
pl[0].set_ylim(-0.03,0.03)
pl[0].grid()
pl[0].legend()

pl[1].plot(EW_names,dY0425,'.',label='0425')
pl[1].plot(EW_names,dY0426,'.',label='0426')
pl[1].plot(EW_names,dYmean0)
pl[1].plot([EW_names[0],EW_names[-1]],[0.005,0.005], '--', color='black',linewidth=0.6)
pl[1].plot([EW_names[0],EW_names[-1]],[-0.005,-0.005], '--', color='black',linewidth=0.6)
pl[1].set_xticklabels(EW_names, rotation=45)
pl[1].set_ylabel('dX [meter]')
pl[1].set_xlabel('baseline')
pl[1].set_ylim(-0.03,0.03)
pl[1].grid()
pl[1].legend()

jsonSaveDict = {}
jsonSaveDict['dX0425'] = dX0425.tolist()
jsonSaveDict['dX0426'] = dX0426.tolist()
jsonSaveDict['dY0425'] = dX0425.tolist()
jsonSaveDict['dY0426'] = dX0426.tolist()

with open('srh1224_dxdy_202204.json','w') as jsonSave:
    json.dump(jsonSaveDict,jsonSave)
