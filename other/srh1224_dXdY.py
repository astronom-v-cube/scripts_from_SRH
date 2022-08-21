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

ap0328 = SrhAntennaPositions('2022-03-28')
bl0328_minus = ap0328.extractBaselines('1224','minus_pi_4')
bl0328_zero = ap0328.extractBaselines('1224','zero')
bl0328_plus = ap0328.extractBaselines('1224','plus_pi_4')

ap0329 = SrhAntennaPositions('2022-03-29')
bl0329_minus = ap0329.extractBaselines('1224','minus_pi_4')
bl0329_zero = ap0329.extractBaselines('1224','zero')
bl0329_plus = ap0329.extractBaselines('1224','plus_pi_4')

ap0401 = SrhAntennaPositions('2022-04-01')
bl0401_minus = ap0401.extractBaselines('1224','minus_pi_4')
bl0401_zero = ap0401.extractBaselines('1224','zero')
bl0401_plus = ap0401.extractBaselines('1224','plus_pi_4')

ap0403 = SrhAntennaPositions('2022-04-03')
bl0403_minus = ap0403.extractBaselines('1224','minus_pi_4')
bl0403_zero = ap0403.extractBaselines('1224','zero')
bl0403_plus = ap0403.extractBaselines('1224','plus_pi_4')

ap0405 = SrhAntennaPositions('2022-04-05')
bl0405_minus = ap0405.extractBaselines('1224','minus_pi_4')
bl0405_zero = ap0405.extractBaselines('1224','zero')
bl0405_plus = ap0405.extractBaselines('1224','plus_pi_4')

ap0406 = SrhAntennaPositions('2022-04-06')
bl0406_minus = ap0406.extractBaselines('1224','minus_pi_4')
bl0406_zero = ap0406.extractBaselines('1224','zero')
bl0406_plus = ap0406.extractBaselines('1224','plus_pi_4')

slopeScale = 1
sl0328_minus = ap0328.calcSlopes(bl0328_minus)*slopeScale
sl0328_zero = ap0328.calcSlopes(bl0328_zero)*slopeScale
sl0328_plus = ap0328.calcSlopes(bl0328_plus)*slopeScale

sl0329_minus = ap0329.calcSlopes(bl0329_minus)*slopeScale
sl0329_zero = ap0329.calcSlopes(bl0329_zero)*slopeScale
sl0329_plus = ap0329.calcSlopes(bl0329_plus)*slopeScale

sl0401_minus = ap0401.calcSlopes(bl0401_minus)*slopeScale
sl0401_zero = ap0401.calcSlopes(bl0401_zero)*slopeScale
sl0401_plus = ap0401.calcSlopes(bl0401_plus)*slopeScale

sl0403_minus = ap0403.calcSlopes(bl0403_minus)*slopeScale
sl0403_zero = ap0403.calcSlopes(bl0403_zero)*slopeScale
sl0403_plus = ap0403.calcSlopes(bl0403_plus)*slopeScale

sl0405_minus = ap0405.calcSlopes(bl0405_minus)*slopeScale
sl0405_plus = ap0405.calcSlopes(bl0405_plus)*slopeScale

sl0406_minus = ap0406.calcSlopes(bl0406_minus)*slopeScale
sl0406_plus = ap0406.calcSlopes(bl0406_plus)*slopeScale

sl0328_minus[19:21] = 0
sl0328_plus[19:21] = 0
sl0329_minus[19:21] = 0
sl0329_plus[19:21] = 0
sl0401_minus[19:21] = 0
sl0401_plus[19:21] = 0
sl0403_minus[19:21] = 0
sl0403_plus[19:21] = 0
sl0405_minus[19:21] = 0
sl0405_plus[19:21] = 0
sl0406_minus[19:21] = 0
sl0406_plus[19:21] = 0

dX0328 = (sl0328_minus - sl0328_plus)*0.707
dX0329 = (sl0329_minus - sl0329_plus)*0.707
dX0401 = (sl0401_minus - sl0401_plus)*0.707
dX0403 = (sl0403_minus - sl0403_plus)*0.707
dX0405 = (sl0405_minus - sl0405_plus)*0.707
dX0406 = (sl0406_minus - sl0406_plus)*0.707
dXmean0 = (dX0328 + dX0329 + dX0401 + dX0403) / 4
dXmean1 = (dX0405 + dX0406) / 2

dY0328 = (sl0328_minus + sl0328_plus)*0.707
dY0329 = (sl0329_minus + sl0329_plus)*0.707
dY0401 = (sl0401_minus + sl0401_plus)*0.707
dY0403 = (sl0403_minus + sl0403_plus)*0.707
dY0405 = (sl0405_minus + sl0405_plus)*0.707
dY0406 = (sl0406_minus + sl0406_plus)*0.707
dYmean0 = (dY0328 + dY0329 + dY0401 + dY0403) / 4
dYmean1 = (dY0405 + dY0406) / 2

dX0328_src = bl0328_minus[1].mean(axis=1) + bl0328_plus[1].mean(axis=1) - bl0328_zero[1].mean(axis=1)
dX0329_src = bl0329_minus[1].mean(axis=1) + bl0329_plus[1].mean(axis=1) - bl0329_zero[1].mean(axis=1)
dX0401_src = bl0401_minus[1].mean(axis=1) + bl0401_plus[1].mean(axis=1) - bl0401_zero[1].mean(axis=1)
dX0403_src = bl0403_minus[1].mean(axis=1) + bl0403_plus[1].mean(axis=1) - bl0403_zero[1].mean(axis=1)
dX_src_mean = (dX0328_src + dX0329_src + dX0401_src + dX0403_src) / 4

dX0405_src = bl0405_minus[1].mean(axis=1) + bl0405_plus[1].mean(axis=1) - bl0405_zero[1].mean(axis=1)
dX0406_src = bl0406_minus[1].mean(axis=1) + bl0406_plus[1].mean(axis=1) - bl0406_zero[1].mean(axis=1)
dX_src_mean1 = (dX0405_src + dX0406_src) / 4

dY0328_src = bl0328_minus[1].mean(axis=1) - bl0328_plus[1].mean(axis=1)
dY0329_src = bl0329_minus[1].mean(axis=1) - bl0329_plus[1].mean(axis=1)
dY0401_src = bl0401_minus[1].mean(axis=1) - bl0401_plus[1].mean(axis=1)
dY0403_src = bl0403_minus[1].mean(axis=1) - bl0403_plus[1].mean(axis=1)
dY_src_mean = (dY0328_src + dY0329_src + dY0401_src + dY0403_src) / 4

dY0405_src = bl0405_minus[1].mean(axis=1) - bl0405_plus[1].mean(axis=1)
dY0406_src = bl0406_minus[1].mean(axis=1) - bl0406_plus[1].mean(axis=1)
dY_src_mean1 = (dY0405_src + dY0406_src) / 2

fig = PL.figure()
pl = fig.subplots(nrows=2,ncols=1)

pl[0].plot(EW_names,dX0328,'.',label='0328')
pl[0].plot(EW_names,dX0329,'.',label='0329')
pl[0].plot(EW_names,dX0401,'.',label='0401')
pl[0].plot(EW_names,dX0403,'.',label='0403')
pl[0].plot(EW_names,dXmean0)
pl[0].plot(EW_names,dXmean1)
pl[0].plot([EW_names[0],EW_names[-1]],[0.005,0.005], '--', color='black',linewidth=0.6)
pl[0].plot([EW_names[0],EW_names[-1]],[-0.005,-0.005], '--', color='black',linewidth=0.6)
pl[0].set_xticklabels(EW_names, rotation=45)
pl[0].set_ylabel('dX [meter]')
pl[0].set_xlabel('baseline')
pl[0].set_ylim(-0.03,0.03)
pl[0].grid()
pl[0].legend()

pl[1].plot(EW_names,dY0328,'.',label='0328')
pl[1].plot(EW_names,dY0329,'.',label='0329')
pl[1].plot(EW_names,dY0401,'.',label='0401')
pl[1].plot(EW_names,dY0403,'.',label='0403')
pl[1].plot(EW_names,dYmean0)
pl[1].plot(EW_names,dYmean1)
pl[1].plot([EW_names[0],EW_names[-1]],[0.005,0.005], '--', color='black',linewidth=0.6)
pl[1].plot([EW_names[0],EW_names[-1]],[-0.005,-0.005], '--', color='black',linewidth=0.6)
pl[1].set_xticklabels(EW_names, rotation=45)
pl[1].set_ylabel('dX [meter]')
pl[1].set_xlabel('baseline')
pl[1].set_ylim(-0.03,0.03)
pl[1].grid()
pl[1].legend()
