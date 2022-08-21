#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 08:03:11 2020

@author: svlesovoi
"""
import numpy as NP
import pylab as PL

freq = 6930 #MHz
ms2arcmin = 1e-3/60*15
L = 4700 #len(micran_corrdata_text)
S = 1/(2000000*49)

#plotName = 'plot_10.03.2020_11.33.18.192.txt'
#corrPlotName = 'plotCorr_10.03.2020_11.33.17.083.txt'
plotName = 'plot_25.02.2020_17.14.02.127.txt'
corrPlotName = 'plotCorr_25.02.2020_17.14.03.174.txt'

micran_data_file = open(plotName)
micran_data_text = micran_data_file.readlines()
micran_data_text = micran_data_text[2:]

micran_corrdata_file = open(corrPlotName)
micran_corrdata_text = micran_corrdata_file.readlines()
micran_corrdata_text = micran_corrdata_text[2:]

micran_corrdata = NP.zeros((10,L))
for i in range(L):
    data_strs = micran_corrdata_text[i].split(';')
    for j in range(10):
        micran_corrdata[j,i] = float(data_strs[j])


micran_data = NP.zeros((14,L))
for i in range(L):
    data_strs = micran_data_text[i].split(';')
    for j in range(14):
        micran_data[j,i] = float(data_strs[j])

PL.figure()
pl0 = PL.subplot(211)
PL.suptitle('%s, %d MHz' % (plotName.split('_')[1], freq))

lcpA = micran_data[1]*S
lcpB = micran_data[2]*S
lcpI = NP.sin(NP.pi/2*micran_corrdata[1]*S) * NP.sqrt(lcpA) * NP.sqrt(lcpB)
lcpQ = NP.sin(NP.pi/2*micran_corrdata[2]*S) * NP.sqrt(lcpA) * NP.sqrt(lcpB)
lcpEnv = NP.sqrt(lcpI**2. + lcpQ**2.)

pl0.plot(micran_corrdata[0]*ms2arcmin, lcpI, label='LCP I')
pl0.plot(micran_corrdata[0]*ms2arcmin, lcpQ, label='LCP Q')
pl0.plot(micran_corrdata[0]*ms2arcmin, lcpEnv, label='LCP Amp')
pl0.xaxis.set_major_locator(PL.MultipleLocator(60))
pl0.grid()
pl0.legend()

lcpA -= lcpA.min()
lcpB -= lcpB.min()
lcpEnv -= lcpEnv.min()
lcpAB = NP.sqrt(lcpA*lcpB)

pl1 = PL.subplot(212)
pl1.plot(micran_data[0]*ms2arcmin, lcpA / lcpA.max(), label='LCP A')
pl1.plot(micran_data[0]*ms2arcmin, lcpB / lcpB.max(), label='LCP B')
pl1.plot(micran_data[0]*ms2arcmin, lcpEnv / lcpEnv.max(), label='LCP Envelope')
pl1.plot(micran_data[0]*ms2arcmin, lcpAB / lcpAB.max(), label='LCP AB')
pl1.xaxis.set_major_locator(PL.MultipleLocator(60))
pl1.set_xlabel('arcmin')
pl1.grid()
pl1.legend()
