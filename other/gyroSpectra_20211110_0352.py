#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 11:49:23 2021

@author: sergey_lesovoi
"""

import numpy as NP
import pylab as PL
import ZirinTb

def freqGHzFormat(f, pos):
    return '%3.2f'%(f * 0.4 + 5.8)

Zi = ZirinTb.ZirinTb()

lcpImages = NP.zeros((16,512,512))
rcpImages = NP.zeros((16,512,512))

for ff in range(16):
    phaseEdit.onFrequencyChannelChanged(ff)
    phaseEdit.onCenterDisk()
    lcpImages[ff] = phaseEdit.lcpData
    skyLevel = lcpImages[ff,470:500,200:300].mean()
    sunLevel = lcpImages[ff,390:420,220:270].mean()
    qSunTb = 0.5*Zi.getTbAtFrequncy(phaseEdit.srhFits.freqList[ff]*1e-6)*1e3
    lcpImages[ff] = (lcpImages[ff] - skyLevel)/(sunLevel - skyLevel) * qSunTb

    rcpImages[ff] = phaseEdit.rcpData
    skyLevel = rcpImages[ff,470:500,200:300].mean()
    sunLevel = rcpImages[ff,390:420,220:270].mean()
    qSunTb = 0.5*Zi.getTbAtFrequncy(phaseEdit.srhFits.freqList[ff]*1e-6)*1e3
    rcpImages[ff] = (rcpImages[ff] - skyLevel)/(sunLevel - skyLevel) * qSunTb

fig = PL.figure()
pl = fig.subplots(4,16)

lcpMaxT1 = []
rcpMaxT1 = []
lcpMaxT2 = []
rcpMaxT2 = []
for ii in range(16):
    pl[0,ii].imshow(lcpImages[ii,300:320,390:410],origin='lower',vmin=0,vmax=2e5)
    pl[0,ii].axis('off')
    pl[1,ii].imshow(rcpImages[ii,300:320,390:410],origin='lower',vmin=0,vmax=2e5)
    pl[1,ii].axis('off')
    pl[2,ii].imshow(lcpImages[ii,140:170,204:234],origin='lower',vmin=0,vmax=2e5)
    pl[2,ii].axis('off')
    pl[3,ii].imshow(rcpImages[ii,140:170,204:234],origin='lower',vmin=0,vmax=2e5)
    pl[3,ii].axis('off')

    lcpMaxT1.append(lcpImages[ii,300:320,390:410].max())
    rcpMaxT1.append(rcpImages[ii,300:320,390:410].max())
    lcpMaxT2.append(lcpImages[ii,140:170,204:234].max())
    rcpMaxT2.append(rcpImages[ii,140:170,204:234].max())
    
fig = PL.figure(figsize=(12,8))
pl = fig.subplots(1)
pl.xaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
pl.xaxis.set_major_locator(PL.MultipleLocator(1.0))

fig.suptitle(phaseEdit.srhFits.dateObs)
pl.plot(rcpMaxT1,'.',color='red', label='RCP NW')
pl.plot(lcpMaxT1,'.',color='blue', label='LCP NW')
pl.plot(rcpMaxT2,color='red', label='RCP SE')
pl.plot(lcpMaxT2,color='blue', label='LCP SE')
pl.set_xlabel('GHz')
pl.set_ylabel('Tb [K]')
pl.legend()
pl.grid()
