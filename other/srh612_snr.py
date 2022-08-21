#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 14:07:52 2021

@author: sergey_lesovoi
"""

import numpy as NP
import pylab as PL
import ZirinTb
import os

def freqGHzFormat(f, pos):
    return '%3.2f'%(f * 0.4 + 5.8)

def hms_format(t, pos):
  hh = int(t / 3600.)
  t -= hh*3600.
  mm = int(t / 60.)
  t -= mm*60.
  ss = int(t)
  return '%02d:%02d:%02d' % (hh,mm,ss)

def closurePhaseLcp(freq, pair):
    if (pair < 127):
        return NP.angle(visLcp[freq,:,pair]) + NP.angle(visLcp[freq,:,pair+1]) - NP.angle(visLcp2[freq,:,pair])
    else:
        return NP.angle(visLcp[freq,:,pair]) + NP.angle(visLcp[freq,:,pair+1]) - NP.angle(visLcp2[freq,:,pair - 1])

def closurePhaseRcp(freq, pair):
    if (pair < 127):
        return NP.angle(visRcp[freq,:,pair]) + NP.angle(visRcp[freq,:,pair+1]) - NP.angle(visRcp2[freq,:,pair])
    else:
        return NP.angle(visRcp[freq,:,pair]) + NP.angle(visRcp[freq,:,pair+1]) - NP.angle(visRcp2[freq,:,pair - 1])

Zi = ZirinTb.ZirinTb()

dP = 15
ANT_NUMBER = 192
VIS_NUMBER = 190

ampsLcp = phaseEdit.srhFits.ampLcp
ampsRcp = phaseEdit.srhFits.ampRcp
visLcp = NP.zeros((ampsRcp.shape[0],ampsRcp.shape[1],63+127),dtype='complex')
visRcp = NP.zeros((ampsRcp.shape[0],ampsRcp.shape[1],63+127),dtype='complex')
visLcp2 = NP.zeros((ampsRcp.shape[0],ampsRcp.shape[1],62+126),dtype='complex')
visRcp2 = NP.zeros((ampsRcp.shape[0],ampsRcp.shape[1],62+126),dtype='complex')

visLcp[:,:,0:127] = phaseEdit.srhFits.visLcp[:,:,10208:10208+127]
visLcp[:,:,127:127+63] = phaseEdit.srhFits.visLcp[:,:,8192:8192+63]
visRcp[:,:,0:127] = phaseEdit.srhFits.visRcp[:,:,10208:10208+127]
visRcp[:,:,127:127+63] = phaseEdit.srhFits.visRcp[:,:,8192:8192+63]

visLcp2[:,:,0:126] = phaseEdit.srhFits.visLcp[:,:,10335:10335+126]
visLcp2[:,:,126:126+62] = phaseEdit.srhFits.visLcp[:,:,8255:8255+62]
visRcp2[:,:,0:126] = phaseEdit.srhFits.visRcp[:,:,10335:10335+126]
visRcp2[:,:,126:126+62] = phaseEdit.srhFits.visRcp[:,:,8255:8255+62]

# fig = PL.figure()
# pl = fig.subplots(2,1)

# pl[0].imshow(ampsRcp[0].T,aspect=1/5)
# pl[1].imshow(ampsLcp[0].T,aspect=1/5)

fig = PL.figure(figsize=(20,4))
pl = fig.subplots(2,16)
PL.tight_layout()
for f in range(1):
    skyStairRcp = ampsRcp[f,-dP:-1].mean(axis=0)
    skyStairLcp = ampsLcp[f,-dP:-1].mean(axis=0)
    sunStairRcp = ampsRcp[f,0:dP] - skyStairRcp
    sunStairLcp = ampsLcp[f,0:dP] - skyStairLcp
    
    fPath = 'ant_test/f%d'%(int(phaseEdit.srhFits.freqList[f]*1e-3))
    os.mkdir(fPath)
    for a in range(ANT_NUMBER):
        fig = PL.figure(figsize=(8,5))
        ax = fig.subplots(1)
        fAmpLcp = ampsLcp[f,:,a] - skyStairLcp[a]
        fAmpRcp = ampsRcp[f,:,a] - skyStairRcp[a]
        ax.set_title(phaseEdit.srhFits.dateObs + ', Antenna %d, %.1f GHz'%(a+1, phaseEdit.srhFits.freqList[f]*1e-6))
        ax.xaxis.set_major_formatter(PL.FuncFormatter(hms_format))
        ax.xaxis.set_major_locator(PL.MultipleLocator(60))
        ax.plot(phaseEdit.srhFits.freqTime[f],fAmpLcp,label='LCP')
        ax.plot(phaseEdit.srhFits.freqTime[f],fAmpRcp,label='RCP')
        ax.set_ylim(-10,4e4)
        ax.set_xlabel('UT')
        ax.set_ylabel('arbitrary')
        ax.grid()
        ax.legend()
        PL.savefig(fPath + '/ant_sun_sky_%.1f_%d'%(phaseEdit.srhFits.freqList[f]*1e-6,a) + '.png')
        PL.close(fig)

    pl[0,f].plot(sunStairRcp.mean(axis=0)*1e-4)
    pl[0,f].plot(sunStairLcp.mean(axis=0)*1e-4)
    pl[0,f].set_title('%.1f GHz'%(phaseEdit.srhFits.freqList[f]*1e-6))
    pl[0,f].set_ylim(0,6)
    if (f > 0):
        pl[0,f].get_yaxis().set_ticklabels([])
    else:
        pl[0,f].set_ylabel('arbitrary')
    pl[0,f].grid()
    
    pl[1,f].plot(sunStairRcp.mean(axis=0)/sunStairRcp.std(axis=0))
    pl[1,f].plot(sunStairLcp.mean(axis=0)/sunStairRcp.std(axis=0))
    pl[1,f].set_ylim(0,1e3)
    if (f > 0):
        pl[1,f].get_yaxis().set_ticklabels([])
    else:
        pl[1,f].set_ylabel('SNR')
    pl[1,f].grid()


outAmpLcp = sunStairLcp.mean(axis=0)
outAmpRcp = sunStairRcp.mean(axis=0)
outSnrLcp = outAmpLcp / sunStairLcp.std(axis=0)
outSnrRcp = outAmpRcp / sunStairRcp.std(axis=0)

outFile = open('srhSNR_' + phaseEdit.srhFits.dateObs + '.txt','w')
outFile.write('Antenna LCP SNR RCP SNR\n')
for f in range(16):
    outFile.write('------------------------------- Frequency %.1f GHz ----------------------------------------\n'%(phaseEdit.srhFits.freqList[f]*1e-6))
    for a in range(192):
        outFile.write('%d %.1f %.1f %.1f %.1f\n'%(a+1,outAmpLcp[a],outSnrLcp[a],outAmpRcp[a],outSnrRcp[a]))

fig = PL.figure(figsize=(20,4))
pl = fig.subplots(2,16)
PL.tight_layout()
for f in range(16):
    skyStairRcp = NP.abs(visRcp[f,-dP:-1])
    skyStairLcp = NP.abs(visLcp[f,-dP:-1])
    sunStairRcp = NP.abs(visRcp[f,0:dP])
    sunStairLcp = NP.abs(visLcp[f,0:dP])
    snrLcp = sunStairLcp.mean(axis=0)/skyStairLcp.std(axis=0)
    snrRcp = sunStairRcp.mean(axis=0)/skyStairRcp.std(axis=0)
    
    fPath = 'vis_test/f%d'%(int(phaseEdit.srhFits.freqList[f]*1e-3))
    os.mkdir(fPath)
    for a in range(VIS_NUMBER):
        fig = PL.figure(figsize=(8,5))
        ax = fig.subplots(1)
        ax.set_title(phaseEdit.srhFits.dateObs + ', Visibility %d, %.1f GHz'%(a+1, phaseEdit.srhFits.freqList[f]*1e-6))
        ax.xaxis.set_major_formatter(PL.FuncFormatter(hms_format))
        ax.xaxis.set_major_locator(PL.MultipleLocator(60))
        ax.plot(phaseEdit.srhFits.freqTime[f],NP.abs(visLcp[f,:,a]),label='LCP SNR %d'%(int(snrLcp[a])))
        ax.plot(phaseEdit.srhFits.freqTime[f],NP.abs(visRcp[f,:,a]),label='RCP SNR %d'%(int(snrRcp[a])))
        ax.set_ylim(0,0.5)
        ax.set_xlabel('UT')
        ax.set_ylabel('arbitrary')
        ax.grid()
        ax.legend()
        PL.savefig(fPath + '/vis_sun_sky_%.1f_%d'%(phaseEdit.srhFits.freqList[f]*1e-6,a) + '.png')
        PL.close(fig)

    pl[0,f].plot(sunStairRcp.mean(axis=0))
    pl[0,f].plot(sunStairLcp.mean(axis=0))
    pl[0,f].set_title('%.1f GHz'%(phaseEdit.srhFits.freqList[f]*1e-6))
    pl[0,f].set_ylim(0,.5)
    pl[0,f].get_xaxis().set_ticklabels([])
    if (f > 0):
        pl[0,f].get_yaxis().set_ticklabels([])
    else:
        pl[0,f].set_ylabel('correlation')
    pl[0,f].grid()
    
    pl[1,f].plot(sunStairRcp.mean(axis=0)/skyStairRcp.std(axis=0))
    pl[1,f].plot(sunStairLcp.mean(axis=0)/skyStairRcp.std(axis=0))
    pl[1,f].set_ylim(0,1e3)
    if (f > 0):
        pl[1,f].get_yaxis().set_ticklabels([])
    else:
        pl[1,f].set_ylabel('SNR')
    pl[1,f].grid()
