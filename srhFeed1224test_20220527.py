#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 27 02:06:14 2022

@author: sergey_lesovoi
"""

rcpAmps0612 = phaseEdit612.srhFits.ampRcp.copy()
lcpAmps0612 = phaseEdit612.srhFits.ampLcp.copy()
freqs = phaseEdit612.srhFits.freqList * 1e3

PL.figure()

ant = 129
f=15
SNR = (rcpAmps0612[f,0:100,ant].mean() - rcpAmps0612[f,170:190,ant].mean())/rcpAmps0612[f,170:190,ant].mean()
PL.plot(rcpAmps0612[f,:,ant],label='antenna %d (1224), SNR %.1f'%(ant,SNR))
ant = 128
SNR = (rcpAmps0612[f,0:100,ant].mean() - rcpAmps0612[f,170:190,ant].mean())/rcpAmps0612[f,170:190,ant].mean()
PL.plot(rcpAmps0612[f,:,ant],label='antenna %d (0612), SNR %.1f'%(ant,SNR))
PL.legend()
#PL.ylim(0,1e5)

PL.figure()
ant = 129
PL.plot(freqs,10*NP.log10(rcpAmps0612[:,0:100,ant].mean(axis=1) - rcpAmps0612[:,170:190,ant].mean(axis=1)),label='antenna %d (1224)'%ant)
ant = 128
PL.plot(freqs,10*NP.log10(rcpAmps0612[:,0:100,ant].mean(axis=1) - rcpAmps0612[:,170:190,ant].mean(axis=1)),label='antenna %d (0612)'%ant)
PL.legend()
ant = 130
PL.plot(freqs,10*NP.log10(rcpAmps0612[:,0:100,ant].mean(axis=1) - rcpAmps0612[:,170:190,ant].mean(axis=1)),label='antenna %d (0612)'%ant)
PL.ylabel('sun a.u. dB')
PL.xlabel('Hz')
PL.grid()
PL.legend()

PL.figure()
ant = 129
PL.plot(freqs,10*NP.log10(rcpAmps0612[:,170:190,ant].mean(axis=1)),label='antenna %d (1224)'%ant)
ant = 128
PL.plot(freqs,10*NP.log10(rcpAmps0612[:,170:190,ant].mean(axis=1)),label='antenna %d (0612)'%ant)
PL.legend()
ant = 130
PL.plot(freqs,10*NP.log10(rcpAmps0612[:,170:190,ant].mean(axis=1)),label='antenna %d (0612)'%ant)
PL.ylabel('sky a.u. dB')
PL.xlabel('Hz')
PL.grid()
PL.legend()

PL.figure()
ant = 129
f=14
PL.plot(rcpAmps0612[f,170:190,:].mean(axis=0),'.', label='f = %.1f GHz'%(freqs[f]*1e-9),color='red')
PL.plot([ant],rcpAmps0612[f,170:190,ant:ant+1].mean(axis=0),'+',color='red')
f=15
PL.plot(rcpAmps0612[f,170:190,:].mean(axis=0),'.', label='f = %.1f GHz'%(freqs[f]*1e-9),color='blue')
PL.plot([ant],rcpAmps0612[f,170:190,ant:ant+1].mean(axis=0),'+',color='blue')

PL.ylabel('sky a.u.')
PL.xlabel('antenna index')
#PL.ylim(0,1e5)
PL.grid()
PL.legend()

freqResp = NP.zeros(16*40)
for ff in range(16):
    freqResp[ff*40:(ff+1)*40] = rcpAmps0612[ff,135:135+40,129]
PL.figure()
PL.plot(freqResp,label='feed 1224')
PL.legend()

t0 = 140
dt = 40
PL.figure()
resp = rcpAmps0612[0,t0:t0+dt,129]
PL.plot((resp - resp.min())/(resp.max() - resp.min()))
for ff in NP.arange(9,16,1):
    resp = rcpAmps0612[ff,t0:t0+dt,129]
    PL.plot((resp - resp.min())/(resp.max() - resp.min()),'.')
