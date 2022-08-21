#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  4 06:00:32 2021

@author: svlesovoi
"""

import numpy as NP
import pylab as PL

N = 256
D = 0.3
Ds = 0.2
sig_arg = NP.linspace(-1,1,N)
sig = NP.exp(-((sig_arg)/D)**2)
sig += NP.exp(-((sig_arg-.2)/Ds)**2)

src_dL = 41
src_arg = NP.linspace(-1,1,src_dL)
srcX, srcY = NP.meshgrid(src_arg, src_arg)
source = NP.exp(-((srcX/.3)**2 + (srcY/.3)**2))
radius = N/4

qSunModel = NP.zeros((N,N))

for i in range(N):
    for j in range(N):
        if NP.sqrt((i-N/2)**2 + (j-N/2)**2) < radius:
            qSunModel[i,j] = 1

#qSunModel[N//3:N//3 + src_dL,N//3:N//3 + src_dL] += source

flux_sig = []
corrplot_sig = []
flux_if_sig = []
sig_trace = []
if_sig_trace = []
#
#for i in range(100):
#    sig = NP.exp(-((sig_arg)/D)**2)
#    sig += 0.01*(100-i)*NP.exp(-((sig_arg-.2)/(Ds*(i+1)/N))**2)
#    sig_trace.append(sig)
#    f_sig = NP.fft.fft(sig)
#    if_sig = NP.fft.ifft(NP.abs(f_sig))
#    flux_f_sig.append(NP.abs(f_sig).sum())
#    flux_sig.append(sig.sum())
#    flux_if_sig.append(sig.sum())
#
#PL.figure()
#PL.imshow(NP.transpose(sig_trace),aspect=.05,origin='lower')
#PL.plot(NP.array(flux_sig)*flux_f_sig[0].real/flux_sig[0].real*16 - flux_f_sig[0].real*16)
#PL.plot(flux_f_sig - flux_f_sig[0])
#

L = 200
srcW = 0.01
for i in range(L):
    if i < L/10:
        source = NP.exp(-((srcX/(srcW*(i+1)))**2 + (srcY/(srcW*(i+1)))**2)) * 9.*(i)
    else:
        source = NP.exp(-((srcX/(srcW*(i+1)))**2 + (srcY/(srcW*(i+1)))**2)) * 1.*(L - i)
    source -= source.min()
    sig = qSunModel.copy()
    sig[N//3:N//3 + src_dL,N//3:N//3 + src_dL] += source
    sig_trace.append(sig[105])
#    sig_trace.append(sig)
    f_sig = NP.fft.fft2(sig)
    if_sig = NP.abs(NP.fft.ifft2(NP.abs(f_sig)))
    if_sig = NP.roll(if_sig,N//2,axis=0)
    if_sig = NP.roll(if_sig,N//2,axis=1)
    corrplot_sig.append(NP.abs(f_sig).sum())
    flux_sig.append(sig.sum())
    flux_if_sig.append(if_sig.sum())

flux_sig = NP.array(flux_sig)
flux_sig -= flux_sig.min()
flux_sig /= flux_sig.max()
corrplot_sig = NP.array(corrplot_sig)**1.1
corrplot_sig -= corrplot_sig.min()
corrplot_sig /= corrplot_sig.max()

PL.figure()
PL.plot(flux_sig,corrplot_sig,'.')
PL.ylabel('flux')
PL.xlabel('corrplot')
PL.grid()
