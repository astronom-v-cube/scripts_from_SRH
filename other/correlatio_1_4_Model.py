#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 09:38:45 2020

@author: sergeylesovoi
"""

import pylab as PL
import numpy as NP
import scipy.signal

N = 50000
L = 500
covar_re = NP.zeros(L)
covar_im = NP.zeros(L)
covar1bit_re = NP.zeros(L)
covar4bit_re = NP.zeros(L)
covar1bit_im = NP.zeros(L)
covar4bit_im = NP.zeros(L)
var_re = NP.zeros(L)
var1bit_re = NP.zeros(L)
var4bit_re = NP.zeros(L)

vanVleck = True

fig = PL.figure(figsize=(10,10))
if (vanVleck):
    fig.suptitle('van Vleck correction is applied')
else:
    fig.suptitle('van Vleck correction is NOT applied')

sub_pl = [PL.subplot(221), PL.subplot(222), PL.subplot(223), PL.subplot(224)]

alphas = NP.array([1., 2., 3., 4.])
#alphas = NP.array([.1, 1., 3., 30.])
#alphas = NP.linspace(0,50,300)
delays = [0, 4, 8, 12]

rhos = NP.zeros(alphas.shape[0])
rhos1bit = NP.zeros(alphas.shape[0])
rhos4bit = NP.zeros(alphas.shape[0])

lf_filter = NP.ones(13)
for j in range(len(alphas)):
    alpha = alphas[j]
    delay = delays[0]
    for i in range(L):
        inputSignal = NP.random.randn(N)
        noise1 = NP.random.randn(N)
        noise2 = NP.random.randn(N)
        
        input1 = alpha*inputSignal + noise1
        input2 = alpha*inputSignal + noise2
        if (delay < 0):
            input1 = NP.convolve(input1,lf_filter)[0:N]
            input2 = NP.convolve(input2,lf_filter)[abs(delay):N + abs(delay)]
        else:
            input1 = NP.convolve(input1,lf_filter)[delay:N + delay]
            input2 = NP.convolve(input2,lf_filter)[0:N]
        input2_a = scipy.signal.hilbert(input2)
        
        oneBitInput1 = NP.sign(input1)
        oneBitInput2_re = NP.sign(input2_a.real)
        oneBitInput2_im = NP.sign(input2_a.imag)
        
        fourScale = .9/input1.max()
        fourBitInput1 = NP.clip(NP.round(input1*fourScale),-7,7)
        fourBitInput2_re = NP.clip(NP.round(input2_a.real*fourScale),-7,7)
        fourBitInput2_im = NP.clip(NP.round(input2_a.imag*fourScale),-7,7)
#        fourBitInput1 = NP.clip(NP.round(input1*fourScale),-6,6)  + NP.sign(input1)
#        fourBitInput2_re = NP.clip(NP.round(input2_a.real*fourScale),-6,6)  + NP.sign(input2_a.real)
#        fourBitInput2_im = NP.clip(NP.round(input2_a.imag*fourScale),-6,6)  + NP.sign(input2_a.imag)
        
        covar_re[i] = NP.mean(input1*input2_a.real)/(NP.std(input1)*NP.std(input2_a.real))
        covar_im[i] = NP.mean(input1*input2_a.imag)/(NP.std(input1)*NP.std(input2_a.imag))

        if (vanVleck):
            covar1bit_re[i] = NP.sin(NP.pi/2*NP.mean(oneBitInput1*oneBitInput2_re))
            covar1bit_im[i] = NP.sin(NP.pi/2*NP.mean(oneBitInput1*oneBitInput2_im))
        else:
            covar1bit_re[i] = NP.mean(oneBitInput1*oneBitInput2_re)
            covar1bit_im[i] = NP.mean(oneBitInput1*oneBitInput2_im)

        var_re[i] = NP.mean(input1*input1)
        var1bit_re[i] = NP.mean(oneBitInput1*oneBitInput1)
        var4bit_re[i] = NP.mean(fourBitInput1*fourBitInput1)
        
        covar4bit_re[i] = NP.mean(fourBitInput1*fourBitInput2_re)/(NP.std(fourBitInput1)*NP.std(fourBitInput2_re))
        covar4bit_im[i] = NP.mean(fourBitInput1*fourBitInput2_im)/(NP.std(fourBitInput1)*NP.std(fourBitInput2_im))
    
    rho = NP.mean(NP.sqrt(covar_re**2 + covar_im**2))
    rhos[j] = rho
    rhos1bit[j] = NP.mean(NP.sqrt(covar1bit_re**2 + covar1bit_im**2))
    rhos4bit[j] = NP.mean(NP.sqrt(covar4bit_re**2 + covar4bit_im**2))
    
    phi = NP.rad2deg(NP.mean(NP.arctan2(covar_im, covar_re)))
    sub_pl[j].set_xlim(-.15,.15)
    sub_pl[j].set_ylim(-.15,.15)
    if (j > 1):
        sub_pl[j].set_xlabel(r'$\delta$Re')
    if (j % 2 == 0):
        sub_pl[j].set_ylabel(r'$\delta$Im')
    sub_pl[j].set_title(r'$\rho$=%1.1f, $\phi$=%3.0f, $\sigma_{re}$=%1.3f, $\sigma_{im}$=%1.3f' % (rho, phi, NP.std(covar4bit_re), NP.std(covar4bit_im)))
    sub_pl[j].plot(covar1bit_re - NP.mean(covar1bit_re), covar1bit_im - NP.mean(covar1bit_im), 'o', markersize=3, label='1 bit')
    sub_pl[j].plot(covar4bit_re - NP.mean(covar4bit_re), covar4bit_im - NP.mean(covar4bit_im), 'o', markersize=3, label='4 bit')
    sub_pl[j].plot(covar_re - NP.mean(covar_re), covar_im - NP.mean(covar_im), 'o', markersize=3, label='float')
    sub_pl[j].legend()
