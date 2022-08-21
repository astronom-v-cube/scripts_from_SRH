#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 09:38:45 2020

@author: sergeylesovoi
"""

import pylab as PL
import numpy as NP
import scipy.signal

N = 1000
L = 500
covar_re = NP.zeros(L)
covar_im = NP.zeros(L)
covar1bit_re = NP.zeros(L)
covar4bit_re = NP.zeros(L)
covar1bit_im = NP.zeros(L)
covar4bit_im = NP.zeros(L)
covarAlpha = NP.zeros(L)

K = 50
deltaK = L//K
inputLevel = NP.zeros(K)
snr = NP.zeros(K)
snr1bit = NP.zeros(K)
snr4bit = NP.zeros(K)

vanVleck = False

fig = PL.figure(figsize=(10,10))
if (vanVleck):
    fig.suptitle('van Vleck correction is applied')
else:
    fig.suptitle('van Vleck correction is NOT applied')

sub_pl = [PL.subplot(221), PL.subplot(222), PL.subplot(223), PL.subplot(224)]

alphas = [.5, .8, 1.5, 30.]
lf_filter = NP.ones(13)
delay = -12
for j in range(4):
    for i in range(L):
        alpha = alphas[j]
        covarAlpha[i] = alpha
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
        
        fourBitInput1 = NP.clip(NP.round(input1*3.5),-7,7)
        fourBitInput2_re = NP.clip(NP.round(input2_a.real*3.5),-7,7)
        fourBitInput2_im = NP.clip(NP.round(input2_a.imag*3.5),-7,7)
        
        covar_re[i] = NP.mean(input1*input2_a.real)/(NP.std(input1)*NP.std(input2_a.real))
        covar_im[i] = NP.mean(input1*input2_a.imag)/(NP.std(input1)*NP.std(input2_a.imag))

        if (vanVleck):
            covar1bit_re[i] = NP.sin(NP.pi/2*NP.mean(oneBitInput1*oneBitInput2_re))
            covar1bit_im[i] = NP.sin(NP.pi/2*NP.mean(oneBitInput1*oneBitInput2_im))
        else:
            covar1bit_re[i] = NP.mean(oneBitInput1*oneBitInput2_re)
            covar1bit_im[i] = NP.mean(oneBitInput1*oneBitInput2_im)

        covar4bit_re[i] = NP.mean(fourBitInput1*fourBitInput2_re)/(NP.std(fourBitInput1)*NP.std(fourBitInput2_re))
        covar4bit_im[i] = NP.mean(fourBitInput1*fourBitInput2_im)/(NP.std(fourBitInput1)*NP.std(fourBitInput2_im))
    
    for i in range(K):
        snr[i] = NP.mean(covar_re[i*deltaK:(i+1)*deltaK])/NP.std(covar_re[i*deltaK:(i+1)*deltaK])
        inputLevel[i] = covarAlpha[i*deltaK]
        snr1bit[i] = NP.mean(covar1bit_re[i*deltaK:(i+1)*deltaK])/NP.std(covar1bit_re[i*deltaK:(i+1)*deltaK])
        snr4bit[i] = NP.mean(covar4bit_re[i*deltaK:(i+1)*deltaK])/NP.std(covar4bit_re[i*deltaK:(i+1)*deltaK])

#fig = PL.figure(figsize=(20,16))
#pl0 = PL.subplot(211)
#pl1 = PL.subplot(212)
#
#pl1.plot(inputLevel, snr1bit*5, label='1 bit')
#pl1.plot(inputLevel, snr4bit*5, label='4 bit')
#
#micranSnrInput = 10**(NP.array([0, 5 , 15, 25, 35])/10)
#micranSnr1bit = [143.3, 606, 2252, 4839, 5470]
#micranSnr4bit = [219, 645, 883, 1146, 1128]
#
#pl0.plot(micranSnrInput, micranSnr1bit, label='1 bit')
#pl0.plot(micranSnrInput, micranSnr4bit, label='4 bit')

#pl0.set_xlabel('input level')
#pl0.set_ylabel('micran snr')
#pl1.set_xlabel('input level')
#pl1.set_ylabel('simulated snr')
#pl0.legend()
#pl1.legend()

    rho = NP.mean(NP.sqrt(covar_re**2 + covar_im**2))
    phi = NP.rad2deg(NP.mean(NP.arctan2(covar_im, covar_re)))
    sub_pl[j].set_xlim(-.3,.3)
    sub_pl[j].set_ylim(-.3,.3)
    if (j > 1):
        sub_pl[j].set_xlabel(r'$\delta$Re')
    if (j % 2 == 0):
        sub_pl[j].set_ylabel(r'$\delta$Im')
    sub_pl[j].set_title(r'$\rho$=%1.1f, $\phi$=%3.0f, $\sigma_{re}$=%1.2f, $\sigma_{im}$=%1.2f' % (rho, phi, NP.std(covar4bit_re), NP.std(covar4bit_im)))
    sub_pl[j].plot(covar1bit_re - NP.mean(covar1bit_re), covar1bit_im - NP.mean(covar1bit_im), 'o', markersize=3, label='1 bit')
    sub_pl[j].plot(covar4bit_re - NP.mean(covar4bit_re), covar4bit_im - NP.mean(covar4bit_im), 'o', markersize=3, label='4 bit')
    sub_pl[j].plot(covar_re - NP.mean(covar_re), covar_im - NP.mean(covar_im), 'o', markersize=3, label='float')
    sub_pl[j].legend()
