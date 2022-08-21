#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 09:38:45 2020

@author: sergeylesovoi
"""

import pylab as PL
import numpy as NP
import scipy.signal

N = 500
L = 500
var_re = NP.zeros(L)
var_im = NP.zeros(L)
var4bit_re = NP.zeros(L)
var4bit_im = NP.zeros(L)

covar_re = NP.zeros(L)
covar_im = NP.zeros(L)
covar1bit_re = NP.zeros(L)
covar4bit_re = NP.zeros(L)
covar1bit_im = NP.zeros(L)
covar4bit_im = NP.zeros(L)

vanVleck = False

fig = PL.figure(figsize=(10,10))
if (vanVleck):
    fig.suptitle('van Vleck correction is applied')
else:
    fig.suptitle('van Vleck correction is NOT applied')

sub_pl = [PL.subplot(221), PL.subplot(222), PL.subplot(223), PL.subplot(224)]

alphas = NP.linspace(0,50,300)
rhos_re = NP.zeros(alphas.shape[0])
rhos1bit_re = NP.zeros(alphas.shape[0])
rhos4bit_re = NP.zeros(alphas.shape[0])
rhos_im = NP.zeros(alphas.shape[0])
rhos1bit_im = NP.zeros(alphas.shape[0])
rhos4bit_im = NP.zeros(alphas.shape[0])

for j in range(len(alphas)):
    alpha = 2.*(alphas[j] // 25) + 1.
    for i in range(L):
        inputSignal = NP.random.randn(N)
        noise1 = NP.random.randn(N)
        noise2 = NP.random.randn(N)
        
        input1 = alpha*inputSignal + noise1
        input2 = alpha*inputSignal + noise2
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
        
        var_re[i] = NP.mean(input2_a.real*input2_a.real)
        var_im[i] = NP.mean(input2_a.imag*input2_a.imag)
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
    
        var4bit_re[i] = NP.mean(fourBitInput2_re*fourBitInput2_re)
        var4bit_im[i] = NP.mean(fourBitInput2_im*fourBitInput2_im)

    rhos_re[j] = NP.mean(covar_re)
    rhos1bit_re[j] = NP.mean(covar1bit_re)
    rhos4bit_re[j] = NP.mean(covar4bit_re)
    rhos_im[j] = NP.mean(covar_im)
    rhos1bit_im[j] = NP.mean(covar1bit_im)
    rhos4bit_im[j] = NP.mean(covar4bit_im)
