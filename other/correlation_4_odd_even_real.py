#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 09:38:45 2020

@author: sergeylesovoi
"""

import pylab as PL
import numpy as NP

N = 1000
L = 500

stdFloat = NP.zeros(L)
covarFloat = NP.zeros(L)
rhoFloat = NP.zeros(L)
varAfloat = NP.zeros(L)
varBfloat = NP.zeros(L)

covar1bit = NP.zeros(L)

std4bit = NP.zeros(L)
covar4bit = NP.zeros(L)
rho4bit = NP.zeros(L)
varA4bit = NP.zeros(L)
varB4bit = NP.zeros(L)


vanVleck = True
oddLevels = True

if (oddLevels):
    f4 = NP.array([1.59304403, 0.68085433, 0.92238376])
else:
    f4 = NP.array([1.28451175, 0.89685013, 0.99309232])

fig, sub_pl = PL.subplots(nrows=1,ncols=1)

title = 'van Vleck correction is'
if (vanVleck): 
    title += ' applied'
else:
    title += ' is NOT applied'
    
if (oddLevels): 
    title += ', odd levels'
else:
    title += ', even levels'

fig.suptitle(title)

alphaMax = 10
alphaN = 100
alphas = NP.linspace(0,alphaMax,alphaN)
betaA = 1
betaB = 1

covar_float = NP.zeros(alphas.shape[0])
std_float = NP.zeros(alphas.shape[0])
rho_float = NP.zeros(alphas.shape[0])
varA_float = NP.zeros(alphas.shape[0])
varB_float = NP.zeros(alphas.shape[0])

covar_1bit = NP.zeros(alphas.shape[0])

covar_4bit = NP.zeros(alphas.shape[0])
std_4bit = NP.zeros(alphas.shape[0])
rho_4bit = NP.zeros(alphas.shape[0])
varA_4bit = NP.zeros(alphas.shape[0])
varB_4bit = NP.zeros(alphas.shape[0])

for j in range(len(alphas)):
    alpha = alphas[j]
    for i in range(L):
        inputSignal = NP.random.randn(N)
        noiseA = NP.random.randn(N)
        noiseB = NP.random.randn(N)
        
        inputA = alpha*inputSignal + betaA*noiseA
        inputB = alpha*inputSignal + betaB*noiseB
        
        oneBitInputA = NP.sign(inputA)
        oneBitInputB = NP.sign(inputB)
        
        fourScale = 6/alphaMax
        if (oddLevels):
            fourBitInputA = NP.clip(NP.round(inputA*fourScale),-7,7)
            fourBitInputB = NP.clip(NP.round(inputB*fourScale),-7,7)
        else:
            fourBitInputA = NP.clip(NP.round(inputA*fourScale),-6,6) + NP.sign(inputA)
            fourBitInputB = NP.clip(NP.round(inputB*fourScale),-6,6) + NP.sign(inputB)
        
        covarFloat[i] = NP.mean(inputA*inputB)
        stdFloat[i] = NP.std(inputA)*NP.std(inputB)
        varAfloat[i] = NP.mean(inputA**2)
        varBfloat[i] = NP.mean(inputB**2)
        rhoFloat[i] = covarFloat[i] / stdFloat[i]

        covar4bit[i] = NP.mean(fourBitInputA*fourBitInputB)
        std4bit[i] = NP.std(fourBitInputA)*NP.std(fourBitInputB)
        varA4bit[i] = NP.mean(fourBitInputA**2)
        varB4bit[i] = NP.mean(fourBitInputB**2)

        if (vanVleck):
            covar1bit[i] = NP.sin(NP.pi/2*NP.mean(oneBitInputA*oneBitInputB))
            rho4bit[i] = f4[0]*NP.sin(f4[1]*(covar4bit[i] / std4bit[i])**f4[2])
        else:
            covar1bit[i] = NP.mean(oneBitInputA*oneBitInputB)
    
    covar_float[j] = NP.mean(covarFloat)
    std_float[j] = NP.mean(stdFloat)
    varA_float[j] = NP.mean(varAfloat)
    varB_float[j] = NP.mean(varBfloat)
    rho_float[j] = NP.mean(rhoFloat)

    covar_1bit[j] = NP.mean(covar1bit)

    covar_4bit[j] = NP.mean(covar4bit)
    std_4bit[j] = NP.mean(std4bit)
    varA_4bit[j] = NP.mean(varA4bit)
    varB_4bit[j] = NP.mean(varB4bit)
    rho_4bit[j] = NP.mean(rho4bit)

sub_pl.set_xlabel('SNR')
sub_pl.set_ylabel('correlation')
sub_pl.plot(alphas,covar_float / std_float, label='float')
if (vanVleck):
    sub_pl.plot(alphas,rho_4bit, label='4 bit')
else:
    sub_pl.plot(alphas,covar_4bit / std_4bit, label='4 bit')
sub_pl.plot(alphas,covar_1bit, label='1 bit')
sub_pl.legend()
sub_pl.grid()
