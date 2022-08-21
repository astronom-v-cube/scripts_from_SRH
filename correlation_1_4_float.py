#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 09:38:45 2020

@author: sergeylesovoi
"""

import pylab as PL
import numpy as NP
import scipy.optimize as opt

def fitFunc(x, A, B, C, D, E):
#    return A + B*x + C*x**2 + D*x**3
#    return A*x + B*x**2 + C*x**3 + D*x**4 + E*x**5
    return A*x + B*x**3 + C*x**5 + D*x**7 + E*x**9

def fitSin(x, A, B, C, D):
    return A*NP.sin(B*x + C*x**2 + D*x**3)

def fitPow(x, A, B, C):
    return A*NP.sin(B*x**C)

sugg = [.6, 0., 0.1, 0., 0.05]
sin_gues = [1.,1.5,0,0]
pow_gues = [1.2,1.,1.]

N = 3000
L = 50
corr = NP.zeros(L)
corr1bit = NP.zeros(L)
corr2bit = NP.zeros(L)
corr3bit = NP.zeros(L)
corr4bit = NP.zeros(L)
corr5bit = NP.zeros(L)
corr6bit = NP.zeros(L)
corr7bit = NP.zeros(L)

vanVleck = True

fig = PL.figure(figsize=(10,10))
if (vanVleck):
    fig.suptitle('van Vleck correction is applied')
else:
    fig.suptitle('van Vleck correction is NOT applied')

sub_pl = [PL.subplot(111)]

alphaMax = 40
alphas = NP.linspace(0,alphaMax,100)

rhos = NP.zeros(alphas.shape[0])
rhos1bit = NP.zeros(alphas.shape[0])
rhos2bit = NP.zeros(alphas.shape[0])
rhos3bit = NP.zeros(alphas.shape[0])
rhos4bit = NP.zeros(alphas.shape[0])
rhos5bit = NP.zeros(alphas.shape[0])
rhos6bit = NP.zeros(alphas.shape[0])
rhos7bit = NP.zeros(alphas.shape[0])

lf_filter = NP.ones(13)
for j in range(len(alphas)):
    alpha = alphas[j]
    for i in range(L):
        inputSignal = NP.random.randn(N)
        noise1 = NP.random.randn(N)
        noise2 = NP.random.randn(N)
        
        signalAndNoise1 = alpha*inputSignal + noise1
        signalAndNoise2 = alpha*inputSignal + noise2
        input1 = NP.convolve(signalAndNoise1,lf_filter)[0:N]
        input2 = NP.convolve(signalAndNoise2,lf_filter)[0:N]
        
        oneBitInput1 = NP.sign(input1)
        oneBitInput2 = NP.sign(input2)
        
        twoBitInput1 = NP.clip(NP.round(input1),-1,1)
        twoBitInput2 = NP.clip(NP.round(input2),-1,1)

        threeBitInput1 = NP.clip(NP.round(2/alphaMax*input1),-3,3)
        threeBitInput2 = NP.clip(NP.round(2/alphaMax*input2),-3,3)

#        fourBitInput1 = NP.clip(NP.round(6/alphaMax*input1),-7,7)
#        fourBitInput2 = NP.clip(NP.round(6/alphaMax*input2),-7,7)
        fourBitInput1 = NP.clip(NP.round(6/alphaMax*input1),-6,6) + NP.sign(input1)
        fourBitInput2 = NP.clip(NP.round(6/alphaMax*input2),-6,6) + NP.sign(input2)
        
        fiveBitInput1 = NP.clip(NP.round(14/alphaMax*input1),-15,15)
        fiveBitInput2 = NP.clip(NP.round(14/alphaMax*input2),-15,15)

        sixBitInput1 = NP.clip(NP.round(input1),-31,31)
        sixBitInput2 = NP.clip(NP.round(input2),-31,31)
        
        sevenBitInput1 = NP.clip(NP.round(input1),-63,63)
        sevenBitInput2 = NP.clip(NP.round(input2),-63,63)
        
        corr[i] = NP.mean(input1*input2)/(NP.std(input1)*NP.std(input2))

        if (vanVleck):
            corr1bit[i] = NP.sin(NP.pi/2*NP.mean(oneBitInput1*oneBitInput2))
        else:
            corr1bit[i] = NP.mean(oneBitInput1*oneBitInput2)
            corr2bit[i] = NP.mean(twoBitInput1*twoBitInput2)/(NP.std(twoBitInput1)*NP.std(twoBitInput2))
            corr3bit[i] = NP.mean(threeBitInput1*threeBitInput2)/(NP.std(threeBitInput1)*NP.std(threeBitInput2))
            corr4bit[i] = NP.mean(fourBitInput1*fourBitInput2)/(NP.std(fourBitInput1)*NP.std(fourBitInput2))
            corr5bit[i] = NP.mean(fiveBitInput1*fiveBitInput2)/(NP.std(fiveBitInput1)*NP.std(fiveBitInput2))
            corr6bit[i] = NP.mean(sixBitInput1*sixBitInput2)/(NP.std(sixBitInput1)*NP.std(sixBitInput2))
            corr7bit[i] = NP.mean(sevenBitInput1*sevenBitInput2)/(NP.std(sevenBitInput1)*NP.std(sevenBitInput2))

    rhos[j] = NP.mean(corr)
    rhos1bit[j] = NP.mean(corr1bit)
    rhos2bit[j] = NP.mean(corr2bit)
    rhos3bit[j] = NP.mean(corr3bit)
    rhos4bit[j] = NP.mean(corr4bit)
    rhos5bit[j] = NP.mean(corr5bit)
    rhos6bit[j] = NP.mean(corr6bit)
    rhos7bit[j] = NP.mean(corr7bit)

fitVector = NP.zeros((7,5))
fit_1_1, _ = opt.curve_fit(fitFunc, rhos1bit, rhos1bit, p0=sugg)
fitVector[0] = fit_1_1
fit_2_1, _ = opt.curve_fit(fitFunc, rhos2bit, rhos1bit, p0=sugg)
fitVector[1] = fit_2_1
fit_3_1, _ = opt.curve_fit(fitFunc, rhos3bit, rhos1bit, p0=sugg)
fitVector[2] = fit_3_1
fit_4_1, _ = opt.curve_fit(fitFunc, rhos4bit, rhos1bit, p0=sugg)
fitVector[3] = fit_4_1
fit_5_1, _ = opt.curve_fit(fitFunc, rhos5bit, rhos1bit, p0=sugg)
fitVector[4] = fit_5_1
#fit_6_1, _ = opt.curve_fit(fitFunc, rhos6bit, rhos1bit, p0=sugg)
#fitVector[5] = fit_6_1
#fit_7_1, _ = opt.curve_fit(fitFunc, rhos7bit, rhos1bit, p0=sugg)
#fitVector[6] = fit_7_1

sub_pl[0].plot(rhos, rhos, 'o', label='float')
sub_pl[0].plot(rhos, rhos1bit, 'o')

#rhosFrom1bit = NP.sin(NP.pi/2*(fitFunc(rhos1bit, fit_1_1[0], fit_1_1[1], fit_1_1[2], fit_1_1[3], fit_1_1[4])))
#sub_pl[0].plot(rhos, rhosFrom1bit,'.', label='from 1 bit')
#rhosFrom2bit = NP.sin(NP.pi/2*(fitFunc(rhos2bit, fit_2_1[0], fit_2_1[1], fit_2_1[2], fit_2_1[3], fit_2_1[4])))
#sub_pl[0].plot(rhos, rhosFrom2bit,'.', label='from 2 bit')
#rhosFrom3bit = NP.sin(NP.pi/2*(fitFunc(rhos3bit, fit_3_1[0], fit_3_1[1], fit_3_1[2], fit_3_1[3], fit_3_1[4])))
#sub_pl[0].plot(rhos, rhosFrom3bit,'.', label='from 3 bit')
#rhosFrom4bit = NP.sin(NP.pi/2*(fitFunc(rhos4bit, fit_4_1[0], fit_4_1[1], fit_4_1[2], fit_4_1[3], fit_4_1[4])))
#sub_pl[0].plot(rhos, rhosFrom4bit,'.', label='from 4 bit')
#rhosFrom5bit = NP.sin(NP.pi/2*(fitFunc(rhos5bit, fit_5_1[0], fit_5_1[1], fit_5_1[2], fit_5_1[3], fit_5_1[4])))
#sub_pl[0].plot(rhos, rhosFrom5bit,'.', label='from 5 bit')
#rhosFrom6bit = NP.sin(NP.pi/2*(fitFunc(rhos6bit, fit_6_1[0], fit_6_1[1], fit_6_1[2], fit_6_1[3], fit_6_1[4])))
#sub_pl[0].plot(rhos, rhosFrom6bit,'.', label='from 6 bit')
#rhosFrom7bit = NP.sin(NP.pi/2*(fitFunc(rhos7bit, fit_7_1[0], fit_7_1[1], fit_7_1[2], fit_7_1[3], fit_7_1[4])))
#sub_pl[0].plot(rhos, rhosFrom7bit,'.', label='from 7 bit')
#

rhosFrom1bit = fitFunc(rhos1bit, fit_1_1[0], fit_1_1[1], fit_1_1[2], fit_1_1[3], fit_1_1[4])
sub_pl[0].plot(rhos, rhosFrom1bit,'.', label='from 1 bit')
rhosFrom2bit = fitFunc(rhos2bit, fit_2_1[0], fit_2_1[1], fit_2_1[2], fit_2_1[3], fit_2_1[4])
sub_pl[0].plot(rhos, rhosFrom2bit,'.', label='from 2 bit')
rhosFrom3bit = fitFunc(rhos3bit, fit_3_1[0], fit_3_1[1], fit_3_1[2], fit_3_1[3], fit_3_1[4])
sub_pl[0].plot(rhos, rhosFrom3bit,'.', label='from 3 bit')
rhosFrom4bit = fitFunc(rhos4bit, fit_4_1[0], fit_4_1[1], fit_4_1[2], fit_4_1[3], fit_4_1[4])
sub_pl[0].plot(rhos, rhosFrom4bit,'.', label='from 4 bit')
#rhosFrom5bit = fitFunc(rhos5bit, fit_5_1[0], fit_5_1[1], fit_5_1[2], fit_5_1[3], fit_5_1[4])
#sub_pl[0].plot(rhos, rhosFrom5bit,'.', label='from 5 bit')
#rhosFrom6bit = fitFunc(rhos6bit, fit_6_1[0], fit_6_1[1], fit_6_1[2], fit_6_1[3], fit_6_1[4])
#sub_pl[0].plot(rhos, rhosFrom6bit,'.', label='from 6 bit')
#rhosFrom7bit = fitFunc(rhos7bit, fit_7_1[0], fit_7_1[1], fit_7_1[2], fit_7_1[3], fit_7_1[4])
#sub_pl[0].plot(rhos, rhosFrom7bit,'.', label='from 7 bit')

sub_pl[0].set_xlabel(r'$\rho$')
sub_pl[0].set_ylabel(r'$\hat\rho$')
sub_pl[0].legend()
sub_pl[0].grid()
sub_pl[0].set_xlim(0,1)
sub_pl[0].set_ylim(0,1)

PL.figure()
PL.plot(fitVector[:,0], label='A')
PL.plot(fitVector[:,1], label='B')
PL.plot(fitVector[:,2], label='C')
PL.plot(fitVector[:,3], label='D')
PL.plot(fitVector[:,4], label='E')
PL.legend()
PL.grid()

#fsin1, _ = opt.curve_fit(fitSin, rhos1bit, rhos, p0 = sin_gues)
#fsin2, _ = opt.curve_fit(fitSin, rhos2bit, rhos, p0 = sin_gues)
#fsin3, _ = opt.curve_fit(fitSin, rhos3bit, rhos, p0 = sin_gues)
#fsin4, _ = opt.curve_fit(fitSin, rhos4bit, rhos, p0 = sin_gues)
#fsin5, _ = opt.curve_fit(fitSin, rhos5bit, rhos, p0 = sin_gues)

fsin1, _ = opt.curve_fit(fitPow, rhos1bit, rhos, p0 = pow_gues)
fsin2, _ = opt.curve_fit(fitPow, rhos2bit, rhos, p0 = pow_gues)
fsin3, _ = opt.curve_fit(fitPow, rhos3bit, rhos, p0 = pow_gues)
fsin4, _ = opt.curve_fit(fitPow, rhos4bit, rhos, p0 = pow_gues)
#fsin5, _ = opt.curve_fit(fitPow, rhos5bit, rhos, p0 = pow_gues)
#fsin6, _ = opt.curve_fit(fitPow, rhos6bit, rhos, p0 = pow_gues)
#fsin7, _ = opt.curve_fit(fitPow, rhos7bit, rhos, p0 = pow_gues)

#fsin_arr = NP.array([fsin1, fsin2, fsin3, fsin4, fsin5, fsin6])
fsin_arr = NP.array([fsin1, fsin2, fsin3, fsin4, fsin5])

PL.figure()
PL.ylim(0,2)
PL.plot(fsin_arr)
PL.grid()


#fsin_arr
#Out[351]: 
#array([[1.00664176, 1.56972987, 1.00124041],
#       [1.00840309, 1.53267617, 1.04988183],
#       [1.01907099, 1.43607287, 1.17617351],
#       [1.09181363, 1.18348301, 1.22178945],
#       [1.89338652, 0.56129703, 1.05600701]])
    