#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 14 11:08:30 2021

@author: maria
"""

import numpy as NP
import pylab as PL
import random
from BadaryRAO import BadaryRAO
from scipy.optimize import least_squares
import base2uvw_36
from astropy.time import Time, TimeDelta

def phaMatrixGenPairsEWN_constrained(pairs, antNumberEW, antNumberN):
    rowsEW = int(((antNumberEW - 1) + (antNumberEW - pairs))/2 * pairs)
    rowsN = int(((antNumberN) + (antNumberN + 1 - pairs))/2 * pairs)
    colsEW = antNumberEW + pairs
    colsN = antNumberN + pairs
    phaMatrix = NP.zeros((rowsEW + rowsN + 2, colsEW + colsN))
    for pair in range(pairs):
        row0 = int(((antNumberEW - 1) + (antNumberEW - pair))/2 * pair)
        row1 = row0 + (antNumberEW - pair - 1)
        phaMatrix[row0:row1,pair] = 1
        for phaPair in range(antNumberEW - pair - 1):
            phaMatrix[phaPair + row0, phaPair + 2*pairs] = 1
            phaMatrix[phaPair + row0, phaPair + 2*pairs + (pair + 1)] = -1
        row0 = int(((antNumberN) + (antNumberN + 1 - pair))/2 * pair)
        row1 = row0 + (antNumberN - pair)
        phaMatrix[row0 + rowsEW:row1 + rowsEW,pairs + pair] = 1
        for phaPair in range(antNumberN - pair):
            if phaPair == 0:
                phaMatrix[rowsEW + row0, 2*pairs + antNumberEW//2] = 1
            else:
                phaMatrix[phaPair + rowsEW + row0, 2*pairs + phaPair + antNumberEW - 1] = 1
            phaMatrix[phaPair + rowsEW + row0, 2*pairs + phaPair + antNumberEW - 1 + (pair + 1)] = -1
    phaMatrix[rowsEW + rowsN, 0] = 1
    phaMatrix[rowsEW + rowsN+1, pairs] = 1
    return phaMatrix.copy()

def phaMatrixGenPairsEWN(pairs, antNumberEW, antNumberN):
    rowsEW = int(((antNumberEW - 1) + (antNumberEW - pairs))/2 * pairs)
    rowsN = int(((antNumberN) + (antNumberN + 1 - pairs))/2 * pairs)
    colsEW = antNumberEW + pairs
    colsN = antNumberN + pairs
    phaMatrix = NP.zeros((rowsEW + rowsN, colsEW + colsN))
    for pair in range(pairs):
        row0 = int(((antNumberEW - 1) + (antNumberEW - pair))/2 * pair)
        row1 = row0 + (antNumberEW - pair - 1)
        phaMatrix[row0:row1,pair] = 1
        for phaPair in range(antNumberEW - pair - 1):
            phaMatrix[phaPair + row0, phaPair + 2*pairs] = 1
            phaMatrix[phaPair + row0, phaPair + 2*pairs + (pair + 1)] = -1
        row0 = int(((antNumberN) + (antNumberN + 1 - pair))/2 * pair)
        row1 = row0 + (antNumberN - pair)
        phaMatrix[row0 + rowsEW:row1 + rowsEW,pairs + pair] = 1
        for phaPair in range(antNumberN - pair):
            if phaPair == 0:
                phaMatrix[rowsEW + row0, 2*pairs + 32] = 1
            else:
                phaMatrix[phaPair + rowsEW + row0, 2*pairs + phaPair + antNumberEW - 1] = 1
            phaMatrix[phaPair + rowsEW + row0, 2*pairs + phaPair + antNumberEW - 1 + (pair + 1)] = -1
    return phaMatrix.copy()

def real_to_complex(z):
    return z[:len(z)//2] + 1j * z[len(z)//2:]
    
def complex_to_real(z):
    return NP.concatenate((NP.real(z), NP.imag(z)))

def allGainsFunc_constrained(x, obsVis, ewAntNumber, nAntNumber, baselineNumber):
        res = NP.zeros_like(obsVis, dtype = complex)
        ewSolarAmp = x[0]
        nSolarAmp = x[1]
        x_complex = real_to_complex(x[2:])
        
        nAntNumber_c = nAntNumber + 1
        
        ewSolarPhase = 0
        nSolarPhase = 0
        
        nGainsNumber = nAntNumber
        ewGainsNumber = ewAntNumber
        nSolVisNumber = baselineNumber - 1
        ewSolVisNumber = baselineNumber - 1
        ewSolVis = NP.append((ewSolarAmp*NP.exp(1j*ewSolarPhase)), x_complex[: ewSolVisNumber])
        nSolVis = NP.append((nSolarAmp*NP.exp(1j*nSolarPhase)), x_complex[ewSolVisNumber : ewSolVisNumber+nSolVisNumber])
        ewGains = x_complex[ewSolVisNumber+nSolVisNumber : ewSolVisNumber+nSolVisNumber+ewGainsNumber]
        nGains = NP.append(ewGains[32], x_complex[ewSolVisNumber+nSolVisNumber+ewGainsNumber :])
        
        solVisArrayN = NP.array(())
        antAGainsN = NP.array(())
        antBGainsN = NP.array(())
        solVisArrayEW = NP.array(())
        antAGainsEW = NP.array(())
        antBGainsEW = NP.array(())
        for baseline in range(1, baselineNumber+1):
            solVisArrayN = NP.append(solVisArrayN, NP.full(nAntNumber_c-baseline, nSolVis[baseline-1]))
            antAGainsN = NP.append(antAGainsN, nGains[:nAntNumber_c-baseline])
            antBGainsN = NP.append(antBGainsN, nGains[baseline:])
            
            solVisArrayEW = NP.append(solVisArrayEW, NP.full(ewAntNumber-baseline, ewSolVis[baseline-1]))
            antAGainsEW = NP.append(antAGainsEW, ewGains[:ewAntNumber-baseline])
            antBGainsEW = NP.append(antBGainsEW, ewGains[baseline:])
            
        res = NP.append(solVisArrayEW, solVisArrayN) * NP.append(antAGainsEW, antAGainsN) * NP.conj(NP.append(antBGainsEW, antBGainsN)) - obsVis
        return complex_to_real(res)  


baselinesNumber = 5
antNumberEW = 97
antNumberN = 31

N = 2048
arcsecPerPix = 2
radPerPix = NP.deg2rad(arcsecPerPix/3600.)
arcsecRadius = 1020
frequency = 3e9
degRadius = NP.deg2rad(arcsecRadius/3600)
radius = int(arcsecRadius/arcsecPerPix +0.5)
model = NP.zeros((N, N))
for i in range(N):
    for j in range(N):
        x=i - N/2
        y=j - N/2
        if (NP.sqrt(x**2 + y**2) < radius):
            model[i, j] = 1.
#model[1024,1024]= 10000
RAO = BadaryRAO('2021-05-14')
gains = NP.ones(128, dtype = 'complex')
for i in range(128):
    gains[i] = random.uniform(1., 2) * NP.exp(1j * random.uniform(-NP.pi/2., NP.pi/2.))
            
ewGains = gains[:97]
nGains = gains[97:]

declination = RAO.declination
noon = RAO.culmination

#uvPlane = NP.zeros((128,128),dtype=complex);
uvPlane2 = NP.zeros((1024,1024),dtype=complex);
#uvPlaneCorrected = NP.zeros((1024,1024),dtype=complex);
uvPlanePSF = NP.zeros((1024,1024),dtype=complex);
   
fitsDate ='2021-05-14T00:00:00';
scan = 0

scanDate = Time(fitsDate, format='isot',scale='utc');
#scanTime = noon
scanTime = 25200.
#scanTime = 28800.
scanDate += TimeDelta(scanTime,format='sec')
hourAngle = NP.deg2rad((scanTime - noon)*15./3600.)
O = 1024//2
FOV = N * radPerPix
x,y = NP.meshgrid(NP.linspace(-.5,.5,N), NP.linspace(-.5,.5,N))
ewSolVis = NP.zeros(baselinesNumber, dtype = 'complex')
nSolVis = NP.zeros(baselinesNumber, dtype = 'complex')
for i in range(baselinesNumber):
    baseline = i+1
    uvw = base2uvw_36.base2uvw(hourAngle,declination, 1, 1+baseline)
    cos_uv = NP.cos(2. * NP.pi * ((uvw[0]*frequency/3e8)*x + (uvw[1]*frequency/3e8)*y) * FOV)
    sin_uv = NP.sin(2. * NP.pi * ((uvw[0]*frequency/3e8)*x + (uvw[1]*frequency/3e8)*y) * FOV)
    real = NP.sum(cos_uv * model)
    imag = NP.sum(sin_uv * model)
    ewSolVis[i] = (real + imag * 1j)/1e6
    
    uvw = base2uvw_36.base2uvw(hourAngle,declination, 98, 98+baseline)
    cos_uv = NP.cos(2. * NP.pi * ((uvw[0]*frequency/3e8)*x + (uvw[1]*frequency/3e8)*y) * FOV)
    sin_uv = NP.sin(2. * NP.pi * ((uvw[0]*frequency/3e8)*x + (uvw[1]*frequency/3e8)*y) * FOV)
    real = NP.sum(cos_uv * model)
    imag = NP.sum(sin_uv * model)
    nSolVis[i] = (real + imag * 1j)/1e6

nGains_C0 = NP.append(ewGains[32], nGains)

ewRedundantVis = NP.array(())
nRedundantVis = NP.array(())

for i in range(baselinesNumber):
    baseline = i+1
    ewRedundantVis = NP.append(ewRedundantVis, ewSolVis[i] * ewGains[:97-baseline] * NP.conj(ewGains[baseline:]))
    nRedundantVis = NP.append(nRedundantVis, nSolVis[i] * nGains_C0[:32-baseline] * NP.conj(nGains_C0[baseline:]))

redundantVisAll = NP.append(ewRedundantVis, nRedundantVis)

# NONLINEAR

x_size = (baselinesNumber-1)*2 + antNumberEW + antNumberN
x_ini = NP.concatenate((NP.ones(x_size+2), NP.zeros(x_size)))
ls_res = least_squares(allGainsFunc_constrained, x_ini, args = (redundantVisAll, antNumberEW, antNumberN, baselinesNumber), max_nfev = 1500)

gains_calc = real_to_complex(ls_res['x'][2:])[(baselinesNumber-1)*2:]
ew_gains_lcp = gains_calc[:antNumberEW]
ewAntPhaLcp= NP.angle(ew_gains_lcp)
n_gains_lcp = gains_calc[antNumberEW:]
nAntPhaLcp = NP.angle(n_gains_lcp)

ewAntAmpLcp = NP.abs(ew_gains_lcp)#/NP.min(NP.abs(ew_gains_lcp))
# ewAntAmpLcp[freqChannel][ewAntAmpLcp[freqChannel]<NP.median(ewAntAmpLcp[freqChannel])*0.6] = 1e6
nAntAmpLcp = NP.abs(n_gains_lcp)#/NP.min(NP.abs(n_gains_lcp))
# nAntAmpLcp[freqChannel][snAntAmpLcp[freqChannel]<NP.median(nAntAmpLcp[freqChannel])*0.6] = 1e6

# LINEAR

allAmp = NP.abs(redundantVisAll)
ampMatrix = NP.abs(phaMatrixGenPairsEWN(baselinesNumber, antNumberEW, antNumberN))

antAmp, c, d, e = NP.linalg.lstsq(ampMatrix,NP.log(allAmp), rcond=None)
antAmp= NP.exp(antAmp[baselinesNumber*2:])
ewAmp = antAmp[:antNumberEW]
nAmp = antAmp[antNumberEW:]
