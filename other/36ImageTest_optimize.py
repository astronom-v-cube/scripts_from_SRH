#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 15:45:09 2021

@author: mariagloba
"""

from astropy.io import fits
import numpy as NP
import pylab as PL
import struct
from astropy import coordinates
from astropy import constants
from BadaryRAO import BadaryRAO
from sunpy import coordinates as sunpy_coordinates
from skimage.transform import warp, AffineTransform
import matrices_conj as matrices
import phaMatrixGen
from scipy.optimize import minimize, least_squares
import os, fnmatch;
import datetime as DT

def getPQScale(size, FOV, time, freq):
    cosP = NP.sin(hAngle) * NP.cos(RAO.declination)
    cosQ = NP.cos(hAngle) * NP.cos(RAO.declination) * NP.sin(RAO.observatory.lat) - NP.sin(RAO.declination) * NP.cos(RAO.observatory.lat)
    FOV_p = 2.*(constants.c / (freq*1e3)) / (RAO.base*NP.sqrt(1. - cosP**2.));
    FOV_q = 2.*(constants.c / (freq*1e3)) / (RAO.base*NP.sqrt(1. - cosQ**2.));
    
    return [int(size*FOV/FOV_p.to_value()), int(size*FOV/FOV_q.to_value())]
    
def getPQ2HDMatrix():
    gP =  NP.arctan(NP.tan(hAngle)*NP.sin(RAO.declination));
    gQ =  NP.arctan(-(NP.sin(RAO.declination) / NP.tan(hAngle) + NP.cos(RAO.declination) / (NP.sin(hAngle)*NP.tan(RAO.observatory.lat))));
    
    if hAngle > 0:
        gQ = NP.pi + gQ;
    g = gP - gQ;4
      
    pqMatrix = NP.zeros((3,3))
    pqMatrix[0, 0] =  NP.cos(gP) - NP.cos(g)*NP.cos(gQ)
    pqMatrix[0, 1] = -NP.cos(g)*NP.cos(gP) + NP.cos(gQ)
    pqMatrix[1, 0] =  NP.sin(gP) - NP.cos(g)*NP.sin(gQ)
    pqMatrix[1, 1] = -NP.cos(g)*NP.sin(gP) + NP.sin(gQ)
    pqMatrix /= NP.sin(g)**2.
    pqMatrix[2, 2] = 1.
    return pqMatrix

def real_to_complex(z):
    return z[:len(z)//2] + 1j * z[len(z)//2:]

def complex_to_real(z):
    return NP.concatenate((NP.real(z), NP.imag(z)))

# phase of first solVis = 0, phase of C0 = 0
def northGainsFunc_constrained(x, obsVis, antNumber, baselineNumber):
    res = NP.zeros_like(obsVis, dtype = complex)
    x_complex = real_to_complex(x)
    solVis = NP.append((1+0j), x_complex[:baselineNumber-1])
    gains = NP.append((1+0j), x_complex[baselineNumber-1:])
    i = 0
    for baseline in range(1, baselineNumber+1):
        for antA in range(antNumber-baseline):
            antB = antA + baseline
            res[i] = solVis[baseline-1] * gains[antA] * NP.conj(gains[antB]) - obsVis[i]
            i+=1
    return complex_to_real(res)

def eastWestGainsFunc_constrained(x, obsVis, antNumber, baselineNumber):
    res = NP.zeros_like(obsVis, dtype = complex)
    x_complex = real_to_complex(x)
    solVis = NP.append((1+0j), x_complex[:baselineNumber-1])
    gains = NP.insert(x_complex[baselineNumber-1:], antNumber//2, (1+0j))
    i = 0
    for baseline in range(1, baselineNumber+1):
        for antA in range(antNumber-baseline):
            antB = antA + baseline
            res[i] = solVis[baseline-1] * gains[antA] * NP.conj(gains[antB]) - obsVis[i]
            i+=1
    return complex_to_real(res)

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;

dateName = DT.datetime.now().strftime("%Y%m%d");
#fitNames = findFits('/home/svlesovoi/SRH_DATA/SRH36/' +  dateName, '*.fit')
#fitNames = findFits('/home/svlesovoi/Documents/Python Scripts/SRH36_temp_20210320_1','*.fit')
#fitNames = findFits('/home/sergeyvlesovoi/SRH36/20210324','*.fit')
#fitNames = findFits('/home/sergeyvlesovoi/SRH36/20210329','*.fit')
#fitNames = findFits('/home/sergeyvlesovoi/SRH36/20210323','*.fit')
fitNames = findFits('/home/sergeyvlesovoi/SRH36/20210404','*.fit')

fitNames.sort()
fitsFile = fits.open(fitNames[-1])

visScale = 1/(2e6*49)
ampScale = visScale / 128

antZeroRow = fitsFile[3].data['ant_zero_row'][:96]
time = fitsFile[1].data['time']
freqList = fitsFile[1].data['frequency']
freqListLength = fitsFile[1].data['frequency'].shape[0]
samplesNumber = time.shape[1]
visibilityNumber = fitsFile[1].data['vis_lcp'].shape[1]//samplesNumber

visLcp = fitsFile[1].data['vis_lcp'].reshape(freqListLength,samplesNumber,visibilityNumber) * visScale
visRcp = fitsFile[1].data['vis_rcp'].reshape(freqListLength,samplesNumber,visibilityNumber) * visScale

antA = fitsFile[4].data['ant_A']
antB = fitsFile[4].data['ant_B']
#antA = fitsFile[4].data['ant_B']
#antB = fitsFile[4].data['ant_A']
antNumberN = 31
antNumberEW = 65
baselinesNumber = 9
scan = 3
freq = 1

flagsN = [17, 22]
flagsEW = [16, 17, 21, 22, 35, 42, 45]

redIndexesN = []
for baseline in range(1, baselinesNumber+1):
    redIndexesN.append(NP.where((antA==98-1+baseline) & (antB==33))[0][0])
    for i in range(antNumberN - baseline):
#        redIndexesN.append(NP.where((antA==98+i) & (antB==98+i+baseline))[0][0])
        redIndexesN.append(NP.where((antB==98+i) & (antA==98+i+baseline))[0][0])

redIndexesEW = []
for baseline in range(1, baselinesNumber+1):
    for i in range(antNumberEW - baseline):
        redIndexesEW.append(NP.where((antA==1+i) & (antB==1+i+baseline))[0][0])
        
        

phaMatrix = phaMatrixGen.phaMatrixGenPairsEWN(baselinesNumber, antNumberEW, antNumberN)
redundantVis = visLcp[freq, scan, NP.append(redIndexesEW, redIndexesN)]
antPha, c, d, e = NP.linalg.lstsq(phaMatrix, NP.angle(redundantVis), rcond=None)
ewAntPha_lin = antPha[baselinesNumber*2:baselinesNumber*2+antNumberEW]
nAntPha_lin = antPha[baselinesNumber*2+antNumberEW:]

redundantVisN = visLcp[freq, scan, redIndexesN]
x_ini = NP.concatenate((NP.ones(baselinesNumber+antNumberN-1), NP.zeros(baselinesNumber+antNumberN-1)))
ls_res = least_squares(northGainsFunc_constrained, x_ini, args = (redundantVisN, antNumberN, baselinesNumber))
n_gains= real_to_complex(ls_res['x'])[baselinesNumber-1:]
nAntPha = NP.angle(n_gains)

redundantVisEW = visLcp[freq, scan, redIndexesEW]
x_ini = NP.concatenate((NP.ones(baselinesNumber+antNumberEW-2), NP.zeros(baselinesNumber+antNumberEW-2)))
ls_res = least_squares(eastWestGainsFunc_constrained, x_ini, args = (redundantVisEW, antNumberEW, baselinesNumber))
gains = real_to_complex(ls_res['x'])[baselinesNumber-1:]
ew_gains= NP.insert(gains, antNumberEW//2, (1+0j))
ewAntPha = NP.angle(ew_gains)

uvSize = 2049
uv = NP.zeros((uvSize, uvSize), dtype = 'complex')

#ewAntPha = ewAntPha_lin
#nAntPha = nAntPha_lin

for i in range(31):
    for j in range(65):
        O = uvSize//2
        uv[O + (i+1)*2, O + (j-32)*2] = visLcp[freq, scan, i*97+j] * NP.exp(1j * (-ewAntPha[j] + nAntPha[i]))
        uv[O - (i+1)*2, O - (j-32)*2] = NP.conj(uv[O + (i+1)*2, O + (j-32)*2])
for i in range(32):
    O = uvSize//2
    uv[O, O + (i-32)*2] = 0.5*(visLcp[freq, scan, antZeroRow[i]] * NP.exp(1j * (-ewAntPha[i] + ewAntPha[32])) + 
                               visLcp[freq, scan, antZeroRow[64 - i]] * NP.exp(1j * (-ewAntPha[32] + ewAntPha[64 - i])))
    uv[O, O + (32-i)*2] = NP.conj(uv[O, O + (i-32)*2])
        
image = NP.fft.fft2(NP.roll(NP.roll(uv,uvSize//2+1,0),uvSize//2+1,1));
image = NP.roll(NP.roll(image,uvSize//2-1,0),uvSize//2-1,1);

PL.figure()
#PL.title('no ZR SRH %s, %s, %d MHz'%(fitsFile[0].header['DATE-OBS'], fitsFile[0].header['TIME-OBS'], freqList[freq]*1e-3))
PL.title('SRH %s, %s, %d MHz'%(fitsFile[0].header['DATE-OBS'], fitsFile[0].header['TIME-OBS'], freqList[freq]*1e-3))
showUV = NP.abs(uv[uvSize//2-65:uvSize//2+65,uvSize//2-65:uvSize//2+65])
#showUV = NP.clip(showUV,0,0.2*NP.max(showUV))
showUV = showUV**.5
showUV = NP.roll(showUV,-1,axis=1)
showUV = NP.roll(showUV,-1,axis=0)
PL.imshow(showUV,cmap='hot')
PL.axis('off')
PL.tight_layout()

PL.figure()
#PL.title('no ZR SRH %s, %s, %d MHz'%(fitsFile[0].header['DATE-OBS'], fitsFile[0].header['TIME-OBS'], freqList[freq]*1e-3))
PL.title('SRH %s, %s, %d MHz'%(fitsFile[0].header['DATE-OBS'], fitsFile[0].header['TIME-OBS'], freqList[freq]*1e-3))
PL.imshow(NP.flip(image.real, 1),cmap='hot')
PL.axis('off')
PL.tight_layout()

#arcsecPerPixel = 4.9/4
#RAO = BadaryRAO(fitsFile[0].header['DATE-OBS'])
#pAngle = NP.deg2rad(sunpy_coordinates.get_sun_P(fitsFile[0].header['DATE-OBS']).to_value())
#omegaEarth = coordinates.earth.OMEGA_EARTH.to_value()
#hAngle = omegaEarth * (time[0,scan] - RAO.culmination)
#scaling = getPQScale(uvSize, NP.deg2rad(arcsecPerPixel * (uvSize - 1)/3600.)*2, time[0,scan], fitsFile[1].data['frequency'][freq])
#scale = AffineTransform(scale=(uvSize/scaling[0], uvSize/scaling[1]))
#shift = AffineTransform(translation=(-uvSize/2,-uvSize/2))
#rotate = AffineTransform(rotation = -pAngle)
#matrix = AffineTransform(matrix = getPQ2HDMatrix())
#back_shift = AffineTransform(translation=(uvSize/2,uvSize/2))
#
#lm_image = NP.flip(image.real, 1)
#dataResult0 = warp(lm_image,(shift + (scale + back_shift)).inverse)
#lcpData = warp(dataResult0,(shift + (matrix + back_shift)).inverse)
#lcpData = warp(NP.flip(lcpData, 0),(shift + (rotate + back_shift)).inverse)
#PL.figure()
#PL.imshow(lcpData, origin = 'lower',cmap='hot')
#PL.title('SRH %s, %s, %d MHz, baselines %d'%(fitsFile[0].header['DATE-OBS'], fitsFile[0].header['TIME-OBS'], freqList[freq]*1e-3, baselinesNumber))
#PL.axis('off')
#PL.tight_layout()
