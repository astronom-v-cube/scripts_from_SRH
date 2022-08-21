#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 15:45:09 2021

@author: mariagloba
"""

from astropy.io import fits
import numpy as NP
import pylab as PL
from astropy import coordinates
from astropy import constants
from BadaryRAO import BadaryRAO
from sunpy import coordinates as sunpy_coordinates
from skimage.transform import warp, AffineTransform
import matrices_conj as matrices
import phaMatrixGen as MG
import scipy.optimize as SPO
import scipy.signal as SPS
from scipy.stats import linregress
import os, fnmatch
import datetime as DT
from scipy.optimize import minimize, least_squares


def antEWFormat(f, pos):
    if f < 64:
        return 'W%02d' % ((64 - f)//2)
    elif f == 64:
        return 'C00' % ((64 - f)//2)
    else:
        return 'E%02d' % ((f - 64)//2)

def antNFormat(f, pos):
    if f < 64:
        return 'N%02d' % ((64 - f)//2)
    elif f == 64:
        return 'C00' % ((64 - f)//2)
    else:
        return 'N%02d' % ((f - 64)//2)

def real_to_complex(z):
    return z[:len(z)//2] + 1j * z[len(z)//2:]

def complex_to_real(z):
    return NP.concatenate((NP.real(z), NP.imag(z)))

def objVisSum(phases, visibilities):
    z = real_to_complex(visibilities)
    res = 0 + 1j*0
    for k in range(len(z)):
        res += NP.exp(1j*(phases[k + 1] - phases[k])) * z[k]
    return -NP.abs(res)

def build1DScan(scanVis, visNumber, antPhases):
    scanSpectrum = NP.zeros(512,dtype='complex')
    pairInd = 0
    for pair in range(visNumber):
        for ant in range(visNumber - pair):
            scanSpectrum[pair + 1] += scanVis[pairInd] * NP.exp(1j*(antPhases[ant + 1 + pair] - antPhases[ant]))
            pairInd += 1
        scanSpectrum[512 - pair - 1] = NP.conj(scanSpectrum[pair + 1])
    return NP.roll(NP.fft.fft(scanSpectrum),256)
    
def objVis1DScanSum(phases, visibilities):
    return -NP.abs(build1DScan(real_to_complex(visibilities), len(phases) - 1, phases).real).sum()

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
    g = gP - gQ;
      
    pqMatrix = NP.zeros((3,3))
    pqMatrix[0, 0] =  NP.cos(gP) - NP.cos(g)*NP.cos(gQ)
    pqMatrix[0, 1] = -NP.cos(g)*NP.cos(gP) + NP.cos(gQ)
    pqMatrix[1, 0] =  NP.sin(gP) - NP.cos(g)*NP.sin(gQ)
    pqMatrix[1, 1] = -NP.cos(g)*NP.sin(gP) + NP.sin(gQ)
    pqMatrix /= NP.sin(g)**2.
    pqMatrix[2, 2] = 1.
    return pqMatrix

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
#fitNames = findFits('/home/svlesovoi/Documents/Python Scripts/SRH36_temp_20210322_2','*.fit')
fitNames = findFits('/home/svlesovoi/Documents/Python Scripts/SRH36_temp_20210410_1','*.fit')

fitNames.sort()
fitsFile = fits.open(fitNames[-1])

visScale = 1/(2e6*49)
ampScale = visScale / 128
westInd0 = 3472
eastInd0 = 3505
northInd0 = 3007
westIndList = []
eastIndList = []
eastWestIndList = []
eastWestInd = 0
westInd = 0

for i in range(31):
    for j in range(31 - i):
        westIndList.append(westInd0 + j + eastWestInd)
    eastWestInd += 96 - i

eastWestInd = 0
for i in range(31):
    for j in range(31 - i):
        eastIndList.append(eastInd0 + j + eastWestInd)
    eastWestInd += 96 - i

eastWestInd = 0
for i in range(64):
    for j in range(64 - i):
        eastWestIndList.append(westInd0 + j + eastWestInd)
    eastWestInd += 96 - i

antZeroRow = fitsFile[3].data['ant_zero_row'][:97]
time = fitsFile[1].data['time']
date_obs = fitsFile[0].header['DATE-OBS']
time_obs = fitsFile[0].header['TIME-OBS']
freqListLength = fitsFile[1].data['frequency'].shape[0]
freqList = fitsFile[1].data['frequency']*1e3
samplesNumber = time.shape[1]
visibilityNumber = fitsFile[1].data['vis_lcp'].shape[1]//samplesNumber

visLcp = fitsFile[1].data['vis_lcp'].reshape(freqListLength,samplesNumber,visibilityNumber) * visScale
visRcp = fitsFile[1].data['vis_rcp'].reshape(freqListLength,samplesNumber,visibilityNumber) * visScale

polarization = 0
calScan = 5
freq = 2
flagsEW = [22, 24, 44] # W10, W8, E12
flagsN = [16, 21, 28]
antA = fitsFile[4].data['ant_A']
antB = fitsFile[4].data['ant_B']
antNumberN = 31
antNumberEW = 65
baselinesNumber = 3

uvSize = 2049
uv = NP.zeros((uvSize, uvSize), dtype = 'complex')
images = NP.zeros((samplesNumber,uvSize,uvSize))

ewnPhaMatrix = MG.phaMatrixGenEWN(65,31)
ewAntPha = NP.zeros(65)
nAntPha = NP.zeros(31)

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
        
redundantVisN = visLcp[freq, calScan, redIndexesN]
x_ini = NP.concatenate((NP.ones(baselinesNumber+antNumberN-1), NP.zeros(baselinesNumber+antNumberN-1)))
ls_res = least_squares(northGainsFunc_constrained, x_ini, args = (redundantVisN, antNumberN, baselinesNumber))
n_gains= real_to_complex(ls_res['x'])[baselinesNumber-1:]
nAntPha = NP.angle(n_gains)

redundantVisEW = visLcp[freq, calScan, redIndexesEW]
x_ini = NP.concatenate((NP.ones(baselinesNumber+antNumberEW-2), NP.zeros(baselinesNumber+antNumberEW-2)))
ls_res = least_squares(eastWestGainsFunc_constrained, x_ini, args = (redundantVisEW, antNumberEW, baselinesNumber))
gains = real_to_complex(ls_res['x'])[baselinesNumber-1:]
ew_gains= NP.insert(gains, antNumberEW//2, (1+0j))
ewAntPha = NP.angle(ew_gains)

for scan in range(samplesNumber):
    eastWestPhase = (NP.angle(visLcp[freq,scan,westInd0:westInd0+64]))
    eastWestPhase[2] = 0
    eastWestPhase[3] = 0
    eastWestPhase[21] = 0
    eastWestPhase[22] = 0
    eastWestPhase[41] = 0
    eastWestPhase[44] = 0
    eastWestPhase[55] = 0
    eastWestPhase[56] = 0
    northRedVis = NP.zeros(31,dtype='complex')
    northRedVis[0]  = visLcp[freq,scan,32]
    northRedVis[1:] = visLcp[freq,scan,northInd0:northInd0+30]
    northRedVis[15] = 0
    northRedVis[16] = 0
    northRedVis[20] = 0
    northRedVis[21] = 0
    northRedVis[27] = 0
    northRedVis[28] = 0
    northPhase = (NP.angle(northRedVis))
    ewnPhase = NP.concatenate((eastWestPhase,northPhase))
    ewnAntPhase, c, d, e = NP.linalg.lstsq(ewnPhaMatrix,ewnPhase)
    ewnAntPhaLcp = ewnAntPhase[2:]
#    ewAntPha += ewnAntPhaLcp[:65]
#    nAntPha += ewnAntPhaLcp[65:]
    
#ewAntPha /= samplesNumber
#nAntPha /= samplesNumber

#        nAntPha[flagsN] = 0
#        ewAntPha[flagsEW] = 0

for scan in range(samplesNumber):
    if polarization == 0:
        for i in range(31):
            for j in range(65):
                O = uvSize//2
                uv[O + (i+1)*2, O + (j-32)*2] = visLcp[freq, scan, i*97+j] * NP.exp(1j * (-ewAntPha[j] + nAntPha[i]))
                uv[O - (i+1)*2, O - (j-32)*2] = NP.conj(uv[O + (i+1)*2, O + (j-32)*2])
        for i in range(32):
            O = uvSize//2
            uv[O, O + (i-32)*2] =  visLcp[freq, scan, antZeroRow[i]] * NP.exp(1j * (-ewAntPha[i] + ewAntPha[32]))
            uv[O, O + (i-32)*2] += visLcp[freq, scan, antZeroRow[64 - i]] * NP.exp(1j * (ewAntPha[64 - i] - ewAntPha[32]))
            uv[O, O + (i-32)*2] /= 2
            uv[O, O + (32-i)*2] = NP.conj(uv[O, O + (i-32)*2])
    else:
        for i in range(31):
            for j in range(65):
                O = uvSize//2
                uv[O + (i+1)*2, O + (j-32)*2] = visRcp[freq, scan, i*97+j] * NP.exp(1j * (-ewAntPha[j] + nAntPha[i]))
                uv[O - (i+1)*2, O - (j-32)*2] = NP.conj(uv[O + (i+1)*2, O + (j-32)*2])
        for i in range(32):
            O = uvSize//2
            uv[O, O + (i-32)*2] = visRcp[freq, scan, antZeroRow[i]] * NP.exp(1j * (-ewAntPha[i] + ewAntPha[32]))
            uv[O, O + (i-32)*2] += visRcp[freq, scan, antZeroRow[64 - i]] * NP.exp(1j * (ewAntPha[64 - i] - ewAntPha[32]))
            uv[O, O + (i-32)*2] /= 2
            uv[O, O + (32-i)*2] = NP.conj(uv[O, O + (i-32)*2])
            
    image = NP.fft.fft2(NP.roll(NP.roll(uv,uvSize//2+1,0),uvSize//2+1,1))
    image = NP.roll(NP.roll(image,uvSize//2-1,0),uvSize//2-1,1);
    images[scan,:,:] = image.real

fig = PL.figure(figsize=(10,9))
fig.subplots_adjust(left=0.01,right=.99,wspace=0.01, hspace=0.01, bottom=0.01)
fig.tight_layout()
ims = fig.subplots(samplesNumber // 5,samplesNumber // 4)
fig.suptitle('SRH lm %s, %s, %s, %d MHz'%(date_obs, time_obs, 'RCP' if polarization  else 'LCP',fitsFile[1].data['frequency'][freq]/1e3))
dW = 1200
for scan in range(samplesNumber):
    ims[scan//5,scan%5].imshow(NP.flip(images[scan,(uvSize - dW)//2:(uvSize + dW)//2,(uvSize - dW)//2:(uvSize + dW)//2], 1),cmap='hot')
    ims[scan//5,scan%5].axis('off')

avgImage = images.real.sum(0)
arcsecPerPixel = 4.9/4
RAO = BadaryRAO(fitsFile[0].header['DATE-OBS'])
pAngle = NP.deg2rad(sunpy_coordinates.get_sun_P(fitsFile[0].header['DATE-OBS']).to_value())
omegaEarth = coordinates.earth.OMEGA_EARTH.to_value()
hAngle = omegaEarth * (time[0,scan] - RAO.culmination)
scaling = getPQScale(uvSize, NP.deg2rad(arcsecPerPixel * (uvSize - 1)/3600.)*2, time[0,scan], fitsFile[1].data['frequency'][freq])
scale = AffineTransform(scale=(uvSize/scaling[0], uvSize/scaling[1]))
shift = AffineTransform(translation=(-uvSize/2,-uvSize/2))
rotate = AffineTransform(rotation = -pAngle)
matrix = AffineTransform(matrix = getPQ2HDMatrix())
back_shift = AffineTransform(translation=(uvSize/2,uvSize/2))

lm_image = NP.flip(avgImage, 1)
dataResult0 = warp(lm_image,(shift + (scale + back_shift)).inverse)
showImage = warp(dataResult0,(shift + (matrix + back_shift)).inverse)
showImage = warp(NP.flip(showImage, 0),(shift + (rotate + back_shift)).inverse)
PL.figure()
PL.imshow(showImage, origin = 'lower',cmap='hot')
PL.title('SRH %s, %s, %d MHz'%(fitsFile[0].header['DATE-OBS'], fitsFile[0].header['TIME-OBS'], freqList[freq]*1e-3))

maxTrace = []

for i in range(samplesNumber):
    maxTrace.append(images[i,:,:].max())
