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
import phaMatrixGen as MG
import scipy.optimize as SPO
import scipy.signal as SPS
import os, fnmatch;


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

fitNames = findFits('/home/svlesovoi/SRH_DATA/SRH36/20210320','*.fit')

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
freqListLength = fitsFile[1].data['frequency'].shape[0]
samplesNumber = time.shape[1]
visibilityNumber = fitsFile[1].data['vis_lcp'].shape[1]//samplesNumber

visLcp = fitsFile[1].data['vis_lcp'].reshape(freqListLength,samplesNumber,visibilityNumber) * visScale
visRcp = fitsFile[1].data['vis_rcp'].reshape(freqListLength,samplesNumber,visibilityNumber) * visScale

antA = fitsFile[4].data['ant_A']
antB = fitsFile[4].data['ant_B']
redIndexesN = []
redIndexesN.append(NP.where((antA==98) & (antB==33))[0][0])
for i in range(30):
    redIndexesN.append(NP.where((antA==98+i) & (antB==98+i+1))[0][0])
redIndexesEW = []
for i in range(31):
    redIndexesEW.append(NP.where((antA==1+i) & (antB==2+i))[0][0])
#    print(i)
for i in range(65):
    redIndexesEW.append(NP.where((antA==32+i) & (antB==33+i))[0][0])
#    print(i)
    
scan = 0
freq = 0
polarization = 0

if polarization == 0:
    redN = visLcp[freq,scan,redIndexesN].copy()
    redEW = visLcp[freq,scan,redIndexesEW].copy()
else:
    redN = visRcp[freq,scan,redIndexesN].copy()
    redEW = visRcp[freq,scan,redIndexesEW].copy()

#
nPha1 = NP.append(NP.angle(redN), 0)
ewPha1 = NP.append(NP.angle(redEW[:64]), 0)
#
#ewAntPha = NP.zeros(65)
nAntPha, c, d, e = NP.linalg.lstsq(matrices.matrixPhaseN, nPha1, rcond=None)
ewAntPha, c, d, e = NP.linalg.lstsq(matrices.matrixPhaseEW, ewPha1, rcond=None)
nAntPha = nAntPha[2:] # nAntPha[1] - central antenna
ewAntPha = ewAntPha[1:]

ewnPhaMatrix = MG.phaMatrixGenEWN(65,31)
eastWestPhase = (NP.angle(visLcp[freq,scan,westInd0:westInd0+64]))
northRedVis = NP.zeros(31,dtype='complex')
northRedVis[0]  = visLcp[freq,scan,32]
northRedVis[1:] = visLcp[freq,scan,northInd0:northInd0+30]
northPhase = (NP.angle(northRedVis))
ewnPhase = NP.concatenate((eastWestPhase,northPhase))
ewnAntPhase, c, d, e = NP.linalg.lstsq(ewnPhaMatrix,ewnPhase)
ewnAntPhaLcp = ewnAntPhase[2:]
northAntPhaLcp = ewnAntPhaLcp[65:]
eastWestAntPhaLcp = ewnAntPhaLcp[:65]
nAntPha = northAntPhaLcp
ewAntPha = eastWestAntPhaLcp

northAntPhaLcpObj = NP.zeros(31)
opt = SPO.minimize(objVisSum,northAntPhaLcpObj,complex_to_real(visLcp[freq,scan,northInd0:northInd0+30]))
northAntPhaLcp = opt.x
northAntPhaRcpObj = NP.zeros(31)
opt = SPO.minimize(objVisSum,northAntPhaRcpObj,complex_to_real(visRcp[freq,scan,northInd0:northInd0+30]))
northAntPhaRcp = opt.x

westAntPhaLcpObj = NP.zeros(32)
opt = SPO.minimize(objVisSum,westAntPhaLcpObj,complex_to_real(visLcp[freq,scan,westIndList[0:31]]))
westAntPhaLcp = opt.x
westAntPhaRcpObj = NP.zeros(32)
opt = SPO.minimize(objVisSum,westAntPhaRcpObj,complex_to_real(visRcp[freq,scan,westIndList[0:31]]))
westAntPhaRcp = opt.x

eastAntPhaLcpObj = NP.zeros(32)
opt = SPO.minimize(objVisSum,eastAntPhaLcpObj,complex_to_real(visLcp[freq,scan,eastIndList[0:31]]))
eastAntPhaLcp = opt.x
eastAntPhaRcpObj = NP.zeros(32)
opt = SPO.minimize(objVisSum,eastAntPhaRcpObj,complex_to_real(visRcp[freq,scan,eastIndList[0:31]]))
eastAntPhaRcp = opt.x

eastWestAntPhaLcpObj = NP.zeros(65)
opt = SPO.minimize(objVisSum,eastWestAntPhaLcpObj,complex_to_real(visLcp[freq,scan,westInd0:westInd0+64]))
eastWestAntPhaLcp = opt.x
eastWestAntPhaRcpObj = NP.zeros(65)
opt = SPO.minimize(objVisSum,eastWestAntPhaRcpObj,complex_to_real(visRcp[freq,scan,westInd0:westInd0+64]))
eastWestAntPhaRcp = opt.x
#------------------------------------------------------------------------------------------------------------------------
northAntPhaLcpObj1D = northAntPhaLcp
opt = SPO.minimize(objVis1DScanSum,northAntPhaLcpObj1D,complex_to_real(visLcp[freq,scan,northInd0:northInd0+int((30+1)/2*30)]),method='Powell')
northAntPhaLcp1D = opt.x
#northAntPhaRcpObj1D = northAntPhaRcp
#opt = SPO.minimize(objVis1DScanSum,northAntPhaRcpObj1D,complex_to_real(visRcp[freq,scan,northInd0:northInd0+int((30+1)/2*30)]),method='Powell')
#northAntPhaRcp1D = opt.x

westAntPhaLcpObj1D = westAntPhaLcp
opt = SPO.minimize(objVis1DScanSum,westAntPhaLcpObj1D,complex_to_real(visLcp[freq,scan,westIndList]),method='Powell')
westAntPhaLcp1D = opt.x
#westAntPhaRcpObj1D = westAntPhaRcp
#opt = SPO.minimize(objVis1DScanSum,westAntPhaRcpObj1D,complex_to_real(visRcp[freq,scan,westIndList]),method='Powell')
#westAntPhaRcp1D = opt.x

eastAntPhaLcpObj1D = eastAntPhaLcp
opt = SPO.minimize(objVis1DScanSum,eastAntPhaLcpObj1D,complex_to_real(visLcp[freq,scan,eastIndList]),method='Powell')
eastAntPhaLcp1D = opt.x
#eastAntPhaRcpObj1D = eastAntPhaRcp
#opt = SPO.minimize(objVis1DScanSum,eastAntPhaRcpObj1D,complex_to_real(visRcp[freq,scan,eastIndList]),method='Powell')
#eastAntPhaRcp1D = opt.x
#------------------------------------------------------------------------------------------------------------------------

#nAntPha = northAntPhaLcp1D
#ewAntPha[:32] =  westAntPhaLcp1D
#ewAntPha[33:] =  eastAntPhaLcp1D

#nAntPha = northAntPhaLcp
#ewAntPha = eastWestAntPhaLcp

uvSize = 2049
uv = NP.zeros((uvSize, uvSize), dtype = 'complex')

for i in range(31):
    for j in range(65):
        O = uvSize//2
        uv[O + (i+1)*2, O + (j-32)*2] = visLcp[freq, scan, i*97+j] * NP.exp(1j * (-ewAntPha[j] + nAntPha[i]))
        uv[O - (i+1)*2, O - (j-32)*2] = NP.conj(uv[O + (i+1)*2, O + (j-32)*2])
for i in range(32):
    O = uvSize//2
    uv[O, O + (i-32)*2] = visLcp[freq, scan, antZeroRow[i]] * NP.exp(1j * (-ewAntPha[i]))
    uv[O, O + (i-32)*2] += visLcp[freq, scan, antZeroRow[64 - i]] * NP.exp(1j * (-ewAntPha[64 - i]))
    uv[O, O + (i-32)*2] /= 2
    uv[O, O + (32-i)*2] = NP.conj(uv[O, O + (i-32)*2])
        
image = NP.fft.fft2(NP.roll(NP.roll(uv,uvSize//2+1,0),uvSize//2+1,1));
image = NP.roll(NP.roll(image,uvSize//2-1,0),uvSize//2-1,1);
PL.clf()
PL.imshow(NP.flip(image.real, 1))
PL.title('f = %d MHz'%(fitsFile[1].data['frequency'][freq]/1e3))

arcsecPerPixel = 4.9/2
RAO = BadaryRAO('2021-02-23')
pAngle = NP.deg2rad(sunpy_coordinates.get_sun_P('2021-02-23').to_value())
omegaEarth = coordinates.earth.OMEGA_EARTH.to_value()
hAngle = omegaEarth * (time[0,scan] - RAO.culmination)
scaling = getPQScale(uvSize, NP.deg2rad(arcsecPerPixel * (uvSize - 1)/3600.)*2, time[0,scan], fitsFile[1].data['frequency'][freq])
scale = AffineTransform(scale=(uvSize/scaling[0], uvSize/scaling[1]))
shift = AffineTransform(translation=(-uvSize/2,-uvSize/2))
rotate = AffineTransform(rotation = -pAngle)
matrix = AffineTransform(matrix = getPQ2HDMatrix())
back_shift = AffineTransform(translation=(uvSize/2,uvSize/2))

lm_image = NP.roll(NP.flip(image.real, 1), -60, axis = 1)
dataResult0 = warp(lm_image,(shift + (scale + back_shift)).inverse)
lcpData = warp(dataResult0,(shift + (matrix + back_shift)).inverse)
lcpData = warp(NP.flip(lcpData, 0),(shift + (rotate + back_shift)).inverse)
PL.figure()
PL.imshow(lcpData, origin = 'lower')
