#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 14:11:02 2021

@author: sergeyvlesovoi
"""
import os, fnmatch;
#from srhFitsFile36_amp import SrhFitsFile
from srhFitsFile36 import SrhFitsFile
import numpy as NP;
import datetime as DT;
import pylab as PL;
import sys;
import matplotlib;
import phaMatrixGen
from BadaryRAO import BadaryRAO
from astropy import coordinates
from sunpy import coordinates as sunpy_coordinates
from skimage.transform import warp, AffineTransform
from astropy import constants

dt_major = 300.;
dt_minor = 30.;
freq_first = 2;
freq_intervals = 5;

def hhmm_format(t, pos):
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  return '%02d:%02d' % (hh,mm);

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
#fitNames = findFits('SRH36_temp_20210416_1/','*.fit')
#fitNames = findFits('SRH36_temp_20210419_1/','*.fit')
#fitNames = findFits('SRH36_temp_20210501_4/','*.fit')
#fitNames = findFits('SRH36_temp_20210508_1/','*.fit')
#fitNames = findFits('SRH36_temp_20210509_2/','*.fit')
#fitNames = findFits('SRH36_temp_20210509_3/','*.fit')
#fitNames = findFits('SRH36_temp_20210522_5/','*.fit')
fitNames = findFits('SRH36_temp_20210527_60/','*.fit')
fitNames.sort();
fitNames = fitNames[1:];

visScale = 1/(2e6*49)
ampScale = visScale / 128
uvSize = 1025

for fName in fitNames:
    if (fName == fitNames[0]):
        sF = SrhFitsFile(fName, uvSize)
    else:
        sF.append(fName)

sF.useNonlinearApproach = True
sF.baselines = 3
freq = 0
sF.getHourAngle(0)
phaseNoiseLevel = .3
amplitudeNoiseLevel = 3

nPhaseLcp = NP.angle(sF.visLcp[freq,:,3007:3007+30+29+28])
ewPhaseLcp = NP.angle(sF.visLcp[freq,:,3472:3472+96+95+94])
nPhaseNoiseLcp = NP.std(nPhaseLcp,axis=0)
ewPhaseNoiseLcp = NP.std(ewPhaseLcp,axis=0)

ewFlaggedVis = NP.array(NP.where(ewPhaseNoiseLcp > phaseNoiseLevel))
nFlaggedVis = NP.array(NP.where(nPhaseNoiseLcp > phaseNoiseLevel))
ewFlaggedVis = ewFlaggedVis[0][::2]
nFlaggedVis = nFlaggedVis[0][::2]

nAmplitudeLcp = NP.abs(sF.visLcp[freq,:,3007:3007+30+29+28])
ewAmplitudeLcp = NP.abs(sF.visLcp[freq,:,3472:3472+96+95+94])
nAmplitudeNoiseLcp = NP.mean(nAmplitudeLcp,axis=0)/NP.std(nAmplitudeLcp,axis=0)
ewAmplitudeNoiseLcp = NP.mean(ewAmplitudeLcp,axis=0)/NP.std(ewAmplitudeLcp,axis=0)

ewFlaggedVis = NP.array(NP.where(ewAmplitudeNoiseLcp < amplitudeNoiseLevel))
nFlaggedVis = NP.array(NP.where(nAmplitudeNoiseLcp < amplitudeNoiseLevel))
ewFlaggedVis = ewFlaggedVis[0][::2]
nFlaggedVis = nFlaggedVis[0][::2]

sF.visLcp[freq,:,3007 + nFlaggedVis] = 0 + 0j
sF.visLcp[freq,:,3472 + ewFlaggedVis] = 0 + 0j
sF.visRcp[freq,:,3007 + nFlaggedVis] = 0 + 0j
sF.visRcp[freq,:,3472 + ewFlaggedVis] = 0 + 0j

flaggedAnts = []

for k in range(ewFlaggedVis.shape[0]):
    antAName = str(sF.antennaA[ewFlaggedVis[k] + 3472])
    flaggedAnts.append(str(sF.hduList[2].data['ant_name'][NP.where(sF.antennaNumbers == antAName)][0]))

for k in range(nFlaggedVis.shape[0]):
    antAName = str(sF.antennaA[nFlaggedVis[k] + 3007])
    flaggedAnts.append(str(sF.hduList[2].data['ant_name'][NP.where(sF.antennaNumbers == antAName)][0]))

#sF.updateAntennaPhase(freq, baselinesNumber=3)
#sF.flag(','.join(flaggedAnts))
sF.calculatePhaseCalibration(baselinesNumber = sF.baselines)
sF.setFrequencyChannel(freq)

arcsecPerPixel = 4.9/2
time = sF.freqTime
RAO = BadaryRAO(sF.hduList[0].header['DATE-OBS'])
pAngle = NP.deg2rad(sunpy_coordinates.get_sun_P(sF.hduList[0].header['DATE-OBS']).to_value())
omegaEarth = coordinates.earth.OMEGA_EARTH.to_value()

imageLcp = []
imageRcp = []
for scan in range(sF.freqTime.shape[1]):
    sF.vis2uv(scan)
    sF.uv2lmImage()
    hAngle = omegaEarth * (time[0,scan] - RAO.culmination)
    scaling = getPQScale(uvSize, NP.deg2rad(arcsecPerPixel * (uvSize - 1)/3600.)*2, time[0,scan], sF.hduList[1].data['frequency'][freq])
    scale = AffineTransform(scale=(uvSize/scaling[0], uvSize/scaling[1]))
    shift = AffineTransform(translation=(-uvSize/2,-uvSize/2))
    rotate = AffineTransform(rotation = -pAngle)
    matrix = AffineTransform(matrix = getPQ2HDMatrix())
    back_shift = AffineTransform(translation=(uvSize/2,uvSize/2))
    
    lmImage = warp(sF.lcp.real,(shift + (scale + back_shift)).inverse)
    hdImage = warp(lmImage,(shift + (matrix + back_shift)).inverse)
    hdImage = warp(NP.flip(hdImage, 0),(shift + (rotate + back_shift)).inverse)
    imageLcp.append(hdImage)

    lmImage = warp(sF.rcp.real,(shift + (scale + back_shift)).inverse)
    hdImage = warp(lmImage,(shift + (matrix + back_shift)).inverse)
    hdImage = warp(NP.flip(hdImage, 0),(shift + (rotate + back_shift)).inverse)
    imageRcp.append(hdImage)

#flare = sF.ampLcp[4,:,:].mean(axis=1) + sF.ampRcp[4,:,:].mean(axis=1)
#flare_max = flare.max()
#flare -= flare.min()
#flare /= flare.max()
#
#sun = sF.ampLcp[4,:,28] + sF.ampRcp[4,:,28]
#sun_max = sun.max()
#sun -= sun.min()
#sun /= sun.max()
#
#SNR_RT32 = sun[50:80].mean()/sun[0:30].std()
#
#fig, pl = PL.subplots(figsize=(8,8))
#fig.suptitle('SRH and RT-32 the Sun response, 20210508, 4900 MHz')
#pl.xaxis.set_major_locator(PL.MultipleLocator(120))
#pl.xaxis.set_minor_locator(PL.MultipleLocator(20))
#pl.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
#pl.set_xlabel('UT')
#pl.set_ylabel('normalized')
#pl.plot(sF.freqTime[4],sun,label='RT-32, SNR %d'%SNR_RT32)
#pl.plot(sF.freqTime[4],flare*flare_max/sun_max,label='SRH')
#pl.legend()
