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
from scipy import signal
import json


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

win = signal.windows.gaussian(11,5)

dateName = DT.datetime.now().strftime("%Y%m%d");
fitNames = findFits('SRH36_temp_20210412_1/','*.fit')
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

freq = 0
sF.setCalibIndex(20)
srhFreqList = sF.hduList[1].data['frequency']
phaseNoiseLevel = 0.1
nPhaseLcp = NP.angle(sF.visLcp[freq,:,3007:3007+30])
ewPhaseLcp = NP.angle(sF.visLcp[freq,:,3472:3472+96])
nPhaseNoiseLcp = NP.std(nPhaseLcp,axis=0)
ewPhaseNoiseLcp = NP.std(ewPhaseLcp,axis=0)

ewFlaggedVis = NP.array(NP.where(ewPhaseNoiseLcp > phaseNoiseLevel))
nFlaggedVis = NP.array(NP.where(nPhaseNoiseLcp > phaseNoiseLevel))
ewFlaggedVis = ewFlaggedVis[0][::2]
nFlaggedVis = nFlaggedVis[0][::2]

flaggedAnts = []
for k in range(ewFlaggedVis.shape[0]):
    antAName = str(sF.antennaA[ewFlaggedVis[k] + 3472])
    flaggedAnts.append(sF.hduList[2].data['ant_name'][NP.where(sF.antennaNumbers == antAName)])

for k in range(nFlaggedVis.shape[0]):
    antAName = str(sF.antennaA[nFlaggedVis[k] + 3007])
    flaggedAnts.append(sF.hduList[2].data['ant_name'][NP.where(sF.antennaNumbers == antAName)])

sF.setFrequencyChannel(freq)

arcsecPerPixel = 4.9/2
time = sF.freqTime
RAO = BadaryRAO(sF.hduList[0].header['DATE-OBS'])
pAngle = NP.deg2rad(sunpy_coordinates.get_sun_P(sF.hduList[0].header['DATE-OBS']).to_value())
omegaEarth = coordinates.earth.OMEGA_EARTH.to_value()

imageLcp = []
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


imArr = NP.array(imageLcp)
cols = 20
rows = imArr.shape[0] // cols
dL = 140
x0 = 400
y0 = 380
srcArr = NP.zeros((rows*dL,cols*dL))

fluxTrace = []
for r in range(rows): 
    for c in range(cols): 
        srcArr[r*dL:(r+1)*dL,c*dL:(c+1)*dL] = imArr[r*cols + c,y0-dL//2:y0+dL//2,x0-dL//2:x0+dL//2]
        fluxTrace.append(srcArr[r*dL:(r+1)*dL,c*dL:(c+1)*dL].mean())

arFlux0 = NP.array(fluxTrace0)*5e-16
arFlux1 = NP.array(fluxTrace1)*5e-16
arFlux2 = NP.array(fluxTrace2)*5e-16
smArFlux0 = signal.convolve(arFlux0,win,mode='same')
smArFlux1 = signal.convolve(arFlux1,win,mode='same')
smArFlux2 = signal.convolve(arFlux2,win,mode='same')

fig = PL.figure()
sub = fig.add_subplot(1,1,1);
sub.set_ylabel('flux');
sub.set_xlabel('UT');
sub.xaxis.set_major_locator(PL.MultipleLocator(900));
sub.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
sub.xaxis.set_minor_locator(PL.MultipleLocator(120));
sub.set_xlim(2.3*3600,2.5*3600)
#sub.set_ylim(0,5e-7)

#sub.plot(xrays_time[t0:],xrays_4_8[t0:],label='GOES X-Ray 0.1-0.8 nm',color='red',markersize=0.2)
#sub.plot(xrays_time[t0:],xrays_005_04[t0:],label='GOES X-Ray 0.05-0.4 nm',color='blue',markersize=0.2)
#for freq in range(srhFreqList.shape[0]):
#    sub.plot(srhTime[freq],srhCorrI[freq]*1e-4,'.',markersize=0.2,color=rainbowColors(100+(srhFreqList.shape[0] - freq)*20),label='SRH %d MHz'%(srhFreqList[freq]*1e-3))
#    sub.plot(srhTime[freq],srhCorrV[freq]*1e-4)
#sub.plot(srhTime[0],(srhMeanFluxI - srhMeanFluxISmoothed)*5e-9 + 1e-8)
sub.plot(sF.freqTime[0],smArFlux0,label='SRH %d MHz'%(srhFreqList[0]*1e-3))
sub.plot(sF.freqTime[0],smArFlux1,label='SRH %d MHz'%(srhFreqList[1]*1e-3))
sub.plot(sF.freqTime[0],smArFlux2,label='SRH %d MHz'%(srhFreqList[2]*1e-3))
sub.grid()
sub.legend()
#sub.set_title('SRH and GOES , %s'%(sF[0].header['DATE-OBS']))

goesFile = open('xrays-1-day.json')
goesData = json.load(goesFile)
goesFile.close()
N = len(goesData)
xrays_time = NP.zeros(N//2)
xrays_4_8 = NP.zeros(N//2)
xrays_005_04 = NP.zeros(N//2)
t0 = 200
for i in range(N):
    if i%2:
        xrays_4_8[i//2] = goesData[i]['flux']
        hhmm = goesData[i]['time_tag'].split('T')[1].split(':')[0:2]
        xrays_time[i//2] = 3600*int(hhmm[0]) + 60*int(hhmm[1])
        if xrays_time[i//2] > 10*3600:
            xrays_time[i//2] -= 24*3600
    else:
        xrays_005_04[i//2] = goesData[i]['flux']
sub.plot(xrays_time[t0:],xrays_4_8[t0:],label='GOES X-Ray 0.1-0.8 nm',color='red',markersize=0.2)
