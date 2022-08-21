#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 09:55:18 2020

@author: mariagloba
"""

import sys;
import numpy as NP
from srhFitsFile_doubleBase import SrhFitsFile
from skimage.transform import warp, AffineTransform
from astropy.io import fits
from astropy.time import Time, TimeDelta
import base2uvw as bl2uvw
import pylab as PL
#from sunpy import coordinates
import os
from zeep import Client
import datetime as DT



def findPhase(currentScan):
    lcpMaxTrace = []
    sPhaseCoefsLcp = NP.zeros(16)
    sLcpPhaseCorrection = NP.zeros(16)
    sRcpPhaseCorrection = NP.zeros(16)
    srhFits.setSizeOfUv(256)
    for sPhaseInd in range(-18,18):
        sLcpPhaseCorrection[:] = NP.deg2rad(sPhaseInd*10)
        srhFits.changeSouthPhase(sLcpPhaseCorrection, sRcpPhaseCorrection)
        srhFits.vis2uv(currentScan, phaseCorrect=True, amplitudeCorrect=False);
        srhFits.uv2lmImage()
        lcpMaxTrace.append(srhFits.lcp.real[128-32:128+32,128-32:128+32].mean())
    srhFits.setSizeOfUv(uvSize)
    
    phaseIndLcp = int(10*(NP.argmax(lcpMaxTrace) - 18) + .5)

  
    sPhaseCoefsLcp[15] = phaseIndLcp
    
    sLcpPhaseCorrection[:] = 0.
    sRcpPhaseCorrection[:] = 0.
    for j in range(16):
        for i in range(16):
            sLcpPhaseCorrection[i] +=  NP.deg2rad(sPhaseCoefsLcp[j] * (-1)**(i // (j + 1)))
    srhFits.changeSouthPhase(sLcpPhaseCorrection,sRcpPhaseCorrection)



with open('/home/serg/SRH/mf_currentFileName.txt') as f:
    fitsName = f.read()

#imSize = 256
imSize = 128
scan = 0
frequencyChannel = 8
#uvSize = 512
uvSize = 256
#arcsecPerPixel = 12#4.91104
arcsecPerPixel = 24

#fitsName = '/home/mariagloba/Work/fits/20200314/mf_20200314_012758.fit'
allImages = NP.zeros((imSize*8, imSize*4))

srhFits = SrhFitsFile(fitsName, uvSize)
srhFits.getHourAngle(scan)
#pAngle = NP.deg2rad(coordinates.get_sun_P(srhFits.dateObs).to_value())
srhFits.setCalibIndex(scan)

pAngle = 0
try:
    client = Client('http://ephemeris.rao.istp.ac.ru/?wsdl')
    result = client.service.Ephemeride('SSRT','sun',srhFits.dateObs)
    pAngle = NP.deg2rad(float(result[0]['PAngle']))
except:
    pass

srhFits.useDoubleBaselinesAmp = True

for frequencyChannel in range(len(srhFits.freqList)):
    srhFits.setFrequencyChannel(frequencyChannel)
    findPhase(scan)
    srhFits.vis2uv(scan, phaseCorrect=True, amplitudeCorrect = True, PSF=False);
    srhFits.uv2lmImage()
    
    scaling = srhFits.getPQScale(uvSize, NP.deg2rad(arcsecPerPixel * (uvSize - 1)/3600.))
    scale = AffineTransform(scale=(uvSize/scaling[0], uvSize/scaling[1]))
    shift = AffineTransform(translation=(-uvSize/2,-uvSize/2))
    rotate = AffineTransform(rotation = pAngle)
    matrix = AffineTransform(matrix = srhFits.getPQ2HDMatrix())
    back_shift = AffineTransform(translation=(uvSize/2,uvSize/2))
    
    lcpData = srhFits.lcp.real
    lcpData = warp(lcpData,(shift + (scale + back_shift)).inverse)
    lcpData = warp(lcpData,(shift + (matrix + back_shift)).inverse)
    lcpData = warp(lcpData,(shift + (rotate + back_shift)).inverse)
    x = frequencyChannel//8
    y = 7 - frequencyChannel%8
    allImages[y*imSize:(y+1)*imSize, x*imSize:(x+1)*imSize] = NP.flip(lcpData[uvSize//2-imSize//2:uvSize//2+imSize//2, uvSize//2-imSize//2:uvSize//2+imSize//2], 0)

allImagesHead = NP.ones((imSize*8+25, imSize*4)) * -1000
allImagesHead[:imSize*8, :] = allImages
PL.figure(figsize=(6,13))
PL.imshow(allImagesHead, origin = 'lower', cmap = 'hot', vmin = -0.0007, vmax = 0.004)
for frequencyChannel in range(len(srhFits.freqList)):
    x = frequencyChannel//8
    y = 7 - frequencyChannel%8
    PL.text(x*imSize, y*imSize+5, str(int(srhFits.freqList[frequencyChannel])), color = 'w', fontsize = 12)
PL.text(10, imSize*8+5, 'SRH ' + srhFits.dateObs.split('T')[0] + ' ' + srhFits.dateObs.split('T')[1]  + '    MHz', color = 'w', fontsize = 15)
PL.xticks([])
PL.yticks([])
dateName = DT.datetime.now().strftime("%Y%m%d")
pngName = '/home/serg/py/srh_'+dateName+'.png'
#pngName = 'srh_'+dateName+'_amp.png'
PL.savefig(pngName, bbox_inches='tight')