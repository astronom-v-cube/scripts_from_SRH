#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 22 01:57:10 2021

@author: svlesovoi
"""

from srhFitsFile36 import SrhFitsFile
import srh36MS2
import numpy as NP
import os, fnmatch, sys
from astropy.io import fits
import skimage
from skimage import measure
from skimage.transform import warp, AffineTransform
from zeep import Client
import scipy.signal
import scipy.constants
from matplotlib.patches import Circle, Ellipse
import matplotlib.colors
import ftplib;
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import scipy.optimize as opt

gXSize = 1
gArcsecPerPixel = 1

def fitFuncTb(f, A, B, C):
    return A + B*f + C*f**-1.8

def arcmin_format(xy, pos):
  return '%2d' % ((xy - gXSize/2) * gArcsecPerPixel / 60);

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;

def prepareCleanImages(fileName, qSunTb):
#    Tb_2800 = 27100
#    Tb_3100 = 24200
    
    fd = fits.open(fileName)
    
    pAngle = 0
    try:
        client = Client('http://ephemeris.rao.istp.ac.ru/?wsdl')
        result = client.service.Ephemeride('SSRT','sun',fd[0].header['DATE-OBS'])
        pAngle = NP.deg2rad(float(result[0]['PAngle']))
    except:
        pass
    
    operatingFrequency = int(fd[0].header['CRVAL3'] * 1e-6 + .5)
    rcpImage = fd[0].data[0,0]
    lcpImage = fd[0].data[1,0]

    srh_x_size = fd[0].header['NAXIS1']
    srh_y_size = fd[0].header['NAXIS2']
    
    arcsecPerPixel = fd[0].header['CDELT1']*3600
    resultArcsecPerPixel = 4.911
    resultScale = arcsecPerPixel / resultArcsecPerPixel
    
    global gXSize
    gXSize = srh_x_size
    global gArcsecPerPixel
    gArcsecPerPixel = resultArcsecPerPixel
    
    iImage = rcpImage + lcpImage
    vImage = rcpImage - lcpImage
    
    iImageHist = NP.histogram(iImage,bins=1000)
    for promFactor in range(5):
        maxInds = scipy.signal.find_peaks(iImageHist[0],prominence = 2000//(promFactor+1))
        if maxInds[0].shape[0] > 1:
            break
    if maxInds[0].shape[0] > 1:
        iImage -= iImageHist[1][maxInds[0][0]]
        iImage /= (iImageHist[1][maxInds[0][1]] - iImageHist[1][maxInds[0][0]])
        vImage /= (iImageHist[1][maxInds[0][1]] - iImageHist[1][maxInds[0][0]])
    else:
        sunLevel = iImage[256-64:256+64,256-64:256+64].mean()
        skyLevel = 0
        iImage -= sunLevel
        iImage /= (sunLevel - skyLevel)
    
        sunLevel = vImage[256-64:256+64,256-64:256+64].mean()
        vImage -= sunLevel
        vImage /= (sunLevel - skyLevel)

    if NP.abs(operatingFrequency - 2.8e9) < NP.abs(operatingFrequency - 3.1e9):
        iImage *= qSunTb
        vImage *= qSunTb
    else:
        iImage *= Tb_3100
        vImage *= Tb_3100

    scale = AffineTransform(scale=(-resultScale,-resultScale))
    shift = AffineTransform(translation=(-srh_y_size/2,-srh_y_size/2))
    rotate = AffineTransform(rotation = -pAngle)
    back_shift = AffineTransform(translation=(srh_y_size/2,srh_y_size/2))
    
    iImage = warp(iImage,(shift + (rotate + back_shift)).inverse)
    vImage = warp(vImage,(shift + (rotate + back_shift)).inverse)
    iImage = warp(iImage,(shift + (scale + back_shift)).inverse)
    vImage = warp(vImage,(shift + (scale + back_shift)).inverse)
    
    pHeader = fits.Header();
    pHeader['DATE-OBS']     = fd[0].header['DATE-OBS']
    pHeader['T-OBS']        = fd[0].header['DATE-OBS']
    pHeader['INSTRUME']     = fd[0].header['INSTRUME']
    pHeader['ORIGIN']       = fd[0].header['ORIGIN']
    pHeader['FREQUENC']     = ('%d') % (fd[0].header['CRVAL3']/1e6 + 0.5)
    pHeader['CDELT1']       = resultArcsecPerPixel
    pHeader['CDELT2']       = resultArcsecPerPixel
    pHeader['CRPIX1']       = srh_x_size // 2
    pHeader['CRPIX2']       = srh_x_size // 2
    pHeader['CTYPE1']       = 'HPLN-TAN'
    pHeader['CTYPE2']       = 'HPLT-TAN'
    pHeader['CUNIT1']       = 'arcsec'
    pHeader['CUNIT2']       = 'arcsec'

    saveFitsIhdu = fits.PrimaryHDU(header=pHeader, data=iImage)
    saveFitsIpath = 'srh_I_%s_%04d.fit'%(fd[0].header['DATE-OBS'].split('.')[0],fd[0].header['CRVAL3']*1e-6 + .5)
    hduList = fits.HDUList(saveFitsIhdu)
    hduList.writeto(saveFitsIpath)
    
    saveFitsVhdu = fits.PrimaryHDU(header=pHeader, data=vImage)
    saveFitsVpath = 'srh_V_%s_%04d.fit'%(fd[0].header['DATE-OBS'].split('.')[0],fd[0].header['CRVAL3']*1e-6 + .5)
    hduList = fits.HDUList(saveFitsVhdu)
    hduList.writeto(saveFitsVpath)
    
    fd.close()
    
#------------------------------------------------------------------------------
fitPath = sys.argv[1]
frequency = int(sys.argv[2])
cleanTresh = sys.argv[3]
scan = sys.argv[4]

fitNames =  findFits(fitPath,'*.fit')
fitNames.sort()

guess  = [1, 1, 1]

ZirinFrequency = [1.4, 1.6, 1.8, 2.0, 2.4, 2.8,  3.2, 3.6,  4.2, 5.0,  5.8, 7.0, 8.2, 9.4, 10.6, 11.8, 13.2, 14.8, 16.4, 18.0]
ZirinTb = [70.5, 63.8, 52.2, 42.9, 32.8, 27.1, 24.2, 21.7, 19.4, 17.6, 15.9, 14.1, 12.9, 12.2, 11.3, 11.0, 10.8, 10.8, 10.7, 10.3]
fitTbParams, _ = opt.curve_fit(fitFuncTb, ZirinFrequency, ZirinTb, p0=guess)

for fileName in fitNames:
    file = SrhFitsFile(fileName, 1025)
    file.useNonlinearApproach = True
    file.getHourAngle(0)
    for frequency in range(4):
        file.solarPhase(frequency)
        file.updateAntennaPhase(frequency, baselinesNumber = 3)
        file.setFrequencyChannel(frequency)
        file.vis2uv(0, average = 20)
        file.centerDisk()
    
        saveName = fileName.split('.')[0] + '.ms'
        ms2Table = srh36MS2.SrhMs2(saveName)
        ms2Table.createMS(file, frequencyChannel = [int(frequency)], phaseCorrect = True, amplitudeCorrect = True)
        
        flags_ew_lcp = NP.where(file.ewAntAmpLcp[frequency] == 1e6)[0] + 1
        flags_ew_rcp = NP.where(file.ewAntAmpRcp[frequency] == 1e6)[0] + 1
        flags_ew = NP.unique(NP.append(flags_ew_lcp, flags_ew_rcp))
        flags_n_lcp = NP.where(file.nAntAmpLcp[frequency] == 1e6)[0]+98
        flags_n_rcp = NP.where(file.nAntAmpRcp[frequency] == 1e6)[0]+98
        flags_n = NP.unique(NP.append(flags_n_lcp, flags_n_rcp))
        flags = ','.join(map(str, NP.append(flags_ew, flags_n)))
        
        command = 'casa -c casaclean.py \'' + saveName + '\'  \'' + flags + '\' \'' +  cleanTresh + '\' \'' + scan  + '\''
        os.system(command)
        
        currentDiskTb = fitFuncTb(file.freqList[frequency]*1e-6,fitTbParams[0],fitTbParams[1],fitTbParams[2])*1e3
        prepareCleanImages(saveName.split('.')[0] + '_clean_image.fit',currentDiskTb)
        
        clnjunks = ['.ms', '.ms.flagversions', '*clean*.fit']
        for clnjunk in clnjunks:
            if os.path.exists(saveName.split('.')[0] + clnjunk):
                os.system('rm -rf ' + saveName.split('.')[0] + clnjunk)
