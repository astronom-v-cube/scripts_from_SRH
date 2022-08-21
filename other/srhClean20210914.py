#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 31 09:41:28 2021

@author: sergeyvlesovoi
"""


from srhFitsFile36 import SrhFitsFile
import srh36MS2
import numpy as NP
import pylab as PL
import os, fnmatch
from astropy.io import fits
from skimage.transform import warp, AffineTransform
from zeep import Client
import scipy.signal
import scipy.constants
import datetime
from ZirinTb import ZirinTb
from optparse import OptionParser
import skimage
from scipy.stats import linregress
import json

gXSize = 1
gArcsecPerPixel = 1
gIImage = []
gVImage = []
gPSF = []
gFrequency = []
gDateObs = []
gImageHist = []
gEll_params = []

def generateEllipse(x0, y0, A, B, theta):
    N = 100
    t = NP.linspace(0,2*NP.pi,100)
    xy = NP.zeros((2,N))
    xy[0,:] = x0 + A*NP.cos(theta)*NP.cos(t) - B*NP.sin(theta)*NP.sin(t)
    xy[1,:] = y0 + A*NP.sin(theta)*NP.cos(t) + B*NP.cos(theta)*NP.sin(t)
    return xy
    
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

def prepareCleanImages(fileName, psfName, qSunTb, saveCleanFitsPath, srhRawFits, frequencyIndex, scan):
    fd = fits.open(fileName)
    
    pAngle = 0
    try:
        client = Client('http://ephemeris.rao.istp.ac.ru/?wsdl')
        result = client.service.Ephemeride('SSRT','sun',fd[0].header['DATE-OBS'])
        pAngle = NP.deg2rad(float(result[0]['PAngle']))
    except:
        pass
    
    rcpImage = fd[0].data[0,0]
    lcpImage = fd[0].data[1,0]

    srh_x_size = fd[0].header['NAXIS1']
    srh_y_size = fd[0].header['NAXIS2']
    
    arcsecPerPixel = NP.abs(fd[0].header['CDELT1']*3600)
    resultArcsecPerPixel = 4.911
    resultScale = arcsecPerPixel / resultArcsecPerPixel
    
    fd_psf = fits.open(psfName)
    contours = (skimage.measure.find_contours(fd_psf[0].data[0,0], 0.5))[0]
    con = NP.zeros_like(contours)
    con[:,1] = contours[:,0]
    con[:,0] = contours[:,1]
    gPSF.append(fd_psf[0].data[0,0])
    sunEll = skimage.measure.EllipseModel()
    sunEll.estimate(con)
    sunEll.params[2] *= resultScale
    sunEll.params[3] *= resultScale
    sunEll.params[4] -= pAngle
    gEll_params.append(NP.array(sunEll.params))

    global gXSize
    gXSize = srh_x_size
    global gArcsecPerPixel
    gArcsecPerPixel = resultArcsecPerPixel
    
#polarization signum
    iImage = rcpImage + lcpImage
    vImage = rcpImage - lcpImage
    
    iImageHist = NP.histogram(iImage,bins=1000)
    for promFactor in range(5):
        maxInds = scipy.signal.find_peaks(iImageHist[0],prominence = 2000//(promFactor+1))
        if maxInds[0].shape[0] > 1:
            break
    global gImageHist
    gImageHist.append(iImageHist)
    if maxInds[0].shape[0] > 1:
        iImage -= iImageHist[1][maxInds[0][0]]
        iImage /= (iImageHist[1][maxInds[0][1]] - iImageHist[1][maxInds[0][0]])
        vImage /= (iImageHist[1][maxInds[0][1]] - iImageHist[1][maxInds[0][0]])
    else:
        sunLevel = iImage[256-64:256+64,256-64:256+64].mean()
        skyLevel = 0
        iImage -= skyLevel
        iImage /= (sunLevel - skyLevel)
        vImage /= (sunLevel - skyLevel)

    iImage *= qSunTb
    vImage *= qSunTb

    rcpImage = iImage + vImage
    lcpImage = iImage - vImage

    imageLength = rcpImage.shape[0]**2
    slope, interc, A, B, C = linregress(NP.clip(rcpImage.reshape(imageLength),0,5e4), NP.clip(lcpImage.reshape(imageLength),0,5e4))
    iImage = 0.5*(slope*rcpImage + interc + lcpImage)
    vImage = 0.5*(slope*rcpImage + interc - lcpImage)

    scale = AffineTransform(scale=(resultScale,resultScale))
    shift = AffineTransform(translation=(-srh_y_size/2,-srh_y_size/2))
    rotate = AffineTransform(rotation = -pAngle)
    back_shift = AffineTransform(translation=(srh_y_size/2,srh_y_size/2))
    
    iImage = warp(iImage,(shift + (rotate + back_shift)).inverse)
    vImage = warp(vImage,(shift + (rotate + back_shift)).inverse)
    iImage = warp(iImage,(shift + (scale + back_shift)).inverse)
    vImage = warp(vImage,(shift + (scale + back_shift)).inverse)
    
    global gIImage
    global gVImage
    
    gIImage.append(iImage)
    gVImage.append(vImage)
    gFrequency.append(fd[0].header['CRVAL3']/1e6 + 0.5)
    
    gDateObs.append(fd[0].header['DATE-OBS'])
    
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
    pHeader['PSF_ELLA']     = NP.ceil(sunEll.params[2] * resultArcsecPerPixel*10)/10 # PSF ellipse A arcsec
    pHeader['PSF_ELLB']     = NP.ceil(sunEll.params[3] * resultArcsecPerPixel*10)/10 # PSF ellipse B arcsec
    pHeader['PSF_ELLT']     = sunEll.params[4] # PSF ellipse theta rad

    saveFitsIhdu = fits.PrimaryHDU(header=pHeader, data=iImage.astype('float32'))
    saveFitsIpath = saveCleanFitsPath + 'srh_I_%s_%04d_%d.fit'%(fd[0].header['DATE-OBS'].split('.')[0],fd[0].header['CRVAL3']*1e-6 + .5,scan)
    ewLcpPhaseColumn = fits.Column(name='ewLcpPhase', format='D', array = srhRawFits.ewAntPhaLcp[frequencyIndex,:] + srhRawFits.ewLcpPhaseCorrection[frequencyIndex,:])
    ewRcpPhaseColumn = fits.Column(name='ewRcpPhase', format='D', array = srhRawFits.ewAntPhaRcp[frequencyIndex,:] + srhRawFits.ewRcpPhaseCorrection[frequencyIndex,:])
    nLcpPhaseColumn = fits.Column(name='nLcpPhase',   format='D', array = srhRawFits.nAntPhaLcp[frequencyIndex,:] + srhRawFits.nLcpPhaseCorrection[frequencyIndex,:])
    nRcpPhaseColumn = fits.Column(name='nRcpPhase',   format='D', array = srhRawFits.nAntPhaRcp[frequencyIndex,:] + srhRawFits.nRcpPhaseCorrection[frequencyIndex,:])
    saveFitsIExtHdu = fits.BinTableHDU.from_columns([ewLcpPhaseColumn, ewRcpPhaseColumn, nLcpPhaseColumn, nRcpPhaseColumn])
    hduList = fits.HDUList([saveFitsIhdu, saveFitsIExtHdu])
    hduList.writeto(saveFitsIpath)
    
    saveFitsVhdu = fits.PrimaryHDU(header=pHeader, data=vImage.astype('float32'))
    saveFitsVpath = saveCleanFitsPath + 'srh_V_%s_%04d.fit'%(fd[0].header['DATE-OBS'].split('.')[0],fd[0].header['CRVAL3']*1e-6 + .5)
    hduList = fits.HDUList(saveFitsVhdu)
    hduList.writeto(saveFitsVpath)
    
    fd.close()
    
#------------------------------------------------------------------------------
currentDate = datetime.datetime.now().date().strftime("%Y%m%d")
parser = OptionParser()
parser.add_option("-f", "--file", dest="fitPath", default = 'SRH36_temp_20210914_1/')
parser.add_option("-t", "--treshold", dest="cleanTresh", default = '0.1mJy')
parser.add_option("-s", "--scan", dest="scan", default = '0')

(clean_options, clean_args) = parser.parse_args()
fitPath = clean_options.fitPath
cleanTresh = clean_options.cleanTresh
scan = clean_options.scan

fitNames =  findFits(fitPath,'*.fit')
fitNames.sort()
fileName = fitNames[3]

ZirinQSunTb = ZirinTb()

file = SrhFitsFile(fileName, 1025)
file.useNonlinearApproach = True
file.getHourAngle(0)
frequencyNumber = 3
if (not os.path.exists('cleanMaps/' + currentDate)):
    os.system('mkdir ' + 'cleanMaps/' + currentDate)
    
frequency = 5
file.solarPhase(frequency)

file.updateAntennaPhase(frequency, baselinesNumber = 5)
file.setFrequencyChannel(frequency)
file.vis2uv(0, average = 20)
file.centerDisk()

saveName = 'MS/' + fileName.split('/')[-1].split('.')[0] + '.ms'
ms2Table = srh36MS2.SrhMs2(saveName)
ms2Table.createMS(file, frequencyChannel = [int(frequency)], phaseCorrect = True, amplitudeCorrect = True)

flags_ew_lcp = NP.where(file.ewAntAmpLcp[frequency] == 1e6)[0] + 1
flags_ew_rcp = NP.where(file.ewAntAmpRcp[frequency] == 1e6)[0] + 1
flags_ew = NP.unique(NP.append(flags_ew_lcp, flags_ew_rcp))
flags_n_lcp = NP.where(file.nAntAmpLcp[frequency] == 1e6)[0]+98
flags_n_rcp = NP.where(file.nAntAmpRcp[frequency] == 1e6)[0]+98
flags_n = NP.unique(NP.append(flags_n_lcp, flags_n_rcp))
flags = ','.join(map(str, NP.append(flags_ew, flags_n)))

currentDiskTb = ZirinQSunTb.getTbAtFrequncy(file.freqList[frequency]*1e-6)*1e3
for s in range(20):
    scan = str(s)
    casaSaveName = saveName.split('.')[0] + '_' + scan + '.ms'
    command = 'casa -c casaclean.py \'' + casaSaveName + '\'  \'' + flags + '\' \'' +  cleanTresh + '\' \'' + scan  + '\''
    os.system(command)
    prepareCleanImages(casaSaveName.split('.')[0] + '_clean_image.fit',casaSaveName.split('.')[0] + '_clean_psf.fit', currentDiskTb, 'cleanMaps/' + currentDate + '/',file, frequency, s)

clnjunks = ['.ms', '.ms.flagversions', '*clean*.fit']
for clnjunk in clnjunks:
    if os.path.exists(saveName.split('.')[0] + clnjunk):
        os.system('rm -rf ' + saveName.split('.')[0] + clnjunk)
os.system('rm -rf MS/*.ms')
os.system('rm -rf MS/*.ms.flagversions')
os.system('rm -rf MS/*clean*.fit')
os.system('rm -rf *casa*.log')
