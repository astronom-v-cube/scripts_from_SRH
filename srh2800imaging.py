#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 22 01:57:10 2021

@author: svlesovoi
"""

from srhFitsFile36 import SrhFitsFile
import pylab as PL
import srh36MS2
import numpy as NP
import os
from astropy.io import fits
import skimage
from skimage import measure
from skimage.transform import warp, AffineTransform
from zeep import Client
import scipy.signal
import scipy.constants
from matplotlib.patches import Circle

gXSize = 1
gArcsecPerPixel = 1
gLcpImage = NP.zeros((100,100))
gRcpImage = NP.zeros((100,100))

def arcmin_format(xy, pos):
  return '%2d' % ((xy - gXSize/2) * gArcsecPerPixel / 60);

def showCleanImages(fileName, psfName):
    filter_dL = 21
    filter_arg = NP.linspace(-1,1,filter_dL)
    filterX, filterY = NP.meshgrid(filter_arg, filter_arg)
    imageFilter = NP.exp(-((filterX/5.7)**2 + (filterY/5.7)**2))
    
    iImageMaxX = []
    iImageMaxY = []
    flux107 = []
    Tb_2800 = 27100
    
    psf_fd = fits.open(psfName)
    psfData = psf_fd[0].data[0][0]
    
    fd = fits.open(fileName)
    
    pAngle = 0
    try:
        client = Client('http://ephemeris.rao.istp.ac.ru/?wsdl')
        result = client.service.Ephemeride('SSRT','sun',fd[0].header['DATE-OBS'])
        pAngle = NP.deg2rad(float(result[0]['PAngle']))
    except:
        pass
    
    lcpImage = fd[0].data[0,0]
    rcpImage = fd[0].data[1,0]
    global gLcpImage
    global gRcpImage
    gLcpImage = lcpImage.copy()
    gRcpImage = rcpImage.copy()
    
    sunMask = NP.zeros_like(lcpImage)
    outRadius = 400
    sunX0 = lcpImage.shape[0]/2
    sunY0 = lcpImage.shape[0]/2
    for i in range(lcpImage.shape[0]):
        for j in range(lcpImage.shape[0]):
            curRadius = NP.sqrt((i - sunX0)**2 + (j - sunY0)**2)
            if curRadius < outRadius:
                sunMask[i,j] = 1.
    
    srh_x_size = fd[0].header['NAXIS1']
    srh_y_size = fd[0].header['NAXIS2']
    
    x0 = int(srh_x_size//2)
    y0 = int(srh_y_size//2)
    arcsecPerPixel = fd[0].header['CDELT1']*3600
    resultArcsecPerPixel = 4.911/2
    resultScale = arcsecPerPixel / resultArcsecPerPixel
    
    gXSize = srh_x_size
    gArcsecPerPixel = resultArcsecPerPixel
    
    sunCir = skimage.measure.CircleModel()
    
    stdV = NP.zeros((10,10))
    dL = 1.2
    for i in range(10):
        for j in range(10):
            shift = AffineTransform(translation=((i-5)*dL,(j-5)*dL))
            stdV[i,j]= ((warp(rcpImage,shift.inverse) - lcpImage)**2).mean()
    
    rlShift = NP.where(stdV==stdV.min())
    rlShift = (rlShift[0]%lcpImage.shape[0], rlShift[1]%lcpImage.shape[1])
    rlShift = [(rlShift[0][0]-5)*dL,(rlShift[1][0]-5)*dL]
    shift = AffineTransform(translation=rlShift)
    rcpImage = warp(rcpImage,shift.inverse)
    iImage = rcpImage + lcpImage
    vImage = rcpImage - lcpImage
    
    iImageHist = NP.histogram(iImage,bins=1000)
    maxInds = scipy.signal.find_peaks(iImageHist[0],prominence=2000)
    iImage -= iImageHist[1][maxInds[0][0]]
    iImage /= (iImageHist[1][maxInds[0][1]] - iImageHist[1][maxInds[0][0]])
    iImage *= Tb_2800
    vImage /= (iImageHist[1][maxInds[0][1]] - iImageHist[1][maxInds[0][0]])
    vImage *= Tb_2800
    
    scale = AffineTransform(scale=(-resultScale,-resultScale))
    shift = AffineTransform(translation=(-srh_y_size/2,-srh_y_size/2))
    rotate = AffineTransform(rotation = -pAngle)
    back_shift = AffineTransform(translation=(srh_y_size/2,srh_y_size/2))
    
    iImage = warp(iImage,(shift + (rotate + back_shift)).inverse)
    vImage = warp(vImage,(shift + (rotate + back_shift)).inverse)
    iImage = warp(iImage,(shift + (scale + back_shift)).inverse)
    vImage = warp(vImage,(shift + (scale + back_shift)).inverse)
    
    conImage = scipy.signal.fftconvolve(iImage,imageFilter) / filter_dL**2
    conImage = conImage[filter_dL//2:filter_dL//2 + srh_x_size,filter_dL//2:filter_dL//2 + srh_y_size]
    contours = measure.find_contours(conImage, 17000)
    contourLength = []
    for n, contour in enumerate(contours):
        contourLength.append(len(contour))
    contourLength = NP.array(contourLength)
    maxInd = NP.argmax(contourLength)
    sunCir.estimate(contours[maxInd])
    i_cx0, i_cy0, i_cR = sunCir.params

    iShift = AffineTransform(translation=(y0 + .5 - i_cy0,x0 + .5 - i_cx0))
    iImage = warp(iImage,iShift.inverse)
    vImage = warp(vImage,iShift.inverse)

    fig, pl = PL.subplots(figsize=(8,8))
    fig.tight_layout()
    cir_patch = Circle((y0,x0), i_cR, edgecolor='red', facecolor='none', linewidth=1.3)
    pl.xaxis.set_major_locator(PL.MultipleLocator(256));
    pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
    pl.xaxis.set_minor_locator(PL.MultipleLocator(64));
    pl.yaxis.set_major_locator(PL.MultipleLocator(256));
    pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
    pl.yaxis.set_minor_locator(PL.MultipleLocator(64));
    pl.text(10,1000,'SRH V, %d MHz, %s' % (fd[0].header['RESTFRQ']*1e-6,fd[0].header['DATE-OBS'].split('.')[0]),color='white',fontsize='14')
    pl.set_xlabel('arcmin')
    pl.set_ylabel('arcmin')
    pl.imshow(vImage,cmap='Greys',vmin=-0.3e6,vmax=0.3e6,origin='lower')
    pl.add_patch(cir_patch)
    fig.savefig('srh_V_%s.png'%fd[0].header['DATE-OBS'].split('.')[0])
    
    maxI = NP.where(iImage==iImage.max())
    iImageMaxX.append(maxI[0][0]%iImage.shape[0])
    iImageMaxY.append(maxI[1][0]%iImage.shape[1])
    flux107.append(2*(iImage*sunMask).sum()*scipy.constants.k/(0.107**2)*(NP.deg2rad(fd[0].header['CDELT1'])**2)/1e-22)

    fig, pl = PL.subplots(figsize=(8,8))
    fig.tight_layout()
    cir_patch = Circle((y0,x0), i_cR, edgecolor='red', facecolor='none', linewidth=.3)
    pl.xaxis.set_major_locator(PL.MultipleLocator(256));
    pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
    pl.xaxis.set_minor_locator(PL.MultipleLocator(64));
    pl.yaxis.set_major_locator(PL.MultipleLocator(256));
    pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
    pl.yaxis.set_minor_locator(PL.MultipleLocator(64));
    pl.text(10,1000,'SRH I, %d MHz, %s, S=%.1f s.f.u.' % (fd[0].header['RESTFRQ']*1e-6,fd[0].header['DATE-OBS'].split('.')[0], flux107[-1]),color='white',fontsize='14')
    pl.set_xlabel('arcmin')
    pl.set_ylabel('arcmin')
#        pl.add_patch(cir_patch)
    pl.imshow(iImage,cmap='hot',vmin=0.,vmax=3e5,origin='lower')
    pl.contour(iImage,cmap='hot',levels=[3e5,5e5,7e5,1e6,1.3e6,1.5e6],linewidths=0.5)
#    pl.contour(vImage,cmap='seismic',levels=[-20e4,-10e4,-5e4,-2e4,-1e4,1e4,2e4,5e4,10e4,20e4],linewidths=0.5)
#    pl.contour(NP.roll(NP.roll(psfData,-srh_x_size//2-100,axis=0),srh_y_size//2+100,axis=1),levels=[0.5],cmap='gray_r')
    pl.contour(NP.roll(NP.roll(psfData,srh_x_size//2+100,axis=0),srh_y_size//2+100,axis=1),levels=[0.5],cmap='gray_r')

    fig.savefig('srh_I_%s.png'%fd[0].header['DATE-OBS'].split('.')[0])

    fd.close()
    
#fileName = 'SRH36_temp_20210522_1/srh_20210522T024143.fit'
#fileName = 'SRH36_temp_20210522_2/srh_20210522T025243.fit'
#fileName = 'SRH36_temp_20210522_3/srh_20210522T033127.fit'
#fileName = 'SRH36_temp_20210521_1/srh_20210521T032752.fit'
#fileName = 'SRH36_temp_20210522_4/srh_20210522T072609.fit'
#fileName = 'SRH36_temp_20210520_1/srh_20210520T032000.fit'
#fileName = 'SRH36_temp_20210522_6/srh_20210522T063533.fit'
#fileName = 'SRH36_temp_20210523_1/srh_20210523T024721.fit'
fileName = 'SRH36_temp_20210523_2/srh_20210523T090705.fit'

file = SrhFitsFile(fileName, 2049)
freq = 1
file.baselines = 3
file.useNonlinearApproach = True
file.getHourAngle(0)
file.solarPhase(freq)
file.updateAntennaPhase(freq, baselinesNumber = 5)
file.setFrequencyChannel(freq)

saveName = fileName.split('.')[0] + '.ms'
ms2Table = srh36MS2.SrhMs2(saveName)
ms2Table.createMS(file, frequencyChannel = [freq], phaseCorrect = True, amplitudeCorrect = True)

flags_ew_lcp = NP.where(file.ewAntAmpLcp[freq] == 1e6)[0] + 1
flags_ew_rcp = NP.where(file.ewAntAmpRcp[freq] == 1e6)[0] + 1
flags_ew = NP.unique(NP.append(flags_ew_lcp, flags_ew_rcp))
flags_n_lcp = NP.where(file.nAntAmpLcp[freq] == 1e6)[0]+98
flags_n_rcp = NP.where(file.nAntAmpRcp[freq] == 1e6)[0]+98
flags_n = NP.unique(NP.append(flags_n_lcp, flags_n_rcp))
flags = ','.join(map(str, NP.append(flags_ew, flags_n)))

command = 'casa -c casaclean.py \'' + saveName + '\'  \'' + flags + '\''
os.system(command)

fileName = saveName.split('.')[0] + '_clean_image.fit'
psfName = saveName.split('.')[0] + '_clean_psf.fit'

showCleanImages(fileName, psfName)

    