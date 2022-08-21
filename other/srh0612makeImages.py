#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 17:27:13 2021

@author: maria
"""


from srhFitsFile612 import SrhFitsFile
# import srh36MS2
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
import ftplib;
import datetime as DT;

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

def prepareCleanImages(fileName, psfName, qSunTb, saveCleanFitsPath, srhRawFits, frequencyIndex):
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
        maxInds = scipy.signal.find_peaks(iImageHist[0],prominence = 1000//(promFactor+1))
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
    saveFitsIpath = saveCleanFitsPath + 'srh_I_%s_%04d.fit'%(fd[0].header['DATE-OBS'].split('.')[0],fd[0].header['CRVAL3']*1e-6 + .5)
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
    
#    T0 = 2e5
#    fig, pl = PL.subplots(figsize=(6,6))
#    fig.tight_layout()
#    pl.clear()
#    pl.xaxis.set_major_locator(PL.MultipleLocator(128));
#    pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
#    pl.xaxis.set_minor_locator(PL.MultipleLocator(32));
#    pl.yaxis.set_major_locator(PL.MultipleLocator(128));
#    pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
#    pl.yaxis.set_minor_locator(PL.MultipleLocator(32));
#    pl.imshow(iImage,cmap='hot',vmax=T0,vmin=0,origin='lower')
#    levels = NP.linspace(T0,NP.max(iImage),10)
#    pl.contour(iImage,cmap='hot',levels=levels,linewidths=0.5,origin='lower')
#    pl.set_title('SRH I %s,  %d MHz'%(fd[0].header['DATE-OBS'].split('.')[0], fd[0].header['CRVAL3']*1e-6 + .5))
#    pl.set_xlabel('arcmin')
#    pl.set_ylabel('arcmin')
#    fig.savefig(saveFitsIpath.split('.')[0] + '.png')
#
#
#    T0 = 1e5
#    fig, pl = PL.subplots(figsize=(6,6))
#    fig.tight_layout()
#    pl.clear()
#    pl.xaxis.set_major_locator(PL.MultipleLocator(128));
#    pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
#    pl.xaxis.set_minor_locator(PL.MultipleLocator(32));
#    pl.yaxis.set_major_locator(PL.MultipleLocator(128));
#    pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
#    pl.yaxis.set_minor_locator(PL.MultipleLocator(32));
#    pl.imshow(vImage,cmap='gray',vmax=T0,vmin=-T0,origin='lower')
#    levels = NP.linspace(T0,NP.max(vImage),3)
#    pl.contour(vImage,cmap='gray',levels=levels,linewidths=0.5,origin='lower')
#    levels = NP.linspace(NP.min(vImage),-T0,3)
#    pl.contour(vImage,cmap='gray',levels=levels,linewidths=0.5,origin='lower')
#    pl.set_title('SRH V %s,  %d MHz'%(fd[0].header['DATE-OBS'].split('.')[0], fd[0].header['CRVAL3']*1e-6 + .5))
#    pl.set_xlabel('arcmin')
#    pl.set_ylabel('arcmin')
#    fig.savefig(saveFitsVpath.split('.')[0] + '.png')
#------------------------------------------------------------------------------
currentDate = datetime.datetime.now().date().strftime("%Y%m%d")
parser = OptionParser()
parser.add_option("-f", "--file", dest="fitPath", default = 'SRH0612/' + currentDate + '/')
parser.add_option("-t", "--treshold", dest="cleanTresh", default = '0.1mJy')
parser.add_option("-s", "--scan", dest="scan", default = '0~19')

(clean_options, clean_args) = parser.parse_args()
fitPath = clean_options.fitPath
cleanTresh = clean_options.cleanTresh
scan = clean_options.scan

fitNames =  findFits(fitPath,'*.fit')
fitNames.sort()
fileName = fitNames[-1]

ZirinQSunTb = ZirinTb()

file = SrhFitsFile(fileName, 1024)
file.useNonlinearApproach = True
file.getHourAngle(0)
freqList = [3,9,15]
frequencyNumber = len(freqList)
gIImage = NP.zeros((frequencyNumber, 512, 512))
gVImage = NP.zeros((frequencyNumber, 512, 512))

# if (not os.path.exists('cleanMaps/' + currentDate)):
#     os.system('mkdir ' + 'cleanMaps/' + currentDate)
    
for frequency in range(frequencyNumber):
    file.calibrate(freqList[frequency], average = 20)
    file.vis2uv(0, average = 20)
    file.uv2lmImage()
    file.lm2Heliocentric()
    gIImage[frequency] = (file.lcp + file.rcp)/2
    gVImage[frequency] = (file.lcp - file.rcp)/2
    gFrequency.append(file.freqList[freqList[frequency]]/1e3)

fig, pl = PL.subplots(nrows=frequencyNumber,ncols=2,figsize=(5.6,7.9))
fig.subplots_adjust(hspace=0.001,wspace=0.001,top=0.95,bottom=0.01,left=0.01,right=0.99)
TI = 2e3
TV = 1e3
imageDate = file.hduList[0].header['DATE-OBS'].replace('-','')
imageTimeShow = file.hduList[0].header['TIME-OBS'][0:5]
imageTime = file.hduList[0].header['TIME-OBS'].replace(':','')[0:4]
combinedImagePath = 'srh_%s_%s.png'%(imageDate,imageTime)
fig.suptitle('%s, %s UT'%(imageDate, imageTimeShow))
for row in range(frequencyNumber):
    pl[row,0].axis('off')
    pl[row,0].imshow(gIImage[row],cmap='hot',vmax=TI,vmin=0,origin='lower')
    levels = NP.linspace(TI,NP.max(gIImage[row]),3)
    pl[row,0].contour(gIImage[row],cmap='hot',levels=levels,linewidths=0.5,origin='lower')
    pl[row,0].text(10,10,'I',color='white')
    # flux = 2*gIImage[row][256-200:256+200,256-200:256+200].sum()*scipy.constants.k/(scipy.constants.c/gFrequency[row]*1e-6)**2 * NP.deg2rad(gArcsecPerPixel/3600)**2 / 1e-22
    # pl[row,0].text(10,gXSize - 30,'F = %.1f sfu'%flux,color='white')

    # ellXY = generateEllipse(50,50,gEll_params[row][2],gEll_params[row][3],gEll_params[row][4])
    # pl[row,0].plot(ellXY[0,:],ellXY[1,:],color='white')

    pl[row,1].axis('off')
    pl[row,1].imshow(gVImage[row],cmap='gray',vmax=TV,vmin=-TV,origin='lower')
    levels = NP.linspace(TV,NP.max(gVImage[row]),3)
    pl[row,1].contour(gVImage[row],cmap='gray',levels=levels,linewidths=0.5,origin='lower')
    levels = NP.linspace(NP.min(gVImage[row]),-TV,3)
    pl[row,1].contour(gVImage[row],cmap='gray',levels=levels,linewidths=0.5,origin='lower')
    pl[row,1].text(10,512 - 30,'%d MHz'%gFrequency[row],color='white')
    pl[row,1].text(10,10,'V',color='white')

fig.savefig(combinedImagePath)
fdSrhImageName = open('srhImageName.txt','w');
fdSrhImageName.write(combinedImagePath);
fdSrhImageName.close();

#fig, pl = PL.subplots()
#pl.imshow(gPSF[2],origin='lower')
#contours = (skimage.measure.find_contours(gPSF[2], 0.3))[0]
#con = NP.zeros_like(contours)
#con[:,1] = contours[:,0]
#con[:,0] = contours[:,1]
#pl.plot(con[:,0],con[:,1],color='red')
#sunEll = skimage.measure.EllipseModel()
#sunEll.estimate(con)
#ellXY = generateEllipse(sunEll.params[0],sunEll.params[1],sunEll.params[2],sunEll.params[3],sunEll.params[4])
#pl.plot(ellXY[0,:],ellXY[1,:])

dateName = DT.datetime.now().strftime("%Y%m%d")
fd = ftplib.FTP('192.168.0.1','sergey','jill21');

fi = open('srhImageName.txt','rb');
fd.storbinary('STOR /Public/sergey/corrPlots/srhImageName.txt',fi);
fi.close();

fi = open('srhImageName.txt');
cleanImageName = fi.readline()
fCleanImage = open(cleanImageName,'rb')
fd.storbinary('STOR /Public/sergey/corrPlots/' + cleanImageName,fCleanImage);
fCleanImage.close()
fi.close()
fd.quit();
