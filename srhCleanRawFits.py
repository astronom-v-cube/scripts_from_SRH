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
import sunpy.coordinates

gXSize = 1
gArcsecPerPixel = 1
gLcpImage = NP.zeros((100,100))
gRcpImage = NP.zeros((100,100))
gSaveIname = ''
gSaveVname = ''

cdict = {'red': ((0.0, 0.0, 0.0),
                 (0.2, 0.0, 0.0),
                 (0.4, 0.2, 0.2),
                 (0.6, 0.0, 0.0),
                 (0.8, 1.0, 1.0),
                 (0.9, 1.0, 1.0),
                 (1.0, 1.0, 1.0)),
        'green':((0.0, 0.0, 0.0),
                 (0.2, 0.0, 0.0),
                 (0.4, 1.0, 1.0),
                 (0.6, 1.0, 1.0),
                 (0.8, 1.0, 1.0),
                 (0.9, 0.0, 0.0),
                 (1.0, 1.0, 1.0)),
        'blue': ((0.0, 0.0, 0.0),
                 (0.2, 1.0, 1.0),
                 (0.4, 1.0, 1.0),
                 (0.6, 0.0, 0.0),
                 (0.8, 0.0, 0.0),
                 (0.9, 0.0, 0.0),
                 (1.0, 1.0, 1.0))}

my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)

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

def showCleanImages(fileName, psfName, qSunLevel, minT, maxT):
    filter_dL = 11
    filter_arg = NP.linspace(-1,1,filter_dL)
    filterX, filterY = NP.meshgrid(filter_arg, filter_arg)
    imageFilter = NP.exp(-((filterX/5.7)**2 + (filterY/5.7)**2))
    
    flux107 = []
    Tb_2800 = 27100
    
    psf_fd = fits.open(psfName)
    psfData = psf_fd[0].data[0][0]
    
    fd = fits.open(fileName)
    
    pAngle = NP.deg2rad(sunpy.coordinates.sun.P(fd[0].header['DATE-OBS']).to_value())
    
    lcpImage = fd[0].data[0,0]
    rcpImage = fd[0].data[1,0]
    global gLcpImage
    global gRcpImage
#    gLcpImage = lcpImage.copy()
#    gRcpImage = rcpImage.copy()
    
    sunMask = NP.zeros_like(lcpImage)
    outRadius = 200
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
    resultArcsecPerPixel = 4.911
    resultScale = arcsecPerPixel / resultArcsecPerPixel
    
    global gXSize
    gXSize = srh_x_size
    global gArcsecPerPixel
    gArcsecPerPixel = resultArcsecPerPixel
    
    sunCir = skimage.measure.CircleModel()
    sunEll = skimage.measure.EllipseModel()
    
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
#    rcpImage = warp(rcpImage,shift.inverse)
    iImage = rcpImage + lcpImage
    vImage = rcpImage - lcpImage
    
    iImageHist = NP.histogram(iImage,bins=1000)
    prominList = NP.linspace(2000,100,10)
    for prominence in prominList:
        maxInds = scipy.signal.find_peaks(iImageHist[0],prominence=2000)
        if maxInds[0].shape[0] > 1:
            iImage -= iImageHist[1][maxInds[0][0]]
            iImage /= (iImageHist[1][maxInds[0][1]] - iImageHist[1][maxInds[0][0]])
            iImage *= Tb_2800
            vImage /= (iImageHist[1][maxInds[0][1]] - iImageHist[1][maxInds[0][0]])
            vImage *= Tb_2800
            break;
#    else:
#        iMean = iImage.mean()
#        iImage /= iMean
#        iImage *= Tb_2800
#        vImage /= iMean
#        vImage *= Tb_2800
        
    gLcpImage = iImage.copy()
    gRcpImage = vImage.copy()

    scale = AffineTransform(scale=(-resultScale,-resultScale))
    shift = AffineTransform(translation=(-srh_y_size/2,-srh_y_size/2))
    rotate = AffineTransform(rotation = -pAngle)
    back_shift = AffineTransform(translation=(srh_y_size/2,srh_y_size/2))
    
    iImage = warp(iImage,(shift + (rotate + back_shift)).inverse)
    vImage = warp(vImage,(shift + (rotate + back_shift)).inverse)
    iImage = warp(iImage,(shift + (scale + back_shift)).inverse)
    vImage = warp(vImage,(shift + (scale + back_shift)).inverse)
    
    psfData = warp(psfData,(shift + (rotate + back_shift)).inverse)

#    conImage = scipy.signal.fftconvolve(iImage,imageFilter) / filter_dL**2
#    conImage = conImage[filter_dL//2:filter_dL//2 + srh_x_size,filter_dL//2:filter_dL//2 + srh_y_size]
#    contours = measure.find_contours(conImage, qSunLevel)
#    contourLength = []
#    for n, contour in enumerate(contours):
#        contourLength.append(len(contour))
#    contourLength = NP.array(contourLength)
#    maxInd = NP.argmax(contourLength)
#    sunCir.estimate(contours[maxInd])
#    i_cx0, i_cy0, i_cR = sunCir.params
#    sunEll.estimate(contours[maxInd])
#    i_ex0, i_ey0, i_eA, i_eB, i_eT = sunEll.params

#    iShift = AffineTransform(translation=(y0 + .5 - i_cy0,x0 + .5 - i_cx0))
#    iImage = warp(iImage,iShift.inverse)
#    vImage = warp(vImage,iShift.inverse)

#    iShift = AffineTransform(translation=(y0 + .5 - i_ey0,x0 + .5 - i_ex0))
#   iImage = warp(iImage,iShift.inverse)
#    vImage = warp(vImage,iShift.inverse)

#    cir_patch = Circle((y0,x0), i_cR, edgecolor='red', facecolor='none', linewidth=1.3)
#    cir_patch = Circle((i_cy0,i_cx0), i_cR, edgecolor='red', facecolor='none', linewidth=1.3)
#    ell_patch = Ellipse((i_ey0,i_ex0), 2*i_eA, 2*i_eB, -i_eT/180*NP.pi, edgecolor='green', facecolor='none', linewidth=1.3)
#    ell_patch = Ellipse((y0,x0), 2*i_eA, 2*i_eB, -i_eT/180*NP.pi, edgecolor='green', facecolor='none', linewidth=1.3)
    
    saveFigPath = '/home/svlesovoi/Pictures/SRH36_pngs/'
    saveFigDate = fd[0].header['DATE-OBS'].split('T')[0].replace('-','')
    fig, pl = PL.subplots(figsize=(6,6))
    fig.tight_layout()
    pl.xaxis.set_major_locator(PL.MultipleLocator(256));
    pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
    pl.xaxis.set_minor_locator(PL.MultipleLocator(64));
    pl.yaxis.set_major_locator(PL.MultipleLocator(256));
    pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
    pl.yaxis.set_minor_locator(PL.MultipleLocator(64));
    pl.text(10,500,'SRH V, %4d MHz, %s' % (fd[0].header['RESTFRQ']*1e-6 + .5,fd[0].header['DATE-OBS'].split('.')[0]),color='white',fontsize='10')
    pl.set_xlabel('arcmin')
    pl.set_ylabel('arcmin')
    pl.imshow(vImage,cmap='Greys',vmin=-0.3e6,vmax=0.3e6,origin='lower')
#    pl.add_patch(cir_patch)
    global gSaveVname
    gSaveVname = 'srh_V_%s_%4d.png'%(fd[0].header['DATE-OBS'].split('.')[0],fd[0].header['RESTFRQ']*1e-6 + .5)
    fig.savefig(saveFigPath + saveFigDate + '/V/' + gSaveVname)
    
    flux107.append(2*(iImage*sunMask).sum()*scipy.constants.k/(0.107**2)*(NP.deg2rad(fd[0].header['CDELT1'])**2)/1e-22)

    PL.close()
    fig, pl = PL.subplots(figsize=(6,6))
    fig.tight_layout()
    pl.xaxis.set_major_locator(PL.MultipleLocator(128));
    pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
    pl.xaxis.set_minor_locator(PL.MultipleLocator(32));
    pl.yaxis.set_major_locator(PL.MultipleLocator(128));
    pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
    pl.yaxis.set_minor_locator(PL.MultipleLocator(32));
    pl.text(10,500,'SRH I, %4d MHz, %s, S=%.1f s.f.u.' % (fd[0].header['RESTFRQ']*1e-6 + .5,fd[0].header['DATE-OBS'].split('.')[0], flux107[-1]),color='white',fontsize='10')
    pl.text(10,480,'Tmax %.0f' % (iImage.max()),color='white',fontsize='10')
    pl.set_xlabel('arcmin')
    pl.set_ylabel('arcmin')
    pl.imshow(iImage,cmap=my_cmap,vmin=minT,vmax=maxT,origin='lower')
    pl.contour(iImage,cmap='hot',levels=[2e5,4e5,7e5,1e6,1.4e6,2.0e6],linewidths=0.5)
    pl.contour(NP.roll(NP.roll(psfData,srh_x_size//2+50,axis=0),srh_y_size//2+50,axis=1),levels=[0.5],cmap='gray_r')

    global gSaveIname
    gSaveIname = 'srh_I_%s_%4d.png'%(fd[0].header['DATE-OBS'].split('.')[0],fd[0].header['RESTFRQ']*1e-6 + .5)
    fig.savefig(saveFigPath + saveFigDate + '/I/' + gSaveIname)
    PL.close()
#
#    store = ftplib.FTP('84.237.21.39','ftpwriter','cn5HeuA4u')
#    fi = open(gSaveIname,'rb')
#    store.storbinary('STOR /SRH/SRH0306/cleanMaps/' + gSaveIname, fi)
#    fi.close()
#    store.close()
    
    fd.close()
    
#------------------------------------------------------------------------------
fitPath = sys.argv[1]
freq = int(sys.argv[2])
cleanTresh = sys.argv[3]#'.13mJy'
scan = sys.argv[4]#'.13mJy'

fitNames =  findFits(fitPath,'*.fit')
fitNames.sort()

for fileName in fitNames:
    file = SrhFitsFile(fileName, 1025)
    file.useNonlinearApproach = True
    file.getHourAngle(0)
    file.solarPhase(freq)
    file.updateAntennaPhase(freq, baselinesNumber = 5)
    file.setFrequencyChannel(freq)
    file.vis2uv(0, average = 20)
    file.centerDisk()
    
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
    
    command = 'casa -c casaclean.py \'' + saveName + '\'  \'' + flags + '\' \'' +  cleanTresh + '\' \'' + scan  + '\''
    os.system(command)
    
    cleanName = saveName.split('.')[0] + '_clean_image.fit'
    psfName = saveName.split('.')[0] + '_clean_psf.fit'
    
    showCleanImages(cleanName, psfName, 13000, 0, 2e5)
    
    clnjunks = ['.ms', '.ms.flagversions']
    for clnjunk in clnjunks:
        if os.path.exists(saveName.split('.')[0] + clnjunk):
            os.system('rm -rf ' + saveName.split('.')[0] + clnjunk)

#    sender_address = 'svlesovoi@gmail.com'
#    sender_pass = 'Jill21_iK'
#    receiver_address = 'svlesovoi@gmail.com'
#    
#    mail_content = '''Your data are ready at ftp://ftp.rao.istp.ac.ru/SRH/SRH0306/cleanMaps/'''  + gSaveIname + ', ' + gSaveVname
#    message = MIMEMultipart()
#    message['From'] = sender_address
#    message['To'] = receiver_address
#    message['Subject'] = 'SRH36 data ready'
#    message.attach(MIMEText(mail_content, 'plain'))
#    
#    session = smtplib.SMTP('smtp.gmail.com', 587)
#    session.starttls()
#    session.login(sender_address, sender_pass)
#    text = message.as_string()
#    session.sendmail(sender_address, receiver_address, text)
#    session.quit()
