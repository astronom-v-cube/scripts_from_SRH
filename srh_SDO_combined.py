#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 22 09:27:48 2022

@author: sergey_lesovoi
"""
import numpy as NP
import pylab as PL
from skimage import measure#, resize
from astropy.io import fits
from scipy import interpolate
import scipy.optimize as opt
from PIL import Image

def arcmin_format(xy, pos):
  return '%2d' % ((xy - 1023/2) * 2.45 / 60);

def create_pl(title,MJ=64):
    fig = PL.figure(figsize=(5,5))
    fig.tight_layout()
    fig.suptitle(title)
    pl = fig.subplots(nrows=1,ncols=1)
    pl.xaxis.set_major_locator(PL.MultipleLocator(MJ))
    pl.xaxis.set_minor_locator(PL.MultipleLocator(MJ//8))
    pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format))
    pl.yaxis.set_major_locator(PL.MultipleLocator(MJ))
    pl.xaxis.set_minor_locator(PL.MultipleLocator(MJ//8))
    pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format))
    pl.set_xlabel('arcmin')
    pl.set_ylabel('arcmin')
    pl.grid(linestyle='--')
    return pl

continuumFit = fits.open('/home/sergey_lesovoi/sunpy/data/hmi_ic_45s_2022_05_15_03_41_15_tai_continuum.fits')
mgFit = fits.open('/home/sergey_lesovoi/sunpy/data/hmi_m_45s_2022_05_15_03_41_15_tai_magnetogram.fits')
sdo171Fit = fits.open('/home/sergey_lesovoi/sunpy/data/aia_lev1_171a_2022_05_15t03_40_45_35z_image_lev1.fits')
sdo304Fit = fits.open('/home/sergey_lesovoi/sunpy/data/aia_lev1_304a_2022_05_15t03_40_41_13z_image_lev1.fits')
srhFit_I3200 = fits.open('../SRH0306/20220515/srh_cm_3200_20220515_03:39:01_I.fit')
srhFit_V3200 = fits.open('../SRH0306/20220515/srh_cm_3200_20220515_03:39:01_V.fit')
srhFit_I6200 = fits.open('../SRH0612/20220515/srh_cm_6200_2_20220515_03:39:05_I.fit')
srhFit_V6200 = fits.open('../SRH0612/20220515/srh_cm_6200_2_20220515_03:39:05_V.fit')

iLevels = [1e4,5e4,2e5,4e5,6e5,8e5,1.0e6,1.4e6,1.8e6,2e6,3e6,4e6]
vLevels = [-5e5,-3e5,-1e5,-2e4,-1e4,1e4,2e4,1e5,3e5,5e5]

hmiData = NP.flipud(NP.fliplr(NP.array(Image.fromarray(mgFit[1].data).resize((1024,1024)))))
hmiData = NP.roll(hmiData,(int(512 - float(mgFit[1].header['CRPIX2'])/4),int(512 - float(mgFit[1].header['CRPIX1'])/4)),axis=(0,1))

continuumData = NP.flipud(NP.fliplr(NP.array(Image.fromarray(continuumFit[1].data).resize((1024,1024)))))
continuumData = NP.roll(continuumData,(int(512 - float(continuumFit[1].header['CRPIX2'])/4),int(512 - float(continuumFit[1].header['CRPIX1'])/4)),axis=(0,1))

FOV_CONTINUUM = float(continuumFit[1].header['CDELT1']*4095)
FOV_SRH = 2.45*1023
dSrh = int(1024*FOV_CONTINUUM/FOV_SRH/2 + .5)

imSrh_I3200 =  NP.array(Image.fromarray(srhFit_I3200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))
imSrh_V3200 =  -NP.array(Image.fromarray(srhFit_V3200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))
imSrh_I6200 =  NP.array(Image.fromarray(srhFit_I6200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))
imSrh_V6200 =  -NP.array(Image.fromarray(srhFit_V6200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))

imSrh_I6200 = NP.roll(imSrh_I6200,(-8,0),axis=(0,1)).copy()
imSrh_V6200 = NP.roll(imSrh_V6200,(-8,0),axis=(0,1)).copy()

pl = create_pl('I-V 20220515 3200')
pl.contourf(imSrh_I3200,cmap='jet',levels=iLevels)
pl.contour(imSrh_V3200,cmap='seismic',levels=vLevels)

pl = create_pl('I-V 20220515 6200')
pl.contourf(imSrh_I6200,cmap='jet',levels=iLevels)
pl.contour(imSrh_V6200,cmap='seismic',levels=vLevels)

iLevels = [4e5,6e5,8e5,1.0e6,1.4e6,1.8e6,2e6]

pl = create_pl('I3200-I6200 20220515')
pl.contourf(imSrh_I3200,cmap='Reds',levels=iLevels)
pl.contourf(imSrh_I6200,cmap='Greens',levels=iLevels)

iLevels3200 = [4e5,6e5,8e5]
iLevels6200 = [2e5,4e5,6e5]
pl = create_pl('continuum-I 20220515 3200, 6200')
pl.imshow(continuumData,origin='lower',cmap='gray',interpolation='bessel')
pl.contourf(imSrh_I3200,cmap='Reds',levels=iLevels3200)
pl.contourf(imSrh_I6200,cmap='Greens',levels=iLevels6200)

FOV_HMI = float(mgFit[1].header['CDELT1']*4095)
FOV_SRH = 2.45*1023
dSrh = int(1024*FOV_HMI/FOV_SRH/2 + .5)

imSrh_I6200 =  NP.array(Image.fromarray(srhFit_I6200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))
imSrh_V6200 =  -NP.array(Image.fromarray(srhFit_V6200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))

pl = create_pl('HMI-V 20220515 6200')
pl.imshow(hmiData,origin='lower',cmap='gray',vmin=-300,vmax=300,interpolation='bessel')
pl.contour(imSrh_V6200,cmap='seismic',levels=vLevels)

pl = create_pl('HMI-I 20220515 6200')
pl.imshow(hmiData,origin='lower',cmap='gray',vmin=-300,vmax=300,interpolation='bessel')
pl.contour(imSrh_I6200,cmap='jet',levels=iLevels)

FOV_AIA = float(sdo171Fit[1].header['CDELT1']*4095)
FOV_SRH = 2.45*1023
dSrh = int(1024*FOV_AIA/FOV_SRH/2 + .5)

imSrh_I6200 =  NP.array(Image.fromarray(srhFit_I6200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))
imSrh_V6200 =  -NP.array(Image.fromarray(srhFit_V6200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))

sdo171Data = NP.array(Image.fromarray(sdo171Fit[1].data).resize((1024,1024)))
sdo171Data = NP.roll(sdo171Data,(int(512 - float(sdo171Fit[1].header['CRPIX2'])/4),int(512 - float(sdo171Fit[1].header['CRPIX1'])/4)),axis=(0,1))
pl = create_pl('20220515 03:34 SDO 171, SRH  I',MJ=128)
pl.imshow(sdo171Data,origin='lower',cmap='gnuplot',vmin=0,vmax=5000)
pl.contour(imSrh_I6200,cmap='jet',levels=iLevels)

FOV_AIA = float(sdo304Fit[1].header['CDELT1']*4095)
FOV_SRH = 2.45*1023
dSrh = int(1024*FOV_AIA/FOV_SRH/2 + .5)

imSrh_I6200 =  NP.array(Image.fromarray(srhFit_I6200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))
imSrh_V6200 =  -NP.array(Image.fromarray(srhFit_V6200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))

sdo304Data = NP.array(Image.fromarray(sdo304Fit[1].data).resize((1024,1024)))
sdo304Data = NP.roll(sdo304Data,(int(512 - float(sdo304Fit[1].header['CRPIX2'])/4),int(512 - float(sdo304Fit[1].header['CRPIX1'])/4)),axis=(0,1))
pl = create_pl('20220515 03:34 SDO 304, SRH  I',MJ=128)
pl.imshow(sdo304Data,origin='lower',cmap='hot',vmin=0,vmax=200)
pl.contour(imSrh_I6200,cmap='jet',levels=iLevels)
