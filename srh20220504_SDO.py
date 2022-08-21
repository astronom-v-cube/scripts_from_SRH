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

mgFit = fits.open('/home/sergey_lesovoi/sunpy/data/hmi_m_45s_2022_05_04_03_35_15_tai_magnetogram.fits')
sdo171Fit = fits.open('/home/sergey_lesovoi/sunpy/data/aia_lev1_171a_2022_05_04t03_34_45_35z_image_lev1.fits')
sdo304Fit = fits.open('/home/sergey_lesovoi/sunpy/data/aia_lev1_304a_2022_05_04t03_34_53_13z_image_lev1.fits')
srhFit_I6200 = fits.open('../SRH0612/20220504/srh_I_2022-05-04T03:40:52_6200.fit')
srhFit_I7000 = fits.open('../SRH0612/20220504/srh_I_2022-05-04T03:40:52_7000.fit')
srhFit_I10200 = fits.open('../SRH0612/20220504/srh_I_2022-05-04T03:40:52_10200.fit')
srhFit_V6200 = fits.open('../SRH0612/20220504/srh_V_2022-05-04T03:40:52_6200.fit')
srhFit_V7000 = fits.open('../SRH0612/20220504/srh_V_2022-05-04T03:40:52_7000.fit')
srhFit_V10200 = fits.open('../SRH0612/20220504/srh_V_2022-05-04T03:40:52_10200.fit')

hmiData = NP.flipud(NP.fliplr(NP.array(Image.fromarray(mgFit[1].data).resize((1024,1024)))))
hmiData = NP.roll(hmiData,(int(512 - float(mgFit[1].header['CRPIX2'])/4),int(512 - float(mgFit[1].header['CRPIX1'])/4)),axis=(0,1))

FOV_HMI = float(mgFit[1].header['CDELT1']*4095)
FOV_SRH = 2.45*1023
dSrh = int(1024*FOV_HMI/FOV_SRH/2 + .5)

imSrh_I6200 =  NP.array(Image.fromarray(srhFit_I6200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))
imSrh_I7000 =  NP.array(Image.fromarray(srhFit_I7000[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))
imSrh_I10200 = NP.array(Image.fromarray(srhFit_I10200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))
imSrh_V6200 =  -NP.array(Image.fromarray(srhFit_V6200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))
imSrh_V7000 =  -NP.array(Image.fromarray(srhFit_V7000[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))
imSrh_V10200 = -NP.array(Image.fromarray(srhFit_V10200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))

vLevels = [-5e5,-4e5,-3e5,-2e5,-1e5,-5e4,-2e4,-1e4,1e4,2e4,5e4,1e5,2e5,3e5,4e5,5e5]
pl = create_pl('V 20220504 6200')
pl.imshow(hmiData,origin='lower',cmap='gray',vmin=-300,vmax=300,interpolation='bessel')
pl.contour(imSrh_V6200,cmap='seismic',levels=vLevels)
pl.set_xlim(500,600)
pl.set_ylim(370,470)

pl = create_pl('V 20220504 7000')
pl.imshow(hmiData,origin='lower',cmap='gray',vmin=-300,vmax=300,interpolation='bessel')
pl.contour(imSrh_V7000,cmap='seismic',levels=vLevels)
pl.set_xlim(500,600)
pl.set_ylim(370,470)

pl = create_pl('V 20220504 10200')
pl.imshow(hmiData,origin='lower',cmap='gray',vmin=-300,vmax=300,interpolation='bessel')
pl.contour(imSrh_V10200,cmap='seismic',levels=vLevels)
pl.set_xlim(500,600)
pl.set_ylim(370,470)

iLevels = [1e5,2e5,4e5,6e5,8e5,1.0e6,1.4e6,1.8e6,2e6]
pl = create_pl('I 20220504 6200')
pl.imshow(hmiData,origin='lower',cmap='gray',vmin=-300,vmax=300,interpolation='bessel')
pl.contour(imSrh_I6200,cmap='hot',levels=iLevels)
# pl.set_xlim(500,600)
# pl.set_ylim(370,470)

pl = create_pl('I 20220504 10200')
pl.imshow(hmiData,origin='lower',cmap='gray',vmin=-300,vmax=300,interpolation='bessel')
pl.contour(imSrh_I10200,cmap='hot',levels=iLevels)
# pl.set_xlim(500,600)
# pl.set_ylim(370,470)

iLevels = [3e4,5e4,1e5,2e5,4e5,6e5,8e5,1.0e6,1.4e6,1.8e6,2e6]
iLevelsHi = [3e4,5e4,1e5,2e5,4e5,6e5,8e5,1.0e6,1.4e6,1.8e6,2e6]

FOV_AIA = float(sdo171Fit[1].header['CDELT1']*4095)
FOV_SRH = 2.45*1023
dSrh = int(1024*FOV_AIA/FOV_SRH/2 + .5)

imSrh_I6200 =  NP.array(Image.fromarray(srhFit_I6200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))
imSrh_I10200 = NP.array(Image.fromarray(srhFit_I10200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))
imSrh_V6200 =  -NP.array(Image.fromarray(srhFit_V6200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))
imSrh_V10200 = -NP.array(Image.fromarray(srhFit_V10200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))

sdo171Data = NP.array(Image.fromarray(sdo171Fit[1].data).resize((1024,1024)))
sdo171Data = NP.roll(sdo171Data,(int(512 - float(sdo171Fit[1].header['CRPIX2'])/4),int(512 - float(sdo171Fit[1].header['CRPIX1'])/4)),axis=(0,1))
pl = create_pl('20220504 03:34 SDO 171, SRH  I',MJ=128)
pl.imshow(sdo171Data,origin='lower',cmap='gnuplot',vmin=0,vmax=5000)
# pl.contourf(imSrh_I6200,cmap='hot',levels=iLevels)
# pl.contourf(imSrh_I10200,cmap='jet',levels=iLevelsHi)
pl.contour(imSrh_V6200,cmap='hot',levels=vLevels)
pl.contour(imSrh_V10200,cmap='jet',levels=vLevels)

FOV_AIA = float(sdo304Fit[1].header['CDELT1']*4095)
FOV_SRH = 2.45*1023
dSrh = int(1024*FOV_AIA/FOV_SRH/2 + .5)

imSrh_I6200 =  NP.array(Image.fromarray(srhFit_I6200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))
imSrh_I10200 = NP.array(Image.fromarray(srhFit_I10200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))
imSrh_V6200 =  -NP.array(Image.fromarray(srhFit_V6200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))
imSrh_V10200 = -NP.array(Image.fromarray(srhFit_V10200[0].data[512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))

sdo304Data = NP.array(Image.fromarray(sdo304Fit[1].data).resize((1024,1024)))
sdo304Data = NP.roll(sdo304Data,(int(512 - float(sdo304Fit[1].header['CRPIX2'])/4),int(512 - float(sdo304Fit[1].header['CRPIX1'])/4)),axis=(0,1))
pl = create_pl('20220504 03:34 SDO 304, SRH  I',MJ=128)
pl.imshow(sdo304Data,origin='lower',cmap='jet',vmin=0,vmax=200)
pl.contour(imSrh_I6200,cmap='gray_r',levels=iLevels)
pl.contourf(imSrh_I10200,cmap='jet',levels=iLevelsHi)

