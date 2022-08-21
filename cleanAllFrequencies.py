#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 10:45:32 2021

@author: sergeyvlesovoi
"""
import pylab as PL
from astropy.io import fits
import os

curDate = '20211227'
freq0 = int(phaseEdit.srhFits.freqList[0] * 1e-3)
dFreq = int((phaseEdit.srhFits.freqList[1] - phaseEdit.srhFits.freqList[0]) * 1e-3)

msPath = '/home/sergeyvlesovoi/SRH_0612/'+curDate+'/'
for ff in range(16):
    print(ff)
    phaseEdit.onFrequencyChannelChanged(ff)
    phaseEdit.onCenterDisk()
    phaseEdit.saveAsMs2(msPath + 'srh%s_%d.ms'%(curDate,int(phaseEdit.srhFits.freqList[ff]*1e-3)))

with open(msPath + 'cleanAllFreqs.py', 'w') as f:
    f.write('freq0 = %d\n'%freq0)
    f.write('curDate = %s\n'%curDate)
    f.write('msPath = \'/home/sergeyvlesovoi/SRH_0612/%s/srh%s_\'\n'%(curDate,curDate))
    f.write('imagePath = \'/home/sergeyvlesovoi/SRH_0612/%s/images/cln_\'\n'%curDate)
    f.write('for curFreq in range(16):\n')
    f.write('\ttclean(vis=msPath+\'%d.ms\'%(freq0 + curFreq*400),imagename=imagePath+\'%d\'%(freq0 + curFreq*400),cell=2.45,imsize=1024,threshold=\'30000Jy\',stokes=\'RRLL\',niter=10000)')
os.system('casa -c ' + msPath + 'cleanAllFreqs.py')

with open(msPath + 'saveFitses.py', 'w') as f:
    f.write('freq0 = %d\n'%freq0)
    f.write('imagePath = \'/home/sergeyvlesovoi/SRH_0612/%s/images/cln_\'\n'%curDate)
    f.write('clnFitsPath = \'/home/sergeyvlesovoi/SRH_0612/%s/images/srh_%s_cln_\'\n'%(curDate,curDate))
    f.write('for curFreq in range(16):\n')
    f.write('\texportfits(imagename=imagePath+\'%d.image\'%(freq0 + curFreq*400),fitsimage=clnFitsPath+\'%d.fit\'%(freq0 + curFreq*400),history=False)')
os.system('casa -c ' + msPath + 'saveFitses.py')

clnFits = []
lcpImages = []
rcpImages = []
for ff in range(16):
    fitName = ('/home/sergeyvlesovoi/SRH_0612/%s/images/srh_%s_cln_%d.fit'%(curDate,curDate,int(phaseEdit.srhFits.freqList[ff]*1e-3)))
    clnFits.append(fits.open(fitName))
    lcpImages.append(clnFits[-1][0].data[0,0])
    rcpImages.append(clnFits[-1][0].data[1,0])

lcpImages = NP.array(lcpImages)
rcpImages = NP.array(rcpImages)

Y0 = 300
X0 = 640
dXY = 128
srcLcp = lcpImages[:,Y0:Y0+dXY,X0:X0+dXY]
srcRcp = rcpImages[:,Y0:Y0+dXY,X0:X0+dXY]

rightSrcLcp = srcLcp[:,75:85,79:87].max(axis=(1,2))
rightSrcRcp = srcRcp[:,75:85,79:87].max(axis=(1,2))
leftSrcLcp = srcLcp[:,55:65,40:60].max(axis=(1,2))
leftSrcRcp = srcRcp[:,55:65,40:60].max(axis=(1,2))


#fig = PL.figure()
#fig.suptitle('LCP')
#pl = fig.subplots(nrows=4, ncols=4)
#for rr in range(4):
#    for cc in range(4):
##        pl[rr,cc].imshow(rcpImages[rr*4 + cc] - 0*rcpImages[rr*4 + cc],origin='lower',vmin=0,vmax=2e5)
#        pl[rr,cc].imshow(srcLcp[rr*4 + cc],origin='lower',vmin=2e4,vmax=1e5)
#        pl[rr,cc].axis('off')
#fig.tight_layout()
#fig = PL.figure()
#fig.suptitle('RCP')
#pl = fig.subplots(nrows=4, ncols=4)
#for rr in range(4):
#    for cc in range(4):
#        pl[rr,cc].imshow(srcRcp[rr*4 + cc],origin='lower',vmin=2e4,vmax=1e5)
#        pl[rr,cc].axis('off')
#fig.tight_layout()
#-------------------------------------------------------------------------------------------------------
arcsecPerPixel = NP.abs(clnFits[15][0].header['CDELT1'] * 3600)
N = 1440
H = int(240/arcsecPerPixel)
XY0 = int(clnFits[15][0].header['CRPIX1'])
R0 = int(930/arcsecPerPixel)
overLimb = NP.zeros((16,N,H))
hfOverLimb = NP.zeros((16,N,H))

alpha = 2*NP.pi*NP.linspace(0,N,N)/N
for ff in range(16):
    for radius in range(H):
        indX = XY0 + NP.array(NP.ceil((radius + R0)*NP.cos(alpha)),dtype='int')
        indY = XY0 + NP.array(NP.ceil((radius + R0)*NP.sin(alpha)),dtype='int')
        overLimb[ff,:,radius] = 0.5*(lcpImages[ff][indY,indX] + rcpImages[ff][indY,indX])

for ff in range(16):
    hfOverLimb[ff] = overLimb[15]
    
PL.figure()
PL.imshow(NP.concatenate(overLimb,axis=1).T,origin='lower',vmin=0,vmax=3e4)
PL.contour(NP.concatenate(hfOverLimb,axis=1).T,origin='lower',levels=[3e4])

saveFitsOverLimbHdu = fits.PrimaryHDU(header=phaseEdit.srhFits.hduList[0].header)
saveFitsOverLimbPath = 'srh_0612_overlimb_' + phaseEdit.srhFits.dateObs + '.fit'

overLimbFreqColumn = fits.Column(name='frequency', format='D', array = phaseEdit.srhFits.freqList)
saveFitsOverLimbFreqExtHdu = fits.BinTableHDU.from_columns([overLimbFreqColumn])

overLimbDataColumn = fits.Column(name='data', format=('%dD'%(N*H)), array = overLimb)
saveFitsOverLimbDataExtHdu = fits.BinTableHDU.from_columns([overLimbDataColumn])

hduList = fits.HDUList([saveFitsOverLimbHdu, saveFitsOverLimbFreqExtHdu, saveFitsOverLimbDataExtHdu])
hduList.writeto(saveFitsOverLimbPath)
#-------------------------------------------------------------------------------------------------------
iSpectrum = NP.zeros((1024,1024,3))
iSpectrum[:,:,0] = (lcpImages+rcpImages)[0:5].mean(axis=0)
iSpectrum[:,:,1] = (lcpImages+rcpImages)[5:10].mean(axis=0)
iSpectrum[:,:,2] = (lcpImages+rcpImages)[10:16].mean(axis=0)
iSpectrum /= iSpectrum.max()

vSpectrum = NP.zeros((1024,1024,3))
vSpectrum[:,:,0] = (lcpImages-rcpImages)[0:5].mean(axis=0)
vSpectrum[:,:,1] = (lcpImages-rcpImages)[5:10].mean(axis=0)
vSpectrum[:,:,2] = (lcpImages-rcpImages)[10:16].mean(axis=0)
vSpectrum /= vSpectrum.max()
vSpectrum /= 2
vSpectrum += .5

fig = PL.figure(figsize=(10,10))
fig.suptitle(phaseEdit.srhFits.dateObs)
pl = fig.subplots(2,2)
PL.tight_layout()
pl[0,0].imshow(iSpectrum[:,:,0],origin='lower',cmap='Reds',vmin=0.0,vmax=0.3)
pl[0,0].axis('off')
pl[0,1].imshow(iSpectrum[:,:,1],origin='lower',cmap='Greens',vmin=0.0,vmax=0.3)
pl[0,1].axis('off')
pl[1,0].imshow(iSpectrum[:,:,2],origin='lower',cmap='Blues',vmin=0.0,vmax=0.3)
pl[1,0].axis('off')
pl[1,1].imshow(iSpectrum + .2,origin='lower',vmin=0.2,vmax=0.01)
pl[1,1].axis('off')


fig = PL.figure(figsize=(10,10))
fig.suptitle(phaseEdit.srhFits.dateObs)
pl = fig.subplots(2,2)
PL.tight_layout()
pl[0,0].imshow(vSpectrum[:,:,0],origin='lower',cmap='Reds',vmin=0.0,vmax=0.6)
pl[0,0].axis('off')
pl[0,1].imshow(vSpectrum[:,:,1],origin='lower',cmap='Greens',vmin=0.0,vmax=0.6)
pl[0,1].axis('off')
pl[1,0].imshow(vSpectrum[:,:,2],origin='lower',cmap='Blues',vmin=0.0,vmax=0.6)
pl[1,0].axis('off')
pl[1,1].imshow(vSpectrum,origin='lower',vmin=0.0,vmax=1.0)
pl[1,1].axis('off')

