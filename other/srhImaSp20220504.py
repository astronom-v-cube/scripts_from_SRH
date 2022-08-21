#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 12 02:24:48 2022

@author: sergeyvlesovoi
"""

from astropy.io import fits
import pylab as PL
import numpy as NP
from ZirinTb import ZirinTb


fit2800I = fits.open('../SRH0306/20220504/srh_cm2800_03:39:42_I.fit')
fit3000I = fits.open('../SRH0306/20220504/srh_cm3000_03:39:42_I.fit')
fit3200I = fits.open('../SRH0306/20220504/srh_cm3200_03:39:42_I.fit')
fit3400I = fits.open('../SRH0306/20220504/srh_cm3400_03:39:42_I.fit')
fit3600I = fits.open('../SRH0306/20220504/srh_cm3600_03:39:42_I.fit')
fit3800I = fits.open('../SRH0306/20220504/srh_cm3800_03:39:42_I.fit')
fit4000I = fits.open('../SRH0306/20220504/srh_cm4000_03:39:42_I.fit')
fit4200I = fits.open('../SRH0306/20220504/srh_cm4200_03:39:42_I.fit')
fit4400I = fits.open('../SRH0306/20220504/srh_cm4400_03:39:42_I.fit')
fit4600I = fits.open('../SRH0306/20220504/srh_cm4600_03:39:42_I.fit')
fit5400I = fits.open('../SRH0306/20220504/srh_cm5400_03:39:42_I.fit')
fit5600I = fits.open('../SRH0306/20220504/srh_cm5600_03:39:42_I.fit')
fit5800I = fits.open('../SRH0306/20220504/srh_cm5800_03:39:42_I.fit')

fit2800V = fits.open('../SRH0306/20220504/srh_cm2800_03:39:42_V.fit')
fit3000V = fits.open('../SRH0306/20220504/srh_cm3000_03:39:42_V.fit')
fit3200V = fits.open('../SRH0306/20220504/srh_cm3200_03:39:42_V.fit')
fit3400V = fits.open('../SRH0306/20220504/srh_cm3400_03:39:42_V.fit')
fit3600V = fits.open('../SRH0306/20220504/srh_cm3600_03:39:42_V.fit')
fit3800V = fits.open('../SRH0306/20220504/srh_cm3800_03:39:42_V.fit')
fit4000V = fits.open('../SRH0306/20220504/srh_cm4000_03:39:42_V.fit')
fit4200V = fits.open('../SRH0306/20220504/srh_cm4200_03:39:42_V.fit')
fit4400V = fits.open('../SRH0306/20220504/srh_cm4400_03:39:42_V.fit')
fit4600V = fits.open('../SRH0306/20220504/srh_cm4600_03:39:42_V.fit')
fit5400V = fits.open('../SRH0306/20220504/srh_cm5400_03:39:42_V.fit')
fit5600V = fits.open('../SRH0306/20220504/srh_cm5600_03:39:42_V.fit')
fit5800V = fits.open('../SRH0306/20220504/srh_cm5800_03:39:42_V.fit')


Zi = ZirinTb()
freqs_0306 = NP.array([2800,3000,3200,3400,3600,3800,4000,4200,4400,4600,5400,5600,5800])

sky2800 = 1900
sky3000 = 0
sky3200 = -1600
sky3400 = -3100
sky3600 = -2900
sky3800 = 3600
sky4000 = 7300
sky4200 = 6200
sky4400 = 30000
sky4600 = 0
sky5400 = 0.000136
sky5600 = 0.00018
sky5800 = 0.00019

sun2800 = 32000
sun3000 = 36000
sun3200 = 37000
sun3400 = 37000
sun3600 = 35000
sun3800 = 28000
sun4000 = 24000
sun4200 = 24000
sun4400 = 24000
sun4600 = 28000
sun5400 = 0.00007
sun5600 = 0.00006
sun5800 = 0.00006

bins = 10000
hist2800 = NP.histogram((fit2800I[0].data - sky2800)/sun2800*Zi.getTbAtFrequency(freqs_0306[0]*.001)*1e3,bins=bins)
hist3000 = NP.histogram((fit3000I[0].data - sky3000)/sun3000*Zi.getTbAtFrequency(freqs_0306[1]*.001)*1e3,bins=bins)
hist3200 = NP.histogram((fit3200I[0].data - sky3200)/sun3200*Zi.getTbAtFrequency(freqs_0306[2]*.001)*1e3,bins=bins)
hist3400 = NP.histogram((fit3400I[0].data - sky3400)/sun3400*Zi.getTbAtFrequency(freqs_0306[3]*.001)*1e3,bins=bins)
hist3600 = NP.histogram((fit3600I[0].data - sky3600)/sun3600*Zi.getTbAtFrequency(freqs_0306[4]*.001)*1e3,bins=bins)
hist3800 = NP.histogram((fit3800I[0].data - sky3800)/sun3800*Zi.getTbAtFrequency(freqs_0306[5]*.001)*1e3,bins=bins)
hist4000 = NP.histogram((fit4000I[0].data - sky4000)/sun4000*Zi.getTbAtFrequency(freqs_0306[6]*.001)*1e3,bins=bins)
hist4200 = NP.histogram((fit4200I[0].data - sky4200)/sun4200*Zi.getTbAtFrequency(freqs_0306[7]*.001)*1e3,bins=bins)
hist4400 = NP.histogram((fit4400I[0].data - sky4400)/sun4400*Zi.getTbAtFrequency(freqs_0306[8]*.001)*1e3,bins=bins)
hist4600 = NP.histogram((fit4600I[0].data - sky4600)/sun4600*Zi.getTbAtFrequency(freqs_0306[9]*.001)*1e3,bins=bins)
hist5400 = NP.histogram((fit5400I[0].data - sky5400)/sun5400*Zi.getTbAtFrequency(freqs_0306[10]*.001)*1e3,bins=bins)
hist5600 = NP.histogram((fit5600I[0].data - sky5600)/sun5600*Zi.getTbAtFrequency(freqs_0306[11]*.001)*1e3,bins=bins)
hist5800 = NP.histogram((fit5800I[0].data - sky5800)/sun5800*Zi.getTbAtFrequency(freqs_0306[12]*.001)*1e3,bins=bins)


fit2800I[0].data = (fit2800I[0].data - sky2800)/sun2800*Zi.getTbAtFrequency(freqs_0306[0]*.001)*1e3
fit3000I[0].data = (fit3000I[0].data - sky3000)/sun3000*Zi.getTbAtFrequency(freqs_0306[1]*.001)*1e3
fit3200I[0].data = (fit3200I[0].data - sky3200)/sun3200*Zi.getTbAtFrequency(freqs_0306[2]*.001)*1e3
fit3400I[0].data = (fit3400I[0].data - sky3400)/sun3400*Zi.getTbAtFrequency(freqs_0306[3]*.001)*1e3
fit3600I[0].data = (fit3600I[0].data - sky3600)/sun3600*Zi.getTbAtFrequency(freqs_0306[4]*.001)*1e3
fit3800I[0].data = (fit3800I[0].data - sky3800)/sun3800*Zi.getTbAtFrequency(freqs_0306[5]*.001)*1e3
fit4000I[0].data = (fit4000I[0].data - sky4000)/sun4000*Zi.getTbAtFrequency(freqs_0306[6]*.001)*1e3
fit4200I[0].data = (fit4200I[0].data - sky4200)/sun4200*Zi.getTbAtFrequency(freqs_0306[7]*.001)*1e3
fit4400I[0].data = (fit4400I[0].data - sky4400)/sun4400*Zi.getTbAtFrequency(freqs_0306[8]*.001)*1e3
fit4600I[0].data = (fit4600I[0].data - sky4600)/sun4600*Zi.getTbAtFrequency(freqs_0306[9]*.001)*1e3
fit5400I[0].data = (fit5400I[0].data - sky5400)/sun5400*Zi.getTbAtFrequency(freqs_0306[10]*.001)*1e3
fit5600I[0].data = (fit5600I[0].data - sky5600)/sun5600*Zi.getTbAtFrequency(freqs_0306[11]*.001)*1e3
fit5800I[0].data = (fit5800I[0].data - sky5800)/sun5800*Zi.getTbAtFrequency(freqs_0306[12]*.001)*1e3

fit2800V[0].data = (fit2800V[0].data)/sun2800*Zi.getTbAtFrequency(freqs_0306[0]*.001)*1e3
fit3000V[0].data = (fit3000V[0].data)/sun3000*Zi.getTbAtFrequency(freqs_0306[1]*.001)*1e3
fit3200V[0].data = (fit3200V[0].data)/sun3200*Zi.getTbAtFrequency(freqs_0306[2]*.001)*1e3
fit3400V[0].data = (fit3400V[0].data)/sun3400*Zi.getTbAtFrequency(freqs_0306[3]*.001)*1e3
fit3600V[0].data = (fit3600V[0].data)/sun3600*Zi.getTbAtFrequency(freqs_0306[4]*.001)*1e3
fit3800V[0].data = (fit3800V[0].data)/sun3800*Zi.getTbAtFrequency(freqs_0306[5]*.001)*1e3
fit4000V[0].data = (fit4000V[0].data)/sun4000*Zi.getTbAtFrequency(freqs_0306[6]*.001)*1e3
fit4200V[0].data = (fit4200V[0].data)/sun4200*Zi.getTbAtFrequency(freqs_0306[7]*.001)*1e3
fit4400V[0].data = (fit4400V[0].data)/sun4400*Zi.getTbAtFrequency(freqs_0306[8]*.001)*1e3
fit4600V[0].data = (fit4600V[0].data)/sun4600*Zi.getTbAtFrequency(freqs_0306[9]*.001)*1e3
fit5400V[0].data = (fit5400V[0].data)/sun5400*Zi.getTbAtFrequency(freqs_0306[10]*.001)*1e3
fit5600V[0].data = (fit5600V[0].data)/sun5600*Zi.getTbAtFrequency(freqs_0306[11]*.001)*1e3
fit5800V[0].data = (fit5800V[0].data)/sun5800*Zi.getTbAtFrequency(freqs_0306[12]*.001)*1e3

Tmin = 1e4
Tmax = 6e6
vLevels = [-7e5,-4e5,-3e5,-2e5,-1e5,-5e4,5e4,1e5,2e5,3e5,4e5,7e5]
iLevels = [5e4, 1e5, 5e5, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6]
vCmap = 'seismic'
iCmap = 'hot'

x0 = 0
dx = 1024
y0 = 0
dy = 1024

x0 = 340
dx = 100
y0 = 330
dy = 150

max2800 = NP.unravel_index(NP.argmax(fit2800I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
max3000 = NP.unravel_index(NP.argmax(fit3000I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
max3200 = NP.unravel_index(NP.argmax(fit3200I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
max3400 = NP.unravel_index(NP.argmax(fit3400I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
max3600 = NP.unravel_index(NP.argmax(fit3600I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
max3800 = NP.unravel_index(NP.argmax(fit3800I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
max4000 = NP.unravel_index(NP.argmax(fit4000I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
max4200 = NP.unravel_index(NP.argmax(fit4200I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
max4400 = NP.unravel_index(NP.argmax(fit4400I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
max4600 = NP.unravel_index(NP.argmax(fit4600I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
max5400 = NP.unravel_index(NP.argmax(fit5400I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
max5600 = NP.unravel_index(NP.argmax(fit5600I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
max5800 = NP.unravel_index(NP.argmax(fit5800I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))

fit3000I[0].data = NP.roll(fit3000I[0].data,NP.array(max2800) - NP.array(max3000),axis=(0,1))
fit3200I[0].data = NP.roll(fit3200I[0].data,NP.array(max2800) - NP.array(max3200),axis=(0,1))
fit3400I[0].data = NP.roll(fit3400I[0].data,NP.array(max2800) - NP.array(max3400),axis=(0,1))
fit3600I[0].data = NP.roll(fit3600I[0].data,NP.array(max2800) - NP.array(max3600),axis=(0,1))
fit3800I[0].data = NP.roll(fit3800I[0].data,NP.array(max2800) - NP.array(max3800),axis=(0,1))
fit4000I[0].data = NP.roll(fit4000I[0].data,NP.array(max2800) - NP.array(max4000),axis=(0,1))
fit4200I[0].data = NP.roll(fit4200I[0].data,NP.array(max2800) - NP.array(max4200),axis=(0,1))
fit4400I[0].data = NP.roll(fit4400I[0].data,NP.array(max2800) - NP.array(max4400),axis=(0,1))
fit4600I[0].data = NP.roll(fit4600I[0].data,NP.array(max2800) - NP.array(max4600),axis=(0,1))
fit5400I[0].data = NP.roll(fit5400I[0].data,NP.array(max2800) - NP.array(max5400),axis=(0,1))
fit5600I[0].data = NP.roll(fit5600I[0].data,NP.array(max2800) - NP.array(max5600),axis=(0,1))
fit5800I[0].data = NP.roll(fit5800I[0].data,NP.array(max2800) - NP.array(max5800),axis=(0,1))

#main sourse
x0 = 495
dx = 100
y0 = 390
dy = 100

#aux sourse 1
#x0 = 530
#dx = 100
#y0 = 630
#dy = 100

#aux sourse 1
#x0 = 310
#dx = 100
#y0 = 360
#dy = 100

fig = PL.figure()
pl = fig.subplots(nrows=2,ncols=7)
pl[0,0].imshow(fit2800I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[0,1].imshow(fit3000I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[0,2].imshow(fit3200I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[0,3].imshow(fit3400I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[0,4].imshow(fit3600I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[0,5].imshow(fit3800I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[0,6].imshow(fit4000I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)

pl[1,0].imshow(fit4200I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[1,1].imshow(fit4200I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[1,2].imshow(fit4400I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[1,3].imshow(fit5400I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[1,4].imshow(fit5600I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[1,5].imshow(fit5800I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)

pl[0,0].contour(fit2800V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[0,1].contour(fit3000V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[0,2].contour(fit3200V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[0,3].contour(fit3400V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[0,4].contour(fit3600V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[0,5].contour(fit3800V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[0,6].contour(fit4000V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)

pl[1,0].contour(fit4200V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[1,1].contour(fit4400V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[1,2].contour(fit4600V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[1,3].contour(fit5400V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[1,4].contour(fit5600V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[1,5].contour(fit5800V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)

iCells_0306 = NP.zeros((13,dy,dx))
iCells_0306[0] = fit2800I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0306[1] = fit3000I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0306[2] = fit3200I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0306[3] = fit3400I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0306[4] = fit3600I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0306[5] = fit3800I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0306[6] = fit4000I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0306[7] = fit4200I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0306[8] = fit4400I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0306[9] = fit4600I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0306[10] = fit5400I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0306[11] = fit5600I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0306[12] = fit5800I[0].data[y0:y0+dy,x0:x0+dx]

vCells_0306 = NP.zeros((13,dy,dx))
vCells_0306[0] = fit2800V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0306[1] = fit3000V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0306[2] = fit3200V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0306[3] = fit3400V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0306[4] = fit3600V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0306[5] = fit3800V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0306[6] = fit4000V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0306[7] = fit4200V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0306[8] = fit4400V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0306[9] = fit4600V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0306[10] = fit5400V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0306[11] = fit5600V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0306[12] = fit5800V[0].data[y0:y0+dy,x0:x0+dx]

rcpCells_0306 = 0.5*(iCells_0306 + vCells_0306)
lcpCells_0306 = 0.5*(iCells_0306 - vCells_0306)

rcpCells_0306[0:-3] /= 1.2
lcpCells_0306[0:-3] /= 1.2

rcpCells_0306[-1] *= 1.2
lcpCells_0306[-1] *= 1.2
rcpCells_0306[-2] *= 1.2
lcpCells_0306[-2] *= 1.2
rcpCells_0306[-3] *= 1.2
lcpCells_0306[-3] *= 1.2

N = 3
S = 3

r0 = 52
c0 = 70
fig = PL.figure(figsize=(6,6))
fig.canvas.set_window_title('0306: %d, %d, %d'%(r0,c0,S))
pl = fig.subplots(nrows=3,ncols=3)
for rr in range(3):
    for cc in range(3):
        pl[rr,cc].plot(freqs_0306,1e-6*rcpCells_0306[:,r0 + S*(rr-N//2)-2:r0 + S*(rr-N//2)+2,c0 + S*(cc-N//2)-2:c0 + S*(cc-N//2)+2].mean(axis=(1,2)),color='red')
        pl[rr,cc].plot(freqs_0306,1e-6*lcpCells_0306[:,r0 + S*(rr-N//2)-2:r0 + S*(rr-N//2)+2,c0 + S*(cc-N//2)-2:c0 + S*(cc-N//2)+2].mean(axis=(1,2)),color='blue')
        pl[rr,cc].plot(freqs_0612,1e-6*rcpCells_0612[:,r0 + S*(rr-N//2)-2:r0 + S*(rr-N//2)+2,c0 + S*(cc-N//2)-2:c0 + S*(cc-N//2)+2].mean(axis=(1,2)),color='red')
        pl[rr,cc].plot(freqs_0612,1e-6*lcpCells_0612[:,r0 + S*(rr-N//2)-2:r0 + S*(rr-N//2)+2,c0 + S*(cc-N//2)-2:c0 + S*(cc-N//2)+2].mean(axis=(1,2)),color='blue')
        pl[rr,cc].set_ylim(0,2)
        pl[rr,cc].grid()
        if (rr == 2):
            pl[rr,cc].set_xlabel('MHz')
        if (cc == 0):
            pl[rr,cc].set_ylabel('10^6 K')

r0 = 52
c0 = 55
fig = PL.figure(figsize=(6,6))
fig.canvas.set_window_title('0306: %d, %d, %d'%(r0,c0,S))
pl = fig.subplots(nrows=3,ncols=3)
for rr in range(3):
    for cc in range(3):
        pl[rr,cc].plot(freqs_0306,1e-6*rcpCells_0306[:,r0 + S*(rr-N//2)-2:r0 + S*(rr-N//2)+2,c0 + S*(cc-N//2)-2:c0 + S*(cc-N//2)+2].mean(axis=(1,2)),color='red')
        pl[rr,cc].plot(freqs_0306,1e-6*lcpCells_0306[:,r0 + S*(rr-N//2)-2:r0 + S*(rr-N//2)+2,c0 + S*(cc-N//2)-2:c0 + S*(cc-N//2)+2].mean(axis=(1,2)),color='blue')
        pl[rr,cc].plot(freqs_0612,1e-6*rcpCells_0612[:,r0 + S*(rr-N//2)-2:r0 + S*(rr-N//2)+2,c0 + S*(cc-N//2)-2:c0 + S*(cc-N//2)+2].mean(axis=(1,2)),color='red')
        pl[rr,cc].plot(freqs_0612,1e-6*lcpCells_0612[:,r0 + S*(rr-N//2)-2:r0 + S*(rr-N//2)+2,c0 + S*(cc-N//2)-2:c0 + S*(cc-N//2)+2].mean(axis=(1,2)),color='blue')
        pl[rr,cc].set_ylim(0,2)
        pl[rr,cc].grid()
        if (rr == 2):
            pl[rr,cc].set_xlabel('MHz')
        if (cc == 0):
            pl[rr,cc].set_ylabel('10^6 K')
r0 = 52
c0 = 44
fig = PL.figure(figsize=(6,6))
fig.canvas.set_window_title('0306: %d, %d, %d'%(r0,c0,S))
pl = fig.subplots(nrows=3,ncols=3)
for rr in range(3):
    for cc in range(3):
        pl[rr,cc].plot(freqs_0306,1e-6*rcpCells_0306[:,r0 + S*(rr-N//2)-2:r0 + S*(rr-N//2)+2,c0 + S*(cc-N//2)-2:c0 + S*(cc-N//2)+2].mean(axis=(1,2)),color='red')
        pl[rr,cc].plot(freqs_0306,1e-6*lcpCells_0306[:,r0 + S*(rr-N//2)-2:r0 + S*(rr-N//2)+2,c0 + S*(cc-N//2)-2:c0 + S*(cc-N//2)+2].mean(axis=(1,2)),color='blue')
        pl[rr,cc].plot(freqs_0612,1e-6*rcpCells_0612[:,r0 + S*(rr-N//2)-2:r0 + S*(rr-N//2)+2,c0 + S*(cc-N//2)-2:c0 + S*(cc-N//2)+2].mean(axis=(1,2)),color='red')
        pl[rr,cc].plot(freqs_0612,1e-6*lcpCells_0612[:,r0 + S*(rr-N//2)-2:r0 + S*(rr-N//2)+2,c0 + S*(cc-N//2)-2:c0 + S*(cc-N//2)+2].mean(axis=(1,2)),color='blue')
        pl[rr,cc].set_ylim(0,2)
        pl[rr,cc].grid()
        if (rr == 2):
            pl[rr,cc].set_xlabel('MHz')
        if (cc == 0):
            pl[rr,cc].set_ylabel('10^6 K')

def arcmin_x_format(x, pos):
  return '%.1f' % ((x - 1024/2 + x0) * 2.45 / 60);

def arcmin_y_format(y, pos):
  return '%.1f' % ((y - 1024/2 + y0) * 2.45 / 60);

fig, pl = PL.subplots(figsize=(4,4))
pl.xaxis.set_major_locator(PL.IndexLocator(24.5,offset=-8));
pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_x_format));
pl.yaxis.set_major_locator(PL.IndexLocator(24.5,offset=0));
pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_y_format));
pl.contourf(iCells_0306.mean(axis=(0)),cmap='hot',levels=[7e4,3e5,1e6,2e6,3e6,4e6])
pl.contour(vCells_0306.mean(axis=(0)),cmap='seismic',levels=[-3e5,-2.7e5,-1.5e5,-5e4,-2e4,-1e4,1e4,2e4,5e4,1.4e5,2.e5])
pl.set_xlabel('arcmin')
pl.set_ylabel('arcmin')
pl.grid(linestyle='--')
