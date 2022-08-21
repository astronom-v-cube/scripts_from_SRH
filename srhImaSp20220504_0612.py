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

def arcmin_format(x, pos):
  return '%.1f' % ((x - 1024/2) * 2.45 / 60);

def arcmin_y_format(y, pos):
  return '%.1f' % ((y - 1024/2 + y0) * 2.45 / 60);

def arcmin_x_format(x, pos):
  return '%.1f' % ((x - 1024/2 + x0) * 2.45 / 60);

def arcmin_y_format(y, pos):
  return '%.1f' % ((y - 1024/2 + y0) * 2.45 / 60);


fit5800I = fits.open('../SRH0612/20220504/srh_I_2022-05-04T03:40:52_5800.fit')
fit6200I = fits.open('../SRH0612/20220504/srh_I_2022-05-04T03:40:52_6200.fit')
fit6600I = fits.open('../SRH0612/20220504/srh_I_2022-05-04T03:40:52_6600.fit')
fit7000I = fits.open('../SRH0612/20220504/srh_I_2022-05-04T03:40:52_7000.fit')
fit7400I = fits.open('../SRH0612/20220504/srh_I_2022-05-04T03:40:52_7400.fit')
fit7800I = fits.open('../SRH0612/20220504/srh_I_2022-05-04T03:40:52_7800.fit')
fit8200I = fits.open('../SRH0612/20220504/srh_I_2022-05-04T03:40:52_8200.fit')
fit8600I = fits.open('../SRH0612/20220504/srh_I_2022-05-04T03:40:52_8600.fit')
fit9000I = fits.open('../SRH0612/20220504/srh_I_2022-05-04T03:40:52_9000.fit')
fit9400I = fits.open('../SRH0612/20220504/srh_I_2022-05-04T03:40:52_9400.fit')
fit9800I = fits.open('../SRH0612/20220504/srh_I_2022-05-04T03:40:52_9800.fit')
fit10200I = fits.open('../SRH0612/20220504/srh_I_2022-05-04T03:40:52_10200.fit')
fit10600I = fits.open('../SRH0612/20220504/srh_I_2022-05-04T03:40:52_10600.fit')
fit11000I = fits.open('../SRH0612/20220504/srh_I_2022-05-04T03:40:52_11000.fit')
fit11400I = fits.open('../SRH0612/20220504/srh_I_2022-05-04T03:40:52_11400.fit')
fit11800I = fits.open('../SRH0612/20220504/srh_I_2022-05-04T03:40:52_11800.fit')

fit5800V = fits.open('../SRH0612/20220504/srh_V_2022-05-04T03:40:52_5800.fit')
fit6200V = fits.open('../SRH0612/20220504/srh_V_2022-05-04T03:40:52_6200.fit')
fit6600V = fits.open('../SRH0612/20220504/srh_V_2022-05-04T03:40:52_6600.fit')
fit7000V = fits.open('../SRH0612/20220504/srh_V_2022-05-04T03:40:52_7000.fit')
fit7400V = fits.open('../SRH0612/20220504/srh_V_2022-05-04T03:40:52_7400.fit')
fit7800V = fits.open('../SRH0612/20220504/srh_V_2022-05-04T03:40:52_7800.fit')
fit8200V = fits.open('../SRH0612/20220504/srh_V_2022-05-04T03:40:52_8200.fit')
fit8600V = fits.open('../SRH0612/20220504/srh_V_2022-05-04T03:40:52_8600.fit')
fit9000V = fits.open('../SRH0612/20220504/srh_V_2022-05-04T03:40:52_9000.fit')
fit9400V = fits.open('../SRH0612/20220504/srh_V_2022-05-04T03:40:52_9400.fit')
fit9800V = fits.open('../SRH0612/20220504/srh_V_2022-05-04T03:40:52_9800.fit')
fit10200V = fits.open('../SRH0612/20220504/srh_V_2022-05-04T03:40:52_10200.fit')
fit10600V = fits.open('../SRH0612/20220504/srh_V_2022-05-04T03:40:52_10600.fit')
fit11000V = fits.open('../SRH0612/20220504/srh_V_2022-05-04T03:40:52_11000.fit')
fit11400V = fits.open('../SRH0612/20220504/srh_V_2022-05-04T03:40:52_11400.fit')
fit11800V = fits.open('../SRH0612/20220504/srh_V_2022-05-04T03:40:52_11800.fit')

Zi = ZirinTb()
freqs_0612 = NP.array([5800,6200,6600,7000,7400,7800,8200,8600,9000,9400,9800,10200,10600,11000,11400,11800])

Tmin = 1e4
Tmax = 1e6
vLevels = [-8e5,-6e5,-4e5,-2e5,-1e5,-5e4,5e4,1e5,2e5,4e5,6e5,8e5]
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

#max2800 = NP.unravel_index(NP.argmax(fit2800I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
#max3000 = NP.unravel_index(NP.argmax(fit3000I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
#max3200 = NP.unravel_index(NP.argmax(fit3200I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
#max3400 = NP.unravel_index(NP.argmax(fit3400I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
#max3600 = NP.unravel_index(NP.argmax(fit3600I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
#max3800 = NP.unravel_index(NP.argmax(fit3800I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
#max4000 = NP.unravel_index(NP.argmax(fit4000I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
#max4200 = NP.unravel_index(NP.argmax(fit4200I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
#max4400 = NP.unravel_index(NP.argmax(fit4400I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
#max4600 = NP.unravel_index(NP.argmax(fit4600I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
#max5400 = NP.unravel_index(NP.argmax(fit5400I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
#max5600 = NP.unravel_index(NP.argmax(fit5600I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
#max5800 = NP.unravel_index(NP.argmax(fit5800I[0].data[y0:y0+dy,x0:x0+dx]),(dy,dx))
#

off_0306_0612 = [5,2]

fit5800I[0].data = NP.roll(fit5800I[0].data,off_0306_0612, axis=(0,1))
fit6200I[0].data = NP.roll(fit6200I[0].data,off_0306_0612, axis=(0,1))
fit6600I[0].data = NP.roll(fit6600I[0].data,off_0306_0612, axis=(0,1))
fit7000I[0].data = NP.roll(fit7000I[0].data,off_0306_0612, axis=(0,1))
fit7400I[0].data = NP.roll(fit7400I[0].data,off_0306_0612, axis=(0,1))
fit7800I[0].data = NP.roll(fit7800I[0].data,off_0306_0612, axis=(0,1))
fit8200I[0].data = NP.roll(fit8200I[0].data,off_0306_0612, axis=(0,1))
fit8600I[0].data = NP.roll(fit8600I[0].data,off_0306_0612, axis=(0,1))
fit9000I[0].data = NP.roll(fit9000I[0].data,off_0306_0612, axis=(0,1))
fit9400I[0].data = NP.roll(fit9400I[0].data,off_0306_0612, axis=(0,1))
fit9800I[0].data = NP.roll(fit9800I[0].data,off_0306_0612, axis=(0,1))
fit10200I[0].data = NP.roll(fit10200I[0].data,off_0306_0612, axis=(0,1))
fit10600I[0].data = NP.roll(fit10600I[0].data,off_0306_0612, axis=(0,1))
fit11000I[0].data = NP.roll(fit11000I[0].data,off_0306_0612, axis=(0,1))
fit11400I[0].data = NP.roll(fit11400I[0].data,off_0306_0612, axis=(0,1))
fit11800I[0].data = NP.roll(fit11800I[0].data,off_0306_0612, axis=(0,1))

fit5800V[0].data = -NP.roll(fit5800V[0].data,off_0306_0612, axis=(0,1))
fit6200V[0].data = -NP.roll(fit6200V[0].data,off_0306_0612, axis=(0,1))
fit6600V[0].data = -NP.roll(fit6600V[0].data,off_0306_0612, axis=(0,1))
fit7000V[0].data = -NP.roll(fit7000V[0].data,off_0306_0612, axis=(0,1))
fit7400V[0].data = -NP.roll(fit7400V[0].data,off_0306_0612, axis=(0,1))
fit7800V[0].data = -NP.roll(fit7800V[0].data,off_0306_0612, axis=(0,1))
fit8200V[0].data = -NP.roll(fit8200V[0].data,off_0306_0612, axis=(0,1))
fit8600V[0].data = -NP.roll(fit8600V[0].data,off_0306_0612, axis=(0,1))
fit9000V[0].data = -NP.roll(fit9000V[0].data,off_0306_0612, axis=(0,1))
fit9400V[0].data = -NP.roll(fit9400V[0].data,off_0306_0612, axis=(0,1))
fit9800V[0].data = -NP.roll(fit9800V[0].data,off_0306_0612, axis=(0,1))
fit10200V[0].data = -NP.roll(fit10200V[0].data,off_0306_0612, axis=(0,1))
fit10600V[0].data = -NP.roll(fit10600V[0].data,off_0306_0612, axis=(0,1))
fit11000V[0].data = -NP.roll(fit11000V[0].data,off_0306_0612, axis=(0,1))
fit11400V[0].data = -NP.roll(fit11400V[0].data,off_0306_0612, axis=(0,1))
fit11800V[0].data = -NP.roll(fit11800V[0].data,off_0306_0612, axis=(0,1))

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

fig = PL.figure(figsize=(16,4))
pl = fig.subplots(nrows=2,ncols=8)
for rr in range(2):
    for cc in range(8):
        pl[rr,cc].xaxis.set_major_locator(PL.IndexLocator(24.5,offset=-8))
        pl[rr,cc].xaxis.set_major_formatter(PL.FuncFormatter(arcmin_x_format))
        pl[rr,cc].yaxis.set_major_locator(PL.IndexLocator(24.5,offset=0))
        pl[rr,cc].yaxis.set_major_formatter(PL.FuncFormatter(arcmin_y_format))
        pl[rr,cc].grid(linestyle='--')

pl[0,0].imshow(fit5800I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[0,1].imshow(fit6200I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[0,2].imshow(fit6600I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[0,3].imshow(fit7000I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[0,4].imshow(fit7400I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[0,5].imshow(fit7800I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[0,6].imshow(fit8200I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[0,7].imshow(fit8600I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)

pl[1,0].imshow(fit9000I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[1,1].imshow(fit9400I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[1,2].imshow(fit9800I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[1,3].imshow(fit10200I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[1,4].imshow(fit10600I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[1,5].imshow(fit11000I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[1,6].imshow(fit11400I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)
pl[1,7].imshow(fit11800I[0].data[y0:y0+dy,x0:x0+dx],origin='lower',vmin=Tmin,vmax=Tmax,cmap=iCmap)

pl[0,0].contour(fit5800V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[0,1].contour(fit6200V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[0,2].contour(fit6600V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[0,3].contour(fit7000V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[0,4].contour(fit7400V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[0,5].contour(fit7800V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[0,6].contour(fit8200V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[0,7].contour(fit8600V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)

pl[1,0].contour(fit9000V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[1,1].contour(fit9400V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[1,2].contour(fit9800V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[1,3].contour(fit10200V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[1,4].contour(fit10600V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[1,5].contour(fit11000V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[1,6].contour(fit11400V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
pl[1,7].contour(fit11800V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)

#pl[0,0].contour(fit2800V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
#pl[0,1].contour(fit3000V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
#pl[0,2].contour(fit3200V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
#pl[0,3].contour(fit3400V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
#pl[0,4].contour(fit3600V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
#pl[0,5].contour(fit3800V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
#pl[0,6].contour(fit4000V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
#
#pl[1,0].contour(fit4200V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
#pl[1,1].contour(fit4400V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
#pl[1,2].contour(fit4600V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
#pl[1,3].contour(fit5400V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
#pl[1,4].contour(fit5600V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)
#pl[1,5].contour(fit5800V[0].data[y0:y0+dy,x0:x0+dx],origin='lower',levels=vLevels,cmap=vCmap)

iCells_0612 = NP.zeros((16,dy,dx))
iCells_0612[0] = fit5800I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0612[1] = fit6200I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0612[2] = fit6600I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0612[3] = fit7000I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0612[4] = fit7400I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0612[5] = fit7800I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0612[6] = fit8200I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0612[7] = fit8600I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0612[8] = fit9000I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0612[9] = fit9400I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0612[10] = fit9800I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0612[11] = fit10200I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0612[12] = fit10600I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0612[13] = fit11000I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0612[14] = fit11400I[0].data[y0:y0+dy,x0:x0+dx]
iCells_0612[15] = fit11800I[0].data[y0:y0+dy,x0:x0+dx]

vCells_0612 = NP.zeros((16,dy,dx))
vCells_0612[0] = fit5800V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0612[1] = fit6200V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0612[2] = fit6600V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0612[3] = fit7000V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0612[4] = fit7400V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0612[5] = fit7800V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0612[6] = fit8200V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0612[7] = fit8600V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0612[8] = fit9000V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0612[9] = fit9400V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0612[10] = fit9800V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0612[11] = fit10200V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0612[12] = fit10600V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0612[13] = fit11000V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0612[14] = fit11400V[0].data[y0:y0+dy,x0:x0+dx]
vCells_0612[15] = fit11800V[0].data[y0:y0+dy,x0:x0+dx]

rcpCells_0612 = 1.2*0.5*(iCells_0612 + vCells_0612)
lcpCells_0612 = 1.2*0.5*(iCells_0612 - vCells_0612)

r0 = 52
c0 = 68
N = 3
S = 3

fig = PL.figure()
fig.suptitle('0612:%d, %d, %d'%(r0,c0,S))
pl = fig.subplots(nrows=3,ncols=3)
for rr in range(3):
    for cc in range(3):
        pl[rr,cc].plot(freqs_0612,rcpCells_0612[:,r0 + S*(rr-N//2),c0 + S*(cc-N//2)])
        pl[rr,cc].plot(freqs_0612,lcpCells_0612[:,r0 + S*(rr-N//2),c0 + S*(cc-N//2)])
        pl[rr,cc].set_ylim(0,3e6)


fig, pl = PL.subplots(figsize=(4,4))
pl.xaxis.set_major_locator(PL.IndexLocator(24.5,offset=-8));
pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_x_format));
pl.yaxis.set_major_locator(PL.IndexLocator(24.5,offset=0));
pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_y_format));
pl.contourf(iCells_0612[-4:].mean(axis=(0)),cmap='hot',levels=[7e4,3e5,1e6,2e6,3e6,4e6])
pl.contour(vCells_0612[-4:].mean(axis=(0)),cmap='seismic',levels=[-3e5,-2.7e5,-1.5e5,-5e4,-2e4,-1e4,1e4,2e4,5e4,1.4e5,2.e5])
pl.set_xlabel('arcmin')
pl.set_ylabel('arcmin')
pl.grid(linestyle='--')


PL.figure()
PL.plot(freqs_0306,rcpCells_0306.mean(axis=(1,2)),'+',color='red')
PL.plot(freqs_0612,rcpCells_0612.mean(axis=(1,2)),'.',color='red')
PL.plot(freqs_0306,lcpCells_0306.mean(axis=(1,2)),'+',color='blue')
PL.plot(freqs_0612,lcpCells_0612.mean(axis=(1,2)),'.',color='blue')

fig, pl = PL.subplots(figsize=(8,8))
pl.xaxis.set_major_locator(PL.MultipleLocator(128));
pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
pl.yaxis.set_major_locator(PL.MultipleLocator(128));
pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
pl.imshow((fit5800I[0].data+fit6200I[0].data+fit6600I[0].data)/3,origin='lower',cmap='jet',vmin=5e3,vmax=1e5)
pl.contourf((fit5800I[0].data+fit6200I[0].data+fit6600I[0].data)/3,origin='lower',cmap='jet',levels=[1e5,1e6,2e6])
