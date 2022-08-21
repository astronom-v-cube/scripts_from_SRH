#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 18 02:48:38 2021

@author: svlesovoi
"""

import numpy as NP
import pylab as PL

def arcmin_format(xy, pos):
  return '%2d' % ((xy - 512/2) * 4.911 / 60);

fig, pl = PL.subplots(figsize=(8,8))
fig.subplots_adjust(left=0.05,right=0.95,top=0.95,bottom=0.1)
pl.xaxis.set_major_locator(PL.MultipleLocator(128));
pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
pl.xaxis.set_minor_locator(PL.MultipleLocator(32));
pl.yaxis.set_major_locator(PL.MultipleLocator(128));
pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
pl.yaxis.set_minor_locator(PL.MultipleLocator(32));
#pl.imshow(phaseEdit.iData,origin='lower',cmap='ocean',vmin=0.5,vmax=4)
#pl.imshow(phaseEdit.iData,origin='lower',cmap='hot',vmin=0.1,vmax=30)
#pl.imshow(phaseEdit.iData,origin='lower',cmap='hot',vmin=10e3,vmax=3.e5)
#pl.imshow(phaseEdit.iData,origin='lower',cmap='gray',vmin=.1,vmax=2)
pl.imshow(phaseEdit.iData,origin='lower',cmap='gnuplot',vmin=3e3,vmax=1e5)
#pl.imshow(phaseEdit.lcpData,origin='lower',cmap='gnuplot',vmin=0.2,vmax=2)
fig.suptitle('SRH %s %.2f GHz'%(phaseEdit.srhFits.dateObs, phaseEdit.srhFits.freqList[phaseEdit.srhFits.frequencyChannel]*1e-6))
pl.set_xlabel('arcmin')
pl.set_ylabel('arcmin')
#pl.imshow(NP.abs(phaseEdit.srhFits.uvLcp),vmax=0.01)
#fig.suptitle('SRH UV 20210918 03:00 UT 11.36 GHz')
#pl.set_xlabel('U')
#pl.set_ylabel('V')
pl.grid(linestyle='--')
