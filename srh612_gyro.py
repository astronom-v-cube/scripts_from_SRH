#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 16:22:25 2021

@author: sergey_lesovoi
"""
#srh_20211109T060252.fit

lcp5 = phaseEdit.lcpData
rcp5 = phaseEdit.rcpData
lcp4 = phaseEdit.lcpData
rcp4 = phaseEdit.rcpData
lcp3 = phaseEdit.lcpData
rcp3 = phaseEdit.rcpData
lcp2 = phaseEdit.lcpData
rcp2 = phaseEdit.rcpData
lcp1 = phaseEdit.lcpData
rcp1 = phaseEdit.rcpData
lcp0 = phaseEdit.lcpData
rcp0 = phaseEdit.rcpData

polpow=[[-35,38],[-35,42],[-50,-58],[-54,64],[-58,76],[-62,74]]

sunLcp0 = lcp0 / lcp0[256-50:256+50,256-50:256+50].mean()*Zi.getTbAtFrequncy(6.1)*1e3
sunRcp0 = rcp0 / rcp0[256-50:256+50,256-50:256+50].mean()*Zi.getTbAtFrequncy(6.1)*1e3
sunLcp1 = lcp1 / lcp1[256-50:256+50,256-50:256+50].mean()*Zi.getTbAtFrequncy(7.25)*1e3
sunRcp1 = rcp1 / rcp1[256-50:256+50,256-50:256+50].mean()*Zi.getTbAtFrequncy(7.25)*1e3
sunLcp2 = lcp2 / lcp2[256-50:256+50,256-50:256+50].mean()*Zi.getTbAtFrequncy(8.72)*1e3
sunRcp2 = rcp2 / rcp2[256-50:256+50,256-50:256+50].mean()*Zi.getTbAtFrequncy(8.72)*1e3
sunLcp3 = lcp3 / lcp3[256-50:256+50,256-50:256+50].mean()*Zi.getTbAtFrequncy(9.00)*1e3
sunRcp3 = rcp3 / rcp3[256-50:256+50,256-50:256+50].mean()*Zi.getTbAtFrequncy(9.00)*1e3
sunLcp4 = lcp4 / lcp4[256-50:256+50,256-50:256+50].mean()*Zi.getTbAtFrequncy(10.76)*1e3
sunRcp4 = rcp4 / rcp4[256-50:256+50,256-50:256+50].mean()*Zi.getTbAtFrequncy(10.76)*1e3
sunLcp5 = lcp5 / lcp5[256-50:256+50,256-50:256+50].mean()*Zi.getTbAtFrequncy(11.36)*1e3
sunRcp5 = rcp5 / rcp5[256-50:256+50,256-50:256+50].mean()*Zi.getTbAtFrequncy(11.36)*1e3

fig = PL.figure()
pl = fig.subplots(2,6)
pl[0,0].imshow(sunLcp0,vmin=0,vmax=2e5,origin='lower')
pl[0,1].imshow(sunLcp1,vmin=0,vmax=2e5,origin='lower')
pl[0,2].imshow(sunLcp2,vmin=0,vmax=2e5,origin='lower')
pl[0,3].imshow(sunLcp3,vmin=0,vmax=2e5,origin='lower')
pl[0,4].imshow(sunLcp4,vmin=0,vmax=2e5,origin='lower')
pl[0,5].imshow(sunLcp5,vmin=0,vmax=2e5,origin='lower')
pl[1,0].imshow(sunRcp0,vmin=0,vmax=2e5,origin='lower')
pl[1,1].imshow(sunRcp1,vmin=0,vmax=2e5,origin='lower')
pl[1,2].imshow(sunRcp2,vmin=0,vmax=2e5,origin='lower')
pl[1,3].imshow(sunRcp3,vmin=0,vmax=2e5,origin='lower')
pl[1,4].imshow(sunRcp4,vmin=0,vmax=2e5,origin='lower')
pl[1,5].imshow(sunRcp5,vmin=0,vmax=2e5,origin='lower')
