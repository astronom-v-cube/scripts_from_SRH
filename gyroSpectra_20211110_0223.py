#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 11:49:23 2021

@author: sergey_lesovoi
"""

lcpImages = NP.zeros((16,512,512))
rcpImages = NP.zeros((16,512,512))

lcpImages[phaseEdit.currentFrequencyChannel] = phaseEdit.lcpData
rcpImages[phaseEdit.currentFrequencyChannel] = phaseEdit.rcpData

for ii in range(16):
    pl[0,ii].imshow(lcpImages[ii,300:320,390:410],origin='lower')
    pl[1,ii].imshow(rcpImages[ii,300:320,390:410],origin='lower')
    pl[0,ii].axis('off')
    pl[1,ii].axis('off')
    
PL.plot(rcpImages[:,290:320,380:410].mean(axis=(1,2)),'.',color='red')
PL.plot(lcpImages[:,290:320,380:410].mean(axis=(1,2)),'.',color='blue')
PL.plot(rcpImages[:,140:170,204:234].mean(axis=(1,2)),color='red')
PL.plot(lcpImages[:,140:170,204:234].mean(axis=(1,2)),color='blue')

