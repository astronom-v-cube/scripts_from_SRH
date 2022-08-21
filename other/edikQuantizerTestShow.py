#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 07:21:01 2022

@author: sergey_lesovoi
"""

import numpy as NP
import pylab as PL

qStep = 2000
qType = 0

# ampRec = phaseEdit.srhFits.ampLcp[:,:,:].mean(axis=(0,2))[1:]
# ampCor = phaseEdit.srhFits.ampLcp_c[:,:,:].mean(axis=(0,2))[1:]
ampRec = phaseEdit612.srhFits.ampLcp[:,:,:].mean(axis=(0,2))[1:]
ampCor = phaseEdit612.srhFits.ampLcp_c[:,:,:].mean(axis=(0,2))[1:]
dateObs = phaseEdit.srhFits.dateObs

PL.figure()
PL.plot(ampCor)

PL.figure()
PL.title(dateObs + ', Step %d, type %d'%(qStep,qType))
PL.plot(ampRec,ampCor,'.')
PL.plot([ampRec[ampRec.argmin()],ampRec[ampRec.argmax()]],[ampCor[ampCor.argmin()],ampCor[ampCor.argmax()]])

PL.xlabel('receiver amp')
PL.ylabel('correlator amp')
PL.grid()