#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 08:16:33 2022

@author: sergey_lesovoi
"""

import numpy as NP
import pylab as PL
from MoonTb import MoonTb
from scipy.stats import linregress
import matplotlib
from ZirinTb import ZirinTb

mt = MoonTb()
zi = ZirinTb()
freqs = phaseEdit.srhFits.freqList * 1e-6

weVis0 = 3472
weVis1 = 3567

nVis0 = 3007
nVis1 = 3036

nVisLcp = phaseEdit.srhFits.visLcp[:,:,nVis0:nVis1]
nVisRcp = phaseEdit.srhFits.visRcp[:,:,nVis0:nVis1]
weVisLcp = phaseEdit.srhFits.visLcp[:,:,weVis0:weVis1]
weVisRcp = phaseEdit.srhFits.visRcp[:,:,weVis0:weVis1]

nLcpPhase = NP.angle(nVisLcp)
nRcpPhase = NP.angle(nVisRcp)
weLcpPhase = NP.angle(weVisLcp)
weRcpPhase = NP.angle(weVisRcp)

PL.figure()
#PL.plot(((nLcpPhase[:,:,:]) - (nRcpPhase[:,:,:])).mean(axis=1),'.')
#PL.plot(((weLcpPhase[:,:,:]) - (weRcpPhase[:,:,:])).mean(axis=1),'.')

# PL.plot(((nLcpPhase[0,:,0]) - (nRcpPhase[0,:,0])),'.')
# PL.plot(((nLcpPhase[0,:,1]) - (nRcpPhase[0,:,1])),'.')
# PL.plot(((nLcpPhase[0,:,2]) - (nRcpPhase[0,:,2])),'.')

PL.plot(((nLcpPhase[0,:,1]) - (nRcpPhase[0,:,1])),'.')
PL.plot(((nLcpPhase[7,:,1]) - (nRcpPhase[7,:,1])),'.')
PL.plot(((nLcpPhase[15,:,1]) - (nRcpPhase[15,:,1])),'.')


