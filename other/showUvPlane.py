#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 11:15:18 2022

@author: mariagloba
"""
import numpy as NP
from srhFitsFile1224 import SrhFitsFile
import pylab as PL
import time

#srhFits = SrhFitsFile('/home/sergey_lesovoi/SRH1224/20220312/srh_1224_20220312T034622.fit', 1025)
#srhFits = SrhFitsFile('/home/sergey_lesovoi/SRH1224/20220312/srh_1224_20220312T065123.fit', 1025)
#srhFits = SrhFitsFile('/home/sergey_lesovoi/SRH1224/20220312/srh_1224_20220312T075704.fit', 1025)
#srhFits = SrhFitsFile('/home/sergey_lesovoi/SRH1224/20220312/srh_1224_20220312T084303.fit', 1025)
srhFits = SrhFitsFile('/home/sergey_lesovoi/SRH1224/20220312/srh_1224_20220312T085612.fit', 1025)
freq = 7
srhFits.setFrequencyChannel(freq)

srhFits.calculatePhaseLcp_nonlinear(freq, baselinesNumber=8)
images = []
for ss in range(20):
    srhFits.vis2uv_fromCoords(ss)
    srhFits.uv2lmImage()
    images.append(srhFits.lcp.real)

images = NP.array(images)
PL.figure()
PL.imshow(images.mean(axis=0))

