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

#srhFits = SrhFitsFile('/home/sergey_lesovoi/SRH1224/20220310/srh_01224_20220310T064459.fit', 1025)
#srhFits = SrhFitsFile('/home/sergey_lesovoi/SRH1224/20220310/srh_01224_20220310T065802.fit', 1025)
#srhFits = SrhFitsFile('/home/sergey_lesovoi/SRH1224/20220310/srh_01224_20220310T071602.fit', 1025)
#srhFits = SrhFitsFile('/home/sergey_lesovoi/SRH1224/20220310/srh_01224_20220310T074850.fit', 1025)
#srhFits = SrhFitsFile('/home/sergey_lesovoi/SRH1224/20220310/srh_01224_20220310T091852.fit', 1025)
#srhFits = SrhFitsFile('/home/sergey_lesovoi/SRH1224/20220310/srh_01224_20220310T094603.fit', 1025)
#srhFits = SrhFitsFile('/home/sergey_lesovoi/SRH1224/20220311/srh_1224_20220311T085935.fit', 1025)
#srhFits = SrhFitsFile('/home/sergey_lesovoi/SRH1224/20220311/srh_1224_20220311T091614.fit', 1025)
#srhFits = SrhFitsFile('/home/sergey_lesovoi/SRH1224/20220311/srh_1224_20220311T092354.fit', 1025)
srhFits = SrhFitsFile('/home/sergey_lesovoi/SRH1224/20220312/srh_1224_20220312T031419.fit', 1025)
srhFits.setFrequencyChannel(1)
srhFits.vis2uv_fromCoords(0)

PL.figure()
PL.imshow(NP.abs(srhFits.uvRcp + srhFits.uvLcp), vmax = 0.01)
# PL.imshow(NP.sqrt(NP.abs(srhFits.uvLcp)), vmax = 0.01)

# antennaNumbers = srhFits.antennaNumbers
# antY = srhFits.antY/srhFits.base

# PL.figure()
# PL.plot(antY[:139], NP.ones(139), '.')