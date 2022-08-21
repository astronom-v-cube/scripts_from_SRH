# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 02:31:16 2018

@author: Sergey
"""

import numpy as NP
import pylab as PL

L = 1000
Q_sun = NP.deg2rad(16/60)
uv = NP.linspace(0,L,L)

vis = NP.sinc(NP.pi*uv*Q_sun/2)
dU = 7
#ivis = NP.array([NP.sum(vis[0:s]) for s in range(0,L - 1)])
ivis = NP.array(NP.diff(vis,dU))

PL.clf()
PL.grid()
PL.xlim(0,1000)
PL.plot(uv,vis)
PL.plot(uv,NP.roll(vis,20) - vis)
#L.plot(uv,vis)
#PL.plot(uv[:-dU],ivis)
