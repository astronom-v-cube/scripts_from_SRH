#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 02:30:53 2020

@author: svlesovoi
"""
import numpy as NP
import pylab as PL

wavelength = 0.05
rms = NP.linspace(0,0.01,100)

PL.figure()
PL.plot(rms, NP.exp(-(4*NP.pi*rms/wavelength)**2))
PL.plot([wavelength/20,wavelength/20],[-10,10], label='$\lambda$ / 20')
PL.ylim(0,1)
PL.xlabel('rms [m]')
PL.ylabel('collecting area / geometric area')
PL.legend()

