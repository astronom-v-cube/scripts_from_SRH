#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 07:23:04 2020

@author: svlesovoi
"""

import numpy as NP
import pylab as PL

N = 1000
x = NP.linspace(0,10,N)
noise = NP.random.randn(N)*.05
y = x % 2*NP.pi
yn = (x + noise) % 2*NP.pi

PL.figure()
PL.plot(x, y)
PL.plot(x, yn)
