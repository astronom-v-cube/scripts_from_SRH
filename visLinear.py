#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 06:39:18 2021

@author: svlesovoi
"""
import numpy as NP
import scipy as SP

M = NP.array([[1, 1, 0,1j,-1j, 0,  1, 0],
             [0, 1, 1, 0, 1j,-1j, 1, 0],
             [1, 0, 1,1j, 0, -1j, 0, 1]])

pM  = NP.linalg.pinv(M)
pM1 = NP.linalg.pinv(M1)

print(M.dot(pM))


