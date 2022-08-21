#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 07:21:11 2022

@author: sergey_lesovoi
"""

import numpy as NP

N = 8
operator = NP.zeros((N-1,2*N))

slopes = NP.zeros((N-1))

h = NP.deg2rad(45)
A = NP.sin(h)
B = NP.cos(h)

for r in range(N-1):
    operator[r,r] = -A
    operator[r,r + 1] = A
    operator[r,r + N] = -B
    operator[r,r + 1 + N] = B
    
slopes = NP.random.randn(N-1)

dR_in = NP.linalg.pinv(operator).dot(slopes)

dR_in = NP.random.randn(2*N)

for b in range(N-1):
    slopes[b] = (dR_in[b + 1] - dR_in[b])*A + (dR_in[b + 1 + N] - dR_in[b + N])*B

dR_out = NP.linalg.pinv(operator).dot(slopes)


operator = NP.zeros((N-1,2*(N - 1)))

for r in range(N-1):
    operator[r,r] = A
    operator[r,r + N - 1] = B

dR_in = NP.random.randn(2*(N - 1))

for b in range(N-1):
    slopes[b] = dR_in[b]*A + dR_in[b + N -  1]*B

dR_out = NP.linalg.pinv(operator).dot(slopes)
    