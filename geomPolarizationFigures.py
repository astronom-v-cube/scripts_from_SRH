#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 15:23:36 2020

@author: sergeylesovoi
"""

import numpy as NP
import pylab as PL
from scipy import signal

fig = PL.figure()
spl0 = fig.add_subplot(221)
spl1 = fig.add_subplot(222)
spl2 = fig.add_subplot(223)
spl3 = fig.add_subplot(224)
spl0.set_xlim(-1.5,1.5)
spl0.set_ylim(-1.5,1.5)
spl0.grid()
spl1.set_xlim(-1.5,1.5)
spl1.set_ylim(-1.5,1.5)
spl1.grid()
spl2.set_xlim(-.2,.2)
spl2.set_ylim(-.2,.2)
spl2.grid()
spl3.set_xlim(-.2,.2)
spl3.set_ylim(-.2,.2)
spl3.grid()

N = 100
N_lp = 50
N_cp = 50
T = 500
t0 = 1000.2
f0 = 1.1
dF = .1

amp_nonp = NP.random.uniform(0, 1, N)
pha_nonp = NP.random.uniform(-NP.pi, NP.pi, N)
ang_nonp = NP.random.uniform(-NP.pi, NP.pi, N)
amp_lp = NP.random.randn(N_lp)
ang_lp = NP.pi/8 + NP.pi/64*NP.random.randn(N_lp)
pha_lp = NP.pi/4*NP.random.randn(N_lp)

E_x_nonp = NP.zeros(T)
E_y_nonp = NP.zeros(T)
E_x_lp = NP.zeros(T)
E_y_lp = NP.zeros(T)

for t in range(T):
    for f in range(N):
        E_t = amp_nonp[f]*NP.cos(2*NP.pi*(f0 + dF/N*f)*(t + t0) + pha_nonp[f])
        E_x_nonp[t] += E_t * NP.cos(ang_nonp[f])
        E_y_nonp[t] += E_t * NP.sin(ang_nonp[f])
        E_x_lp[t] += E_t * NP.cos(ang_nonp[f])
        E_y_lp[t] += E_t * NP.sin(ang_nonp[f])
    for f in range(N_lp):
        E_t = amp_lp[f]*NP.cos(2*NP.pi*(f0 + dF/N*f)*(t + t0) + pha_lp[f])
        E_x_lp[t] += E_t * NP.cos(ang_lp[f])
        E_y_lp[t] += E_t * NP.sin(ang_lp[f])

headLength = 0.05
headWidth = 0.03

for pol in range(N):
    spl0.arrow(0, 0, amp_nonp[pol]*NP.cos(ang_nonp[pol])/2, amp_nonp[pol]*NP.sin(ang_nonp[pol])/2,
             color='red', shape='full', head_length=headLength, head_width=headWidth)
    spl0.arrow(0, 0, -amp_nonp[pol]*NP.cos(ang_nonp[pol])/2, -amp_nonp[pol]*NP.sin(ang_nonp[pol])/2,
             color='red', shape='full', head_length=headLength, head_width=headWidth)
    
for pol in range(N):
    spl1.arrow(0, 0, amp_nonp[pol]*NP.cos(ang_nonp[pol])/2, amp_nonp[pol]*NP.sin(ang_nonp[pol])/2,
             color='red', shape='full', head_length=headLength, head_width=headWidth)
    spl1.arrow(0, 0, -amp_nonp[pol]*NP.cos(ang_nonp[pol])/2, -amp_nonp[pol]*NP.sin(ang_nonp[pol])/2,
             color='red', shape='full', head_length=headLength, head_width=headWidth)

for pol in range(N_lp):
    spl1.arrow(0, 0, amp_lp[pol]*NP.cos(ang_lp[pol])/2, amp_lp[pol]*NP.sin(ang_lp[pol])/2,
             color='green', shape='full', head_length=headLength, head_width=headWidth)
    spl1.arrow(0, 0, -amp_lp[pol]*NP.cos(ang_lp[pol])/2, -amp_lp[pol]*NP.sin(ang_lp[pol])/2,
             color='green', shape='full', head_length=headLength, head_width=headWidth)

spl2.plot((E_x_nonp)/N, (E_y_nonp)/N, '.')
spl3.plot((E_x_lp)/N, (E_y_lp)/N, '.')

spl2.set_title(r'$R_{xy} = $' + ('%.1f'%(NP.corrcoef(E_x_nonp,E_y_nonp)[0,1])))
spl3.set_title(r'$R_{xy} = $' + ('%.1f'%(NP.corrcoef(E_x_lp,E_y_lp)[0,1])))

spl0.set_xlabel('X')
spl0.set_ylabel('Y')
spl1.set_xlabel('X')
spl1.set_ylabel('Y')
spl2.set_xlabel(r'$E_x$')
spl2.set_ylabel(r'$E_y$')
spl3.set_xlabel(r'$E_x$')
spl3.set_ylabel(r'$E_y$')

fig = PL.figure()
spl0 = fig.add_subplot(221)
spl1 = fig.add_subplot(222)
spl2 = fig.add_subplot(223)
spl3 = fig.add_subplot(224)
spl0.set_xlim(-2,2)
spl0.set_ylim(-2,2)
spl0.grid()
spl1.set_xlim(-2,2)
spl1.set_ylim(-2,2)
spl1.grid()
spl2.set_xlim(-.2,.2)
spl2.set_ylim(-.2,.2)
spl2.grid()
spl3.set_xlim(-.2,.2)
spl3.set_ylim(-.2,.2)
spl3.grid()

amp_nonp = NP.random.uniform(0, 1, N)
pha_nonp = NP.random.uniform(-NP.pi, NP.pi, N)
ang_nonp = NP.random.uniform(-NP.pi, NP.pi, N)
amp_cp = NP.random.randn(N)
pha_cp = NP.pi/4*NP.random.randn(N)

E_x_nonp = NP.zeros(T)
E_y_nonp = NP.zeros(T)
E_x_cp = NP.zeros(T)
E_y_cp = NP.zeros(T)

for t in range(T):
    for f in range(N):
        E_t = amp_nonp[f]*NP.cos(2*NP.pi*(f0 + dF/N*f)*(t + t0) + pha_nonp[f])
        E_x_nonp[t] += E_t * NP.cos(ang_nonp[f])
        E_y_nonp[t] += E_t * NP.sin(ang_nonp[f])
        E_x_cp[t] += E_t * NP.cos(ang_nonp[f])
        E_y_cp[t] += E_t * NP.sin(ang_nonp[f])
    for f in range(N_cp):
        arg_cp = 2*NP.pi*(f0 + dF/N*f)*(t + t0) + pha_cp[f]
        E_x_cp[t] += amp_cp[f]*NP.cos(arg_cp)
        E_y_cp[t] -= amp_cp[f]*NP.sin(arg_cp)
E_y_hcp = signal.hilbert(E_y_cp).imag

spl2.plot(E_x_cp/N/2, E_y_cp/N/2, '.')
spl3.plot(E_x_cp/N/2, E_y_hcp/N/2, '.')
spl2.set_title(r'$R_{xy} = $' + ('%.1f'%(NP.corrcoef(E_x_cp,E_y_cp)[0,1])))
spl3.set_title(r'$R_{xy} = $' + ('%.1f'%(NP.corrcoef(E_x_cp,E_y_hcp)[0,1])))
