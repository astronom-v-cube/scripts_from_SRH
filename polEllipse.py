# -*- coding: utf-8 -*-
"""
Created on Sat Mar  2 00:40:33 2019

@author: Sergey
"""

import numpy as NP;
import pylab as PL;

alpha = .8
couple = 0.01
 
RCP = 1
LCP = 1
 
def V_vs_alpha(LCP, RCP, alpha, couple):
    out_L = LCP*(1 + alpha)/2 + RCP*(1 - alpha)/2
    out_R = RCP*(1 + alpha)/2 + RCP*(1 - alpha)/2
    
    out_LL = out_L + couple*out_R
    out_RR = out_R + couple*out_L
     
    return (out_LL - out_RR)/(out_LL + out_RR)*100

polDeg = []
RCP = []

for i in NP.linspace(0,1,100):
    RCP.append(i)
    polDeg.append(V_vs_alpha(1, i, alpha, couple))


PL.plot(RCP, polDeg)
PL.grid()