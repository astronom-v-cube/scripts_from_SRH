#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 06:58:27 2022

@author: sergey_lesovoi
"""
import numpy as NP
import pylab as PL

EW_names = []

for ant in range(24):
    EW_names.append('W%02d'%(24-ant))

EW_names.append('C1')
    
for ant in range(24):
    EW_names.append('E%02d'%(ant + 1))

EW_base = NP.array([
    2461,\
    2438,\
    2458,\
    2455,\
    2438,\
    2452,\
    2460,\
    2460,\
    
    2434,\
    2458,\
    2431,\
    2466,\
    2448,\
    2448,\
    2471,\
    2422,\

    2454,\
    2464,\
    2444,\
    2438,\
    2451,\
    2463,\
    2438,\
    2468,\
#-------------------------        
    2450,\
    2448,\
    2435,\
    2461,\
    2460,\
    2435,\
    2466,\
    2437,\
    
    2462,\
    2456,\
    2427,\
    2451,\
    2458,\
    2444,\
    2448,\
    2450,\
        
    2454,\
    2467,\
    2428,\
    2446,\
    2453,\
    2466,\
    2435,\
    2437   
    ])
    
antsWE = NP.zeros(49)

for ant in range(24):
    antsWE[23-ant] = -EW_base[23-ant:24].sum(0)
    antsWE[25+ant] =  EW_base[24:24+ant+1].sum(0)

idealAnts = NP.linspace(-24,24,49)*2450

PL.figure(figsize=(14,6))
PL.ylim(-30,30)
PL.grid()
PL.title('SRH 1224 antenna positions 20220402')
PL.ylabel('mm')
PL.plot(EW_names,antsWE - idealAnts)
PL.plot(EW_names,antsWE - idealAnts,'.')
PL.xticks(rotation=45)


b_69_70 = NP.array([  3.63764128,  -1.34388628,   6.16583072, -10.71523426])
b_70_71 = NP.array([ -0.99791828,  -4.78472572,  -7.2994653 , -10.24980725])
b_71_72 = NP.array([ 2.7970959 , -7.81632994,  1.88758247, -4.10357888])
b_72_73 = NP.array([  2.70606091,  -5.98780197,  11.08579546,-18.83314023])
b_73_74 = NP.array([ 1.02046384,  2.74735392, -5.38064774, 13.50784481])
