#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 07:06:05 2021

@author: svlesovoi
"""
import numpy as NP

def phaMatrixGen(antNumber):
    phaMatrix = NP.zeros((antNumber - 1, antNumber + 1))
    phaMatrix[:,0] = 1
    for pair in range(antNumber - 1):
        phaMatrix[pair, pair + 1] = 1
        phaMatrix[pair, pair + 2] = -1
    return phaMatrix.copy()

def phaMatrixGenEWN(antNumberEW, antNumberN):
    phaMatrix = NP.zeros((antNumberEW - 1 + antNumberN, antNumberEW + 1 + antNumberN + 1))
    phaMatrix[:antNumberEW-1,0] = 1
    phaMatrix[antNumberEW-1:antNumberEW-1+antNumberN,1] = 1
    for pair in range(antNumberEW - 1):
        phaMatrix[pair, pair+1 + 1] = 1
        phaMatrix[pair, pair+1 + 2] = -1
    for pair in range(antNumberN):
        if pair == 0:
            phaMatrix[antNumberEW - 1, 2 + antNumberEW//2] = 1
        else:
            phaMatrix[pair + antNumberEW - 1, pair + 1 + antNumberEW] = 1
        phaMatrix[pair + antNumberEW - 1, pair + 2 + antNumberEW] = -1
    return phaMatrix.copy()

def phaMatrixGenPairsEWN(pairs, antNumberEW, antNumberN):
    rowsEW = int(((antNumberEW - 1) + (antNumberEW - pairs))/2 * pairs)
    rowsN = int(((antNumberN) + (antNumberN + 1 - pairs))/2 * pairs)
    colsEW = antNumberEW + pairs
    colsN = antNumberN + pairs
    phaMatrix = NP.zeros((rowsEW + rowsN, colsEW + colsN))
    for pair in range(pairs):
        row0 = int(((antNumberEW - 1) + (antNumberEW - pair))/2 * pair)
        row1 = row0 + (antNumberEW - pair - 1)
        phaMatrix[row0:row1,pair] = 1
        for phaPair in range(antNumberEW - pair - 1):
            phaMatrix[phaPair + row0, phaPair + 2*pairs] = 1
            phaMatrix[phaPair + row0, phaPair + 2*pairs + (pair + 1)] = -1
        row0 = int(((antNumberN) + (antNumberN + 1 - pair))/2 * pair)
        row1 = row0 + (antNumberN - pair)
        phaMatrix[row0 + rowsEW:row1 + rowsEW,pairs + pair] = 1
        for phaPair in range(antNumberN - pair):
            if phaPair == 0:
                phaMatrix[rowsEW + row0, 2*pairs + antNumberEW//2] = 1
            else:
                phaMatrix[phaPair + rowsEW + row0, 2*pairs + phaPair + antNumberEW - 1] = 1
            phaMatrix[phaPair + rowsEW + row0, 2*pairs + phaPair + antNumberEW - 1 + (pair + 1)] = -1
    return phaMatrix.copy()
