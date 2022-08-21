#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 15:40:18 2022

@author: maria
"""
import json

ew_phase_dif = NP.unwrap(phaseEdit36.srhFits.ewAntPhaLcp - phaseEdit36.srhFits.ewAntPhaRcp)
n_phase_dif = NP.unwrap(phaseEdit36.srhFits.nAntPhaLcp - phaseEdit36.srhFits.nAntPhaRcp)
ew_amp_dif = phaseEdit36.srhFits.ewAntAmpLcp/phaseEdit36.srhFits.ewAntAmpRcp
n_amp_dif = phaseEdit36.srhFits.nAntAmpLcp/phaseEdit36.srhFits.nAntAmpRcp

RLdif_Dict = {}
RLdif_Dict['ewPhaseDif'] = ew_phase_dif.tolist()
RLdif_Dict['nPhaseDif'] = n_phase_dif.tolist()
RLdif_Dict['ewAmpDif'] = ew_amp_dif.tolist()
RLdif_Dict['nAmpDif'] = n_amp_dif.tolist()
filename = 'RL_dif.json'
with open(filename, 'w') as saveGainFile:
    json.dump(RLdif_Dict, saveGainFile)