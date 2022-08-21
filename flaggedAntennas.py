#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 01:43:50 2021

@author: svlesovoi
"""

import numpy as NP
import srh36Utils


import requests

def telegram_bot_sendtext(bot_message):
   bot_token = '1807793449:AAFSfsfxlZjM1Du-MrnIkwnEeITVnmlpcy0'
   bot_chatID = '-1001496182817'
   send_text = 'https://api.telegram.org/bot' + bot_token + '/sendMessage?chat_id=' + bot_chatID + '&parse_mode=Markdown&text=' + bot_message

   response = requests.get(send_text)

   return response.json()



flagDataFile = open('flagdata.last')
flagDataText = flagDataFile.readlines()
flagDataFile.close()
flaggedAntennaText = flagDataText[-1].split(',antenna=')[1].split(',uvrange=')[0]
flaggedAntennaNumber = NP.array(flaggedAntennaText.replace('"',"").split(','),dtype='int')
n2nDict = srh36Utils.number2name()
flaggedMsg = 'Flagged antennas are '
for ant in flaggedAntennaNumber:
    flaggedMsg += '%s '%n2nDict[ant]

test = telegram_bot_sendtext(flaggedMsg)
