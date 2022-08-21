# -*- coding: utf-8 -*-
"""
Created on Sun Dec 11 02:35:24 2016

@author: serg
"""

from sys import stdin
line = ''
while line[0:4] != 'exit':
    wrongLine = False
    line = stdin.readline() # stdin
    line = line[0:len(line) - 1]
    length = len(line)
    if length % 2 == 0:
        print('even')
        for i in range (length):
            if not(line[i] =='(' or line[i] == ')'):
                wrongLine = True
        if wrongLine != True:
           correctSequence = 0
           for i in range (length):
               if line[i] == '(':
                   correctSequence += 1
               else:
                   line[i] == ')'
                   correctSequence -= 1
        print(correctSequence)
    
            
             
        
    else:
        print('odd')
    
   