# -*- coding: utf-8 -*-
"""
Created on Sat Sep 11 17:56:21 2021

@author: HUAWEI
"""

import numpy as np
import matplotlib.pyplot as plt

def myround(x,base,Max_Lev):
    X=np.zeros(len(x));
    for i in range(0,len(x)):
        Z=base*round(x[i]/base)
        if (Z>Max_Lev):
            X[i]=Max_Lev
        elif (Z<-Max_Lev):
            X[i]=-Max_Lev
        else:
            X[i]=Z
    return (X)

# Задаем сигнал с шумовой составляющей звезды и временную шкалу
amp=10;
dots=1000;
quant=3 #4-битное квантование

def correlator(mu1,mu2,alpha):
    S=np.random.randn(dots)
    AmpNoise1=3;
    AmpNoise2=2;
    signal1= alpha*S+np.random.randn(dots)*AmpNoise1
    signal2= alpha*S+np.random.randn(dots)*AmpNoise2
    #Выбираем максимальный уровень (по модулю), который может фиксировать антенна
    Max_Lev=7
    # Выбираем шаг квантования
    base=Max_Lev/2**(quant-1)
    #Производим квантование шумового сигнала
    signal1_q=myround(signal1,base,Max_Lev)
    signal2_q=myround(signal2,base,Max_Lev)
    # Построение графика сигнала 1 до и после квантования 
    if alpha==1.:
        plt.figure(1)
        plt.subplot(211)
        plt.plot(signal1)
        plt.xlabel('Time,count')
        plt.ylabel('Signal')
        plt.subplot(212)
        plt.plot(signal1_q)
        plt.title( str(quant)+'-битовое квантование с максимальным уровнем '+str(Max_Lev))
        plt.ylabel('Quantized signal')
        plt.xlabel('Time,count')
    # Импульсная характеристика характеризует полосовой фильтр антенн
    h1=np.zeros(100)
    h2=np.zeros(100)
    indeces = np.arange(0,100,1)
    # Отличие фильтров сделано путем разных параметров смещения mu1, mu2
    h1[:] = np.sin(0.5*indeces) *1/np.sqrt(2*np.pi*2000)*np.exp(-1/2000*(indeces-mu1)**2)
    h2[:] = np.sin(0.5*indeces) *1/np.sqrt(2*np.pi*2000)*np.exp(-1/2000*(indeces-mu2)**2)
    # Пропускаем сигнал через фильтр
    signal1_filt=np.convolve(signal1_q,h1,"same")
    signal2_filt=np.convolve(signal2_q,h2,"same")
    # Нахождение коэффициента корреляции
    CrossCorr=np.mean(signal1_filt*signal2_filt)/(np.std(signal1_filt)*np.std(signal2_filt))
    return(CrossCorr)


alpha=np.arange(0,100,0.1)
# Задаем параметры смещения
mu1=60;
mu2=20;
CrossCorr=[]
for i in range(0,len(alpha)):
    CrossCorr.append(correlator(mu1,mu2,alpha[i]))
#CrossCorr.append(correlator(mu1,mu2,alpha))
    
plt.figure(2)
plt.xlabel('alpha');
plt.ylabel('Correlation Coefficient')
plt.plot(CrossCorr,linestyle='',marker='*')
