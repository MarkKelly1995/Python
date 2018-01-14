# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 21:56:32 2018

@author: Mark
"""

import numpy as np
import scipy as sp
from scipy import constants as const
from numpy.random import rand
import matplotlib.pyplot as plt


#Create random matrix of 1 and -1. High = 2 as high is never chosen

rows = 20
cols = 20
high = 2
low = -1
step = 2


#Define this initial state, t=0. Call it state
def initialstate(rows,cols): 
    global state

    state = np.random.choice([x for x in range(low,high,step)],rows*cols)
    state.resize(rows,cols)
    return state

state=initialstate(rows, cols)
print (state)

#temp=2.25     Uncomment when geting contours
dt = 80      #no. of temp points
sweep = 800       #no. of sweeps
equil = 800       #no. of sweeps for equilibrium

#Create zeros for properties evaluating
Energy       = np.zeros(dt)
Magnetisation  = np.zeros(dt)
Magnetic Susceptibility = np.zeros(dt)
Specific Heat = np.zeros(dt)

#monte carlo using metro
def mont(point, T):
    global p
    global n_p
    for i in range(rows):
        for j in range(rows):
            a = np.random.randint(0, rows)      #row co-ord
            b = np.random.randint(0, rows)      #column co-ord
            p =  point[a, b]
            n_p = point[(a+1)%rows,b] + point[a,(b+1)%rows] + point[(a-1)%rows,b] + point[a,(b-1)%rows]  #Calc neightbours spin
            e = 2*p*n_p     #difference 
               
            if e <= 0:
                p *= -1        #if diff <= 0, flip spin
            elif np.exp(-e/T) > rand():
                p *= -1
              
            point[a, b] = p
    
    return point

e1 = 0
m1 = 0
temp=np.random.normal(2, .5, dt)   #creating gaus distrib for temp
dt = np.size(temp)     #equating dt to size of T

#Statistical Information
for m in range(len(temp)):
    point = initialstate(rows,cols)
    T1=temp[m] * 0.1 #0.1
    #Let system reach equilibrium before taking results
    for l in range(equil):         
        mont(point, T1)  
    point = initialstate(rows,cols)     
      #Define energy
    def configenergy(point):    
           energy = 0
           for i in range(len(point)):
               for j in range(len(point)):
                   energy += -n_p*p
           return energy/(2**2)
    configmag = np.abs(np.sum(point))
    #run to find properties  
    for k in range(sweep):
        mont(point, T1)           
        e1 = (e1 + configenergy(point) ) * 1/(sweep*(rows**2))
        m1 = (m1 + configmag(point) ) * 1/(sweep*(rows**2))
        ms = (T1 * 1/const.Boltzmann)     #For mag sus
        ms2 = ms * T1
        Energy[m]         = e1
        Magnetisation[m]  = m1
        
        Magnetic Susceptibility [m] = (ms * (((m1**2)*(1/(sweep)))-(mf**2)))
        Specific Heat [m] =  (ms2 * (((e1**2)*(1/(sweep))-(ef**2)))  
            	

plt.plot(T, Energy, 'b.')
plt.xlabel("Temperature (K)")
plt.ylabel("Energy (J) ")
plt.title("Energy per site vs Temperture ")
plt.axvline(2.25, linestyle='--')  #plots vertical line at critical temperature
plt.show

plt.plot(temp, Magnetisation, 'b.')
plt.xlabel("Temperature (K)")
plt.ylabel("Magnetisation per spin (A) ")
plt.title("Magnetisation per site vs Temperture ")
plt.grid()
plt.axvline(2.25, linestyle='--')  #plots vertical line at critical temperature
plt.show

plt.plot(T, Magnetic Susceptibility, 'b.')
plt.xlabel("Temperature (K)")
plt.ylabel("Magnetic Susceptibility ")
plt.ylabel("Magnetic Susceptibility per spin (A) ")
plt.axvline(2.25, linestyle='--')  #plots vertical line at critical temperature
plt.grid()
plt.show

plt.plot(T, Specific Heat, 'b.')
plt.xlabel("Temperature (K)")
plt.ylabel("Specific Heat ")
plt.ylabel("Specific heat per spin (A) ")
plt.axvline(2.25, linestyle='--')  #plots vertical line at critical temperature
plt.grid()
plt.show

'''uncomment when getting contours
print (mont(point, temp)) 
plt.imshow(mont(point, temp),cmap="hot")'''
