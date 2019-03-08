#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 15:44:10 2019

@author: ty

This script is where I developed and verified my box-muller sampling
"""
import numpy as np
import matplotlib.pyplot as plt

eps = 1.67e-21 #joules
sig = 3.40e-10 #ang
mass = 6.63e-26 #kg
kb = 1.38064852e-23 #J/K
j2ev = 6.242e18 #J to eV

val = 20 #tempertature to sample

num = 2**14 #number of elements
xarr = np.arange(-350,350,.1)
mb = (mass/kb/val/2.0)**.5 * np.exp( -xarr**2*mass/2.0/kb/val )
mb = mb/np.trapz(mb,xarr)

if num%2 != 0:
    num = num+1

loop = 100

k1 = np.zeros(loop) #KE for each loop
k2 = np.zeros(loop) #KE for each loop, zero'd
h1 = np.zeros((num/2,loop)) #histogram for each loop
h2 = np.zeros((num/2,loop)) #histogram for each loop, zero'd
t1 = np.zeros(loop) #tmp for each loop
t2 = np.zeros(loop) #temp for each loop, zero'd
for i in range(loop):
            
    x1 = np.random.uniform(0,1,(num*3/2))
    x2 = np.random.uniform(0,1,(num*3/2))
    y1 = (kb*val/mass)**.5*(-2*np.log(x1))**.5*np.cos(2*np.pi*x2)
    y2 = (kb*val/mass)**.5*(-2*np.log(x1))**.5*np.sin(2*np.pi*x2)
    vels = np.append(y1,y2)
    
    h1[:,i] = np.histogram(vels,num/2)[0]
    if i == 0:
        hx = np.histogram(vels,num/2)[1]
    h1[:,i] = h1[:,i]/np.trapz(h1[:,i],hx[0:-1])
    ke = (mass*vels**2).sum()/2.0
    t1[i] = 2/3.0/num/kb*ke #equipartition
    k1[i] = ke*j2ev #in eV
    
    ##### ZERO THE CENTER OF MASS VELOCITY VECTOR #####
    vels = np.reshape(vels,(num,3)) 
    cmv = (mass*vels).sum(axis=0)/(mass*num)
    vels = vels-cmv
    
    ##### CHECK VALIDITY OF ZEROING #####
    cmv2 = (mass*vels).sum(axis=0)/(mass*num)
    vels = np.reshape(vels,num*3)
    h2[:,i] = np.histogram(vels,num/2)[0]
    h2[:,i] = h2[:,i]/np.trapz(h2[:,i],hx[0:-1])
    ke2 = (mass*vels**2).sum()/2.0 
    t2[i] = 2/3.0/num/kb*ke2 #equipartition
    k2[i] = ke2*j2ev #in eV

fig1, ax1 = plt.subplots()
for i in range(loop):
    ax1.plot(hx[0:-1],h1[:,i])
ax1.plot(xarr,mb,'k')
plt.title('Velocity Distribution')
plt.xlabel('vel.')
plt.ylabel('Rel. Amplitude (arb. units)')
fig1.show()
#fig1.savefig('vdist.pdf',format='pdf',dpi=300)

fig2, ax2 = plt.subplots()
for i in range(loop):
    ax2.plot(hx[0:-1],h2[:,i])
ax2.plot(xarr,mb,'k')
plt.title('Zero\'d Center of Mass Velocity Dist.')
plt.xlabel('vel.')
plt.ylabel('Rel. Amplitude (arb. units)')
fig2.show()
#fig2.savefig('zerod.pdf',format='pdf',dpi=300)

xarr = np.arange(0,loop,1)
t1u = np.mean(t1)
t2u = np.mean(t2)
t1sd = np.std(t1,ddof=1)
t2sd = np.std(t2,ddof=1)

fig3, ax3 = plt.subplots()
ax3.plot(xarr,t1,'r')
plt.axis([0,100,10,30])
plt.title('Avg. Temp over 100 runs')
plt.xlabel('run')
plt.ylabel('Temp')
fig3.show()
#fig3.savefig('temps.pdf',format='pdf',dpi=300)