#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 10:41:01 2019
@author: Ty Sterling

This is the main function for my MD code for problem 3:
This function calls other functions from 'moduleMD.py' to do the work.
"""
import moduleMD as md

## Input parameters ##
infile = 'ArgonCluster32.xyz' #positions
dump = 100
thermo = 100
elem = 'argon'
nondim = 'argon' #enter the atomic type for dimensionless LJ, otherwise no
dist = 'constant' #velocity distribution
val = 0 #value of constant velocity distibution
dt = 0.003 #0.002 #timestep in ps ?
tTot = 8 #ps ?

###############################################################
## You shouldn't have to change anything below here to run a
## general simulation for an FCC solid (as of 2.22.2019)
###############################################################

## Initial conditions for simulation
steps = int(tTot/dt) #number of steps to run simulation

# Call function to read initial positions
num, pos, types = md.readXYZ(infile,nondim) #non dimensional lj units
del infile #clean up variables 

# Call function to initialize velocities (trivial here)
vels = md.vInit(pos,dist,val,nondim) #get initial velocites

# Call function to initialize forces
fij, vij, vTot = md.fLJ(pos,num,nondim,elem)
###############################################################
md.tic()
print('\nRunning MD simulation!\n')
for k in range(steps): #run the MD simulation
    if (k+1)%1000 == 0:
        print('\n\tNow on step:\t'+str(k+1)+' out of '+str(steps))
    pos, vels, fij, vTot = md.vVerlet(num,pos,vels,fij,dt,nondim,elem)
    
    try: dump
    except NameError: dump = 'no'
    if type(dump) == int and (k+1)%dump == 0:
        md.dump(k,dump,num,pos,types)
        
    try: thermo
    except NameError: thermo = 'no'
    if type(thermo) == int and (k+1)%thermo == 0:
        md.thermo(k,thermo,vels,vTot,nondim,elem)
        
print('\tALL DONE!')
md.toc()

                






