#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 10:41:01 2019
Updated on Tues. Mar 5 2019

@author: Ty Sterling
ty.sterling@colorado.edu

This is the main function for my MD code for the Midterm Exam:
This function calls other functions from 'moduleMD.py' to do the work.

Note: I simplified the code quite a bit by minimizing the number of arguments
passed to each function. 
"""
import moduleMD as md

##### SIMULATION INPUTS #####
infile = 'ArgonCluster32.xyz' #positions
dump = 50
thermo = 50
dist = 'mb' #velocity distribution
val = 20 #Kelvin #argument for velocity sampling, see docstring
dt = 0.002 #timestep in ps, nondimensionalized below
tTot = 10 #ps 
##############################


###############################################################
## You shouldn't have to change anything below here to run a
## generic simulation for a FCC argon (as of 3.5.2019)
###############################################################

## Initial conditions for simulation
steps = int(tTot/dt) #number of steps to run simulation
dt = md.ndTime(dt) #nondimesnional timestep

# Call function to read initial positions
num, pos, types = md.readXYZ(infile) #non dimensional lj units
del infile #clean up variables 

# Call function to initialize velocities (trivial here)
vels = md.vInit(pos,dist,val) #get initial velocites

# Call function to initialize forces
fij, vij, vTot = md.fLJ(pos,num)

###############################################################
md.tic()
print('\nRunning MD simulation!\n')
for k in range(steps): #run the MD simulation
    if (k+1)%500 == 0:
        print('\tNow on step:\t'+str(k+1)+' out of '+str(steps))
    pos, vels, fij, vTot = md.vVerlet(num,pos,vels,fij,dt)
    
    try: dump
    except NameError: dump = 'no'
    if type(dump) == int and (k+1)%dump == 0:
        md.dump(k,dump,num,pos,types)
        
    try: thermo
    except NameError: thermo = 'no'
    if type(thermo) == int and (k+1)%thermo == 0:
        md.thermo(k,thermo,vels,vTot)
        
print('\tALL DONE!')
md.toc()

                






