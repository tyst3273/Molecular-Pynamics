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
import mdFunctions as md
import numpy as np

##### SIMULATION INPUTS #####
infile = 'cube.xyz' #positions
dump = 50
thermo = 10
dist = 'mb' #velocity distribution
val = 20 #Kelvin #argument for velocity sampling, see docstring
dt = 0.002 #timestep in ps, nondimensionalized below
tTot = 0.5 #ps 

pbc = 'yes'
rcut = 1.5 #cut off distance in units of sigma, neighbor sphere
#Set cutoff distance to a large number (larger than box size) to avoid 
#generating verlet lists.
skin = 0.5 #verlet skin distance
##############################


###############################################################
## You shouldn't have to change anything below here to run a
## generic simulation for a FCC argon (as of 3.5.2019)
###############################################################

## Initial conditions for simulation
steps = int(tTot/dt) #number of steps to run simulation
dtMD = md.ndTime(dt) #nondimesnional timestep

# Call function to read initial positions
num, pos, types, box = md.readXYZ(infile) #non dimensional lj units
del infile #clean up variables 

# Call function to initialize velocities (trivial here)
vels = md.vInit(pos,dist,val) #get initial velocites
   
# Initializer verlet lists 
vlist, vcoord = md.verletList(num,pos,rcut,pbc,box) #initialize verlet lists



vlist, vcoord = md.checkVerlet(num,pos,rcut,skin,vlist,vcoord,pbc='no',box=0)
# Check validity of Verlet list, put this part in the MD loop

    
    

## Call function to initialize forces
#fij, vij, vTot = md.fLJ(pos,num)


################################################################
#md.tic()
#print('\n\tRunning MD simulation!\n')
#print('\tNo. of atoms:\t\t'+str(num)+'\t--')
#print('\tNo. of steps:\t\t'+str(steps)+'\t--')
#print('\tTime step:\t\t'+str(dt)+'\tps')
#print('\tTotal duration:\t\t'+str(tTot)+'\tps\n')
#print('\t------------------------------------\n')
#
#for k in range(steps): #run the MD simulation
#    if (k+1)%500 == 0:
#        print('\tNow on step:\t'+str(k+1)+' out of '+str(steps))
#    pos, vels, fij, vTot = md.vVerlet(num,pos,vels,fij,dtMD)
#    
#    try: dump
#    except NameError: dump = 'no'
#    if type(dump) == int and (k+1)%dump == 0:
#        md.dump(k,dump,num,pos,types)
#        
#    try: thermo
#    except NameError: thermo = 'no'
#    if type(thermo) == int and (k+1)%thermo == 0:
#        md.thermo(k,thermo,vels,vTot)
#        
#print('\n\tALL DONE!')
#md.toc()

                






