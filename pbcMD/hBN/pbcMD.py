#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 10:41:01 2019
Updated on Tues. Mar 5 2019
Updated on Tues. Mar 10 2019 to include periodic boundaries, cut off distance
    for the force calculation, and Verlet lists.

@author: Ty Sterling
ty.sterling@colorado.edu
"""
import functionsLJSF as md

##### SIMULATION INPUTS #####
infile = 'test.xyz' #positions
dump = 10
thermo = 10
dist = 'mb' #velocity distribution
val = 300 #Kelvin #argument for velocity sampling, see docstring
dt = 0.0005 #timestep in ps, nondimensionalized below
tTot = 1 #200 #ps 

rcut = 5 #cut off distance in units of sigma, neighbor sphere
skin = .5 #verlet skin distance, units of sigma

###############################################################
## You shouldn't have to change anything below here to run a
## Simulation for a FCC argon using periodic boundary conditions,
## Verlet lists, and cutoff LJ forces (NOT DONE WITH CUTOFF)
###############################################################
## Initial conditions for simulation

# Some Paramters
steps = int(tTot/dt) #number of steps to run simulation

# Call function to read initial positions
num, pos, types, box, bounds = md.readXYZ(infile) #non dimensional lj units
# Call function to initialize velocities k
vels = md.vInit(pos,dist,val) #get initial velocites
# Initialize verlet lists 
vlist, vcoord = md.verletList(num,pos,rcut,box) #initialize verlet lists
# Call function to initialize forces
##fij, vij, vTot = md.fLJ(pos,num,vlist,box)
fij, vij, vTot = md.fLJ_SF(pos,num,rcut,vlist,box)
#
################################################################
#md.tic()
#md.printParams(num,steps,dt,tTot)
#
#for k in range(steps): #run the MD simulation
#    if (k+1)%500 == 0:
#        print(('\tNow on step:\t'+str(k+1)+' out of '+str(steps)))
#        md.toc()
#    # Check verlet list eachs step. Update when necessary,
#    vlist, vcoord = md.checkVerlet(num,pos,rcut,skin,vlist,vcoord,box)
#    # Velocity Verlet algorithm.
#    pos, vels, fij, vTot = md.vVerlet(num,pos,rcut,vels,fij,vlist,dtMD,box)
#    # Impose periodic boundary conditions to particle positions.
#    pos = md.pbcCoords(num,pos,bounds,box)
#    # Write trajectory to file
#    try: dump
#    except NameError: dump = 'no'
#    if type(dump) == int and (k+1)%dump == 0:
#        md.dump(k,dump,num,pos,types)
#    # Write thermodynamic outputs to file
#    try: thermo
#    except NameError: thermo = 'no'
#    if type(thermo) == int and (k+1)%thermo == 0:
#        md.thermo(k,thermo,vels,vTot)
#        
#print('\tALL DONE!')
#

                






