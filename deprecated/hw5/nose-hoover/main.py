#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 10:41:01 2019
Updated on Tues. Mar 5 2019
Updated on Tues. Mar 10 2019 to include periodic boundaries, cut off distance
    for the force calculation, and Verlet lists.
Updated on Wed. Apr 24 2019 to include Nose-Hoover thermostating and time-
    integration based on the Nose-Hooover equations of motion

This is GNU licensed so use it however you want but don't claim it as your own!

@author: Ty Sterling
ty.sterling@colorado.edu
"""
import functionsLJSF as md

##### SIMULATION INPUTS #####
infile = 'liquid256.xyz' #positions
dump = 5
thermo = 5
dist = 'mb' #velocity distribution
val = 100 #Kelvin #argument for velocity sampling, see docstring
dt = 0.002 #timestep in ps, nondimensionalized below
tTot = 200 #200 ps -> 100,000 steps since dt and tTot are scaled the same
tdamp = 1.0711 #Nose-Hoover damping constant in ps = 0.05 in LJ 

rcut = 2.5 #cut off distance in units of sigma, neighbor sphere
skin = 0.1 #verlet skin distance, units of sigma

###############################################################
## You shouldn't have to change anything below here to run MD
###############################################################
## Initial conditions for simulation

# Some Paramters
steps = int(tTot/dt) #number of steps to run simulation
dtMD = md.ndTime(dt) #nondimesnional timestep
tdampMD = md.ndTime(tdamp) #nondimensional damping parameter

# Call function to read initial positions
num, pos, types, box, bounds = md.readXYZ(infile) #non dimensional lj units
# Call function to initialize velocities 
vels = md.vInit(pos,dist,val) #get initial velocites
eta = 0 #Initial value for eta in Nose-Hoover equations
# Initialize verlet lists 
vlist, vcoord = md.verletList(num,pos,rcut,box) #initialize verlet lists
# Call function to initialize forces
fij, vij, vTot, vir = md.fLJ_SF(pos,num,rcut,vlist,box)

###############################################################
md.tic()
md.printParams(num,steps,dt,tTot)

for k in range(steps): #run the MD simulation
    if (k+1)%500 == 0:
        print(('\tNow on step:\t'+str(k+1)+' out of '+str(steps)))
        md.toc()
    # Check verlet list eachs step. Update when necessary,
    vlist, vcoord = md.checkVerlet(num,pos,rcut,skin,vlist,vcoord,box)
    # Nose-Hoover using Velocity Verlet algorithm.
    pos, vels, fij, eta, vTot, vir = md.vvNoseHoover(num,pos,rcut,vels,fij,eta,
                                                     val,vlist,dtMD,tdampMD,box)
    # Impose periodic boundary conditions to particle positions.
    pos = md.pbcCoords(num,pos,bounds,box)
    # Write trajectory to file
    try: dump
    except NameError: dump = 'no'
    if type(dump) == int and (k+1)%dump == 0:
        md.dump(k,dump,num,pos,types)
    # Write thermodynamic outputs to file
    try: thermo
    except NameError: thermo = 'no'
    if type(thermo) == int and (k+1)%thermo == 0:
        md.thermo(k,thermo,vels,vTot,vir,box)
        
print('\tALL DONE!')


                






