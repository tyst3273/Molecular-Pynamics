#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 10:41:01 2019
Updated on Tues. Mar 5 2019

@author: Ty Sterling
ty.sterling@colorado.edu

module containing functions called by mainMDfx.py
"""
import numpy as np
import copy as cp
import sys

eps = 1.67e-21 #joules
sig = 3.40 #ang
mass = 6.63e-26 #kg
kb = 1.38064852e-23 #J/K

####################################
def findNN(pos,num):
    """
    This function takes in a position array with elements [id, type, x,y,z]
    and find nearest neighbors for each atom in pos. NOTE: the ids in pos 
    should be sequencial from 1 to num where num is the total number of 
    particles. 
    
    This function returns nl: the neighbor list. Each row in nl corresponds
    to each particle in pos. You can see this by looking at the FIRST column
    in nl; it should be the ids in pos. The columns in nl are the ids of the 
    unsorted nearest neighbors to each particle.
    
    As a check, the diagonal of the nd matrix should be zeros
    
    nd, the 'neighbors distances' are the distance values corresponding to 
    each neighbor in nl.
    
    NOTE: if the input pos array is nondimensional, then the distances 
    computed here are also nondimensional
    """
    nl = np.zeros((num,num)) #neighbor list that contains unique IDs for 
    # each neigbor; each row is an atom in pos, the columns in each row are
    # it's neighbors: i.e. column 1 should identically be the ids in pos
    nd = np.zeros((num,num)) #neighbor distances, each element corresponds
    # to neighbor in nl
    rij = np.zeros((num,num,3)) #vectors from each particle to all neighbors
    # needed to compute forces from LJ
    # the first index corresponds to each particle in pos, the next is
    # each other particle in pos, and the last are the x,y,z vectors
    
    for i in range(num): #loop over particles
        tmp = np.zeros((2,num)) # ids, distance
        for j in range(num): #loop over neighbors
            rx = pos[j,2]-pos[i,2] # vector from partice i to j
            ry = pos[j,3]-pos[i,3] # vector from partice i to j
            rz = pos[j,4]-pos[i,4] # vector from partice i to j
            dist = np.sqrt(rx**2+ry**2+rz**2)
            tmp[0,j] = pos[j,0] #id
            tmp[1,j] = dist #distance
            
            rij[i,j,0] = -rx
            rij[i,j,1] = -ry
            rij[i,j,2] = -rz
        
        nl[i,:] = tmp[0,:]
        nd[i,:] = tmp[1,:]
        
    #nondimensional    
    return [nl, nd, rij] #return the neighbor lists
####################################
    


#####################################
def readXYZ(infile):
    """
    This function reads in positions from an xyz file and return 3 objects:
    the total number of atoms 'num', the types and positions of atoms with rows
    format [id,type, x, y, z] in 'pos', and the names of the types in 'types'
    
    It is assumed that .xyz coords are in angstrom
    """
    with open(infile,'r') as fid:
        num = int(fid.readline()) #number of atoms
        pos = np.zeros((num,3)) #positions 
        types = np.zeros((num,1)).astype(str) #types of atoms, doesn't 
        # matter now but will be useful when we have compounds
        
        tmp = fid.readline() #skip comments
        for i in range(num): #read in positions
            tmp = fid.readline().strip().split()
            types[i,0] = tmp[0]
            pos[i,:] = tmp[1:4]
            
    utypes = np.unique(types) 
    n = len(utypes) #number of species
    tmp = np.zeros((num,1))
    for i in range(n): #change types to numbers
        tmp[np.argwhere(types[:,0] == utypes[i])] = i+1
        
    pos = np.append(tmp,pos,axis=1) #append types to pos
    types = utypes.reshape(n,1)
    pos = np.append(np.arange(1,num+1).reshape(num,1),pos,axis=1) 
    #particle ids

    #nondimensionalize    
    pos[:,2:5] = pos[:,2:5]/sig #assumed to be angstroms
    
    #nondimensional    
    return [num ,pos, types]
###################################
    


###################################
def vInit(pos,dist='constant',val=0):
    """
    This function takes in a pos array (easiest way to track types and ids)
    and return an initial velocity array based on the specified distribution 
    and initial value. 
    
    If constant, val should be the value of the velocity in ang/ps. 
    
    If maxwell boltzmann, dist should be 'mb' and val should be the 
    the temperature in kelvin
    
    Returns nondimensional velocities
    
    """   
    vels = cp.deepcopy(pos) #same elemets as pos except x->vx, y->vy, z->vz
    
    if dist == 'constant': #if constant vel across whole structure
        vels[:,2:5] = val #initialize velocities to given value
    
    elif dist == 'mb':
        num = len(pos[:,0])
        if num%2 != 0:
            num = num+1 #need an even number, will get rid of extra later
            
        x1 = np.random.uniform(0,1,(num*3/2))
        x2 = np.random.uniform(0,1,(num*3/2))
        y1 = (kb*val/mass)**.5*(-2*np.log(x1))**.5*np.cos(2*np.pi*x2)
        y2 = (kb*val/mass)**.5*(-2*np.log(x1))**.5*np.sin(2*np.pi*x2)
        num = len(pos[:,0])
        vels[:,2:5] = np.reshape(np.append(y1,y2),(num,3))[0:num,:]
        
        cmv = (mass*vels[:,2:5]).sum(axis=0)/(mass*num)
        vels[:,2:5] = vels[:,2:5]-cmv #zero the center of mass velocity
        
    else:
        sys.exit('\n\tUSAGE ERROR: Give a valid velocity distribution\n')
    
    
    #nondimensionalize the velocities
    vels[:,2:5] = vels[:,2:5]*(mass/eps)**.5 
     
    #nondimensional
    return vels
##################



####################
def fLJ(pos,num):
    """
    This function computes the interparticle pairwise force from the LJ 
    potential as well as the interparticle and total potential energy.
    
    How it works:
    it calls findNN to calculate the neighbor distances and the rij 
    vectors needed to compute the interparticle forces. it then loops over 
    all particles, and all other particles for each, to calculate the forces
    and potential
    """
    nl, nd, rij = findNN(pos,num)
    
    fij = cp.deepcopy(pos) #to copy ids, types
    vij = cp.deepcopy(pos[:,0:3]) #to copy ids, types
    fij[:,2:5] = 0 #forces, same layout as pos and vels
    vij[:,2] = 0 #potential, same layout as pos and vels
    for i in range(num):
        tmpf = np.zeros((num,3)) #tmp array to store forces
        tmpv = np.zeros((1,num)) #tmp array to store PE
        for j in range(num):
            if i == j:
                tmpf[j,:] = 0
                tmpv[0,j] = 0
            else:
                rv = rij[i,j,:] #vector from i to j
                rd = nd[i,j]
                tmpf[j,:] = (48/rd**2 * ((1/rd)**12 - 0.5*(1/rd)**6))*rv
                # nondimensional
                tmpv[0,j] = ( 4* ((1/rd)**12 - (1/rd)**6))
        fij[i,2:5] = tmpf.sum(axis=0) 
        vij[i,2] = tmpv.sum() 
    vTot = vij[:,2].sum()/2.0 #total potential, dont double count
            
    #nondimensional
    return [fij, vij, vTot]
#####################################
    
    

###################################
def vVerlet(num,pos,vels,fij,dt):
    """
    Use velocity verlet algorithm to update positions and velocities
    given initial positions, velocities, forces, and the time step.
    
    returns positions, velocities, and forces at the next time step.
    Also returns total potential energy each time step.
    
    Call within a loop to run MD simulation.
    """
    
        
    ## Compute half step velocity
    for i in range(num): #loop over all atoms
        for j in range(3): #loop over x y z
            vels[i,j+2] = vels[i,j+2]+dt*fij[i,j+2]/2.0
            #compute the velocity at the next half step for each atom in
            #x, y, and z
        
    ## Compute the next positions at full step from vhalf
    for i in range(num): #loop over all atoms
        for j in range(3): #loop over x y z
            pos[i,j+2] = pos[i,j+2]+dt*vels[i,j+2]
            #compute the position at the next full step for each atom in
            #x, y, and z based on velocity at half step
            
    ## Update the forces for the new positions
    fij, vij, vTot = fLJ(pos,num)
    
    ## Compute full step velocity
    for i in range(num): #loop over all atoms
        for j in range(3): #loop over x y z
            vels[i,j+2] = vels[i,j+2]+dt*fij[i,j+2]/2.0
            #compute the velocity at the next half step for each atom in
            #x, y, and z
            
    #nondimensional
    return [pos, vels, fij, vTot]
##################################
    


#######################################
def ndTime(dt):
    """
    Nondimensionalize timestep
    """
    dt = dt*1e-12/(sig*1e-10)*(eps/mass)**.5
    
    return dt
#############################################
    


####################################
def dump(k,dump,num,pos,types):
    """
    dump is the frequency you want to print trajectory to the file
    'traj.xyz'
    
    The rest is behind the scenes.
    """
    cpos = cp.deepcopy(pos)
    cpos[:,2:5] = cpos[:,2:5]*sig
    
    if k+1 == dump: #if the first step, erase the file. 
        with open('traj.xyz','w') as fid: #write pos trajectory to file
            fid.write(str(num)+'\n')
            fid.write('Time Step: '+str(k+1)+'\n')
            for l in range(num-1):
                spec = str(types[int(cpos[l,1])-1,0])
                fid.write(spec+'\t'+str(cpos[l,2])+'\t'+
                          str(cpos[l,3])+'\t'+str(cpos[l,4])+'\n')
            spec = str(types[int(cpos[-1,1])-1,0])
            fid.write(spec+'\t'+str(cpos[-1,2])+'\t'+
                      str(cpos[-1,3])+'\t'+str(cpos[-1,4]))
    else:
        with open('traj.xyz','a') as fid: #write pos trajectory to file
            fid.write('\n'+str(num)+'\n')
            fid.write('Time Step: '+str(k+1)+'\n')
            for l in range(num-1):
                spec = str(types[int(cpos[l,1])-1,0])
                fid.write(spec+'\t'+str(cpos[l,2])+'\t'+
                          str(cpos[l,3])+'\t'+str(cpos[l,4])+'\n')
            spec = str(types[int(cpos[-1,1])-1,0])
            fid.write(spec+'\t'+str(cpos[-1,2])+'\t'+
                      str(cpos[-1,3])+'\t'+str(cpos[-1,4]))
#####################################################################
            
            
###############################################################
def thermo(k,thermo,vels,vTot):
    """
    Computes Temp, KE, PE, Etot
        
    Writes to file 'log.MD' with frequency thermo
    """
    j2ev = 6.242e18
    num = len(vels[:,0])
    ke = np.multiply(vels[:,2:5],vels[:,2:5])/2.0 #unit mass 
    ke = ke.sum(axis=1).sum()
    temp = 2.0/3/num*ke #unit kb = 1
    temp = np.round(temp*eps/kb,decimals=3)
    ke = np.round(ke*eps*j2ev,decimals=3)
    pe = np.round(vTot*eps*j2ev,decimals=3)
    etot = np.round(ke+pe,decimals=3)
    
    mx = np.round(vels[:,2].sum()*mass*(eps/mass)**.5*1e12*1e10*1e12,decimals=3)
    my = np.round(vels[:,3].sum()*mass*(eps/mass)**.5*1e12*1e10*1e12,decimals=3)
    mz = np.round(vels[:,4].sum()*mass*(eps/mass)**.5*1e12*1e10*1e12,decimals=3)
    # pg*ang/ps
    
    if (k+1)%thermo == 0:
        if k+1 == thermo: #if the first step, erase the file. 
            with open('log.MD','w') as fid: #write energy to file
                fid.write('STEP\tT\tKE\tPE\tE\tPx\tPy\tPz\n')
                fid.write('----\tK\teV\teV\teV\tpgA/ps\tpgA/ps'
                          '\tpgA/ps\n')
                fid.write(str(k+1)+'\t'+str(temp)+'\t'+str(ke)+'\t'+
                          str(pe)+'\t'+str(etot)+'\t'+str(mx)+'\t'+
                          str(my)+'\t'+str(mz)+'\n')
        else:
            with open('log.MD','a') as fid: #append energy to file
                fid.write(str(k+1)+'\t'+str(temp)+'\t'+str(ke)+'\t'+
                          str(pe)+'\t'+str(etot)+'\t'+str(mx)+'\t'+
                          str(my)+'\t'+str(mz)+'\n')
                
#######################################################################
def tic():
    """
    Same as MATLAB tic and toc functions. Use ty.tic() at the beginning of
    code you want to time and ty.toc() at the end. Once ty.toc() is reached,
    elapsted time will be printed to screen and optionally (by default) written
    to 'log.txt' file.
    """
    import time
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    """
    Same as MATLAB tic and toc functions. Use ty.tic() at the beginning of
    code you want to time and ty.toc() at the end. Once ty.toc() is reached,
    elapsted time will be printed to screen and optionally (by default) written
    to 'log.txt' file.
    """
    import numpy as np
    import time
    if 'startTime_for_tictoc' in globals():
        print("\n\tElapsed time is "+
              str(np.round(time.time()-
                       startTime_for_tictoc,decimals=3))+" seconds.")
    else:
        print("\n\t\tToc: start time not set")


