#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Module containing python functions that are called by my main MD code.
author: Tyler Sterling
"""
import numpy as np
import copy as cp
import sys

######################################################
def ljND(elem):
    """
    This function takes an element type as an argument and return LJ 
    parameters
    Only supports el = 'argon' right now!
    """
    if elem == 'argon':
        eps = 1.67e-21 #joules
        sig = 3.40 #ang
        mass = 6.63e-26 #kg
    else:
        sys.exit('\n\tArgument of ljND must be a valid element name!\n')
    
    return [eps, sig, mass]
#####################################
    


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
        
#        tmp = tmp[:,np.argsort(tmp[1,:])] #sort by distances
        #dont sort, will waste time later computing forces
        
        nl[i,:] = tmp[0,:]
        nd[i,:] = tmp[1,:]
        
    return [nl, nd, rij] #return the neighbor lists
####################################
    


#####################################
def readXYZ(infile,nondim='no'):
    """
    This function reads in positions from an xyz file and return 3 objects:
    the total number of atoms 'num', the types and positions of atoms with rows
    format [id,type, x, y, z] in 'pos', and the names of the types in 'types'
    
    It is assumed that .xyz coords are in angstrom
    
    set nondim to the atom type to nondimensionalize
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
    
    if nondim != 'no':
        eps, sig, mass = ljND(nondim)
        pos[:,2:5] = pos[:,2:5]/sig #assumed to be angstroms
        
    return [num ,pos, types]
###################################
    


###################################
def vInit(pos,dist='constant',val=0,nondim='no'):
    """
    This function takes in a pos array (easiest way to track types and ids)
    and return an initial velocity array based on the specified distribution 
    and initial value. If constant, val should be the value of the velocity
    in ang/ps. Support for other distributions will come later. 
    
    nondim is an optional argument which should be the type of particle;
    the velocity data will be returned in non-dimensionaliozed length/ps units
    
    """   
    vels = cp.deepcopy(pos) #same elemets as pos except x->vx, y->vy, z->vz
    
    if dist == 'constant': #if constant vel across whole structure
        vels[:,2:5] = val #initialize velocities to given value
    
    if nondim != 'no': #nondimensionalize length component of velocity
        eps, sig, mass = ljND(nondim)
        vels[:,2:5] = vels[:,2:5]/sig #assumed to be ang/ps
        #nondimensionalize the length component of velocity
        
    return vels
##################



####################
def fLJ(pos,num,nondim='no',elem='argon'):
    """
    This function computes the interparticle pairwise force from the LJ 
    potential as well as the interparticle and total potential energy.
    
    How it works:
        it calls findNN to calculate the neighbor distances and the rij 
        vectors needed to compute the interparticle forces. 
        
        it then loops over all particles, and all other particles for each, 
        to calculate the 
        
    set nondim to the atom type to nondimensionalize, set the elem
    to the atom type to use real parameters for otherwise. If nondim is not
    specified, default is no and default parameters for argon are used.
    Change elem to the species you want to use if not nondimensionalized;
    if not specified, the default is argon.
    """
    nl, nd, rij = findNN(pos,num)
    if nondim != 'no': #nondimensional
        eps = 1
        sig = 1
        mass = 1
    else:
        eps, sig, mass = ljND(elem)
    
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
                tmpf[j,:] = (48*eps/rd**2 * ((sig/rd)**12-0.5*(sig/rd)**6))*rv
                # for general expression, see written part of HW assignment
                tmpv[0,j] = ( 4*eps * ((sig/rd)**12 - (sig/rd)**6))
        fij[i,2:5] = tmpf.sum(axis=0) 
        vij[i,2] = tmpv.sum() 
        vTot = vij[:,2].sum()/2.0 #total potential, dont double count
            
    return [fij, vij, vTot]
#####################################
    


###################################
def vVerlet(num,pos,vels,fij,dt,nondim='no',elem='argon'):
    """
    Use velocity verlet algorithm to update positions and velocities
    given initial positions, velocities, forces, and the time step.
    
    returns positions, velocities, and forces at the next time step.
    Also returns total potential energy each time step.
    
    Call within a loop to run MD simulation.
    
    set nondim to the atom type to nondimensionalize, set the elem
    to the atom type to use real parameters for otherwise. If nondim is not
    specified, default is no and default parameters for argon are used.
    Change elem to the species you want to use if not nondimensionalized;
    if not specified, the default is argon.
    """
    
    ### nondimensionalize 
    if nondim != 'no': #nondimensional
        eps = 1
        sig = 1
        mass = 1
    else:
        eps, sig, mass = ljND(elem)
        
    ## Compute half step velocity
    for i in range(num): #loop over all atoms
        for j in range(3): #loop over x y z
            vels[i,j+2] = vels[i,j+2]+dt*fij[i,j+2]/2.0/mass
            #compute the velocity at the next half step for each atom in
            #x, y, and z
        
    ## Compute the next positions at full step from vhalf
    for i in range(num): #loop over all atoms
        for j in range(3): #loop over x y z
            pos[i,j+2] = pos[i,j+2]+dt*vels[i,j+2]
            #compute the position at the next full step for each atom in
            #x, y, and z based on velocity at half step
            
    ## Update the forces for the new positions
    fij, vij, vTot = fLJ(pos,num,nondim,elem)
    
    ## Compute full step velocity
    for i in range(num): #loop over all atoms
        for j in range(3): #loop over x y z
            vels[i,j+2] = vels[i,j+2]+dt*fij[i,j+2]/2.0/mass
            #compute the velocity at the next half step for each atom in
            #x, y, and z
            
    return [pos, vels, fij, vTot]
##################################
    


####################################
def dump(k,dump,num,pos,types):
    """
    dump is the frequency you want to print trajectory to the file
    'traj.xyz'
    
    The rest is behind the scenes.
    """
    if k+1 == dump: #if the first step, erase the file. 
        with open('traj.xyz','w') as fid: #write pos trajectory to file
            fid.write(str(num)+'\n')
            fid.write('Time Step: '+str(k+1)+'\n')
            for l in range(num-1):
                spec = str(types[int(pos[l,1])-1,0])
                fid.write(spec+'\t'+str(pos[l,2])+'\t'+
                          str(pos[l,3])+'\t'+str(pos[l,4])+'\n')
            spec = str(types[int(pos[-1,1])-1,0])
            fid.write(spec+'\t'+str(pos[-1,2])+'\t'+
                      str(pos[-1,3])+'\t'+str(pos[-1,4]))
    else:
        with open('traj.xyz','a') as fid: #write pos trajectory to file
            fid.write('\n'+str(num)+'\n')
            fid.write('Time Step: '+str(k+1)+'\n')
            for l in range(num-1):
                spec = str(types[int(pos[l,1])-1,0])
                fid.write(spec+'\t'+str(pos[l,2])+'\t'+
                          str(pos[l,3])+'\t'+str(pos[l,4])+'\n')
            spec = str(types[int(pos[-1,1])-1,0])
            fid.write(spec+'\t'+str(pos[-1,2])+'\t'+
                      str(pos[-1,3])+'\t'+str(pos[-1,4]))
#####################################################################
            
            
###############################################################
def thermo(k,thermo,vels,vTot,nondim,elem):
    """
    Computes KE, PE, Etot
        
    Writes to file 'log.MD' with frequency thermo
    """
        
    ### nondimensionalize 
    if nondim != 'no': #nondimensional
        eps = 1
        sig = 1
        mass = 1
    else:
        eps, sig, mass = ljND(elem) 
        
    ke = np.multiply(vels[:,2:5],vels[:,2:5])*mass/2.0 
    #assumes only 1 species
    ke = np.round((ke.sum(axis=1)).sum(),decimals=3) #total KE
    pe = np.round(vTot,decimals=3)
    etot = np.round(ke+pe,decimals=3)
    
    if (k+1)%thermo == 0:
        if k+1 == thermo: #if the first step, erase the file. 
            with open('log.MD','w') as fid: #write energy to file
                fid.write('STEP\t\tKE\t\tPE\t\tE\n')
                fid.write(str(k+1)+'\t\t'+str(ke)+'\t\t'+str(pe)+'\t\t'+
                          str(etot)+'\n')
        else:
            with open('log.MD','a') as fid: #append energy to file
                fid.write(str(k+1)+'\t\t'+str(ke)+'\t\t'+str(pe)+'\t\t'+
                          str(etot)+'\n')
                
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

def toc(logFlag='yes'):
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


