#!/usr/bin/env python3
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

#eps = 1.67e-21 #joules
#sig = 3.40 #ang
#mass = 6.63e-26 #kg
#kb = 1.38064852e-23 #J/K

eB = 0.002733 #eV
eN = 0.001723 #eV
eBN = np.sqrt(eB*eN) #eV
sB = 3.4003 #A
sN = 3.2177 #A 
sBN = np.sqrt(sB*sN) #A
kb = 8.6173303e-5 #eV/K
amu2kg = 1.66054e-27 
mB = 10.811*amu2kg #AMU to kg
mN = 14.0067*amu2kg #AMU to kg
kbSI = 1.38064852e-23 #J/K
j2eV = 6.241509e18

def printParams(num,steps,dt,tTot):
    """
    Simple function that writes input parameters to terminal
    """
    print('\n\tRunning MD simulation!\n')
    print('\tNo. of atoms:\t\t'+str(num)+'\t--')
    print('\tNo. of steps:\t\t'+str(steps)+'\t--')
    print('\tTime step:\t\t'+str(dt)+'\tps')
    print('\tTotal duration:\t\t'+str(tTot)+'\tps\n')
    print('\t------------------------------------\n')
    
###################################################################
def doPBC(rx,ry,rz,lx,ly,lz):
    """
    shift particle distance to impose periodic boundary conds
    """
    if rx >= lx:
        rx = rx-lx*2
    elif rx <= -lx:
        rx = rx+lx*2
    
    if ry >= ly:
        ry = ry-ly*2
    elif ry <= -ly:
        ry = ry+ly*2
            
    if rz >= lz:
        rz = rz-lz*2
    elif rz <= -lz:
        rz = rz+lz*2
        
    return [rx ,ry, rz]
######################################################################
    
############################################################
def pbcCoords(num,pos,bounds,box):
   """
   Impose periodic boundary conditions to atom coordinates
   """
   lx = box[0] #distance to shift for pbc
   ly = box[1]
   lz = box[2]
   
   for i in range(num):
      
      if pos[i,2] <= bounds[0,0]:
         pos[i,2] = pos[i,2]+lx
      elif pos[i,2] > bounds[0,1]:
         pos[i,2] = pos[i,2]-lx
         
      if pos[i,3] <= bounds[1,0]:
         pos[i,3] = pos[i,3]+ly
      elif pos[i,3] > bounds[1,1]:
         pos[i,3] = pos[i,3]-ly
         
      if pos[i,4] <= bounds[2,0]:
         pos[i,4] = pos[i,4]+lz
      elif pos[i,4] > bounds[2,1]:
         pos[i,4] = pos[i,4]-lz
         
   return pos
   
#################################################################

####################################
def findNN(pos,num,box):
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

    lx = box[0]/2 #distance to shift for pbc
    ly = box[1]/2
    lz = box[2]/2
    
    for i in range(num): #loop over particles
        tmp = np.zeros((2,num)) # ids, distance
        for j in range(num): #loop over neighbors
            rx = pos[j,2]-pos[i,2] # vector from partice i to j
            ry = pos[j,3]-pos[i,3] # vector from partice i to j
            rz = pos[j,4]-pos[i,4] # vector from partice i to j
            rx, ry, rz = doPBC(rx,ry,rz,lx,ly,lz)
            dist = np.sqrt(rx**2+ry**2+rz**2)
            tmp[0,j] = pos[j,0]-1 #index
            tmp[1,j] = dist #distance
            
            rij[i,j,0] = -rx
            rij[i,j,1] = -ry
            rij[i,j,2] = -rz
        
        nl[i,:] = tmp[0,:]
        nd[i,:] = tmp[1,:]
        
    return [nl, nd, rij] #return the neighbor lists
####################################

#######################################
def findRij(pos,vlist,box):
    """
    This function calculates rij for all atoms in the verlet list
    """
    lx = box[0]/2 #distance to shift for pbc
    ly = box[1]/2
    lz = box[2]/2
        
    numNN = len(vlist)
    rij = np.zeros((numNN,3))
    nd = np.zeros(numNN)
    i = int(vlist[0,0])
    for j in range(numNN): #loop over neighbors
        if j == 0:
            rij[j,:] = 0
            dist = 0
        else:
            k = int(vlist[j,0])
            rx = pos[k,2]-pos[i,2] # vector from partice i to j
            ry = pos[k,3]-pos[i,3] # vector from partice i to j
            rz = pos[k,4]-pos[i,4] # vector from partice i to j
            rx, ry, rz = doPBC(rx,ry,rz,lx,ly,lz)
            dist = np.sqrt(rx**2+ry**2+rz**2)
            
            rij[j,0] = -rx
            rij[j,1] = -ry
            rij[j,2] = -rz
    
        nd[j] = dist
        
    #nondimensional    
    return [nd, rij] #return the neighbor lists
#######################################

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
            
        x1 = np.random.uniform(0,1,int(num*3/2))
        x2 = np.random.uniform(0,1,int(num*3/2))
        y1 = (-2*np.log(x1))**.5*np.cos(2*np.pi*x2)
        y2 = (-2*np.log(x1))**.5*np.sin(2*np.pi*x2)
        num = len(pos[:,0])
        vels[:,2:5] = np.reshape(np.append(y1,y2),(num,3))[0:num,:]
        
        cmv = np.zeros(3)
        mass = 0
        for i in range(num):
            if pos[i,1] == 1:
                vels[i,2:5] = vels[i,2:5]*np.sqrt(kbSI*val/mB) #m/s
                cmv = cmv+vels[i,2:5]*mB
                mass = mass+mB
            else:
                vels[i,2:5] = vels[i,2:5]*np.sqrt(kbSI*val/mN) #m/s
                cmv = cmv+vels[i,2:5]*mN
                mass = mass+mN
                
        cmv = cmv/mass
                
        cmv = (mass*vels[:,2:5]).sum(axis=0)/(mass*num)
        vels[:,2:5] = vels[:,2:5]-cmv #zero the center of mass velocity
        
    else:
        sys.exit('\n\tUSAGE ERROR: Give a valid velocity distribution\n')
    
    vels[:,2:5] = vels[:,2:5]/100 #angstom/ps
    
    
    ### Print temp for debugging
    num = len(vels[:,0])
    vels2 = cp.deepcopy(vels)*100 #m/s
    ke = 0
    for i in range(num):
       if vels2[i,1] == 1:
           mass = mB
       else:
           mass = mN
       ke = ke+(vels2[i,2]**2+vels2[i,3]**2+vels2[i,4]**2)*mass/2
       
    temp = 3*ke/num/kbSI #unit kb = 1
    tempR = np.round(temp,decimals=3)
    keR = np.round(j2eV*ke,decimals=3)
    print('\t\tKE = '+str(keR)+'\t\tT = '+str(tempR))    
    return vels
##################

####################
def fLJ_SF(pos,num,rcut,vlist,box):
    """
    This function computes the shifted-force version of the interparticle 
    pairwise force from the LJ potential as well as the interparticle and 
    total potential energy.
    
    How it works:
    it calls findRij to calculate the neighbor distances and the rij 
    vectors needed to compute the interparticle forces for all atoms in the 
    verlet list. it then loops over all particles, and all other particles f
    or each, to calculate the forces and potential.
    """
    vcBB = 4*eB*((sB/rcut)**12-(sB/rcut)**6) #potential at cutoff
    dvdrBB = -48*eB*((sB/rcut)**13-0.5*(sB/rcut**7)) #slope of potential at cutoff
    fcBB = (48*eB/rcut**2*((sB/rcut)**12-0.5*(sB/rcut)**6)) #force at cutoff
    
    vcNN = 4*eN*((sN/rcut)**12-(sN/rcut)**6) #potential at cutoff
    dvdrNN = -48*eN*((sN/rcut)**13-0.5*(sN/rcut**7)) #slope of potential at cutoff
    fcNN = (48*eN/rcut**2*((sN/rcut)**12-0.5*(sN/rcut)**6)) #force at cutoff
    
    vcBN = 4*eBN*((sBN/rcut)**12-(sBN/rcut)**6) #potential at cutoff
    dvdrBN = -48*eBN*((sBN/rcut)**13-0.5*(sBN/rcut**7)) #slope of potential at cutoff
    fcBN = (48*eBN/rcut**2*((sBN/rcut)**12-0.5*(sBN/rcut)**6)) #force at cutoff

    fij = cp.deepcopy(pos) #to copy ids, types
    vij = cp.deepcopy(pos[:,0:3]) #to copy ids, types
    fij[:,2:5] = 0 #forces, same layout as pos and vels
    vij[:,2] = 0 #potential, same layout as pos and vels
    
    for i in range(num):
        nd, rij = findRij(pos,vlist[i],box)
        numNN = len(vlist[i])
        tmpf = np.zeros((numNN,3)) #tmp array to store forces
        tmpv = np.zeros((1,numNN)) #tmp array to store PE
        for j in range(numNN):
            
            if j == 0:
                tmpf[j,:] = 0
                tmpv[0,j] = 0
                
            elif pos[i,1] == 1 and pos[j,1] == 1:
                rv = rij[j,:]
                rd = nd[j]
                tmpf[j,:] = ((48*eB/rd**2*((sB/rd)**12-0.5*(sB/rd)**6))-fcBB)*rv
                tmpv[0,j] = (4*eB*((sB/rd)**12-(sB/rd)**6))-vcBB-(rd-rcut)*dvdrBB
                
            elif pos[i,1] == 1 and pos[j,1] == 2:
                rv = rij[j,:]
                rd = nd[j]
                tmpf[j,:] = ((48*eBN/rd**2*((sBN/rd)**12-0.5*(sBN/rd)**6))-fcBN)*rv
                tmpv[0,j] = (4*eBN*((sBN/rd)**12-(sBN/rd)**6))-vcBN-(rd-rcut)*dvdrBN
                
            elif pos[i,1] == 2 and pos[j,1] == 1:
                rv = rij[j,:]
                rd = nd[j]
                tmpf[j,:] = ((48*eBN/rd**2*((sBN/rd)**12-0.5*(sBN/rd)**6))-fcBN)*rv
                tmpv[0,j] = (4*eBN*((sBN/rd)**12-(sBN/rd)**6))-vcBN-(rd-rcut)*dvdrBN
                
            else:
                rv = rij[j,:]
                rd = nd[j]
                tmpf[j,:] = ((48*eN/rd**2*((sN/rd)**12-0.5*(sN/rd)**6))-fcNN)*rv
                tmpv[0,j] = (4*eN*((sN/rd)**12-(sN/rd)**6))-vcNN-(rd-rcut)*dvdrNN
                
        fij[i,2:5] = tmpf.sum(axis=0) 
        vij[i,2] = tmpv.sum() 
    vTot = vij[:,2].sum()/2.0 #total potential, dont double count
            
    #nondimensional
    return [fij, vij, vTot]
#####################################

###################################
def vVerlet(num,pos,rcut,vels,fij,vlist,dt,box):
    """
    Use velocity verlet algorithm to update positions and velocities
    given initial positions, velocities, forces, and the time step.
    
    returns positions, velocities, and forces at the next time step.
    Also returns total potential energy each time step.
    
    Call within a loop to run MD simulation.
    """
    c = 0.0001*1.60218e-19 #conver ps*eV/ang/kg to ang/ps

    ## Compute half step velocity
    for i in range(num): #loop over all atoms
        if pos[i,1] == 1:
           mass = mB
        else:
           mass = mN
        for j in range(3): #loop over x y z
            vels[i,j+2] = vels[i,j+2]+dt*c*fij[i,j+2]/2.0/mass
            #compute the velocity at the next half step for each atom in
            #x, y, and z
        
    ## Compute the next positions at full step from vhalf
    for i in range(num): #loop over all atoms
        for j in range(3): #loop over x y z
            pos[i,j+2] = pos[i,j+2]+dt*vels[i,j+2]
            #compute the position at the next full step for each atom in
            #x, y, and z based on velocity at half step
            
    ## Update the forces for the new positions
    fij, vij, vTot = fLJ_SF(pos,num,rcut,vlist,box)
    
    ## Compute full step velocity
    for i in range(num): #loop over all atoms
        if pos[i,1] == 1:
           mass = mB
        else:
           mass = mN
        for j in range(3): #loop over x y z
            vels[i,j+2] = vels[i,j+2]+dt*c*fij[i,j+2]/2.0/mass
            #compute the velocity at the next half step for each atom in
            #x, y, and z
            
    #nondimensional
    return [pos, vels, fij, vTot]
##############################################
    
##################################################
def verletList(num,pos,rcut,box):     
    """
    Generate verlet lists for force cutoffs. Impose periodic boundary 
    conditions if necessary.
    """
    nl, nd, rij = findNN(pos,num,box) #get NN to make neighborlist
    vlist = [0]*num #verlet list
    vcoord = [0]*num #verlet list
    for i in range(num): #find neighbor lists
        tmp = np.zeros((2,num))
        tmp[0,:] = nl[i,:] #look up nn list
        tmp[1,:] = nd[i,:] 
        
        tmp[:,:] = tmp[:,np.argsort(tmp[1,:])] #sort nn list
        nns = tmp[0,np.argwhere(tmp[1,:] <= rcut)] #check cutoff against distance
        
        vlist[i] = nns
        vcoord[i] = pos[nns[:,0].astype(int),2:5] #coord of atoms in verlet list
                
    return [vlist, vcoord]
###########################################################

#########################################################
def checkVerlet(num,pos,rcut,skin,vlist,vcoord,box):
    """
    Check if the current verlet list is still valid. If not, make new verlet
    lists.
    """
    for i in range(num):
        coords = pos[vlist[i][:,0].astype(int),2:5] #current coord of atoms
        maxdisp = np.abs(vcoord[i]-coords).max() #max displacement of any atom in list
        if maxdisp >= skin: #check if it's moved beyond skin distance
            vlist, vcoord = verletList(num,pos,rcut,box)
            break 
            
    return [vlist, vcoord] #valid verlet list
######################################
    
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
        
        box = np.array(fid.readline().strip().split()).astype(float) #box boundaries
        box = box.reshape(3,2)
        bounds = cp.deepcopy(box)
        box = box[:,1]-box[:,0]
        for i in range(num): #read in positions
            tmp = fid.readline().strip().split()
            types[i,0] = tmp[0]
            pos[i,:] = tmp[1:4]
            
    utypes = np.unique(types) 
    n = len(utypes) #number of species
    tmp = np.zeros((num,1))
    
    tmp[np.argwhere(types[:,0] == 'B')] = 1
    tmp[np.argwhere(types[:,0] == 'N')] = 2
      
    pos = np.append(tmp,pos,axis=1) #append types to pos
    types = utypes.reshape(n,1)
    pos = np.append(np.arange(1,num+1).reshape(num,1),pos,axis=1) 
    #particle ids
    
    return [num ,pos, types, box, bounds]
###################################

####################################
def dump(k,dump,num,pos,types):
    """
    dump is the frequency you want to print trajectory to the file
    'traj.xyz'
    
    The rest is behind the scenes.
    """
    cpos = cp.deepcopy(pos)
    
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
    num = len(vels[:,0])
    vels2 = cp.deepcopy(vels)*100 #m/s
    ke = 0
    mx = 0
    my = 0
    mz = 0
    for i in range(num):
       if vels2[i,1] == 1:
           mass = mB
       else:
           mass = mN
       ke = ke+(vels2[i,2]**2+vels2[i,3]**2+vels2[i,4]**2)*mass
       mx = mx+vels2[i,2]*mass
       my = my+vels2[i,3]*mass
       mz = mz+vels2[i,4]*mass
       
#    temp = eV2J*ke/3/num/kbSI 
    temp = ke/3/num/kbSI
    tempR = np.round(temp,decimals=3)
    keR = np.round(j2eV*ke,decimals=3)
    pe = np.round(vTot,decimals=3)
    etot = keR+pe
    mxR = np.round(mx,decimals=3)
    myR = np.round(my,decimals=3)
    mzR = np.round(mz,decimals=3)
    
    if (k+1)%thermo == 0:
        if k+1 == thermo: #if the first step, erase the file. 
            with open('log.MD','w') as fid: #write energy to file
                fid.write('STEP\tT\tKE\tPE\tE\tPx\tPy\tPz\n')
                fid.write('----\tK\teV\teV\teV\tkg.m/s\tkg.m/s'
                          '\tkg.m/s\n')
                fid.write(str(k+1)+'\t'+str(tempR)+'\t'+str(keR)+'\t'+
                          str(pe)+'\t'+str(etot)+'\t'+str(mxR)+'\t'+
                          str(myR)+'\t'+str(mzR)+'\n')
        else:
            with open('log.MD','a') as fid: #append energy to file
                fid.write(str(k+1)+'\t'+str(tempR)+'\t'+str(keR)+'\t'+
                          str(pe)+'\t'+str(etot)+'\t'+str(mxR)+'\t'+
                          str(myR)+'\t'+str(mzR)+'\n')
    return ke, temp, mx, my, mz
########################################################################
                
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
        print(("\t\tElapsed time is "+
              str(np.round(time.time()-
                       startTime_for_tictoc,decimals=3))+" seconds.\n"))
    else:
        print("\n\t\tToc: start time not set")
################################################
              
###################################################
def readLog(logfile):
    """
    Read the output log file from run into data structures
    """
    num = sum(1 for line in open(logfile,'r'))-2
    
    with open(logfile,'r') as fid:
        col = fid.readline().strip().split()
        ncol = len(col)
        data =  np.zeros((num,ncol))
        fid.readline()
        
        for i in range(num):
            data[i,:] = fid.readline().strip().split()
            
    return data
            
            
