#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 15:38:41 2019

@author: ty
"""

def flinear(xi):
    """
    Force only depends on position
    """
    f = -xi
    return f

def fNonlinear(xi):
    """
    Force from the nonlinear spring
    """
    f = -4*xi**3+4*xi
    return f

def vhalf(vi,fi,dt):
    """
    Compute velocity for next half step
    """
    vh = vi+dt*fi/2.0
    return vh

def rplust(ri,vi,dt):
    """
    compute position at full timestep
    """
    rp = ri+dt*vi
    return rp

def vplust(vi,fi,dt):
    """
    compute velocities at the next full step
    Note: have to compute new forces before calling this
    """
    vp = vi+dt*fi/2.0
    return vp

def runMD(x0,v0,dt,tTot,ftype):
    """
    Run the MD simulation for a single particle on a 
    spring
    
    ALL VALUES ARE UNITLESS
    
    x0 and v0 are the initial position and velocity respectively
    
    dt is the times step and tTot is the total simulation time
    
    ftype is the force type; ftype = 'linear' for the linear spring
    model and ftype = 'nonlinear' for the other
    """
    import numpy as np
    import sys 
    
    if ftype != 'linear' and ftype != 'nonlinear':
        sys.exit('\n\tftype must be \'linear\' or \'nonlinear\' \n')
    
    n = int(tTot/dt) #total numbe of steps
    t = np.arange(0,tTot,dt) #time array 
    vh = np.zeros(n) #half-step velocities
    v = np.zeros(n) #velocity as fx of time
    x = np.zeros(n) #pos as fx of time
    f = np.zeros(n) #force as fx of time
    if ftype == 'linear':
        f0 = flinear(x0) #initial force
    else:
        f0 = fNonlinear(x0)
    for i in range(n):
        if i == 0:
            vh[i] = vhalf(v0,f0,dt) #compute half step velocity 
            x[i] = rplust(x0,v0,dt) #compute next position
            if ftype == 'linear': 
                f[i] = flinear(x[0])
            else:
                f[i] = fNonlinear(x[0])
#            f[i] = flinear(x[0])
            v[i] = vplust(vh[0],f[0],dt) #compute next full step velocity
        else:
            vh[i] = vhalf(v[i-1],f[i-1],dt)
            x[i] = rplust(x[i-1],vh[i],dt)
            if ftype == 'linear': 
                f[i] = flinear(x[i])
            else:
                f[i] = fNonlinear(x[i])
#            f[i] = flinear(x[i])
            v[i] = vplust(vh[i],f[i],dt)
    return [x, v, t]