#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 10:21:12 2019

@author: ty
"""

import numpy as np
import copy as cp

nx = int(input('Enter number of unit cells in x:\n\t>'))
ny = int(input('Enter numer of unit cells in y:\n\t>'))
xyz = str(input('Enter the name of the xyz data file or "none" if '
                       'you dont want one\n\t>'))

ax = 2.510 #2*1.255
ay = 4.348 # 6*0.724574

basis = np.array([[1,   0,   0,   0],
                  [2,   1,   1,   0],
                  [1,   1,   3,   0],
                  [2,   0,   4,   0]])

pos = cp.deepcopy(basis)
tmp = cp.deepcopy(basis)
for i in range(nx-1):
    tmp[:,1] = tmp[:,1]+2
    pos = np.append(pos,tmp,axis=0)
    
tmp = cp.deepcopy(pos)
for i in range(ny-1):
    tmp[:,2] = tmp[:,2]+6
    pos = np.append(pos,tmp,axis=0)
    
size = nx*ny*4
pos = pos[np.lexsort((pos[:,3],pos[:,2],pos[:,1]))]
pos = np.append(np.arange(1,size+1,1).reshape(size,1),pos,axis=1) #ids

pos = pos.astype(float)
pos[:,2] = pos[:,2]*ax/2 
pos[:,3] = pos[:,3]*ay/6

bx = ax/2.0
by = ay/6
xmin = pos[:,2].min()-bx
xmax = pos[:,2].max()+bx
ymin = pos[:,3].min()-by
ymax = pos[:,3].max()+by
zmin = pos[:,4].min()-1000
zmax = pos[:,4].max()+1000
       
if xyz != 'none':
    with open(xyz,'w') as fid:
        fid.write(str(size)+'\n')
        fid.write(str(xmin)+' '+str(xmax)+' '+str(ymin)+' '+str(ymax)+' '+
                  str(zmin)+' '+str(zmax)+'\n')
        for i in range(size-1):
            if pos[i,1] == 1:
                fid.write('B ')
            else: 
                fid.write('N ')
            fid.write(str(pos[i,2])+' '+str(pos[i,3])+' '+str(pos[i,4])+'\n')
        if pos[-1,1] == 1:
            fid.write('B ')
        else: 
            fid.write('N ')
        fid.write(str(pos[-1,2])+' '+str(pos[-1,3])+' '+str(pos[-1,4]))

