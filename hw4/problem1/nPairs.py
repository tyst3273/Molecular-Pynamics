#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 14:59:19 2019

@author: ty
"""

import numpy as np

r = 1.0 #nm
rd = 10.0 #nm
r0 = rd-r #NM

rho = 33.4 #molecules per nm**3

c1 = 2*np.pi**2*rho
c2 = np.pi/3
p1 = 4/3*np.pi*r**3*r
p2 = r0**2*(r0+3*r)*r
p3 = (r0+r)*(rd**3-r0**3)
p4 = 3/2*r0*(r0+2*r)*(rd**2-r0**2)
p5 = 1/4*(rd**4-r0**4)

n = c1*(p1-c2*(p2+p3-p4-p5))

ntot = 4/3*np.pi*r0**3*rho+n

pairs = (9*ntot**2-3*ntot)/2