#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Just for plotting the potential and force. 

Only for visualization purposes
"""

import numpy as np
import matplotlib.pyplot as plt
import problem2fx as p

#### Problem 2 MD ###
#initial conditions
x0 = 0 #initial position
v0 = 2**.5 #initial velocity

tTot = 20
dt = 0.01

#run the MD simulation
x, v, t = p.runMD(x0,v0,dt,tTot,ftype='linear')

### Plot pos and vel ###
fig1, ax1 = plt.subplots()
ax1.plot(t,x,'r',t,v,'b')
ax1.axis([0,20,-2,2])
plt.xlabel('t')
plt.ylabel('x(t) and v(t)')
fig1.legend(('x(t)-MD','v(t)-MD'))
plt.title('MD Trajectory - Linear dt = 0.01 (Problem 2)')
plt.show()
#fig1.savefig('MDTrajectoryLinear_problem2.pdf',format='pdf',dpi=300)

## Plot pos and pot ###
ux = np.arange(-2,2,0.0001)
u1 = ux**2/2.0
e1 = np.ones(len(x))*v0**2/2.0 #energy for x0 = 0 

fig2, ax2 = plt.subplots()
ax2.plot(ux,u1,'g',x,e1,'r')
ax2.axis([-2,2,0,2])
plt.xlabel('x')
plt.ylabel('Pot, Pos')
fig2.legend(('U(x)','x(t)'))
plt.title('Linear Spring (Problem 2)')
plt.show()
#fig2.savefig('pot_vs_posProblem2-linear.pdf',format='pdf',dpi=300)

     
## Constant energy trajectory for linear spring ###
#eMan = np.zeros((len(x),len(x))) #energy manifold
#for i in range(len(x)):
#    for j in range(len(x)): #constant energy = 1
#        tmp = v[i]**2/2.0+x[j]**2/2.0
#        if tmp <= 1.05 and tmp >= 0.95:
#            eMan[i,j] = 1
#        else:
#            eMan[i,j] = 0
            
fig3, ax3 = plt.subplots()
ax3.scatter(x,v)
plt.xlabel('x(t)')
plt.ylabel('v(t)')
plt.title('Constant Energy Contour, SHO')
plt.show()
#fig3.savefig('energyContourSHO.pdf',format='pdf',dpi=300)
#####################################################

### Check the convergence of dt for linear spring using the MD code
x1, v1, t1 = p.runMD(x0,v0,5,tTot,ftype='linear')
x2, v2, t2 = p.runMD(x0,v0,1,tTot,ftype='linear')
x3, v3, t3 = p.runMD(x0,v0,0.5,tTot,ftype='linear')
x4, v4, t4 = p.runMD(x0,v0,0.1,tTot,ftype='linear')
x5, v5, t5 = p.runMD(x0,v0,0.05,tTot,ftype='linear')
x6, v6, t6 = p.runMD(x0,v0,0.01,tTot,ftype='linear')
x7, v7, t7 = p.runMD(x0,v0,0.005,tTot,ftype='linear')

fig4, ax4 = plt.subplots()
ax4.plot(t1,x1,color=(1,0,0))
ax4.plot(t2,x2,color=(.7,0,.3))
ax4.plot(t3,x3,color=(.5,0,.5))
ax4.plot(t4,x4,color=(.3,0,.7))
ax4.plot(t5,x5,color=(0,0,1))
ax4.plot(t6,x6,color=(0,.3,.7))
ax4.plot(t7,x7,color=(0,.5,.5))
ax4.axis([0,20,-2,2])
plt.xlabel('t')
plt.ylabel('x(t)')
fig4.legend(('5','1','.5','.1','.05','.01','.005'))
plt.title('dt convergence (Problem 2)')
plt.show()
#fig4.savefig('dtConvergence_problem2.pdf',format='pdf',dpi=300)
###############################################################

### Plot the analytical trajectories ###
x0 = 0
v0 = 2**.5
t = np.arange(0,20,0.0001)
x = x0*np.cos(t)+v0*np.sin(t)
v = -x0*np.sin(t)+v0*np.cos(t)

fig5, ax5 = plt.subplots()
ax5.plot(t,x,'r',t,v,'b')
ax5.axis([0,20,-2,2])
plt.xlabel('t')
plt.ylabel('x(t) and v(t)')
fig5.legend(('x(t)','v(t)'))
plt.title('Trajectory of Harmonic System (Problem 2)')
plt.show()
##fig5.savefig('trajectory_problem2b.pdf',format='pdf',dpi=300)
####
