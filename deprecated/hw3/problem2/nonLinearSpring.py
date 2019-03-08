#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 16:59:15 2019

@author: ty
"""
import numpy as np
import matplotlib.pyplot as plt
import problem2fx as p

#### Problem 2 nonlinear spring MD ###
tTot = 20
dt = 0.01

x = np.zeros((15,int(tTot/dt)))
v = np.zeros((15,int(tTot/dt)))

######### Total energy = .25 ###
# KE = PE = 0.125 
#initial conditions 1
x0 = .804019 #initial position
v0 = 0.5 #initial velocity
x[0,:], v[0,:], t = p.runMD(x0,v0,dt,tTot,ftype='nonlinear')

#initial conditions 2
x0 = -.804019 #initial position
v0 = 0.5 #initial velocity
x[1,:], v[1,:], t = p.runMD(x0,v0,dt,tTot,ftype='nonlinear')

###################
###################
# KE = 0.25, PE = 0
#initial conditions 3
x0 = -1 #initial position
v0 = 2**.5/2.0 #initial velocity 
x[2,:], v[2,:], t = p.runMD(x0,v0,dt,tTot,ftype='nonlinear')

#####################
######################
# KE = 0, PE = 0.25
#initial conditions 4
x0 = 2**.5/2.0 #initial position
v0 = 0 #initial velocity
x[3,:], v[3,:], t = p.runMD(x0,v0,dt,tTot,ftype='nonlinear')

#initial conditions 5
x0 = -2**.5/2.0 #initial position
v0 = 0 #initial velocity
x[4,:], v[4,:], t = p.runMD(x0,v0,dt,tTot,ftype='nonlinear')

#fig2, ax2 = plt.subplots()
#ax2.plot(t,x[0,:],'r',t,v[0,:],'r:')
#ax2.plot(t,x[2,:],'b',t,v[2,:],'b:')
#ax2.plot(t,x[4,:],'g',t,v[4,:],'g:')
#plt.xlabel('t')
#plt.ylabel('v(t) and x(t)')
#plt.title('E-tot = 0.25, dt = 0.01')
#fig2.legend(('x: x0 = 0.8, v0 = 0.5','v: x0 = -1, v0 = 0.707',
#             'x: x0 = -.707, v0 = 0','v: x0 = -.707, v0 = 0'))
#plt.show()
#fig2.savefig('eTot_dot25.pdf',format='pdf',dpi=300)
########################


######### Total energy = 1 ###
# KE = PE = 0.5
#initial conditions 6
x0 = 0.541196 #initial position
v0 = 1 #initial velocity
x[5,:], v[5,:], t = p.runMD(x0,v0,dt,tTot,ftype='nonlinear')

#initial conditions 7
x0 = -0.541196 #initial position
v0 = 1 #initial velocity
x[6,:], v[6,:], t = p.runMD(x0,v0,dt,tTot,ftype='nonlinear')

###################
###################
# KE = 1, PE = 0
#initial conditions 8
x0 = 1 #initial position
v0 = 2**.5 #initial velocity 
x[7,:], v[7,:], t = p.runMD(x0,v0,dt,tTot,ftype='nonlinear')

#####################
######################
# KE = 0, PE = 1
#initial conditions 9
x0 = 0 #initial position
v0 = 0 #initial velocity
x[8,:], v[8,:], t = p.runMD(x0,v0,dt,tTot,ftype='nonlinear')

#initial conditions 10
x0 = 2**.5 #initial position
v0 = 0 #initial velocity
x[9,:], v[9,:], t = p.runMD(x0,v0,dt,tTot,ftype='nonlinear')

#initial conditions 11
x0 = -2**.5 #initial position
v0 = 0 #initial velocity
x[10,:], v[10,:], t = p.runMD(x0,v0,dt,tTot,ftype='nonlinear')

#fig3, ax3 = plt.subplots()
#ax3.plot(t,x[5,:],'r',t,v[5,:],'r:')
#ax3.plot(t,x[7,:],'b',t,v[7,:],'b:')
#ax3.plot(t,x[8,:],'k',t,v[8,:],'k:')
#ax3.plot(t,x[10,:],'g',t,v[10,:],'g:')
#plt.xlabel('t')
#plt.ylabel('v(t) and x(t)')
#plt.title('E-Tot = 1, dt = 0.01')
#fig3.legend(('x: x0 = 0.542, v0 = 1','v: x0 = 1, v0 = 1.414',
#             'x: x0 = 1, v0 = 1.414','v: x0 = 1, v0 = 1.414',
#             'x: x0 = 0, v0 = 0','v: x0 = 0, v0 = 0',
#             'x: x0 = -1.141, v0 = 0','v: x0 = -1.414, v0 = 0'))
#plt.show()
#fig3.savefig('eTot_1.pdf',format='pdf',dpi=300)
########################


######### Total energy = 2 ###
# KE = PE = 1
#initial conditions 12
x0 = 2**.5 #initial position
v0 = 2**.5 #initial velocity
x[11,:], v[11,:], t = p.runMD(x0,v0,dt,tTot,ftype='nonlinear')

#initial conditions 13
x0 = -2**.5 #initial position
v0 = 2**.5 #initial velocity
x[12,:], v[12,:], t = p.runMD(x0,v0,dt,tTot,ftype='nonlinear')

###################
###################
# KE = 2, PE = 0
# 1 solutions for x
#initial conditions 14
x0 = -1 #initial position
v0 = 2 #initial velocity 
x[13,:], v[13,:], t = p.runMD(x0,v0,dt,tTot,ftype='nonlinear')

######################
# KE = 0, PE = 2
# 1 solutions for x
#initial conditions 15
x0 = (1+2**.5)**.5 #initial position
v0 = 0 #initial velocity
x[14,:], v[14,:], t = p.runMD(x0,v0,dt,tTot,ftype='nonlinear')

#fig4, ax4 = plt.subplots()
#ax4.plot(t,x[11,:],'r',t,v[11,:],'r:')
#ax4.plot(t,x[13,:],'b',t,v[13,:],'b:')
#ax4.plot(t,x[14,:],'g',t,v[14,:],'g:')
#plt.xlabel('t')
#plt.ylabel('v(t) and x(t)')
#plt.title('E-tot = 2, dt = 0.01')
#fig4.legend(('x: x0 = 1.414, v0 = 1.414','v: x0 = 1.414, v0 = 1.414',
#             'x: x0 = -1, v0 = 2','v: x0 = -1, v0 = 2',
#             'x: x0 = 1.554, v0 = 0','v: x0 = 1.554, v0 = 0'))
#plt.show()
#fig4.savefig('eTot_2.pdf',format='pdf',dpi=300)
########################


### Plot x vs Pot surface ###
n = len(x[0,:])
ux = np.arange(-2,2,0.0001)
u2 = ux**4-2*ux**2+1
e1 = np.ones(n)*.25
e2 = np.ones(n)*1
e3 = np.ones(n)*2

#fig5, ax5 = plt.subplots()
#ax5.plot(x[0,:],e1,'b',x[4,:],e1,'b')
#ax5.plot(x[5,:],e2,'r',x[10,:],e2,'r:',x[8,:],e2,'k*')
#ax5.plot(x[11,:],e3,'g')
#ax5.plot(ux,u2,'m')
#ax5.axis([-2,2,0,2.5])
#plt.xlabel('x')
#plt.ylabel('U(x) and x(t)')
#plt.title('Pos vs. Nonlinear Pot')
#fig5.legend(('E=0.25','E=0.25','E=1','E=1','E=1','E=2','U(x)'))
#plt.show()
#fig5.savefig('x_vs_nonlinearPot.pdf',format='pdf',dpi=300)


#### Plot cont. energy surfaces ####

        
fig6, ax6 = plt.subplots()
ax6.scatter(x[0,:],v[0,:],color='blue',s=1)
ax6.scatter(x[4,:],v[4,:],color='blue',s=1)
ax6.scatter(x[5,:],v[5,:],color='red',s=1)
ax6.scatter(x[8,:],v[8,:],color='black',s=3)
ax6.scatter(x[10,:],v[10,:],color='green',s=1)

plt.xlabel('x(t)')
plt.ylabel('v(t)')
plt.title('Const. E contours, nonlinear')
fig6.legend(('E=0.25','E=0.25','E=1','E=1','E=2'))
ax6.axis([-0.5,0.5,-0.5,0.5])
fig6.savefig('nonLinear_contoursZOOM.pdf',format='pdf',dpi=300)
ax6.axis([-1.5,1.5,-1.5,1.5])
plt.show()
fig6.savefig('nonLinear_contours.pdf',format='pdf',dpi=300)



#### Check the convergence of dt for nonlinear spring using the MD codes
x4, v4, t4 = p.runMD(x0,v0,0.1,tTot,ftype='nonlinear')
x5, v5, t5 = p.runMD(x0,v0,0.05,tTot,ftype='nonlinear')
x6, v6, t6 = p.runMD(x0,v0,0.01,tTot,ftype='nonlinear')
x7, v7, t7 = p.runMD(x0,v0,0.005,tTot,ftype='nonlinear')
x8, v8, t8 = p.runMD(x0,v0,0.001,tTot,ftype='nonlinear')
x9, v9, t9 = p.runMD(x0,v0,0.0005,tTot,ftype='nonlinear')

#fig1, ax1 = plt.subplots()
#ax1.plot(t4,x4,color=(.3,0,.7))
#ax1.plot(t5,x5,color=(0,0,1))
#ax1.plot(t6,x6,color=(0,.3,.7))
#ax1.plot(t7,x7,color=(0,.5,.5))
#ax1.plot(t8,x8,color=(.7,0,.3))
#ax1.plot(t9,x9,color=(.5,0,.5))
#ax1.axis([0,20,-2,2])
#plt.xlabel('t')
#plt.ylabel('x(t)')
#fig1.legend(('.1','.05','.01','.005','.001','.0005'))
#plt.title('dt convergence - Nonlinear(Problem 2)')
#plt.show()
##fig1.savefig('dtConvNonlinear_problem2.pdf',format='pdf',dpi=300)
##############################################################