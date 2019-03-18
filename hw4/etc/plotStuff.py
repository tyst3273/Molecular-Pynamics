#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 11:28:21 2019

@author: ty
"""

import functionsLJSF as md
import matplotlib.pyplot as plt

#dat1 = md.readLog('fine.log.MD') #print thermo info every 2 steps for 100
#dat1 = md.readLog('full.log.MD') #full run
dat1 = md.readLog('equil.log.MD') #full run

fig1, ax1 = plt.subplots()
ax1.plot(dat1[:,0],dat1[:,1],'r')
plt.xlabel('step')
plt.ylabel('T (K)')
plt.title('Constant Temperature')
plt.axis([0,1000,50,150])
plt.show()
#plt.savefig('equilibratedTemp.pdf',format='pdf',dpi=300)

fig2, ax2 = plt.subplots()
ax2.plot(dat1[:,0],dat1[:,2],'r',dat1[:,0],dat1[:,3],'b',dat1[:,0],dat1[:,4],'m')
plt.xlabel('step')
plt.ylabel('E (eV)')
plt.legend(('Kinetic Energy','Potential Energy','Total Energy'))
plt.axis([0,1000,-15,10])
plt.title('Energy Conservation')
plt.show()
#plt.savefig('equilibratedEnergy.pdf',format='pdf',dpi=300)

fig3, ax3 = plt.subplots()
ax3.plot(dat1[:,0],dat1[:,5],'r',dat1[:,0],dat1[:,6],'b',dat1[:,0],dat1[:,7],'m')
plt.xlabel('step')
plt.ylabel('P (pg*A/ps)')
plt.legend(('Px','Py','Pz'))
plt.axis([0,1000,-0.5,0.5])
plt.title('Total Momentum')
plt.show()
#plt.savefig('equilibratedMomentum.pdf',format='pdf',dpi=300)