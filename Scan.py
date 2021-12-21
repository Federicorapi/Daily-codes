#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 19:54:18 2021

@author: federicorapisarda
"""

## Useful commands to use in console:
# pdb.set_trace()
# %matplotlib qt to show plots in separate windows
# %matplotlib inline to show plots in separate windows




import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import math
import pdb
#Open h5 file
f=h5.File('Prova per GUI.h5','r')
wvl=np.array(f['/WL'])
eV=1240/wvl
data=np.array(f['/Data'])
xy=np.array(f['/Consigne'])
pixel=len(wvl);
step_volt=xy[1,0]-xy[0,0];
scansizex=((xy[-1,0]-xy[0,0])/step_volt);       # Insert number of steps along x
scansizey=((xy[-1,1]-xy[0,1])/step_volt);
scansizex=math.ceil(scansizex)
scansizey=math.ceil(scansizey)
scansizex = scansizex+1
scansizey = scansizey+1
## Making the map
intensities=np.ones((scansizex,scansizey))
background=np.zeros((scansizex,scansizey))

for y in range(scansizey):
        # pdb.set_trace()
 for x in range(scansizex):   
     intensities[x,y] = np.amax(data[x+y*scansizex,:]) 
     
     ## Muon killing
treshold=1500;
for ss in range(scansizey):
    for kk in range(scansizex):
        if np.amax(intensities[kk,ss])>treshold:
            if kk != 1 and kk != scansizex and ss != 1 and ss != scansizey:
                intensities[kk,ss]=np.mean([intensities[kk-1,ss-1],intensities[kk,ss-1],intensities[kk+1,ss-1],intensities[kk-1,ss],intensities[kk+1,ss],intensities[kk-1,ss+1],intensities[kk,ss+1],intensities[kk+1,ss+1]])
            elif ss != 1 and ss != scansizey and (kk == 1 or kk == scansizex):
              intensities[kk,ss] = np.mean([intensities[kk,ss-1],intensities[kk,ss+1]])
            elif kk != 1 and kk != scansizex and (ss == 1 or ss == scansizey):
              intensities[kk,ss] = np.mean([intensities[kk-1,ss],intensities[kk+1,ss]])

fig1 = plt.figure(1)    
ax1 = fig1.gca()
p=plt.imshow(intensities,cmap = plt.cm.ocean) # "ocean" map is interesting too
fig1.colorbar(p)

    