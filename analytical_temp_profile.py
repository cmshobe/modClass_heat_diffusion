# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 08:30:00 2016

@author: Charlie

Week 1 assignment, part 1:
steady state linear geotherm, analytical sol'n only
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

##initialize
tsurf = -12 #degrees c
qm = 0.045 #w/m2
k = 2.5 #w/mk
z_min = 0 #m
z_max = 800 #m
dz = 0.01 #m
z = np.arange(z_min,z_max+dz,dz)

#for variable k:
k = np.zeros(len(z))
k[0:30/dz] = 1.2
k[30/dz:] = 2.5

##run
temp = np.zeros(len(z))
temp[0:30/dz] = tsurf + (qm / k[0:30/dz]) * z[0:30/dz]
tsurf_lower = temp[30/dz - dz]
z_lower = z[30/dz:] - z[30/dz-dz]
temp[30/dz:] = tsurf_lower + (qm / k[30/dz:]) * z_lower

##finalize
temp_prof = plt.figure(figsize=(3,3))
t = plt.subplot()
t.plot(temp, z, color='k', linewidth=3)
t.plot([0, 0], [z_min, z_max], color='b', linestyle='--', linewidth=3)
plt.gca().invert_yaxis()
plt.xlabel('Temperature [C]')
plt.ylabel('Depth [m]')
plt.show()
temp_prof.savefig('lin_prof_part1b.png', dpi=300, bbox_inches='tight')