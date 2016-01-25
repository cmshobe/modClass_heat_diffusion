# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 08:57:29 2016

@author: Charlie

Week 1 assignment, part 2:
transient sinusoidal geotherm, analytical sol'n only
with animations and such
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
plt.cla()
## initialize
plot_freq = 100
Ts_bar = -12 #celsius
amp_annual = 15 #amplitude of annual fluc
amp_day = 10 #amp of daily fluc
qm = 0.045 #w/m2
k = 2.5 #w/mk
kappa = 1e-6 #m2/s
z_min = 0
z_max = 1 #meters
dz = 0.01 #m
z = np.arange(z_min, z_max+dz, dz)
dt = (1/4) * 3600 * 24 #days * secs/day
dt_day = 10 #seconds
period = 365 * 3600 * 24 #secs/yr
period_day = 3600 * 24;
t_min = 0
t_max = period 
t = np.arange(t_min, t_max+dt, dt)

z_star_annual = np.sqrt(kappa * period / np.pi)
z_star_daily = np.sqrt(kappa * 3600 * 24 / np.pi)

Tenv_hot=Ts_bar+amp_annual*np.exp(-z/z_star_annual)+amp_day*np.exp(-z/z_star_daily)
Tenv_cold=Ts_bar-amp_annual*np.exp(-z/z_star_annual)-amp_day*np.exp(-z/z_star_daily)

temp_memory_1 = np.zeros(len(t))
temp_memory_2 = np.zeros(len(t))
temp_memory_3 = np.zeros(len(t))
temp_memory_4 = np.zeros(len(t))
temp_memory_5 = np.zeros(len(t))
temp_memory_6 = np.zeros(len(t))
temp_memory_7 = np.zeros(len(t))
temp_memory_8 = np.zeros(len(t))
temp_memory_9 = np.zeros(len(t))
temp_memory_10 = np.zeros(len(t))

##run
temp_prof = plt.figure(figsize=(3,3))
tplot = plt.subplot()
tplot.plot(Tenv_hot, z, color='r', linewidth=2)
tplot.plot(Tenv_cold, z, color='b', linewidth=2)
tplot.plot([0, 0], [z_min, z_max], color='k', linestyle='--', linewidth=2)
plt.gca().invert_yaxis()
plt.ylim([1,0])
plt.xlabel('Temperature [C]')
plt.ylabel('Depth [m]')
iter_count = 0
for i in range(len(t)):
    iter_count += 1
    temp_annual = Ts_bar + amp_annual * np.exp(-z/z_star_annual) * np.sin((2 * np.pi * t[i]) / period - (z/z_star_annual))
    temp = temp_annual + amp_day* np.exp(-z/z_star_daily) * np.sin((2 * np.pi * t[i]/period_day) - (z/z_star_daily))
    temp_memory_1[i] = temp[.1/dz]
    temp_memory_2[i] = temp[.2/dz]
    temp_memory_3[i] = temp[.3/dz]
    temp_memory_4[i] = temp[.4/dz]
    temp_memory_5[i] = temp[.5/dz]
    temp_memory_6[i] = temp[.6/dz]
    temp_memory_7[i] = temp[.7/dz]
    temp_memory_8[i] = temp[.8/dz]
    temp_memory_9[i] = temp[.9/dz]
    temp_memory_10[i] = temp[1/dz]  
    if iter_count == 1 or np.remainder(iter_count, plot_freq) == 0:
        tplot.plot(temp, z, linewidth=0.5)   
    else: 
        pass
plt.show()
temp_prof.savefig('sin_prof_part2_withdaily.png', dpi=300, bbox_inches='tight')
t_years = t / 365 / 24 / 3600
#finalize, plot temperatures at depth over time
temps_at_depth = plt.figure(figsize=(10,10))
temps_at_depth.subplots_adjust(bottom=0.1)
t1 = plt.subplot(431)
plt.ylim([-25,0])
plt.ylabel('Temperature [C]')
plt.title('0.1 m depth')
plt.setp( t1.get_xticklabels(), visible=False)
#temp_at_2 = plt.figure(figsize=(10,10))
t2 = plt.subplot(432)
plt.ylim([-25,0])
plt.setp( t2.get_yticklabels(), visible=False)
plt.setp( t2.get_xticklabels(), visible=False)
plt.title('0.2 m depth')
#temp_at_3 = plt.figure(figsize=(10,10))
t3 = plt.subplot(433)
plt.ylim([-25,0])
plt.setp( t3.get_yticklabels(), visible=False)
plt.setp( t3.get_xticklabels(), visible=False)
plt.title('0.3 m depth')
#temp_at_4 = plt.figure(figsize=(10,10))
t4 = plt.subplot(434)
plt.ylim([-25,0])
plt.ylabel('Temperature [C]')
plt.title('0.4 m depth')
plt.setp( t4.get_xticklabels(), visible=False)
#temp_at_5 = plt.figure(figsize=(10,10))
t5 = plt.subplot(435)
plt.ylim([-25,0])
plt.setp( t5.get_yticklabels(), visible=False)
plt.setp( t5.get_xticklabels(), visible=False)
plt.title('0.5 m depth')
#temp_at_6 = plt.figure(figsize=(10,10))
t6 = plt.subplot(436)
plt.ylim([-25,0])
plt.setp( t6.get_yticklabels(), visible=False)
plt.setp( t6.get_xticklabels(), visible=False)
plt.title('0.6 m depth')
#temp_at_7 = plt.figure(figsize=(10,10))
t7 = plt.subplot(437)
plt.ylim([-25,0])
plt.ylabel('Temperature [C]')
plt.title('0.7 m depth')
plt.setp( t7.get_xticklabels(), visible=False)
#temp_at_8 = plt.figure(figsize=(10,10))
t8 = plt.subplot(438)
plt.ylim([-25,0])
plt.setp( t8.get_yticklabels(), visible=False)
plt.xlabel('Time [yr]')
plt.title('0.8 m depth')
#temp_at_9 = plt.figure(figsize=(10,10))
t9 = plt.subplot(439)
plt.ylim([-25,0])
plt.xlabel('Time [yr]')
plt.setp( t9.get_yticklabels(), visible=False)
plt.title('0.9 m depth')
#temp_at_10 = plt.figure(figsize=(10,10))
t10 = plt.subplot(4,3,10)
plt.ylim([-25,0])
plt.xlabel('Time [yr]')
plt.ylabel('Temperature [C]')
plt.title('1 m depth')

t1.plot(t_years, temp_memory_1, color='r')
t2.plot(t_years, temp_memory_2, color='r')
t3.plot(t_years, temp_memory_3, color='r')
t4.plot(t_years, temp_memory_4, color='r')
t5.plot(t_years, temp_memory_5, color='r')
t6.plot(t_years, temp_memory_6, color='r')
t7.plot(t_years, temp_memory_7, color='r')
t8.plot(t_years, temp_memory_8, color='r')
t9.plot(t_years, temp_memory_9, color='r')
t10.plot(t_years, temp_memory_10, color='r')

plt.show()
temps_at_depth.savefig('sin_prof_part2_depths_withdaily.png', dpi=300, bbox_inches='tight')