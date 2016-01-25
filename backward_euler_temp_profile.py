# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 16:38:47 2016

@author: Charlie Shobe

This is a finite difference solution to a thermal temperature profile.
Euler backward solution
"""
from __future__ import division
import numpy as np
from pyThomas import thomas #Thomas algorithm efficiently solves linear system
import matplotlib.pyplot as plt
import matplotlib

#plotting parameters:
matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams['xtick.major.pad']='10'

#INITIALIZE################################################
#domain setup:
ts_bar = -5 #mean surface temp [C]
temp_amplitude = 14 #annual temperature variation [C]
t_min = 0 #starting time [years]
dt_years = .01 #timestep [years]
t_max_years = 10 #end time [years]
z_min = 0 #minimum depth [m]
z_max = 150 #maximum depth [m]
dz = 0.01 #spatial step in depth [m]
rho = 2650 #density of subsurface material [kg/m3]
c = 920 #heat capacity of subsurface material [J/kg/K], value is for ~sandstone
k = 2.0 #thermal conductivity of subsurface material []
should_plot = 1 #0 for no plots, 1 for plots
t_plot_years = 0.1 #plot every __ years

#unit conversions:
period_year = 24 * 3600 * 365
t_min = t_min * 365 * 24 * 3600
dt = dt_years * 365 * 24 * 3600
t_plot_seconds = t_plot_years * 365 * 24 * 3600
t_max = t_max_years * 365 * 24 * 3600; #seconds

#instantiate spatial and temporal domain, other important vars
z = np.arange(z_min,z_max, dz)
times = np.arange(t_min,t_max + dt,dt)
surface_condition = ts_bar + (temp_amplitude*np.sin(2*np.pi*times/period_year))
t_benchmark_years = t_max_years - dt_years
t_benchmark_seconds = t_max - dt
kappa = k / (rho * c)
mu = (kappa * dt) / np.power(dz, 2)
z_star_annual = np.sqrt(kappa * period_year / np.pi)

#boundary condition constants for Jacobian matrix
c1 = 1
c2 = 0
c3 = 0
c4 = 1

#record keeping array:
temp_record = np.zeros((len(times), len(z)))

#assemble vectors making up tridiagonal Jacobian matrix:
a_guts = np.repeat(-mu, len(z)-2)
a_add_elem1 = np.insert(a_guts, 0, 0)
a_final = np.append(a_add_elem1, -c4/dz)
b_guts = np.repeat(1 + 2*mu, len(z)-2)
b_add_elem1 = np.insert(b_guts, 0, (c1 - c2 / dz))
b_final = np.append(b_add_elem1, (c3 + c4 / dz))
c_guts = np.repeat(-mu, len(z)-2)
c_add_elem1 = np.insert(c_guts, 0, c2 / dz)
c_final = np.append(c_add_elem1, 0)

#RUN######################################################

t_old = ts_bar + temp_amplitude * np.exp(-z/z_star_annual) * np.sin((2 * np.pi * 0) / period_year - (z/z_star_annual))
temp_record[0, :] = t_old
it = 0
if should_plot == 1:
    profiles = plt.figure(figsize=(4, 8))
    temp_prof = plt.subplot()
    plt.xlabel('Temp [C]')
    plt.ylabel('Depth [m]')
    temp_prof.set_xlim(ts_bar-temp_amplitude-1, ts_bar+temp_amplitude+1)
    plt.gca().invert_yaxis()
    text = plt.text(ts_bar-temp_amplitude, z_max-1, 'Time [yrs]:   %.1f' % 0)
    plt.ion()
    plt.show()
else:
    pass

for t in range(len(times)-1):
    current_time = times[t]
    current_time_years = current_time / 365 / 24 / 3600
    it += 1
    f1 = surface_condition[it]
    f2 = c4 * ((t_old[-2]-t_old[-1]) / dz)
    interior_nodes = t_old[1:-1]
    rtemp = np.insert(interior_nodes, 0, f1) #apply BC's
    r = np.append(rtemp, f2)    
    t_new = thomas(a_final,b_final,c_final,r)
    temp_record[it, :] = t_new
    t_old = t_new
    if current_time == t_benchmark_seconds:
        benchmark = t_new
        bench_index = t
    else:
        pass
    if (should_plot == 1) & (current_time % t_plot_seconds == 0):
        temp_prof.plot(t_new, z)
        text.set_text('Time [yrs]:   %.1f' % current_time_years)
        plt.draw()
        plt.pause(0.001)
    else:
        pass

#benchmarking model against analytical solution
temp_annual = ts_bar + temp_amplitude * np.exp(-z/z_star_annual) * np.sin((2 * np.pi * times[bench_index]) / period_year - (z/z_star_annual))
model_error = max(abs(temp_annual-benchmark))
print '-------------------------------'
print 'Model compared to analytical solution at %d years.' % t_benchmark_years
print 'Maximum node error: %.2f C.' % model_error
print '-------------------------------'