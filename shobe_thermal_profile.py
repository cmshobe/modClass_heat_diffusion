# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 16:38:47 2016
Last updated Jan 25 2016
@author: Charlie Shobe

This is a finite difference solution to a thermal temperature profile. This
model uses a Backward Euler implicit solving scheme, with a matrix solution
solved by the Thomas algorithm.

Dependencies:
-Numpy
-Matplotlib

Inputs:
-'cape_thompson.dat': a record of borehole temperatures at depth

Outputs:
-'borehole_data_fitting.png': a linear fit to the lower portion of borehole temp profile
-'thermal_profiles.png': thermal profiles over time
-'data_model_comparison.png': modeled thermal profile compared with field data
-the error between data and model (sum of squares) is printed to the screen

User guide:
***To run the model and test against borehole data, simply run with defaults***
-Set the value of sinusoidal_temp_or_borehole_evolution to tell the script 
    whether you'd like to explore a sinusoidal surface condition, or replicate 
    and fit the borehole data.
-If you choose the sinusoidal option, z_max values of tens of meters are ideal
    so that you can see the profile evolve.
-If you choose to replicate the borehole data, z_max should be set to 400 m to
    match the borehole data from 'cape_thompson.dat'. The script will fit the
    linear portion of the borehole data to find its initial condition,
    and will implement the step change from Lachenbruch and Marshall. If 
    you choose this option, t_max should be set to 100 years.
"""
#import modules:
from __future__ import division #true division
import numpy as np #numpy
import matplotlib.pyplot as plt #plotting
import matplotlib #plotting

#define matrix equation solver (Thomas algorithm)
def thomas(a,b,c,r):
    '''
    Solves tridiagonal system of equations, in this case for the purpose of
    implicit solutions for numerical models.
    '''
    J = len(b) #number of nodes
    #j = np.arange(1,J) #node vector

    f = np.zeros(len(b))
    g = np.zeros(len(b))
    f[0] = b[0]
    g[0] = c[0] / f[0]
    
    e = a #assumes zero has already been added by code that calls Thomas

    for i in range(1,J):
        f[i] = b[i] - e[i] * g[i-1]
        g[i] = c[i] / f[i]
    y = np.zeros(len(b))
    y[0] = r[0] / f[0]
    for i2 in range(1, J):
        y[i2] = (r[i2] - e[i2] * y[i2-1]) / f[i2]
    
    x = np.zeros(len(b))
    x[J-1] = y[J-1]
    for i3 in range(J-1, 0, -1): #iterate backwards through nodes
        x[i3-1] = y[i3-1] - g[i3-1] * x[i3]

    return x
    
#user defined parameters:#####################################################
sinusoidal_temp_or_borehole_evolution = 1 #0 for sinusoidal temp, 1 for borehole evolution
ts_bar = -5 #mean surface temp [C]
temp_amplitude = 14 #annual temperature variation [C]
t_min = 0 #starting time [years]
dt_years = 1 #timestep [years]
t_max_years = 100 #end time [years]. SET TO 100 FOR BOREHOLE DATA COMPARISON.
z_min = 0 #minimum depth [m]
z_max = 400 #maximum depth [m] SET TO 400 FOR BOREHOLE DATA COMPARISON, <100 FOR SINUSOIDAL.
dz = 0.01 #spatial step in depth [m]
rho = 2000 #density of subsurface material [kg/m3]
c = 920 #heat capacity of subsurface material [J/kg/K], value is for ~sandstone
k = 2.0 #thermal conductivity of subsurface material []
should_plot = 1 #0 for no plots, 1 for plots
t_plot_years = 10 #plot every __ years
##############################################################################

#this section imports and fits the borehole data from cape_thompson.dat.'#####
#plotting parameters:
matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams['xtick.major.pad']='10'

if sinusoidal_temp_or_borehole_evolution == 1:
    #import borehole data:
    cape_thompson = np.fromfile('cape_thompson.dat', dtype=float, sep = " ")
    depth_indices = np.arange(0, 46, 2)
    temp_indices = np.arange(1, 46, 2)
    depths = cape_thompson[depth_indices]
    temps = cape_thompson[temp_indices]
    
    temp_data = plt.figure(figsize=(6, 8))
    thompson_data = plt.subplot()
    thompson_data.plot(temps, depths, marker='o', markersize=8, linewidth=2, label='Borehole Data')
    plt.gca().invert_yaxis()
    plt.xlabel('Temp [C]')
    plt.ylabel('Depth [m]')
    plt.title('Borehole Data')
    
    linfit = np.polyfit(temps[9:], depths[9:], 1)
    line_plotting_x_values = np.arange(-7, 0.1, 0.1)
    line_values = linfit[0] * line_plotting_x_values + linfit[1]
    thompson_data.plot(line_plotting_x_values, line_values, color='r', linewidth=2, linestyle='--', label = 'Linear Fit')
    plt.legend(loc='lower left')
    temp_data.savefig('borehole_data_fitting.png')
elif sinusoidal_temp_or_borehole_evolution == 0:
    print 'Sinusoidal surface signal, skipping borehole data analysis.'
else:
    print 'CHOOSE EITHER 0 OR 1 IN LINE 48!'
##############################################################################

#INITIALIZE################################################
#unit conversions:
period_year = 24 * 3600 * 365
t_min = t_min * 365 * 24 * 3600
dt = dt_years * 365 * 24 * 3600
t_plot_seconds = t_plot_years * 365 * 24 * 3600
t_max = t_max_years * 365 * 24 * 3600; #seconds
 
#instantiate spatial and temporal domain, other important vars
z = np.arange(z_min,z_max, dz)
times = np.arange(t_min,t_max + dt,dt)
if sinusoidal_temp_or_borehole_evolution == 0:
    surface_condition = ts_bar + (temp_amplitude*np.sin(2*np.pi*times/period_year))
elif sinusoidal_temp_or_borehole_evolution == 1:
    surface_condition = np.concatenate((np.repeat(-7, 59/dt_years), np.repeat(-5, 41/dt_years)))
else:
    print 'CHOOSE EITHER 0 OR 1 IN LINE 48!'
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

#initial conditions:
if sinusoidal_temp_or_borehole_evolution == 1:
    t_initial = (z - linfit[1]) / linfit[0]
    t_old = t_initial
elif sinusoidal_temp_or_borehole_evolution == 0:
    t_old = ts_bar + temp_amplitude * np.exp(-z/z_star_annual) * np.sin((2 * np.pi * 0) / period_year - (z/z_star_annual))

#RUN######################################################
temp_record[0, :] = t_old
it = 0
if should_plot == 1:
    profiles = plt.figure(figsize=(6, 8))
    temp_prof = plt.subplot()
    plt.xlabel('Temp [C]')
    plt.ylabel('Depth [m]')
    plt.title('Temp Profiles')
    temp_prof.set_xlim(ts_bar-temp_amplitude-1, ts_bar+temp_amplitude+1)
    plt.gca().invert_yaxis()
    text = plt.text(ts_bar-temp_amplitude, z_max-.1*z_max, 'Time [yrs]:   %.1f' % 0)
    plt.ion()
    plt.show()
else:
    pass

for t in range(len(times)-1):
    current_time = times[t] + dt
    current_time_years = current_time / 365 / 24 / 3600
    it += 1
    f1 = surface_condition[it-1]
    if sinusoidal_temp_or_borehole_evolution == 0:
        f2 = c4 * ((t_old[-2]-t_old[-1]) / dz)
    elif sinusoidal_temp_or_borehole_evolution == 1:
        f2 = c4 * .0195
    else:
        print 'CHOOSE EITHER 0 OR 1 IN LINE 48!'
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
        text.set_text('Time [yrs]: %.1f' % current_time_years)
        plt.draw()
        plt.pause(0.001)
    else:
        pass

if sinusoidal_temp_or_borehole_evolution == 1:
    #comparing model output to borehole data
    comparison = plt.figure(figsize=(6, 8))
    compare = plt.subplot()
    plt.gca().invert_yaxis()
    compare.plot(t_initial, z, color='r', linestyle='--', label='Initial Condition')
    #compare.plot(line_plotting_x_values, line_values, linestyle='--', linewidth=2)
    compare.plot(temps, depths, marker='o', linewidth=0, label='Borehole Data')
    compare.plot(temp_record[-1, :], z, label='Model Prediction')
    plt.xlabel('Temp [C]')
    plt.ylabel('Depth [m]')
    plt.title('Model, meet data')
    plt.legend(loc='lower left')
    comparison.savefig('data_model_comparison.png')
    #calculating error by sum of the squares
    rounded_borehole_depths = np.round(depths, decimals=2)*1/dz
    borehole_indices = rounded_borehole_depths.astype(int)
    model_temps_to_compare = temp_record[-1, borehole_indices]
    error = sum(np.power(model_temps_to_compare - temps, 2))
    print '------------------------------------------------'
    print 'Exporting figures...'
    print 'Simulation complete.'
    print 'Error (sum of squares): %.2f' % error
    print '------------------------------------------------'
    
elif sinusoidal_temp_or_borehole_evolution == 0:
    print 'Exporting figures...'
    print 'Simulation complete.'
    
#save model figure
if should_plot == 1:
    profiles.savefig('thermal_profiles.png', bbox_inches='tight')
else:
    pass