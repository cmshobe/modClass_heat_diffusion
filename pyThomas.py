# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 12:47:05 2016

@author: Charlie

PYTHON VERSION OF THE THOMAS ALGORITHM
"""
import numpy as np
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