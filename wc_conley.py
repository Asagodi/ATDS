#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 12:55:35 2019

@author: abel
"""

from scipy import *
import networkx as nx
import numpy as np
import pickle
import subprocess
import conley_functions as confun
import matplotlib.pyplot as plt

runfile('/home/abel/Downloads/Py_XPPCALL-master/Py_XPPCALL_Example.py', wdir='/home/abel/Downloads/Py_XPPCALL-master')
path_to_ode = '/home/abel/Downloads/xppaut8.0ubuntu/ode/'

#RCs
xlo=-.125;ylo=-.125;xhi=1;yhi=1
grid=[]
delta = .5
for i in range(int((xhi-xlo)/delta)):
    for j in range(int((yhi-ylo)/delta)):
        gc=(xlo+i*delta, ylo+j*delta)
        grid.append(gc)
        
time = 100
data = []
for gc in grid:
    npa, vn = xpprun(path_to_ode+'lif.ode', inits={'x':gc[0],'y':gc[1]}, parameters={'total':time}, clean_after=True)
    data.append(npa[:,1:])        


arrdata = np.array(data)
data = np.reshape(arrdata, (arrdata.shape[0]*arrdata.shape[1], arrdata.shape[2]))

delta = .001
cubes, time_cubes, cube_ind = confun.get_cubes_from_data(data, delta)
    
G = nx.DiGraph()
number_of_cubes = len(cubes)
conn_graph = zeros((number_of_cubes,number_of_cubes))

G = confun.make_graph(cube_ind, time)
    
scc=nx.strongly_connected_component_subgraphs(G)
RCs = []
for cc in scc:
    if len(cc.nodes())>1:
        RCs.append(cc.nodes())
        

fig, ax = plt.subplots(dpi=141)

x = np.array(cubes)[RCs[0]][:,0]
y = np.array(cubes)[RCs[0]][:,1]


ax.set_xlim(xlo, xhi)
ax.set_ylim(ylo, yhi)

ax.set_aspect(1)
fig.canvas.draw()
s = max(1,((ax.get_window_extent().width  / (xhi-xlo+1.) * delta*72./fig.dpi) ** 2))  

for i in range(len(RCs)):
    plt.scatter(np.array(cubes)[RCs[i]][:,0], np.array(cubes)[RCs[i]][:,1], s=s, marker='s')

plt.show()
        

G_mcs = confun.maximal_closed_subgraph(G)
    
attractors = []
for component in RCs:
    A = set(component)

    for i in range(1000):

        Aprime = set()
        for zeta in A:
            for newzeta in G_mcs.in_edges(zeta):
                Aprime.add(newzeta[0])

        if A == Aprime:
            #check whether really invariant:
#            if confun.calc_F(A, G) == A:
            attractors.append(A)
            break
        else:
            for r in A:
                Aprime.add(r)
            A = Aprime.copy()
    
    
repellers = []
for component in RCs:
    R = set(component)

    for i in range(1000):

        Rprime = set()
        for zeta in R:
            for newzeta in G_mcs.in_edges(zeta):
                Rprime.add(newzeta[0])

        if R == Rprime:
            #check whether really invariant:
#            if calc_F_inv(R, G) == R:
            repellers.append(R)
        
            break
        else:
            for r in R:
                Rprime.add(r)
            R = Rprime.copy()
    
        
        
    
morse_sets = set()
for A in attractors:
    for R in repellers:
        morse_sets.add(frozenset(A.intersection(R)))
        
        
for i in range(len(morse_sets)):
    plt.scatter(np.array(cubes)[list(list(morse_sets)[i])][:,0], np.array(cubes)[list(list(morse_sets)[i])][:,1], s=s, marker='s')
    plt.xlim(xlo, xhi)
    plt.ylim(ylo, yhi)
    plt.show()
    
    

