#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 18:00:56 2019

@author: abel
"""
from scipy import *
import networkx as nx
import numpy as np
import pickle
import subprocess

runfile('/home/abel/Downloads/Py_XPPCALL-master/Py_XPPCALL_Example.py', wdir='/home/abel/Downloads/Py_XPPCALL-master')

path_to_ode = '/home/abel/Downloads/xppaut8.0ubuntu/ode/'

def get_cubes_from_data(data, delta):
    cubes = []
    time_cubes = []
    cube_ind = [-1]
    for i,t in enumerate(data[:-1, 0]):
        cube = (np.sign(data[i, 0])*delta*math.floor(abs(data[i, 0])/delta)+np.sign(data[i, 0])*delta/2.,
                  np.sign(data[i, 1])*delta*math.floor(abs(data[i, 1])/delta)+np.sign(data[i, 1])*delta/2.)
        time_cubes.append(cube)
        if cube in cubes:
            ind = cubes.index(cube)
            cube_ind.append(ind)
        else:
            cube_ind.append(max(cube_ind)+1)
            cubes.append(cube)
    cube_ind.remove(-1)
    return cubes, time_cubes, cube_ind

path_to_ode = '/home/abel/Downloads/xppaut8.0ubuntu/ode/'
#RCs
xlo=0;xhi=74;ylo=0;yhi=52
grid=[]
delta = 1.
for i in range(int((xhi-xlo)/delta)):
    for j in range(int((yhi-ylo)/delta)):
        gc=(xlo+i*delta, ylo+j*delta)
        grid.append(gc)
        
time = 100
data = []
for gc in grid:
    npa, vn = xpprun(path_to_ode+'leslie.ode', inits={'x':gc[0],'y':gc[1]}, parameters={'total':time}, clean_after=True)
    data.append(npa[:,1:])        


data = np.array(data)
data = np.reshape(data, (data.shape[0]*data.shape[1], data.shape[2]))

delta = .25
cubes, time_cubes, cube_ind = get_cubes_from_data(data, delta)
    
G = nx.DiGraph()
number_of_cubes = len(cubes)
conn_graph = zeros((number_of_cubes,number_of_cubes))

for t,ci2 in enumerate(cube_ind[:]):
    if t % (time+1) == 0:
        ci1=cube_ind[t]
        continue

    conn_graph[ci1,ci2] = 1.
    G.add_edge(ci1, ci2)
    ci1 = ci2
    
scc=nx.strongly_connected_component_subgraphs(G)
RCs = []
for cc in scc:
    if len(cc.nodes())>1:
        RCs.append(cc.nodes())
        

fig, ax = plt.subplots(dpi=141)

xmin = 0
xmax = 74
ymin = 0
ymax = 52

x = np.array(cubes)[RCs[0]][:,0]
y = np.array(cubes)[RCs[0]][:,1]


ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

ax.set_aspect(1)
fig.canvas.draw()
s = max(1,((ax.get_window_extent().width  / (xmax-xmin+1.) * delta*72./fig.dpi) ** 2))  

for i in range(len(RCs)):
    plt.scatter(np.array(cubes)[RCs[i]][:,0], np.array(cubes)[RCs[i]][:,1], s=s, marker='s')

plt.show()
        
        
#A = list(G.out_edges(RCs[0][0]))
#for i in range(10):
#    
#    Aprime = A
#    for zeta in A:
#        for newzeta in G.out_edges(zeta):
#            Aprime.append(newzeta)
#    A = list(set(Aprime))
#    
#    
#zeta1 = np.where(np.abs(data[:,0]-np.array(cubes)[RCs[0][0]][0])<delta/2.)
#zeta2 = np.where(np.abs(data[:,1]-np.array(cubes)[RCs[0][0]][1])<delta/2.)
#plt.scatter(data[np.array(list(set(list(zeta1[0])).intersection(set(list(zeta2[0]))))),:][0,:,0],data[np.array(list(set(list(zeta1[0])).intersection(set(list(zeta2[0]))))),:][0,:,1])
#    


candidates = set()

for component in RCs:
    for zeta in component:
        candidates.add(zeta)


attractors = []
while len(candidates) > 0:
    A = set([list(candidates)[0]])
    for i in range(1000):

        Aprime = set()
        for zeta in A:
            for newzeta in G.out_edges(zeta):
                Aprime.add(newzeta[1])
        
        if A == Aprime:
            #check whether invariant
            if calc_F_inv(A, G):
                attractors.append(A)
            break
        else:
            for a in A:
                Aprime.add(a)
            A = Aprime.copy()
            
    candidates = candidates - A#this is not a good thing to do
    
        
    
candidates = set()

for component in RCs:
    for zeta in component:
        candidates.add(zeta)


repellers = []
while len(candidates) > 0:
    R = set([list(candidates)[0]])
    for i in range(1000):

        Rprime = R.copy()
        for zeta in R:
            for newzeta in G.in_edges(zeta):
                Rprime.add(newzeta[0])

        if R == Rprime:
            repellers.append(R)
            break
        else:
            R = Rprime.copy()
    candidates = candidates - R#this is not a good thing to do
    
        
        
for i in range(len(attractors)):
    plt.scatter(np.array(cubes)[list(attractors[i])][:,0], np.array(cubes)[list(attractors[i])][:,1], s=s, marker='s')
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.show()
        
    
morse_sets = set()
for A in attractors:
    for R in repellers:
        morse_sets.add(frozenset(A.intersection(R)))
        
        
for i in range(len(morse_sets)):
    plt.scatter(np.array(cubes)[list(list(morse_sets)[i])][:,0], np.array(cubes)[list(list(morse_sets)[i])][:,1], s=s, marker='s')
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.show()
    
    
def calc_F(U, graph):
    Uprime = set()
    for zeta in U:
        for edge in graph.out_edges(zeta):
            Uprime.add(edge[1])
    return Uprime


def calc_F_inv(U, graph):
    Uprime = set()
    for zeta in U:
        for edge in graph.in_edges(zeta):
            Uprime.add(edge[0])
    return Uprime


def maximal_closed_subgraph(G):
    nnodes = len(G.nodes())
    while True:
        for i in G.nodes():
            if len(G.in_edges(i)) == 0:
                G.remove_node(i)
            elif len(G.out_edges(i)) == 0:
                G.remove_node(i)
        if len(G.nodes())==nnodes:
            break
        else:
            nnodes = len(G.nodes())
    return G