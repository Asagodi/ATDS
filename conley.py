#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 09:32:10 2019

@author: abel
"""

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


def index_pair(N, graph, delta):
    A = N.copy()
    B = N.copy()
    while True:
        Aprime = calc_F(A, graph)
        Bprime = calc_F_inv(B, graph)
        if A == Aprime and B == Bprime:
            break
        A = Aprime.copy()
        B = Bprime.copy()
    C = A.intersection(B)
    r = delta
    intN = 0
    if B(C,R) in (intN):
        P1 = A
        P0 = A - B
        return (P1, P0)
    else:
        print("Failure")
        
def make_graph(cube_ind, time):
    for t,ci2 in enumerate(cube_ind[:]):
        if t % (time+1) == 0:
            ci1=cube_ind[t]
            continue
    
        conn_graph[ci1,ci2] = 1.
        G.add_edge(ci1, ci2)
        ci1 = ci2
    return G