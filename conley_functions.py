#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 09:32:10 2019

@author: abel
"""
import numpy as np 
import math
import networkx as nx
from scipy.sparse import csr_matrix, lil_matrix
from copy import deepcopy
import itertools


def heaviside(x1):
    if x1<0:
        return -1
    if x1>=0:
        return 1


class Dynamical_System(object):
    def __init__(self, data, delta):
        sefl.data = data
        self.delta = delta



def get_cubes_from_data(data, delta, cubes = [], time_cubes = [], cube_ind = [-1]):
    dim = data.shape[1]
    sigdelta = str(delta)[::-1].find('.')+1
    deltadiv2 = round(delta/2., sigdelta)
    cubes = []
    time_cubes = []
    cube_ind = [-1]
    cube_ind_dict = {}
    for i,t in enumerate(data[:-1, 0]):
        cube = []
        for j in range(dim):
            coord = heaviside(data[i, j])*(round(delta*np.floor_divide(data[i,j],delta), sigdelta)+deltadiv2)
            cube.append(coord)
        time_cubes.append(cube)
        if cube in cubes:
            ind = cubes.index(cube)
            cube_ind.append(ind)
        else:
            cube_ind.append(max(cube_ind)+1)
            cubes.append(cube)
            cube_ind_dict[tuple(cube)] = max(cube_ind)+1
    cube_ind.remove(-1)
    return cubes, time_cubes, cube_ind, cube_ind_dict

def get_cubes_from_data_intervals(data, delta, cubes = [], time_cubes = [], cube_ind = [-1]):
    sigdelta = str(delta)[::-1].find('.')+1
    dim = data.shape[1]
    cubes = []
    time_cubes = []
    cube_ind = [-1]
    for i,t in enumerate(data[:-1, 0]):
        cube = []
        for j in range(dim):
            interval = tuple([heaviside(data[i, j])*round(delta*np.floor_divide(data[i,j],delta), sigdelta),
                              heaviside(data[i, j])*round(delta*np.floor_divide(data[i,j],delta), sigdelta) + delta])
            cube.append(interval)
        time_cubes.append(cube)
        if cube in cubes:
            ind = cubes.index(cube)
            cube_ind.append(ind)
        else:
            cube_ind.append(max(cube_ind)+1)
            cubes.append(cube)
    cube_ind.remove(-1)
    return cubes, time_cubes, cube_ind

def get_cubes_from_datalist(data_list, delta):
    cubes_list = []
    time_cubes_list = []
    cube_ind_list = []
    for data in data_list:
        cubes, time_cubes, cube_ind = get_cubes_from_data(data, delta)
        cubes_list.append(cubes)
        time_cubes_list.append(time_cubes)
        cube_ind_list.append(cube_ind)
    return cubes_list, time_cubes_list, cube_ind_list



def get_recurrent_components(G, includeselfedges=False):
    
    scc=nx.strongly_connected_components(G)
    RCs = []
    for cc in scc:
        if len(list(G.subgraph(cc).nodes()))>1:
            RCs.append(list(G.subgraph(cc).nodes()))
    if includeselfedges == True:
        for edge in G.edges():
            if edge[0]==edge[1]:
                RCs.append([edge[0]])
    return RCs

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

def invariantPart(N, G):
    '''Combinatorial invariant set S for set N and graph (function) G'''
    H = restricted_map(G, N)
    S = set(N)
    while True:  
        Sp = S
        Spp = S.intersection(calc_F(S, H))
        S = Spp.intersection(calc_F_inv(S, H))
        if S == Sp:
            break
    return S

def index_pair(N, G, delta):
    S = invariantPart(N,G)
    M = cubical_wrap(S, delta) 
#     print(M,N)
    if M.issubset(N):
        F = restricted_map(G, M)
        
        C = get_neighours_cubicalset(S, delta)#collar(S)
        P0 = calc_F(S, F).intersection(C)
        while True:
            lastP0 = P0
            P0 = calc_F(P0, F).intersection(C)
            P0 = P0.union(lastP0)
            if P0 == lastP0:
                break
        P1 = S.union(P0)
        Pbar1 = calc_F(P1, F)
        Pbar0 = Pbar1 - S
        return P1, P0, Pbar1, Pbar0
    else:
        return "Failure"
        
def get_neighours_cube(cube, delta, kmax):
    '''takes the neighbours of a cube'''
    L = set()
    for some in itertools.product([-1,0,1], repeat=kmax):
        face = list(deepcopy(cube))
        for i,s in enumerate(some):
            face[i] = cube[i]+s*delta
        L.add(tuple(face))
    L.remove(cube)
    
    return L

def get_neighours_cubicalset(S, delta):
    '''takes the neighbours of a cubical set,
    same as collar'''
    N = cubical_wrap(S, delta)
    return N-set(S)

def cubical_wrap(S, delta):
    '''takes the neighbours of a cubical set,
    same as collar'''
        #now only works for cubes of dim kmax
    maxdim = 0
    for cube in S:
        maxdim = max(get_dim(cube, delta), maxdim)
    print(maxdim)
    
    N = set()
    for cube in S:
        L = get_neighours_cube(cube, delta, maxdim)
        for l in L:
            N.add(l)
    return N
        
# def make_graph(cubes, cube_ind, time):
#     ncubes = len(cubes)
#     A=csr_matrix((ncubes, ncubes), dtype=np.int8)
#     G = nx.DiGraph()
#     for t,ci2 in enumerate(cube_ind[:]):
#         if t % (time+1) == 0:
#             ci1=cube_ind[t]
#             continue
    
#         G.add_edge(ci1, ci2)
#         ci1 = ci2
#     return G, A

def make_graph_fl(cube_ind_list):
    G = nx.DiGraph()
    for cube_ind in cube_ind_list:
        for t,ci2 in enumerate(cube_ind[:]):
            ci1=cube_ind[t]
            continue
    
            G.add_edge(ci1, ci2)
            ci1 = ci2
    return G

def make_graph(cubes, cube_ind, data_length_list):
    ncubes = len(cubes)
    A=lil_matrix((ncubes, ncubes), dtype=np.int8)
    G = nx.DiGraph()
    i=0
    ci1=3
    s=-1
    for t,ci2 in enumerate(cube_ind):
        s+=1
        if s % data_length_list[i] == 0:
            ci1=cube_ind[t+1]
            if s == data_length_list[i]:
                i+=1
                s=0
            continue
        if s % data_length_list[i] == 1:
            continue
        A[ci1, ci2] = 1.
        G.add_edge(ci1, ci2)
        ci1 = ci2
    return G, A

def get_recurrent_components(G, includeselfedges=False):
    
    scc=nx.strongly_connected_components(G)
    RCs = []
    for cc in scc:
        if len(list(G.subgraph(cc).nodes()))>1:
            RCs.append(list(G.subgraph(cc).nodes()))
    if includeselfedges == True:
        for edge in G.edges():
            if edge[0]==edge[1]:
                RCs.append([edge[0]])
    return RCs


def convert_to_chomp_format(cubical_set, delta):
    filetxt = ""
    for cube in cubical_set:
        filetxt += '['
        for coord in cube:
            filetxt+=str(tuple((np.array(coord)/delta).astype(int)))+ ' '
        filetxt = filetxt[:-1]+']\n'

    return filetxt




def get_dim(cube, delta):
    "returns dimension of cube"
    return np.where( np.abs(np.mod(cube, delta)-delta/2.) < 0.00001 )[0].shape[0]

def primaryFaces(cube, delta):
    "cube = list(coordinate of centre)"
    "return primary faces of cube"
    relcoord = np.where( np.abs(np.mod(cube, delta)-delta/2.) < 0.00001 )[0]
    k = relcoord.shape[0]
    L = set()
    
    for i in range(k):
        for j in range(2):
            face = list(deepcopy(cube))
            face[relcoord[i]] = cube[relcoord[i]]+(-1)**j*delta/2.
            L.add(tuple(face))
    return L


def boundaryOperator(cube, delta):
    sgn = 1
    chain = {}
    relcoord = np.where( np.abs(np.mod(cube, delta)-delta/2.) < 0.00001 )[0]
    k = relcoord.shape[0]
    
    for i in range(k):
        for j in range(2):
            face = list(deepcopy(cube))
            face[relcoord[i]] = cube[relcoord[i]]+(-1)**j*delta/2.
            chain[tuple(face)] = (-1)**j*sgn
        sgn = -sgn
    return chain


def cubicalChainGroups(K, delta):
    "K: cubical set (list of its maximal faces)"
    "return the groups of cubical chains of a cubical set"
    E = []
    
    while K != set():
        Q = list(K)[0]
        K = set(list(K)[1:])
        k = get_dim(Q, delta)
        L = primaryFaces(Q, delta)
        if L != set():
            K = K.union(L)
        try:
            E[k] = E[k].union([Q])
        except:
            maxk = len(E)
            for ktok in range(k-maxk+1):
                E.append(set())
            E[k] = E[k].union([Q])
        if L != set():
            E[k-1] = E[k-1].union(L)
    return E

def canonicalCoordinates(chain, K):
    #K is list of all elementary cubes (can be calc by unrolledE)
    v = np.zeros(len(K))
    for i in range(len(K)):
#         print(chain[tuple(list(K)[i])])
        try:
            v[i] = chain[tuple(list(K)[i])]
        except:
            0
    return v
def boundaryOperatorMatrix(E, delta):
    D = [np.array([]) for k in range(len(E))]
    for k in range(1, len(E)): #range(0, ... ??
        m = len(E[k-1])
        l = len(E[k])
        D[k] = np.zeros((m, l))
        for j in range(l):
            chain = boundaryOperator(list(E[k])[j], delta)
#             print(chain)
            v = canonicalCoordinates(chain, E[k-1])
#             print(v)
            D[k][:,j] = v
    D[0] = np.array([[0]*len(E[0])])
    return D


def rowSwap(A, i, j):
   temp = np.copy(A[i, :])
   A[i, :] = A[j, :]
   A[j, :] = temp
 
def colSwap(A, i, j):
   temp = np.copy(A[:, i])
   A[:, i] = A[:, j]
   A[:, j] = temp
 
def scaleCol(A, i, c):
   A[:, i] *= c*np.ones(A.shape[0])
 
def scaleRow(A, i, c):
   A[i, :] *= c*np.ones(A.shape[1])
 
def colCombine(A, addTo, scaleCol, scaleAmt):
   A[:, addTo] += scaleAmt * A[:, scaleCol]
 
def rowCombine(A, addTo, scaleRow, scaleAmt):
   A[addTo, :] += scaleAmt * A[scaleRow, :]

def numPivotCols(A):
   z = np.zeros(A.shape[0])
   return [np.all(A[:, j] == z) for j in range(A.shape[1])].count(False)
 
def numPivotRows(A):
   z = np.zeros(A.shape[1])
   return [np.all(A[i, :] == z) for i in range(A.shape[0])].count(False)

def simultaneousReduce(A, B):
   if A.shape[1] != B.shape[0]:
      raise Exception("Matrices have the wrong shape.")
 
   numRows, numCols = A.shape # col reduce A
 
   i,j = 0,0
   while True:
      if i >= numRows or j >= numCols:
         break
 
      if A[i][j] == 0:
         nonzeroCol = j
         while nonzeroCol < numCols and A[i,nonzeroCol] == 0:
            nonzeroCol += 1
 
         if nonzeroCol == numCols:
            i += 1
            continue
 
         colSwap(A, j, nonzeroCol)
         rowSwap(B, j, nonzeroCol)
 
      pivot = A[i,j]
      scaleCol(A, j, 1.0 / pivot)
      scaleRow(B, j, 1.0 / pivot)
 
      for otherCol in range(0, numCols):
         if otherCol == j:
            continue
         if A[i, otherCol] != 0:
            scaleAmt = -A[i, otherCol]
            colCombine(A, otherCol, j, scaleAmt)
            rowCombine(B, j, otherCol, -scaleAmt)
 
      i += 1; j+= 1
 
   return A,B

def singleReduce(A):
 
   numRows, numCols = A.shape 
 
   i,j = 0,0
   while True:
      if i >= numRows or j >= numCols:
         break
 
      if A[i][j] == 0:
         nonzeroCol = j
         while nonzeroCol < numCols and A[i,nonzeroCol] == 0:
            nonzeroCol += 1
 
         if nonzeroCol == numCols:
            i += 1
            continue
 
         colSwap(A, j, nonzeroCol)
 
      pivot = A[i,j]
      scaleCol(A, j, 1.0 / pivot)
 
      for otherCol in range(0, numCols):
         if otherCol == j:
            continue
         if A[i, otherCol] != 0:
            scaleAmt = -A[i, otherCol]
            colCombine(A, otherCol, j, scaleAmt)
 
      i += 1; j+= 1
 
   return A

def bettiNumber(k, d_k, d_kplus1):
    A, B = np.copy(d_k), np.copy(d_kplus1)
    if k == 0:
        B = singleReduce(B.T).T

    simultaneousReduce(A, B)

    dimKChains = A.shape[1]
    kernelDim = dimKChains - numPivotCols(A)
    imageDim = min(numPivotRows(B), B.shape[1]) #is this the right fix?
    
    print(dimKChains, kernelDim, imageDim)

    return kernelDim - imageDim


def restricted_map(G, N):
    '''for a graph G takes subgraph with nodes in N'''
    return nx.subgraph(G, N)