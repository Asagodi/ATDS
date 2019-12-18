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


def heaviside(x):
    if x<0:
        return -1
    if x>=0:
        return 1


class Combinatorial_Dynamical_System(object):
    "Combinatorial representation of datapoints from samples (from a dynamical system)"
    def __init__(self, delta):
        """
        delta: diameter of cubes
        self.cubes are the cubes in the grid that have a datapoint in them
        a cube is determined by the coordinates of its centre
        self.tuplecubes cubes as tuples
        self.cube_ind: the indices for the cubes in the order as they come in the data
        self.cube_ind_dict: dictionary between the cubes and the indices
        self.cube_ind_dict: dictionary between the indices and the cubes
        the multivalued map \mathcal{F}=graph G
        """
#         self.data = data
#         self.data_length_list = data_length_list
        self.delta = delta
        self.cubes = []
        self.tuplecubes = []
        self.cube_ind = [-1]
        self.cube_ind_dict = {}
        self.index_cube_dict = {}
        ncubes = 0
        self.A=lil_matrix((ncubes, ncubes), dtype=np.int8)
        self.G = nx.DiGraph()
        self.G.add_nodes_from([i for i in range(ncubes)])
        
    def initialise_with_data(self, data, data_length_list):
        cube_ind = self.get_cubes(data)
        self.update_graph(cube_ind, data_length_list)
        
        
    def get_cubes(self, data):
        """Creates cubes for all the data with elementary cubes with diameter self.delta
        data has shape (n, d) with n: number of samples, d:dimension of space
        also updates"""
        #should data be translated here? (data is supposed to be positive in homcubes).
        max_cubeind = max(self.cube_ind)
        cube_ind = [max_cubeind]
        dim = data.shape[1]
        sigdelta = str(self.delta)[::-1].find('.')+1
        deltadiv2 = round(self.delta/2., sigdelta)
        for i,t in enumerate(data[:-1, 0]):
            cube = []
            for j in range(dim):
                coord = round(heaviside(data[i, j])*(self.delta*np.floor_divide(data[i,j],self.delta)+deltadiv2), sigdelta)
                cube.append(coord)
            if cube in self.cubes:
                ind = self.cubes.index(cube)
                cube_ind.append(ind)
            else:
                max_cubeind += 1
                cube_ind.append(max_cubeind)
                self.cubes.append(cube)
                self.cube_ind_dict[tuple(cube)] = max_cubeind
        try:
            cube_ind.remove(max(self.cube_ind))
        except:
            0
        self.cube_ind.extend(cube_ind)
        self.index_cube_dict = {v: k for k, v in self.cube_ind_dict.items()}
        self.tuplecubes = []
        for cube in self.cubes:
            self.tuplecubes.append(tuple(cube))
        return cube_ind
    
    def update_graph(self, cube_ind, data_length_list, calc_matrix=False):
        """
        Updates graph G and transition matrix A from the cube_ind
        takes data_length_list as the list of the lengths 
        of the samples of the dynamics
        """
        i=0
        ci1=0
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
            #how to change size A?
            if calc_matrix:
                self.A[ci1, ci2] = 1.
            self.G.add_edge(ci1, ci2)
            ci1 = ci2
        
        

    def update_cubesandgraph(self, data, data_length_list):
        "Update cube set and graph for new data set"
        cube_ind = self.get_cubes(data)
        self.update_graph(cube_ind, data_length_list)

    def get_recurrent_components(self, includeselfedges=False):
        """Calculates recurrent components/strongly connected components
        these correnspond to the Morse sets
        """
        scc=nx.strongly_connected_components(self.G)
        RCs = []
        for cc in scc:
            if len(list(self.G.subgraph(cc).nodes()))>1:
                RCs.append(list(self.G.subgraph(cc).nodes()))
        if includeselfedges == True:
            for edge in self.G.edges():
                if edge[0]==edge[1]:
                    RCs.append([edge[0]])
        self.RCs = RCs
        return RCs

    def calc_F(self, U):
        "Calculates \mathcal{F}(U) for graph \mathcal{F}=self.G and subset U"
        Uprime = set()
        for zeta in U:
            for edge in self.G.out_edges(zeta):
                Uprime.add(edge[1])
        return Uprime

    def calc_F_inv(self, U):
        "Calculates \mathcal{F}^{-1}(U) for graph \mathcal{F}=self.G and subset U"
        Uprime = set()
        for zeta in U:
            for edge in self.G.in_edges(zeta):
                Uprime.add(edge[0])
        return Uprime
    
    def get_forward_invariant_subset(self, U, limit=1000):
        "Iterates until an invariant set is reached"
        A = set(U)
        #what should be the limit?
        for i in range(limit):

            Aprime = set()
            for zeta in A:
                for newzeta in self.G.in_edges(zeta):
                    Aprime.add(newzeta[0])

            if A == Aprime:
                break
            else:
                for r in A:
                    Aprime.add(r)
                A = Aprime.copy()
        return A
    
    def get_backward_invariant_subset(self, U, limit=1000):
        "Iterates back in time until an invariant set is reached"
        R = set(U)
        for i in range(limit):
            Rprime = set()
            for zeta in R:
                for newzeta in self.G.in_edges(zeta):
                    Rprime.add(newzeta[0])

            if R == Rprime:
                #check whether really invariant:
    #            if calc_F_inv(R, G) == R:
#                 repellers.append(R)
                break
            else:
                for r in R:
                    Rprime.add(r)
                R = Rprime.copy()
        return R
    

    def evaluate(self, F, U):
        "Calculates F(U) for give grap F and subset U"
        Uprime = set()
        for zeta in U:
            for edge in F.out_edges(zeta):
                Uprime.add(edge[1])
        return Uprime
    
    def evaluate_inv(self, F, U):
        "Calculates F^-1(U) for give grap F and subset U"
        Uprime = set()
        for zeta in U:
            for edge in F.in_edges(zeta):
                Uprime.add(edge[0])
        return Uprime

    def maximal_closed_subgraph(self):
        "Takes subgraph of all the nodes that have connetions both forwards and backwards"
        nnodes = len(self.G.nodes())
        while True:
            for i in self.G.nodes():
                if len(self.G.in_edges(i)) == 0:
                    self.G.remove_node(i)
                elif len(self.G.out_edges(i)) == 0:
                    self.G.remove_node(i)
            if len(self.G.nodes())==nnodes:
                break
            else:
                nnodes = len(self.G.nodes())
        self.max_closed = self.G
        return self.G

    def invariantPart(self, N):
        """Combinatorial invariant set S inside set N and graph G=\mathcal{F}
        returns the set Inv(N , F)"""
        H = deepcopy(nx.subgraph(self.G, N)) #self.restricted_map(N)
        S = set(N).copy()
        while True:  
            Sp = S
            Spp = S.intersection(self.evaluate(H, S))
            S = Spp.intersection(self.evaluate_inv(H, S))
            if S == Sp:
                break
        return S
    
    def convert_cubes_to_indices(self, cubicalset):
        "Gets indices of set of cubes"
        indexset = set()
        for cube in cubicalset:
            indexset.add(self.cube_ind_dict[cube]) 
        return indexset
            
    def convert_indices_to_cubes(self, indexset):
        "Gets coordinates of the centres of cubes of set of cube indices (inverse of convert_cubes_to_indices)"
        cubicalset = set()
        for cube in indexset:
            cubicalset.add(self.index_cube_dict[cube])
        return cubicalset
    
    def index_pair(self, N):
        """N is set of coordinates of the centres of cubes
        returns index pair of N, if it is an isolating neighbourhood,
        otherwise returns Failure
        """
        N_ind = []
        for cube in N:
            N_ind.append(self.cube_ind_dict[cube])

        S_ind = self.invariantPart(N_ind)
        S_cubes = []
        for cube in S_ind:
            S_cubes.append(self.index_cube_dict[cube])
        
        M = self.cubical_wrap(S_cubes).intersection(set(self.tuplecubes))
        if M.issubset(N):
            M_ind = []
            for cube in N:
                M_ind.append(self.cube_ind_dict[cube])
            F = self.restricted_map(M_ind)
            C = self.get_neighours_cubicalset(S_cubes)#collar(S)
            C = C.intersection(set(self.tuplecubes))
            C_ind = set()
            for cube in C:
                C_ind.add(self.cube_ind_dict[cube])
            P0 = set(self.evaluate(F, S_ind)).intersection(C_ind)
            while True:
                lastP0 = P0
                P0 = self.evaluate(F, P0).intersection(C)
                P0 = P0.union(lastP0)
                if P0 == lastP0:
                    break
            P1 = S_ind.union(P0)
            Pbar1 = self.evaluate(F, P1)
            Pbar0 = Pbar1 - S_ind
            return P1, P0, Pbar1, Pbar0
        else:
            return "Failure"
        
    def get_neighours_cube(self, cube, kmax):
        """Takes the neighbours of a cube"""
        L = set()
        for some in itertools.product([-1,0,1], repeat=kmax):
            face = list(deepcopy(cube))
            for i,s in enumerate(some):
                face[i] = round(cube[i]+s*self.delta, 5)
            L.add(tuple(face))
        L.remove(cube)

        return L
    
    def cubicalwrap_cube(self, cube, kmax):
        """takes the neighbours of a cube"""
        L = set()
        for some in itertools.product([-1,0,1], repeat=kmax):
            face = list(deepcopy(cube))
            for i,s in enumerate(some):
                face[i] = round(cube[i]+s*self.delta, 5)
            L.add(tuple(face))
        return L

    def get_neighours_cubicalset(self, S):
        """takes the neighbours of a cubical set,
        same as collar in Computational Homology (Kaczynski, 2006)"""
        N = self.cubical_wrap(S)
        return N-set(S)

    def cubical_wrap(self, S):
        """takes the neighbours of a cubical set,
        same as collar"""
        #now only works for cubes of dim kmax
        maxdim = 0
        for cube in S:
            maxdim = max(get_dim(cube, self.delta), maxdim)


        N = set()
        for cube in S:
            L = self.cubicalwrap_cube(cube, maxdim)
#             print("CUBE", cube)
            for l in L:
#                 print("l", l)
                N.add(l)
        return N
        


    def restricted_map(self, N):
        """For the graph self.G takes the subgraph with nodes in list of nodes N"""
        return nx.subgraph(self.G, N)
    
    def convert_to_invertal_representation(self, cubes):
        sigdelta = str(self.delta)[::-1].find('.')+3
        intervalcubes = set()
        dim = len(list(cubes)[0])
#         print(dim)
        for cube in cubes:
            newcube = []
            for j in range(dim):
                interval = tuple([round(heaviside(cube[j])*self.delta*np.floor_divide(cube[j],self.delta), sigdelta),
                                  round(heaviside(cube[j])*self.delta*np.floor_divide(cube[j],self.delta) + self.delta, sigdelta)])
                newcube.append(interval)
            intervalcubes.add(tuple(newcube))
        return intervalcubes


def convert_to_chomp_format(cubical_set, delta):
    """Takes the cubical_set in coordindate representation to the representation for chomp
    (0.5, 1.5, 5.5)->[(0,1) (1,2) (5,6)] if delta=1"""
    filetxt = ""
    for cube in cubical_set:
        
        filetxt += '['
        for coord in cube:
            new = np.array(np.round((np.array(coord)/delta)),dtype='int')

                
            filetxt+=str(tuple(np.array(np.round((np.array(coord)/delta)),dtype='int')))+ ' '
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
    "calculates the boundary operator for an elementary cube"
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
    "returns E, the groups of cubical chains of a cubical set"
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
    """returns the vector representation of chain
    K is list of all elementary cubes for a certain dimension (can be calc by unrolledE)"""
    v = np.zeros(len(K))
    for i in range(len(K)):
        try:
            v[i] = chain[tuple(list(K)[i])]
        except:
            0
    return v


def boundaryOperatorMatrix(E, delta):
    "Calculates the boundary operator for a cubical chain E"
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

def get_bettiNumbers_of_cubicalset(tuplecubes, delta):
    """
    returns betti numbers up to maximal dimension of cubes in tuplecubes
    tuplecubes is the set of cubes as tuples of the coordinates of the center
    """
    E = cubicalChainGroups(tuplecubes, delta)
    D = boundaryOperatorMatrix(E, delta)
    bettinums = []
    for i in range(len(D)-1):
        bettinums.append(bettiNumber(i, D[i], D[i+1]))
    return bettinums
