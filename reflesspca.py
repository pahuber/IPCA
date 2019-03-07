# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 13:43:21 2019

@author: philipp
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import random as rd


'''create t x n x n data cube of t images consisting of star signal, planet signal and noise'''
def airy(x, y, x_offset, y_offset, i0):
    r = np.sqrt((x+x_offset)**2+(y+y_offset)**2)
    return i0*(2*sp.special.j1(r)/r)**2

frames = 10 #corresponds to t
boundary = 30 #corresponds to n
steps = 80 #resolution

cube_ideal = []
cube_real = []

for i in range(frames):

    star = []
    planet= []
    noise = []

    for x in np.linspace(boundary, -boundary, steps):
        temp = []
        for y in np.linspace(-boundary, boundary, steps):
            temp.append(airy(x, y, 0.00000001, 0.00000001, 10)) #if 0, then nan
        star.append(temp)
        
    for x in np.linspace(boundary, -boundary, steps):
        temp = []
        for y in np.linspace(-boundary, boundary, steps):
            temp.append(airy(x, y, 15, -14, 1))
        planet.append(temp)
    
    for x in np.linspace(boundary, -boundary, steps):
        temp = []
        for y in np.linspace(-boundary, boundary, steps):
            temp.append(rd.uniform(0, 0.3))
        noise.append(temp)
    
    frame_ideal = np.asarray(planet)    
    frame_real = np.asarray(star) + np.asarray(planet) + np.asarray(noise)

    cube_ideal.append(frame_ideal)
    cube_real.append(frame_real)    
    
cube_ideal = np.asarray(cube_ideal)
cube_real = np.asarray(cube_real)


'''reshape data cube into t x n^2 matrix'''
Y_ideal = []
Y_real = []

for i in range(len(cube_ideal)):
    temp = []
    for j in cube_ideal[i]:
        for k in j:
            temp.append(k)
    Y_ideal.append(temp)

for i in range(len(cube_real)):
    temp = []
    for j in cube_real[i]:
        for k in j:
            temp.append(k)
    Y_real.append(temp)
    
Y_ideal = np.asarray(Y_ideal)   
Y_real = np.asarray(Y_real)   


'''PCA'''
def trimatmul(A, B, C):
    return np.matmul(A, np.matmul(B, C))
     
def SVD(A):
    U, sigma, Vh = sp.linalg.svd(A)
    Sigma = np.diag(sigma) #only works that simple for QUADRATIC matrices!
    return U, Sigma, Vh

def LRA(A, rank):
    U, Sigma, Vh = SVD(A)
    #calculate Sigma_r
    Sigma_r = []
    for i in range(rank):
        temp = []
        for j in range(rank):
            temp.append(Sigma[i][j])
        Sigma_r.append(temp)
    Sigma_r = np.asarray(Sigma_r)
    #calculate U_r
    U_r = []
    m = len(U)
    for i in range(m):
        temp = []
        for j in range(rank):
            temp.append(U[i][j])
        U_r.append(temp)
    U_r = np.asarray(U_r)
    #calculate Vh_r
    Vh_r = []
    n = len(Vh[0])
    for i in range(rank):
        temp = []
        for j in range(n):
            temp.append(Vh[i][j])
        Vh_r.append(temp)
    Vh_r = np.asarray(Vh_r)
    
    L = trimatmul(U_r, Sigma_r, Vh_r)
    
    #set negative parts of L to zero
    for i in range(len(L)):
        for j in range(len(L[0])):
            if L[i][j] < 0:
                L[i][j] = 0
    return L
    
    
'''reference-less PCA'''    
def RLPCA(Y, p):   
    def S_r(rank):
        if rank == 0:
            return Y - LRA(Y, 1)
        return Y - LRA(Y - S_r(rank-1), rank)
    return S_r(p)
    
def RLPCA2(Y, p):
    S = Y - LRA(Y, 1)
    for i in range(1, p+1):
        S = Y - LRA(S, i)
    return S

S_pca = Y_real - LRA(Y_real, 2)
S_p = RLPCA(Y_real, 3)


'''make plot'''
plt.suptitle("Simulated Exoplanet")
plt.subplots_adjust(hspace=0.5)

plt.subplot(2, 2, 1)
plt.title("Ideal")
plt.imshow(Y_ideal, vmin=0, vmax=1)

plt.subplot(2, 2, 2)
plt.title("w/ Noise")
plt.imshow(Y_real, vmin=0, vmax=1)

plt.subplot(2, 2, 3)
plt.title("PCA (Rank 2)")
plt.imshow(S_pca, vmin=0, vmax=1)

plt.subplot(2, 2, 4)
plt.title("RLPCA")
plt.imshow(S_p, vmin=0, vmax=1)

plt.savefig("fig.png", dpi=400)