# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 13:43:21 2019

@author: philipp
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import random as rd


'''create data consisting of star signal, planet signal and noise'''
def airy(x, y, x_offset, y_offset, i0):
    r = np.sqrt((x+x_offset)**2+(y+y_offset)**2)
    return i0*(2*sp.special.j1(r)/r)**2

star = []
planet= []
noise = []

boundary = 30
steps = 80

for x in np.linspace(boundary, -boundary, steps):
    temp = []
    for y in np.linspace(-boundary, boundary, steps):
        temp.append(airy(x, y, 0, 0, 10))
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

Y_ideal = np.asarray(planet)

Y_real = np.asarray(star) + np.asarray(planet) + np.asarray(noise)


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
    S_p = S_r(p)
    return S_p
    
#def S_r(Y, rank):
#    if rank == 0:
#        return Y - LRA(Y, 1)
#    return Y - LRA(Y - S_r(Y, rank-1), rank)
    
    
def tester(Y, p):
    S = Y - LRA(Y, 1)
    for i in range(1, p+1):
        S = Y - LRA(S, i)
    return S

test = np.asarray([[1, 2, 1], [2, 2, 3], [3, 1, 2]])
print(RLPCA(test, 2))
print(tester(test, 2))

S_pca = Y_real - LRA(Y_real, 2)
S_p = tester(Y_real, 3)
S_p2 = tester(Y_real, 73)


#print(Y_real - LRA(Y_real, 1))
#print("\n\n\n")
#print(RLPCA(Y_real, 1))
#print("\n\n\n")
#print(RLPCA(Y_real, 10))
#print("\n\n\n")
#print(RLPCA(Y_real, 30))
#print("\n\n\n")
#print(RLPCA(Y_real, 50))

                        
            

'''make plot'''
plt.suptitle("Simulated Exoplanet")
plt.subplot(2, 2, 1)
plt.title("Ideal")
plt.imshow(Y_ideal, vmin=0, vmax=1)

plt.subplot(2, 2, 2)
plt.title("w/ Noise")
plt.imshow(Y_real, vmin=0, vmax=1)

plt.subplot(2, 2, 3)
plt.title("PCA (Rank 2)")
plt.imshow(S_p2, vmin=0, vmax=1)

plt.subplot(2, 2, 4)
plt.title("RL PCA")
plt.imshow(S_p, vmin=0, vmax=1)

plt.savefig("fig.png", dpi=400)