# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 13:43:21 2019

@author: philipp
"""

'''TO DO: use gaussian noise in 2d, rotate planet, V komponenten visualisieren'''


import numpy as np
import scipy as sp
from scipy import ndimage
import matplotlib.pyplot as plt
import random as rd




'''create t x n x n data cube of t images consisting of star signal, planet signal and noise'''
def airy(x, y, x_offset, y_offset, i0):
    r = np.sqrt((x+x_offset)**2+(y+y_offset)**2)
    return i0*(2*sp.special.j1(r)/r)**2

frames = 25 #corresponds to t
boundary = 30 #resolution
steps = 80 #corresponds to n
max_ang = 270

cube_ideal = []
cube_real = []

for i in range(frames):

    star = []
    planet= []
    noise = []
    #current_ang = 10
    current_ang = np.linspace(0, max_ang, frames)[i]
    #print(current_ang)
    
    for x in np.linspace(boundary, -boundary, steps):
        temp = []
        for y in np.linspace(-boundary, boundary, steps):
            temp.append(airy(x, y, 0.00000001, 0.00000001, 20)) #if 0, then nan
        star.append(temp)
        
    for x in np.linspace(boundary, -boundary, steps):
        temp = []
        for y in np.linspace(-boundary, boundary, steps):
            temp.append(airy(x, y, 15, -14, 1))
        planet.append(temp)
    
    for x in np.linspace(boundary, -boundary, steps):
        temp = []
        for y in np.linspace(-boundary, boundary, steps):
            temp.append(rd.uniform(0, 0.8))
        noise.append(temp)
    
    frame_ideal = ndimage.rotate(np.array(planet), current_ang, reshape = False)
    frame_real = np.array(star) + ndimage.rotate(np.array(planet), current_ang, reshape = False) + np.array(noise)
    #frame_real = ndimage.rotate(frame_real, current_ang, reshape = False)
    
    cube_ideal.append(frame_ideal)
    cube_real.append(frame_real)    
    
cube_ideal = np.array(cube_ideal)
cube_real = np.array(cube_real)




'''functions'''
def cube2mat(A): #take t x n x n data cube and reshape it to t x n^2 matrix
   return A.reshape((len(A), len(A[0])**2))

def mat2frame(A): #take t x n^2 matrix and reshape it to t x n x n data cube, then take the mean to create final frame
    A_cube = np.reshape(A, (len(A), int(np.sqrt(len(A[0]))), int(np.sqrt(len(A[0])))))
    return np.mean(A_cube, axis = 0)

def trimatmul(A, B, C):
    return np.matmul(A, np.matmul(B, C))
     
def SVD(A):
    U, sigma, Vh = sp.linalg.svd(A)
    #create corresponding matrix Sigma from list sigma
    Sigma = np.zeros((len(A), len(A[0])))
    for i in range(len(sigma)):
        Sigma[i][i] = sigma[i]
    return U, Sigma, Vh

def LRA(A, rank):
    U, Sigma, Vh = SVD(A)    
    U_r = np.compress(np.ones(rank), U, axis = 1)
    Sigma_r = np.compress(np.ones(rank), np.compress(np.ones(rank), Sigma, axis = 0), axis = 1)
    Vh_r = np.compress(np.ones(rank), Vh, axis = 0)
    L = trimatmul(U_r, Sigma_r, Vh_r)
    #only keep positive parts of L, i.e. set negative parts of L to zero
    L = L.clip(min = 0)
    return L

def PCA(Y, rank):
    return Y - LRA(Y, rank)
    
def RLPCA(Y, p):   
    def S_r(rank):
        if rank == 0:
            return Y - LRA(Y, 1)
        return Y - LRA(Y - S_r(rank-1), rank) #implement rotation
    return S_r(p)

def RLPCA2(Y, p): #alternative calculation, yields different result ?!?!?
    S = Y - LRA(Y, 1)
    for i in range(1, p+1):
        S = Y - LRA(S, i)
    return S




'''algorithm process'''
#reshape data cubes to matrices Y
Y_ideal = cube2mat(cube_ideal)
Y_real = cube2mat(cube_real)

#apply RLPCA algorithm to calculate S_p
S_PCA = PCA(Y_real, 2)
S_p = RLPCA(Y_real, 2)

#reshape matrix to final time averaged frame
processed_frame = mat2frame(S_p)




'''make plot'''
plt.suptitle("Simulated Exoplanet")
plt.subplots_adjust(hspace=0.5)

plt.subplot(2, 2, 1)
plt.title("Ideal")
plt.imshow(cube_ideal[0], vmin=0, vmax=1)

plt.subplot(2, 2, 2)
plt.title("Real")
plt.imshow(cube_real[0], vmin=0, vmax=1)

plt.subplot(2, 2, 3)
plt.title("PCA (Rank 2)")
plt.imshow(mat2frame(S_PCA))

plt.subplot(2, 2, 4)
plt.title("RLPCA ($p=2$)")
plt.imshow(processed_frame)

plt.savefig("figure.png", dpi=400)