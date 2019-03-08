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

frames = 25 #corresponds to t
boundary = 30 #resolution
steps = 80 #corresponds to n

cube_ideal = []
cube_real = []

for i in range(frames):

    star = []
    planet= []
    noise = []

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
    
    frame_ideal = np.asarray(planet)    
    frame_real = np.asarray(star) + np.asarray(planet) + np.asarray(noise)

    cube_ideal.append(frame_ideal)
    cube_real.append(frame_real)    
    
cube_ideal = np.asarray(cube_ideal)
cube_real = np.asarray(cube_real)


'''reshape data cube into t x n^2 matrix'''
def cube2mat(A):
    A_mat = []
    
    for i in range(len(A)):
        temp = []
        for j in A[i]:
            for k in j:
                temp.append(k)
        A_mat.append(temp)
    A_mat = np.asarray(A_mat)  
    return A_mat


'''PCA'''
def trimatmul(A, B, C):
    return np.matmul(A, np.matmul(B, C))
     
def SVD(A):
    U, sigma, Vh = sp.linalg.svd(A)
    #create matrix with same dimensions as A, sigma on diagonal and fill additional cols/rows with zeros if A is not quadratic
    rows = len(A)
    cols = len(A[0])   
    Sigma = np.diag(sigma)
    if rows < cols:
        Sigma_new = []
        for i in Sigma:
            i = i.tolist()
            for j in range(cols-rows):
                i.append(0)
            Sigma_new.append(i)
        Sigma = Sigma_new
    elif rows > cols:
        Sigma = Sigma.tolist()
        for i in range(rows-cols):
            temp = []
            for j in range(len(Sigma[0])):
                temp.append(0)
            Sigma.append(temp)
    Sigma = np.asarray(Sigma)
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
    #calculate L
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

#alternative calculation, yields different result ?!?!?
def RLPCA2(Y, p):
    S = Y - LRA(Y, 1)
    for i in range(1, p+1):
        S = Y - LRA(S, i)
    return S


'''reshape t x n^2 matrix back to time averaged n x n matrix, i. e. final processed frame'''
def mat2frame(A):
    #create matrix with dimensions of original image containing only zeros
    summed = []
    for i in range(steps):
        temp = []
        for j in range(steps):
            temp.append(0)
        summed.append(temp)
    summed = np.asarray(summed)
    #create a n x n matrix out of every row in the t x n^2 matrix and sum it to matrix summed
    temp_row = []
    temp_matrix = []
    counter = 0
    for row in A:
        for i in row:
            temp_row.append(i)
            counter += 1
            #print(temp_row)
            if counter == steps:
                temp_matrix.append(temp_row)
                temp_row = []
                counter = 0
        summed = np.add(summed, np.asarray(temp_matrix))
        temp_matrix = []
        
    processed_frame = 1/frames*summed
    return processed_frame


'''algorithm process'''
#reshape data cubes to matrices Y
Y_ideal = cube2mat(cube_ideal)
Y_real = cube2mat(cube_real)

#apply algorithm to calculate S_p
S_pca = Y_real - LRA(Y_real, 2)
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
plt.title("w/ Noise")
plt.imshow(cube_real[0], vmin=0, vmax=1)

plt.subplot(2, 2, 3)
plt.title("PCA (Rank 2)")
plt.imshow(mat2frame(S_pca))

plt.subplot(2, 2, 4)
plt.title("RLPCA ($p=2$)")
plt.imshow(processed_frame)

plt.savefig("fig.png", dpi=400)