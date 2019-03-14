# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 13:43:21 2019

@author: philipp
"""

'''TO DO: use gaussian noise in 2d, V komponenten visualisieren'''


import numpy as np
import scipy as sp
from scipy import ndimage
import matplotlib.pyplot as plt
import random as rd
from astropy.io import fits


a = None
print(a)
print(type(a))

'''create t x n x n data cube of t images consisting of star signal, planet signal and noise'''
def airy(x, y, x_offset, y_offset, i0):
    r = np.sqrt((x+x_offset)**2+(y+y_offset)**2)
    return i0*(2*sp.special.j1(r)/r)**2

frames = 25 #corresponds to t
boundary = 30 #resolution
steps = 80 #corresponds to n
angles = np.linspace(0, 70, frames)
#noise2 = random.multivariate_normal(np.array([0, 0]), np.array([[1, 0], [0, 10]]), [(steps, steps)])

cube_ideal = []
cube_real = []

for i in range(frames):

    star = []
    planet= []
    noise = []
    current_ang = angles[i]
    
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
    #frame_noise = np.array(noise)    
    
    cube_ideal.append(frame_ideal)
    cube_real.append(frame_real)    
    
cube_ideal = np.array(cube_ideal)
cube_real = np.array(cube_real)




'''functions'''
def cube2mat(A): #take t x n x n data cube and reshape it to t x n^2 matrix
   return A.reshape((len(A), len(A[0])**2))

def mat2frame(A, angles = None): #take t x n^2 matrix and reshape it to t x n x n data cube, then take the mean to create final frame
    A_cube = np.reshape(A, (len(A), int(np.sqrt(len(A[0]))), int(np.sqrt(len(A[0])))))
    #derotate each frame if angles list given
    if angles is not None:
        for i in range(len(A_cube)):
            A_cube[i] = ndimage.rotate(A_cube[i], -1*angles[i], reshape = False)
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
    S = Y - LRA(Y, 1)
    for i in range(1, p+1):
        S = Y - LRA(S, i)
    return S

def RLPCA2(Y, p): #alternative calculation, yields worse results   
    def S_r(rank):
        if rank == 0:
            return Y - LRA(Y, 1)
        return Y - LRA(Y - S_r(rank-1), rank)
    return S_r(p)




'''algorithm process for simulated data'''
#reshape data cubes to matrices Y
Y_ideal = cube2mat(cube_ideal)
Y_real = cube2mat(cube_real)

#apply RLPCA algorithm to calculate S_p
rank = 10
S_PCA = PCA(Y_real, rank)
S_p = RLPCA(Y_real, rank)

#reshape matrix to final time averaged frame
processed_frame = mat2frame(S_p, angles)




'''make plot for simulated data'''
plt.suptitle("Simulated Exoplanet")
plt.subplots_adjust(hspace=0.5)

plt.subplot(2, 2, 1)
plt.title("Ideal")
plt.imshow(cube_ideal[0], vmin=0, vmax=1)
plt.colorbar()

plt.subplot(2, 2, 2)
plt.title("Real")
plt.imshow(mat2frame(cube2mat(cube_real), angles), vmin=0, vmax=1)
plt.colorbar()

plt.subplot(2, 2, 3)
plt.title("PCA (Rank " + str(rank) + ")")
plt.imshow(mat2frame(S_PCA, angles))
plt.colorbar()

plt.subplot(2, 2, 4)
plt.title("RLPCA (Rank " + str(rank) + ")", fontsize = 10)
plt.imshow(processed_frame)
plt.colorbar()

plt.savefig("simulated.png", dpi=400)
plt.show()




'''get real data'''
images_1 = fits.open('data/PSFsub_mask0.05.fits')
images_2 = fits.open('data/stack100_rad1.6as.fits')
#images.info()
data_cube_1 = images_1[0].data
data_cube_2 = images_2[0].data

#import angles list
parangs = []
with open("data/parang.txt") as file:
    counter = 0    
    for line in file.readlines():
        if counter != 0:
            parangs.append(float(line.replace("\n", "")))
        counter += 1




'''algorithm process for real data'''
#reshape data cubes to matrices Y
Y_1 = cube2mat(data_cube_1)
Y_2 = cube2mat(data_cube_2)

#apply RLPCA algorithm to calculate S_p
rank = 10
S_PCA_1 = PCA(Y_1, rank)
S_p_1 = RLPCA(Y_1, rank)
S_PCA_2 = PCA(Y_2, rank)
S_p_2 = RLPCA(Y_2, rank)

#reshape matrix to final time averaged frame
processed_frame_1 = mat2frame(S_p_1)
processed_frame_2 = mat2frame(S_p_2)




'''make plot for real data'''
plt.suptitle("Real Exoplanet")
plt.subplots_adjust(hspace=0.5)

plt.subplot(2, 2, 1)
plt.title("PSFsub Mean")
plt.imshow(np.mean(data_cube_1, axis = 0))
plt.colorbar()

plt.subplot(2, 2, 2)
plt.title("PSFsub RLPCA")
plt.imshow(processed_frame_1)
plt.colorbar()

plt.subplot(2, 2, 3)
plt.title("Stack Mean")
plt.imshow(np.mean(data_cube_2, axis = 0))
plt.colorbar()

plt.subplot(2, 2, 4)
plt.title("Stack RLPCA")
plt.imshow(processed_frame_2)
plt.colorbar()

plt.savefig("real.png", dpi=400)
plt.show()