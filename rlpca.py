# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 11:49:40 2019

@author: philipp
"""

import numpy as np
import scipy as sp
from scipy import ndimage


'''TO DO: cli [:r, :], implement rotation theta, implement red'''


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
    U_r = U[:, :rank]
    Sigma_r = Sigma[:rank, :rank]
    Vh_r = Vh[:rank, :]
    #U_r = np.compress(np.ones(rank), U, axis = 1)
    #Sigma_r = np.compress(np.ones(rank), np.compress(np.ones(rank), Sigma, axis = 0), axis = 1)
    #Vh_r = np.compress(np.ones(rank), Vh, axis = 0)
    L = trimatmul(U_r, Sigma_r, Vh_r)
    #only keep positive parts of L, i.e. set negative parts of L to zero
    #L = L.clip(min = 0)
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