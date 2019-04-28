'''
Contains the functions used to perform iterative principal components
analysis (IPCA).
'''


import numpy as np
import scipy as sp
from scipy import ndimage
import sys
from copy import copy


def cube2mat(A): #take t x n x n data cube and reshape it to t x n^2 matrix

   return A.reshape((len(A), len(A[0])**2))


def mat2cube(A): #take t x n^2 matrix and reshape it to t x n x n cube

    A_cube = np.reshape(A, (len(A), int(np.sqrt(len(A[0]))), int(np.sqrt(len(A[0])))))
    return A_cube #np.mean(A_cube, axis = 0)

    
def frame2cube(frame,
               original_cube): #takes a (PCA processed) n x n frame and returns a t x n x n cube with t copies of it, gets value of t from length of original cube

    t = len(original_cube)
    return np.stack([frame]*t, axis = 0)


def trimatmul(A,
              B,
              C):

    return np.matmul(A, np.matmul(B, C))

     
def SVD(A):

    U, sigma, Vh = sp.linalg.svd(A)
    #create corresponding matrix Sigma from list sigma
    Sigma = np.zeros((len(A), len(A[0])))
    for i in range(len(sigma)):
        Sigma[i][i] = sigma[i]
    return U, Sigma, Vh


def LRA(A,
        rank):
            
    U, Sigma, Vh = SVD(A)
    U_r = U[:, :rank]
    Sigma_r = Sigma[:rank, :rank]
    Vh_r = Vh[:rank, :]
    L = trimatmul(U_r, Sigma_r, Vh_r)
    return L


def red(S,
        angles = None): #takes t x n^2 matrix S, reshapes it to cube S_cube and rotates each frame if angles list is given and returns mean of cube, i.e. processed frame
    
    S_cube = mat2cube(S)    
    if angles is not None:
        for i in range(len(S_cube)):
            S_cube[i] = ndimage.rotate(S_cube[i], -1*angles[i], reshape = False)
    return np.mean(S_cube, axis = 0)

    
def theta(frame,
          original_cube,
          angles = None): #takes a (PCA processed) frame, sets negative parts of it to zero, reshapes it into t x n x n cube, rotates frames according to list and returns t x n^2 matrix
          
    d = frame.clip(min = 0)
    d_cube = frame2cube(d, original_cube)
    if angles is not None:
        for i in range(len(d_cube)):
            d_cube[i] = ndimage.rotate(d_cube[i], angles[i], reshape = False)
    D = cube2mat(d_cube)
    return D


def PCA(A_cube,
        rank,
        angles = None): #takes an unprocessed data cube and an angles list and returns PCA processed frame

    sys.stdout.write("Sarting PCA ... ")
    Y = cube2mat(A_cube)
    S = Y - LRA(Y, rank)
    sys.stdout.write("[DONE] ")
    return red(S, angles)


def IPCA(A_cube,
         pmax,
         pmin = 1,
         angles = None): #takes an unprocessed data cube, a max rank and an angles list and returns IPCA processed frame

    sys.stdout.write("Starting IPCA ... ")    
    Y = cube2mat(A_cube)
    S = Y - LRA(Y, pmin) #S_0
    for counter, i in enumerate(range(pmin, pmax+1)):
        S = Y - LRA(Y-theta(red(S, angles), A_cube, angles), i)
        sys.stdout.write(str(100*(counter+1)//len(range(pmin, pmax+1))) + "% Progress ")
        sys.stdout.flush()
    sys.stdout.write("[DONE] ")
    return red(S, angles)

    
def IPCA_savesteps(output_path,
                   prefix,
                   interval,
                   A_cube,
                   pmax,
                   pmin = 1,
                   angles = None): #takes an unprocessed data cube, a max rank and an angles list and returns IPCA processed frame

    sys.stdout.write("Starting IPCA ... ") 
    Y = cube2mat(A_cube)
    S = Y - LRA(Y, pmin) #S_0
    for counter, i in enumerate(range(pmin, pmax+1)):
        S = Y - LRA(Y-theta(red(S, angles), A_cube, angles), i)
        A = copy(S) #copy is necesarry! putting A=S will create a reference and yield change S when red(S, angles) is called for savetext
        if i > pmin:
            if (i == 1) or ((i % interval) == 0 ):
                np.savetxt(output_path + prefix + "_ipca_" + str(pmin) + "_" + str(i) + ".txt", red(A, angles))
        sys.stdout.write(str(100*(counter+1)//len(range(pmin, pmax+1))) + "% Progress ")
        sys.stdout.flush()
    sys.stdout.write("[DONE] ")
    return red(S, angles)

    
def printV(A_cube,
           p,
           angles = None): #IPCA funtion to visualize frame just containing the star

    Y = cube2mat(A_cube)
    S = Y - LRA(Y, 1) #S_0
    for i in range(1, p+1):
        V = Y - theta(red(S, angles), A_cube, angles)
    return V