# -*- coding: utf-8 -*-
"""
Created on Sat Mar  2 22:29:37 2019

@author: philipp
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import random as rd


'''define matrix Y'''
A = np.array([[1, 0, 3, 10], [0, 6, 7, 8], [0, 10, 0, 12]])


'''compute SVD'''
U, sigma, Vh = sp.linalg.svd(A)

Sigma = np.array([[sigma[0], 0, 0, 0], [0, sigma[1], 0, 0], [0, 0, sigma[2], 0]])

print(A)


'''compute low rank approximation'''
#choose 2x2 matrix as Sigma_1
Sigma_1 = np.array([[Sigma[0][0], 0], [0, Sigma[1][1]]])

U_1 = np.array([[U[0][0], U[0][1]], [U[1][0], U[1][1]], [U[2][0], U[2][1]]])

Vh_1 = np.array([[Vh[0][0], Vh[0][1], Vh[0][2], Vh[0][3]], [Vh[1][0], Vh[1][1], Vh[1][2], Vh[1][3]]])

B = np.matmul(U_1, np.matmul(Sigma_1, Vh_1))

print(A-B)

plt.subplot(1, 3, 1)
plt.title("A")
plt.imshow(A, vmin=0, vmax=12)

plt.subplot(1, 3, 2)
plt.title("B")
plt.imshow(B, vmin=0, vmax=12)

plt.subplot(1, 3, 3)
plt.title("A-B")
plt.imshow(A-B, vmin=0, vmax=12)
plt.colorbar()

#plt.savefig("principle.png", dpi=300)