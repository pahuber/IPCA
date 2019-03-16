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
from rlpca import *



#'''create t x n x n data cube of t images consisting of star signal, planet signal and noise'''
#def airy(x, y, x_offset, y_offset, i0):
#    r = np.sqrt((x+x_offset)**2+(y+y_offset)**2)
#    return i0*(2*sp.special.j1(r)/r)**2
#
#frames = 25 #corresponds to t
#boundary = 30 #resolution
#steps = 80 #corresponds to n
#angles = np.linspace(0, 70, frames)

#cube_ideal = []
#cube_real = []
#
#for i in range(frames):
#
#    star = []
#    planet= []
#    noise = []
#    current_ang = angles[i]
#    
#    for x in np.linspace(boundary, -boundary, steps):
#        temp = []
#        for y in np.linspace(-boundary, boundary, steps):
#            temp.append(airy(x, y, 0.00000001, 0.00000001, 20)) #if 0, then nan
#        star.append(temp)
#        
#    for x in np.linspace(boundary, -boundary, steps):
#        temp = []
#        for y in np.linspace(-boundary, boundary, steps):
#            temp.append(airy(x, y, 15, -14, 1))
#        planet.append(temp)
#    
#    for x in np.linspace(boundary, -boundary, steps):
#        temp = []
#        for y in np.linspace(-boundary, boundary, steps):
#            temp.append(rd.uniform(0, 0.8))
#        noise.append(temp)
#        
#    frame_ideal = ndimage.rotate(np.array(planet), current_ang, reshape = False)
#    frame_real = np.array(star) + ndimage.rotate(np.array(planet), current_ang, reshape = False) + np.array(noise)
#    #frame_noise = np.array(noise)    
#    
#    cube_ideal.append(frame_ideal)
#    cube_real.append(frame_real)    
#    
#cube_ideal = np.array(cube_ideal)
#cube_real = np.array(cube_real)



#
#'''algorithm process for simulated data'''
##reshape data cubes to matrices Y
#Y_ideal = cube2mat(cube_ideal)
#Y_real = cube2mat(cube_real)
#
##apply RLPCA algorithm to calculate S_p
#rank = 10
#S_PCA = PCA(Y_real, rank)
#S_p = RLPCA(Y_real, rank)
#
##reshape matrix to final time averaged frame
#processed_frame = mat2frame(S_p, angles)
#
#
#
#
#'''make plot for simulated data'''
#plt.suptitle("Simulated Exoplanet")
#plt.subplots_adjust(hspace=0.5)
#
#plt.subplot(2, 2, 1)
#plt.title("Ideal")
#plt.imshow(cube_ideal[0], vmin=0, vmax=1)
#plt.colorbar()
#
#plt.subplot(2, 2, 2)
#plt.title("Real")
#plt.imshow(mat2frame(cube2mat(cube_real), angles), vmin=0, vmax=1)
#plt.colorbar()
#
#plt.subplot(2, 2, 3)
#plt.title("PCA (Rank " + str(rank) + ")")
#plt.imshow(mat2frame(S_PCA, angles))
#plt.colorbar()
#
#plt.subplot(2, 2, 4)
#plt.title("RLPCA (Rank " + str(rank) + ")", fontsize = 10)
#plt.imshow(processed_frame)
#plt.colorbar()
#
#plt.savefig("output/simulated.png", dpi=400)
#plt.show()
#



'''get real data'''
images = fits.open('data/stack100_rad1.6as.fits')
#images.info()
data_cube = images[0].data

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
Y = cube2mat(data_cube)

#apply RLPCA algorithm to calculate S_p
rank = 10
S_PCA = PCA(Y, rank)
S_p = RLPCA(Y, rank)

#reshape matrix to final time averaged frame
processed_frame_pca = mat2frame(S_PCA, parangs)
processed_frame = mat2frame(S_p, parangs)




'''make plot for real data'''
plt.suptitle("NACO Exoplanet Imaging Data")
plt.subplots_adjust(hspace=0.5)

plt.subplot(2, 2, 1)
plt.title("Unprocessed")
plt.imshow(np.mean(data_cube, axis = 0), origin = "lower")
plt.colorbar()

plt.subplot(2, 2, 2)
plt.title("PCA (Rank " + str(rank) + ")")
plt.imshow(processed_frame_pca, origin = "lower")
plt.colorbar()

plt.subplot(2, 2, 3)
plt.title("RLPCA (Rank " + str(rank) + ")")
plt.imshow(processed_frame, origin = "lower")
plt.colorbar()

plt.subplot(2, 2, 4)
plt.title("Nothing yet")
plt.imshow(processed_frame, origin = "lower")
plt.colorbar()

plt.savefig("output/real.png", dpi=400)
plt.show()