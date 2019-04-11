'''
Performs iterative principal component analysis (IPCA) with simulated ADI
data of an exoplanet.
'''


import numpy as np
import scipy as sp
from scipy import ndimage
import matplotlib.pyplot as plt
import random as rd
from ipca import PCA, IPCA, red, cube2mat
import time


'''declare input and output paths'''
#input_path = "/home/philipp/Documents/BachelorProjectInput/planets/"
output_path = "output/"


'''get start time'''
start = time.time()


'''create t x n x n data cube of t images consisting of star signal, planet signal and noise'''
def airy(x, y, x_offset, y_offset, i0):
    r = np.sqrt((x+x_offset)**2+(y+y_offset)**2)
    return i0*(2*sp.special.j1(r)/r)**2

frames = 50 #corresponds to t
boundary = 30 #resolution
steps = 80 #corresponds to n
angles = np.linspace(0, 70, frames)

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


'''process data'''
rank_pca = 10

rank_ipca_init = 5
rank_ipca_end = 20

frame_pca = PCA(cube_real, rank_pca, angles)
np.savetxt(output_path + "arrays/simulation_pca_" + str(rank_pca) + ".txt", frame_pca)
frame_ipca = np.loadtxt(output_path + "arrays/simulation_pca_" + str(rank_pca) + ".txt")

frame_ipca = IPCA(cube_real, rank_ipca_end, rank_ipca_init, angles)
np.savetxt(output_path + "arrays/simulation_ipca_" + str(rank_ipca_init) + "_" + str(rank_ipca_end) + ".txt", frame_ipca)
frame_ipca = np.loadtxt(output_path + "arrays/simulation_ipca_" + str(rank_ipca_init) + "_" + str(rank_ipca_end) + ".txt")


'''plot data'''
plt.suptitle("Simulated Exoplanet Imaging Data", y = 0.8)
plt.subplots_adjust(wspace=0.6)
font = 8
font_title = 10

plt.subplot(1, 3, 1)
plt.title("Unprocessed", fontsize = font_title)
plt.imshow(np.mean(cube_real, axis = 0), origin = "lower", extent = [-0.8235, 0.8235, -0.8235, 0.8235])
plt.ylabel("Dec. offset [arcsec]", fontsize = font)
plt.xlabel("R. A. offset [arcsec]", fontsize = font)
plt.tick_params(labelsize=font)
plt.colorbar(fraction = 0.0455).ax.tick_params(labelsize=font)

plt.subplot(1, 3, 2)
plt.title("PCA (Rank " + str(rank_pca) + ")", fontsize = font_title)
plt.imshow(frame_pca, origin = "lower", extent = [-0.8235, 0.8235, -0.8235, 0.8235])
plt.xlabel("R. A. offset [arcsec]", fontsize = font)
plt.tick_params(labelsize=font)
plt.colorbar(fraction = 0.0455).ax.tick_params(labelsize=font)

plt.subplot(1, 3, 3)
plt.title("IPCA [" + str(rank_ipca_init) + ", " + str(rank_ipca_end) + "]", fontsize = font_title)
plt.imshow(frame_ipca, origin = "lower", extent = [-0.8235, 0.8235, -0.8235, 0.8235])
plt.xlabel("R. A. offset [arcsec]", fontsize = font)
plt.tick_params(labelsize=font)
plt.colorbar(fraction = 0.0455).ax.tick_params(labelsize=font)

plt.savefig("output/simulation_" + str(rank_pca) + "_" + str(rank_ipca_init) + "_" + str(rank_ipca_end) + ".png", dpi=400)
plt.show()


'''print runtime'''
print(time.time() - start)