''' Performs iterative principal component analysis (IPCA) with simulated ADI
    data of an exoplanet.
'''


import numpy as np
import scipy as sp
from scipy import ndimage
import matplotlib.pyplot as plt
import random as rd
from ipca import PCA, IPCA, red, cube2mat
import time


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




'''algorithm process for simulated data'''
rank_pca = 10
rank_ipca = 10

frame_pca = PCA(cube_real, rank_pca, angles)
frame_ipca = IPCA(cube_real, rank_ipca, angles)



'''plot'''
plt.suptitle("Simulated Exoplanet")
plt.subplots_adjust(hspace=0.5)

plt.subplot(2, 2, 1)
plt.title("Ideal")
plt.imshow(cube_ideal[0], vmin=0, vmax=1)
plt.colorbar()

plt.subplot(2, 2, 2)
plt.title("Real")
plt.imshow(red(cube2mat(cube_real), angles), vmin=0, vmax=1)
plt.colorbar()

plt.subplot(2, 2, 3)
plt.title("PCA (Rank " + str(rank_pca) + ")")
plt.imshow(frame_pca)
plt.colorbar()

plt.subplot(2, 2, 4)
plt.title("IPCA (Rank " + str(rank_ipca) + ")", fontsize = 10)
plt.imshow(frame_ipca)
plt.colorbar()

plt.savefig("output/simulation.png", dpi=400)
plt.show()

print(time.time() - start)