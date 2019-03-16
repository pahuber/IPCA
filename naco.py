# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 13:43:21 2019

@author: philipp
"""


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from rpca import PCA, RPCA
import time


start = time.time()


'''get data'''
images = fits.open('input/stack100_rad1.6as.fits')
#images.info()
data_cube = images[0].data

#import angles list
parangs = []
with open("input/parang.txt") as file:
    counter = 0    
    for line in file.readlines():
        if counter != 0:
            parangs.append(float(line.replace("\n", "")))
        counter += 1




'''process data'''
rank_pca = 10
rank_rpca = 10

frame_pca = PCA(data_cube, rank_pca, parangs)
frame_rpca = RPCA(data_cube, rank_rpca, parangs)



'''plot data'''
plt.suptitle("NACO Exoplanet Imaging Data")
plt.subplots_adjust(hspace=0.5)

plt.subplot(2, 2, 1)
plt.title("Unprocessed")
plt.imshow(np.mean(data_cube, axis = 0), origin = "lower")
plt.colorbar()

plt.subplot(2, 2, 2)
plt.title("PCA (Rank " + str(rank_pca) + ")")
plt.imshow(frame_pca, origin = "lower")
plt.colorbar()

plt.subplot(2, 2, 3)
plt.title("RPCA (Rank " + str(rank_rpca) + ")")
plt.imshow(frame_rpca, origin = "lower", vmax=7)
plt.colorbar()

#plt.subplot(2, 2, 4)
#plt.title("Nothing yet")
#plt.imshow(processed_frame, origin = "lower")
#plt.colorbar()

plt.savefig("output/naco.png", dpi=400)
plt.show()

print(time.time() - start)