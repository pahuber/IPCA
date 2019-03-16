# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 13:43:21 2019

@author: philipp
"""


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from rlpca import cube2mat, PCA, RLPCA, mat2frame




'''get data'''
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




'''algorithm process'''
#reshape data cubes to matrices Y
Y = cube2mat(data_cube)

#apply RLPCA algorithm to calculate S_p
rank = 15
S_PCA = PCA(Y, rank)
S_p = RLPCA(Y, rank)

#reshape matrix to final time averaged frame
processed_frame_pca = mat2frame(S_PCA, parangs)
processed_frame = mat2frame(S_p, parangs)




'''plot'''
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

plt.savefig("output/naco.png", dpi=400)
plt.show()