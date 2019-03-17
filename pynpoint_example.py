# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 13:43:21 2019

@author: philipp
"""


import numpy as np
import matplotlib.pyplot as plt
from ipca import PCA, IPCA
import time
import h5py

start = time.time()


'''get data'''
images = h5py.File("input/betapic_naco_mp.hdf5", "r")
#for key in images.keys():
#    print(key)

header = images["header_stack"]
stack = images["stack"]
#for key in stack.keys():
#    print(key)

parangs = header["PARANG"].value
data_cube = stack.value




'''process data'''
rank_pca = 20
rank_ipca = 263

frame_pca = PCA(data_cube, rank_pca, parangs)
frame_ipca = IPCA(data_cube, rank_ipca, parangs)



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
plt.title("IPCA (Rank " + str(rank_ipca) + ")")
plt.imshow(frame_ipca, origin = "lower", vmax=7)
plt.colorbar()

#plt.subplot(2, 2, 4)
#plt.title("Nothing yet")
#plt.imshow(processed_frame, origin = "lower")
#plt.colorbar()

plt.savefig("output/pynpoint_example.png", dpi=400)
plt.show()

print(time.time() - start)