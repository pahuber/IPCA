''' Performs iterative principal component analysis (IPCA) with ADI data of
    Beta Pic from the NACO instrument at VLT.
'''


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
#frame_ipca = IPCA(data_cube, rank_ipca, parangs)
#np.savetxt("output/arrays/pynpoint_example_" + str(rank_ipca) + ".txt", frame_ipca)

frame_ipca = np.loadtxt("output/arrays/pynpoint_example_" + str(rank_ipca) + ".txt")




'''plot data'''
plt.suptitle("Beta Pic b - NACO Exoplanet Imaging Data", y = 0.8)
plt.subplots_adjust(wspace=0.6)
font = 8
font_title = 10

plt.subplot(1, 3, 1)
plt.title("Unprocessed", fontsize = font_title)
plt.imshow(np.mean(data_cube, axis = 0), origin = "lower", extent = [-0.8235, 0.8235, -0.8235, 0.8235])
plt.ylabel("R. A. offset [arcsec]", fontsize = font)
plt.xlabel("Dec. offset [arcsec]", fontsize = font)
plt.tick_params(labelsize=font)
plt.colorbar(fraction = 0.0455).ax.tick_params(labelsize=font)

plt.subplot(1, 3, 2)
plt.title("PCA (Rank " + str(rank_pca) + ")", fontsize = font_title)
plt.imshow(frame_pca, origin = "lower", extent = [-0.8235, 0.8235, -0.8235, 0.8235])
plt.xlabel("Dec. offset [arcsec]", fontsize = font)
plt.tick_params(labelsize=font)
plt.colorbar(fraction = 0.0455).ax.tick_params(labelsize=font)

plt.subplot(1, 3, 3)
plt.title("IPCA (Rank " + str(rank_ipca) + ")", fontsize = font_title)
plt.imshow(frame_ipca, origin = "lower", vmax=13, extent = [-0.8235, 0.8235, -0.8235, 0.8235])
plt.xlabel("Dec. offset [arcsec]", fontsize = font)
plt.tick_params(labelsize=font)
plt.colorbar(fraction = 0.0455).ax.tick_params(labelsize=font)

plt.savefig("output/pynpoint_example_" + str(rank_ipca) + ".png", dpi=400)
plt.show()

print(time.time() - start)