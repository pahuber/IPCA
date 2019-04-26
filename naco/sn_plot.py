'''
Performs iterative principal component analysis (IPCA) with ADI data.
'''


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os


'''stack, pca rank and ipca ranks'''
stack = 100
rank_pca = 5
rank_ipca_init = 1
rank_ipca_end = 3


'''input and output paths'''
input_path = "/home/philipp/Documents/BachelorProjectInput/pynpoint/pc1/"
output_path = "output/stack" + str(stack) + "/"


'''get data'''
ranks = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80]
snr_inner = np.loadtxt(input_path + "test_inner.txt")
snr_outer = np.loadtxt(input_path + "test_outer.txt")

lst_inner = []
lst_outer = []

for i in snr_inner:
    lst_inner.append(i[4])

for i in snr_outer:
    lst_outer.append(i[4])

'''plot data'''
plt.plot(ranks, lst_inner, "b.")
plt.plot(ranks, lst_outer, "r.")
plt.plot(ranks, lst_inner, "b--")
plt.plot(ranks, lst_outer, "r--")
plt.show()

#plt.suptitle("Planetary System", y = 0.8)
#plt.subplots_adjust(wspace=0.6)
#font = 8
#font_title = 10
#
#plt.subplot(1, 3, 1)
#plt.title("Unprocessed", fontsize = font_title)
#plt.imshow(np.mean(data_cube, axis = 0), origin = "lower", extent = [-0.8235, 0.8235, -0.8235, 0.8235])
#plt.ylabel("Dec offset [arcsec]", fontsize = font)
#plt.xlabel("R. A. offset [arcsec]", fontsize = font)
#plt.tick_params(labelsize=font)
#plt.colorbar(fraction = 0.045).ax.tick_params(labelsize=font)
#
#plt.subplot(1, 3, 2)
#plt.title("PCA (Rank " + str(rank_pca) + ")", fontsize = font_title)
#plt.imshow(frame_pca, origin = "lower", extent = [-0.8235, 0.8235, -0.8235, 0.8235])
#plt.xlabel("R. A. offset [arcsec]", fontsize = font)
#plt.tick_params(labelsize=font)
#plt.colorbar(fraction = 0.045).ax.tick_params(labelsize=font)
#
#plt.subplot(1, 3, 3)
#plt.title("IPCA [" + str(rank_ipca_init) + ", " + str(rank_ipca_end) + "]", fontsize = font_title)
#plt.imshow(frame_ipca, origin = "lower", extent = [-0.8235, 0.8235, -0.8235, 0.8235])
#plt.xlabel("R. A. offset [arcsec]", fontsize = font)
#plt.tick_params(labelsize=font)
#plt.colorbar(fraction = 0.045).ax.tick_params(labelsize=font)
#
#plt.savefig(output_path + "plots/planets_" + str(rank_pca) + "_" + str(rank_ipca_init) + "_" + str(rank_ipca_end) + ".png", dpi=400)
#plt.show()