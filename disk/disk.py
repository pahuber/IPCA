'''
Performs iterative principal component analysis (IPCA) with ADI data.
'''


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import time
import os
os.chdir("..")
from ipca import PCA, IPCA
os.chdir("disk")


'''decalre stack, pca rank and ipca ranks'''
stack = 10
rank_pca = 10
rank_ipca_init = 1
rank_ipca_end = 2


'''declare input and output paths'''
input_path = "/home/philipp/Documents/BachelorProjectInput/disk/"
output_path = "output/stack" + str(stack) + "/"


'''get start time'''
start = time.time()


'''get data'''
images = fits.open(input_path + "stack" + str(stack) + ".fits")
#images.info()
data_cube = images[0].data

#import angles list
parangs = []
with open(input_path + "parang_stack" + str(stack) + ".txt") as file:
    counter = 0    
    for line in file.readlines():
        if counter != 0:
            parangs.append(float(line.replace("\n", "")))
        counter += 1


'''process data'''
frame_pca = PCA(data_cube, rank_pca, parangs)
np.savetxt(output_path + "arrays/disk_pca_" + str(rank_pca) + ".txt", frame_pca)
#frame_ipca = np.loadtxt(output_path + "arrays/disk_pca_" + str(rank_pca) + ".txt")

frame_ipca = IPCA(data_cube, rank_ipca_end, rank_ipca_init, parangs)
np.savetxt(output_path + "arrays/disk_ipca_" + str(rank_ipca_init) + "_" + str(rank_ipca_end) + ".txt", frame_ipca)
#frame_ipca = np.loadtxt(output_path + "arrays/disk_ipca_" + str(rank_ipca_init) + "_" + str(rank_ipca_end) + ".txt")


##'''plot V from Y vs. V from Y-D'''
#image_number = 10
#
##old V
#V_old = SVD(cube2mat(data_cube))[2].T
#V_old_cube = mat2cube(V_old)
#V_old_frame = V_old_cube[image_number]
#
##increased V
#V = SVD(IPCA_V(data_cube, rank_ipca, parangs))[2].T # from SVD(Y - theta(...))
#V_cube = mat2cube(V)
#V_frame = V_cube[image_number]
#
##plot
#plt.suptitle("Comparison of V built from different Matrices")
#plt.subplot(1, 2, 1)
#plt.title("V from Y")
#plt.imshow(V_frame)
#plt.colorbar()
#plt.subplot(1, 2, 2)
#plt.title("V from Y - D")
#plt.imshow(V_old_frame)
#plt.colorbar()
#plt.savefig(output_path + "disk_comparison.png", dpi = 400)
#plt.show()


'''plot data'''
plt.suptitle("Disk System", y = 0.8)
plt.subplots_adjust(wspace=0.6)
font = 8
font_title = 10

plt.subplot(1, 3, 1)
plt.title("Unprocessed", fontsize = font_title)
plt.imshow(np.mean(data_cube, axis = 0), origin = "lower", extent = [-0.8235, 0.8235, -0.8235, 0.8235])
plt.ylabel("Dec. offset [arcsec]", fontsize = font)
plt.xlabel("R. A. offset [arcsec]", fontsize = font)
plt.tick_params(labelsize=font)
plt.colorbar(fraction = 0.045).ax.tick_params(labelsize=font)

plt.subplot(1, 3, 2)
plt.title("PCA (Rank " + str(rank_pca) + ")", fontsize = font_title)
plt.imshow(frame_pca, origin = "lower", extent = [-0.8235, 0.8235, -0.8235, 0.8235])
plt.xlabel("R. A. offset [arcsec]", fontsize = font)
plt.tick_params(labelsize=font)
plt.colorbar(fraction = 0.045).ax.tick_params(labelsize=font)

plt.subplot(1, 3, 3)
plt.title("IPCA [" + str(rank_ipca_init) + ", " + str(rank_ipca_end) + "]", fontsize = font_title)
plt.imshow(frame_ipca, origin = "lower", vmax = 7, extent = [-0.8235, 0.8235, -0.8235, 0.8235])
plt.xlabel("R. A. offset [arcsec]", fontsize = font)
plt.tick_params(labelsize=font)
plt.colorbar(fraction = 0.045).ax.tick_params(labelsize=font)

plt.savefig(output_path + "plots/disk_" + str(rank_pca) + "_" + str(rank_ipca_init) + "_" + str(rank_ipca_end) + ".png", dpi=400)
plt.show()


'''print runtime'''
print(time.time() - start)