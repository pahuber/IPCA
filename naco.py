'''
Performs iterative principal component analysis (IPCA) with ADI data from
the NACO instrument at VLT.
'''


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from ipca import PCA, IPCA
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
rank_ipca = 230

frame_pca = PCA(data_cube, rank_pca, parangs)
#frame_ipca = IPCA(data_cube, rank_ipca, parangs)
#np.savetxt("output/arrays/naco_" + str(rank_ipca) + ".txt", frame_ipca)
frame_ipca = np.loadtxt("output/arrays/naco_230.txt")





#'''plot V from Y vs. V from Y-D'''
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
#plt.savefig("output/naco_comparison.png", dpi = 400)
#plt.show()





'''plot data'''
plt.suptitle("NACO Exoplanet Imaging Data", y = 0.8)
plt.subplots_adjust(wspace=0.6)
font = 8
font_title = 10

plt.subplot(1, 3, 1)
plt.title("Unprocessed", fontsize = font_title)
plt.imshow(np.mean(data_cube, axis = 0), origin = "lower", extent = [-0.8235, 0.8235, -0.8235, 0.8235])
plt.ylabel("R. A. offset [arcsec]", fontsize = font)
plt.xlabel("Dec. offset [arcsec]", fontsize = font)
plt.tick_params(labelsize=font)
plt.colorbar(fraction = 0.045).ax.tick_params(labelsize=font)

plt.subplot(1, 3, 2)
plt.title("PCA (Rank " + str(rank_pca) + ")", fontsize = font_title)
plt.imshow(frame_pca, origin = "lower", extent = [-0.8235, 0.8235, -0.8235, 0.8235])
plt.xlabel("Dec. offset [arcsec]", fontsize = font)
plt.tick_params(labelsize=font)
plt.colorbar(fraction = 0.045).ax.tick_params(labelsize=font)

plt.subplot(1, 3, 3)
plt.title("IPCA (Rank " + str(rank_ipca) + ")", fontsize = font_title)
plt.imshow(frame_ipca, origin = "lower", vmax=7, extent = [-0.8235, 0.8235, -0.8235, 0.8235])
plt.xlabel("Dec. offset [arcsec]", fontsize = font)
plt.tick_params(labelsize=font)
plt.colorbar(fraction = 0.045).ax.tick_params(labelsize=font)

plt.savefig("output/naco_" + str(rank_ipca) + ".png", dpi=400)
plt.show()

print(time.time() - start)