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
stack = 2
rank_pca = 80 #maximum pca value
rank_ipca_init_list = [1, 3, 6, 10, 15] #different inital values
rank_ipca_end = 80 #max ipca value, must be larger than all values in init_list
interval = 5 #plot every 'interval' rank


'''declare input and output paths'''
input_path = "/home/philipp/Documents/BachelorProjectInput/disk/"
output_path = "output/stack" + str(stack) + "/"


'''get start time'''
start = time.time()


'''process data'''
fig, axes = plt.subplots(nrows=len(rank_ipca_init_list)+1, ncols=int(rank_pca/interval)+1)
plt.suptitle("PCA/IPCA Evolution Matrix (Stack " + stack + ")", fontsize = 8)
#plt.xlabel("R. A. Offset [arcsec]")
#plt.subplots_adjust(wspace=0.2)
subplot_counter = 0
size = 1
vmin = None
vmax = None
fontsize = 2

for i in range(1, rank_pca+1):
    if (i == 1) or ((i % interval) == 0):
        im = axes.flat[subplot_counter].imshow(np.loadtxt(output_path + "arrays_matrix/disk_pca_" + str(i) + ".txt"), origin='lower', vmin=vmin, vmax=vmax, extent=[size, -size, -size, size])
        axes.flat[subplot_counter].set_title("PCA (Rank " + str(i) + ")", fontsize=fontsize)    
        axes.flat[subplot_counter].tick_params(labelsize=fontsize, axis='both', left=False, bottom=False, colors=(1, 1, 1, 0))
        subplot_counter += 1
    
for init in rank_ipca_init_list:
    for i in range(subplot_counter, subplot_counter+(init // interval)+1):
        fig.delaxes(axes.flatten()[i])
        subplot_counter += 1
    for end in range(init+1, rank_ipca_end+1):
        if (end == 1) or ((end % interval) == 0):
            im = axes.flat[subplot_counter].imshow(np.loadtxt(output_path + "arrays_matrix/disk_ipca_" + str(init) + "_" + str(end) + ".txt"), origin='lower', vmin=vmin, vmax=vmax, extent=[size, -size, -size, size])
            axes.flat[subplot_counter].set_title("IPCA [" + str(init) + ", " + str(end) + "]", fontsize=fontsize)
            if init == 6 and end == 30:
                axes.flat[subplot_counter].set_title("IPCA [" + str(init) + ", " + str(end) + "]", fontsize=fontsize, color="green")    
            axes.flat[subplot_counter].tick_params(labelsize=fontsize, axis='both', left=False, bottom=False, colors=(1, 1, 1, 0))
            subplot_counter += 1

fig.tight_layout()
#plt.subplot_tool()
fig.colorbar(im, ax=axes.ravel().tolist())
plt.savefig(output_path +"plots/matrix.png", bbox_inches='tight', orientation="landscape", dpi=1000)



'''print runtime'''
print(time.time() - start)