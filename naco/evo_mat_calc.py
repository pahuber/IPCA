'''
Calculates arrays to create evolution matrix of IPCA processed data.
'''


import numpy as np
from astropy.io import fits
import os
os.chdir("..")
from ipca import PCA, IPCA, IPCA_savesteps
os.chdir("planets")


'''stack, pca rank and ipca ranks'''
stack = 100
rank_pca = 40 #maximum pca value
rank_ipca_init_list = [3, 4, 5, 6] #different inital values
rank_ipca_end = 40 #max ipca value, must be larger than all values in init_list
interval = 5 #plot every 'interval' rank


'''input and output paths'''
input_path = "/home/philipp/Documents/BachelorProjectInput/planets/"
output_path = "output/stack" + str(stack) + "/arrays_matrix/"


'''get data'''
images = fits.open(input_path + "stack" + str(stack) + "_rad1.6as.fits")
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

cube_firsthalf = data_cube[:len(data_cube)//2, :]
parangs_firsthalf = parangs[:len(parangs)//2]

cube_everysecond = data_cube[0::2]
parangs_everysecond = parangs[0::2]


'''process data'''
for i in range(1, rank_pca+1):       
    if (i == 1) or ((i % interval) == 0):
        frame_pca = PCA(data_cube, i, parangs)
        np.savetxt(output_path + "planets_pca_" + str(i) + ".txt", frame_pca)
        

for init in rank_ipca_init_list:
    frame_ipca = IPCA_savesteps(output_path, "planets", interval, data_cube, rank_ipca_end, init, parangs)