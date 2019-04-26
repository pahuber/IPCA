'''
Creates fits file from (multiple) IPCA processed array(s).
'''


import numpy as np
from astropy.io import fits
import os


'''stack, pca rank and ipca ranks'''
stack = 100
rank_pca = 80 #maximum pca value
rank_ipca_init_list = [1, 3, 6, 10, 15, 20, 25] #different inital values
rank_ipca_end = rank_pca #max ipca value, must be larger than all values in init_list
interval = 5 #plot every 'interval' rank
arrays_folder = "pc1"


'''input and output paths'''
input_path = "output/stack" + str(stack) + "/" + "arrays_matrix/" + arrays_folder + "/"
output_path = "/home/philipp/Documents/BachelorProjectInput/pynpoint/" + arrays_folder + "/"


'''process data'''
for j in rank_ipca_init_list:
    im_shape = np.shape(np.loadtxt(input_path + "planets_ipca_" + str(rank_ipca_init_list[0]) + "_" + str(rank_ipca_end) + ".txt"))
    
    rank_list = []
    for i in range(j, rank_pca+1):       
        if (j != i) and ((i % interval) == 0):
            rank_list.append(i)
    
    data = np.zeros((len(rank_list), im_shape[0], im_shape[1])) #np.shape(np.loadtxt(output_path + "planets_ipca_" + str(rank_ipca_init) + "_" + str(rank_ipca_end) + ".txt")[0], np.shape(output_path + "planets_ipca_" + str(rank_ipca_init) + "_" + str(rank_ipca_end) + ".txt")[1]))
    
    for i, item in enumerate(rank_list):
        data[i] = np.loadtxt(input_path + "planets_ipca_" + str(j) + "_" + str(int(item)) + ".txt")
        
    hdu = fits.PrimaryHDU(data=data)
    hdu.writeto(os.path.join(output_path, "naco_ipca_" + str(j) + "_" + str(rank_ipca_end) + "_" + str(interval) + ".fits"))
   
   

im_shape = np.shape(np.loadtxt(input_path + "planets_pca_" + str(rank_ipca_end) + ".txt"))

rank_list = []
for i in range(1, rank_pca+1):       
    if (i == 1 ) or ((i % interval) == 0):
        rank_list.append(i)

data = np.zeros((len(rank_list), im_shape[0], im_shape[1])) #np.shape(np.loadtxt(output_path + "planets_ipca_" + str(rank_ipca_init) + "_" + str(rank_ipca_end) + ".txt")[0], np.shape(output_path + "planets_ipca_" + str(rank_ipca_init) + "_" + str(rank_ipca_end) + ".txt")[1]))

for i, item in enumerate(rank_list):
    data[i] = np.loadtxt(input_path + "planets_pca_" + str(int(item)) + ".txt")
   
hdu = fits.PrimaryHDU(data=data)
hdu.writeto(os.path.join(output_path, "naco_pca_" + str(rank_ipca_end) + "_" + str(interval) + ".fits"))