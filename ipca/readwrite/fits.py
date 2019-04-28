'''
Contains functions to read fits files or create fits files from array data.
'''


import numpy as np
from astropy.io import fits
import os


def read_fits(input_path,
              name_fits,
              name_parangs):
    
    images = fits.open(input_path + name_fits)
    #images.info()
    data_cube = images[0].data
    
    parangs = []
    with open(input_path + name_parangs) as file:
        counter = 0    
        for line in file.readlines():
            if counter != 0:
                parangs.append(float(line.replace("\n", "")))
            counter += 1
    return data_cube, parangs


def create_multi_frame_fits(ipca_init_list,
                pca_ipca_end,
                interval,
                input_path,
                output_path,
                prefix,
                name_out="out"):
        
    im_shape = np.shape(np.loadtxt(input_path + prefix + "_ipca_" + str(ipca_init_list[0]) + "_" + str(pca_ipca_end) + ".txt"))    
            
    for j in ipca_init_list:
        rank_list = []
        for i in range(j, pca_ipca_end+1):       
            if (j != i) and ((i % interval) == 0):
                rank_list.append(i)
        
        data = np.zeros((len(rank_list), im_shape[0], im_shape[1]))
        
        for i, item in enumerate(rank_list):
            data[i] = np.loadtxt(input_path + prefix + "_ipca_" + str(j) + "_" + str(int(item)) + ".txt")
            
        hdu = fits.PrimaryHDU(data=data)
        hdu.writeto(os.path.join(output_path, prefix + "_ipca_" + str(j) + "_" + str(pca_ipca_end) + "_" + str(interval) + ".fits"))
    
    rank_list = []
    for i in range(1, pca_ipca_end+1):       
        if (i == 1 ) or ((i % interval) == 0):
            rank_list.append(i)
    
    data = np.zeros((len(rank_list), im_shape[0], im_shape[1]))    
    for i, item in enumerate(rank_list):
        data[i] = np.loadtxt(input_path + prefix + "_pca_" + str(int(item)) + ".txt")
       
    hdu = fits.PrimaryHDU(data=data)
    hdu.writeto(os.path.join(output_path, prefix + "_pca_" + str(pca_ipca_end) + "_" + str(interval) + ".fits"))
    

def create_single_frame_fits(ipca_init_list,
                pca_ipca_end,
                interval,
                input_path,
                output_path,
                prefix,
                name_out="out"):
            
    for j in ipca_init_list:
        rank_list = []
        for i in range(j, pca_ipca_end+1):       
            if (j != i) and ((i % interval) == 0):
                rank_list.append(i)
        
        for item in rank_list:
            data = np.loadtxt(input_path + prefix + "_ipca_" + str(j) + "_" + str(int(item)) + ".txt")
            hdu = fits.PrimaryHDU(data = data)
            hdu.writeto(os.path.join(output_path, prefix + "_ipca_" + str(j) + "_" + str(int(item)) + "_single.fits"))    
    
    rank_list = []
    for i in range(1, pca_ipca_end+1):       
        if (i == 1 ) or ((i % interval) == 0):
            rank_list.append(i)
            
    for item in rank_list:
        data = np.loadtxt(input_path + prefix + "_pca_" + str(int(item)) + ".txt")
        hdu = fits.PrimaryHDU(data = data)
        hdu.writeto(os.path.join(output_path, prefix + "_pca_" + str(pca_ipca_end) + "_single.fits"))    
    
    
def create_one_single_frame_fits(ipca_init,
                pca_ipca_end,
                interval,
                input_path,
                output_path,
                prefix,
                name_out="out"):
                    

    hdu = fits.PrimaryHDU(data=np.loadtxt(input_path + prefix + "_ipca_" + str(ipca_init) + "_" + str(pca_ipca_end) + ".txt"))
    hdu.writeto(os.path.join(output_path, prefix + "_ipca_" + str(ipca_init) + "_" + str(pca_ipca_end) + "_" + str(interval) + ".fits"))