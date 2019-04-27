# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 10:59:33 2019

@author: philipp
"""


from astropy.io import fits
from ipca.processing.main import IPCA, PCA
from ipca.readwrite.fits import create_fits, read_fits, create_single_fits
from ipca.processing.processing import evaluate_multiple
from ipca.tool.plot import plot_matrix, plot_snr


'''input and output paths'''
input_path = "/home/philipp/Documents/BA_In_out/processed/test/stack100/arrays/pc1/"
output_path = "/home/philipp/Documents/BA_In_out/processed/test/stack100/fits/pc1/"


'''set variables'''
pca_ipca_end = 80 #maximum pca/ipca value, must be larger than all values in init_list
ipca_init_list = [6] #different inital values
interval = 5 #plot every 'interval' rank
stack = 100

name_fits = "stack" + str(stack) + "_rad1.6as.fits"
name_parangs = "parang_stack" + str(stack) + ".txt"


'''get data'''
#data_cube, parangs = read_fits(input_path,
#                          name_fits,
#                          name_parangs)
#                         
#frame_pca = PCA(data_cube, 10, parangs)
#frame_ipca = IPCA(data_cube, 10, 5, parangs)


'''calculate arrays'''
#mat_calc(data_cube, 
#         parangs,
#         ipca_init_list,
#         pca_ipca_end,
#         interval,
#         output_path,
#         prefix="out")


'''create evolution matrix plot'''
#plot_matrix(ipca_init_list,                         
#         pca_ipca_end,
#         interval,
#         input_path,
#         output_path,
#         prefix="challenge",
#         name_out="challenge",
#         title="PCA/IPCA Evolution Matrix",
#         vmin=None,
#         vmax=None)


'''create fits files'''
#create_fits(ipca_init_list,
#            pca_ipca_end,
#            interval,
#            input_path,
#            output_path,
#            prefix="challenge",
#            name_out="challenge")


'''create single fits file'''
create_single_fits(6,
            80,
            interval,
            input_path,
            output_path,
            prefix="planets",
            name_out="planets")



'''create SNR plots'''
#plot_snr(ipca_init_list,
#         pca_ipca_end,
#         interval,
#         input_path,
#         output_path,
#         array_file_name="outer.txt",
#         name_out="snr_outer",
#         title="SNR vs. Rank")