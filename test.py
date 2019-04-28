# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 10:59:33 2019

@author: philipp
"""


from astropy.io import fits
from ipca.processing.main import IPCA, PCA
from ipca.readwrite.fits import read_fits, create_multi_frame_fits, create_single_frame_fits, create_one_single_frame_fits
from ipca.processing.processing import evaluate_multiple
from ipca.tool.plot import plot_matrix, plot_snr


'''input and output paths'''
input_path = "/home/philipp/Documents/BA_In_out/processed/planets/stack100/snr/pc1/"
output_path = "/home/philipp/Documents/BA_In_out/processed/planets/stack100/plots/pc1/"


'''set variables'''
pca_ipca_end = 80 #maximum pca/ipca value, must be larger than all values in init_list
ipca_init_list = [1, 3, 6, 10, 15, 20, 25] #different inital values
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


'''create multiple multi-frame fits files'''
#create_multi_frame_fits(ipca_init_list,
#            pca_ipca_end,
#            interval,
#            input_path,
#            output_path,
#            prefix="challenge",
#            name_out="challenge")


'''create multiple single-frame fits files (use this for snr plots)'''
#create_single_frame_fits(ipca_init_list,
#            pca_ipca_end,
#            interval,
#            input_path,
#            output_path,
#            prefix="planets",
#            name_out="planets")


'''create one single-frame fits file'''
#create_one_single_frame_fits(6,
#            80,
#            interval,
#            input_path,
#            output_path,
#            prefix="planets",
#            name_out="planets")


'''create SNR plots'''
plot_snr(ipca_init_list,
         pca_ipca_end,
         interval,
         input_path,
         output_path,
         array_prefix="inner",
         name_out="snr_inner",
         title="SNR vs. Rank")