'''
Contains functions to process data.
'''


import numpy as np
import matplotlib.pyplot as plt
from ipca.processing.main import PCA, IPCA_savesteps


def evaluate_multiple(data,
             parangs,
             ipca_init_list,
             pca_ipca_end,
             interval,
             output_path,
             prefix):
                 
    for i in range(1, pca_ipca_end+1):       
        if (i == 1) or ((i % interval) == 0):
            frame_pca = PCA(data, i, parangs)
            np.savetxt(output_path + prefix + "_pca_" + str(i) + ".txt", frame_pca)
        
    for i in ipca_init_list:
        IPCA_savesteps(output_path, prefix, interval, data, pca_ipca_end, i, parangs)