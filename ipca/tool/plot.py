'''
Contains functions to create various kinds of plots.
'''


import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler


def plot_matrix(ipca_init_list,                         
             pca_ipca_end,
             interval,
             input_path,
             output_path,
             prefix,
             name_out="out",
             title="PCA/IPCA Evolution Matrix",
             vmin=None,
             vmax=None):
                 
    fig, axes = plt.subplots(nrows=len(ipca_init_list)+1, ncols=int(pca_ipca_end/interval)+1)
    plt.suptitle(title, fontsize = 8)
    #plt.xlabel("R. A. Offset [arcsec]")
    #plt.subplots_adjust(wspace=0.2)
    subplot_counter = 0
    size = 1
    fontsize = 2
    
    for i in range(1, pca_ipca_end+1):
        if (i == 1) or ((i % interval) == 0):
            im = axes.flat[subplot_counter].imshow(np.loadtxt(input_path + prefix + "_pca_" + str(i) + ".txt"), origin='lower', vmin=vmin, vmax=vmax, extent=[size, -size, -size, size])
            axes.flat[subplot_counter].set_title("PCA [" + str(i) + "]", fontsize=fontsize)    
            axes.flat[subplot_counter].tick_params(labelsize=fontsize, axis='both', left=False, bottom=False, colors=(1, 1, 1, 0))
            subplot_counter += 1
        
    for init in ipca_init_list:
        for i in range(subplot_counter, subplot_counter+(init // interval)+1):
            fig.delaxes(axes.flatten()[i])
            subplot_counter += 1
        for end in range(init+1, pca_ipca_end+1):
            if (end == 1) or ((end % interval) == 0):
                im = axes.flat[subplot_counter].imshow(np.loadtxt(input_path + prefix + "_ipca_" + str(init) + "_" + str(end) + ".txt"), origin='lower', vmin=vmin, vmax=vmax, extent=[size, -size, -size, size])
                axes.flat[subplot_counter].set_title("IPCA [" + str(init) + ", " + str(end) + "]", fontsize=fontsize)
                #if init == 6 and end == 30:
                #    axes.flat[subplot_counter].set_title("IPCA [" + str(init) + ", " + str(end) + "]", fontsize=fontsize, color="green")    
                axes.flat[subplot_counter].tick_params(labelsize=fontsize, axis='both', left=False, bottom=False, colors=(1, 1, 1, 0))
                subplot_counter += 1
    
    #fig.tight_layout()
    #plt.subplot_tool()
    fig.colorbar(im, ax=axes.ravel().tolist())
    plt.savefig(output_path + name_out + "_" + str(ipca_init_list[0]) + "_" + str(pca_ipca_end) + "_" + str(interval) + ".png", bbox_inches='tight', orientation="landscape", dpi=1000)


def plot_snr(ipca_init_list,
             pca_ipca_end,
             interval,
             input_path,
             output_path,
             array_file_name,
             name_out="out",
             title="SNR vs. Rank"):
    
    snr_arrays = np.loadtxt(input_path + array_file_name)

    pca_lst_inner = []
    list_fraction_inner = int(pca_ipca_end/interval)+1
    for i in range(list_fraction_inner):
        pca_lst_inner.append(snr_arrays[i][4])
    
    x_axis = np.concatenate((np.array([1]), np.linspace(5, 80, len(pca_lst_inner)-1)), axis=0)
    
    plt.rc("axes", prop_cycle=(cycler("color", ["k", "r", "y", "g", "c", "b", "m"])))
    plt.plot(x_axis, pca_lst_inner, marker=".", linestyle="--", label="PCA")
    
    ipca_dic = {}
    extra = list_fraction_inner
    
    for j in ipca_init_list:
        ipca_dic["{0}".format(j)]= []
        list_fraction_inner = int((pca_ipca_end-j-1)//interval)+1
        for i in range(list_fraction_inner):
            ipca_dic["{0}".format(j)].append(snr_arrays[i+extra][4])
        extra += list_fraction_inner
        if len(ipca_dic["{0}".format(j)]) != len(pca_lst_inner):
            ipca_dic["{0}".format(j)] = [None] * (len(pca_lst_inner)-len(ipca_dic["{0}".format(j)])) + ipca_dic["{0}".format(j)]
        #print(ipca_dic["ipca_{0}_lst_inner".format(j)])
        plt.plot(x_axis, ipca_dic["{0}".format(j)], marker=".", linestyle="--", label= "IPCA " + str(j))
    
    plt.title(title)
    plt.xlabel("Rank")
    plt.ylabel("SNR")
    plt.legend()
    plt.grid()
    #plt.show()
    plt.savefig(output_path + name_out + ".png", dpi=1000)