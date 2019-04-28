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
             array_prefix,
             name_out="out",
             title="SNR vs. Rank"):
                 
    x_axis = []
    for counter, i in enumerate(range(ipca_init_list[0], pca_ipca_end+1)):       
        if counter == 0 or ((i % interval) == 0):
            x_axis.append(i)
    
    pca_lst = np.loadtxt(input_path + "pca_" + array_prefix + "_.txt")
    
    plt.rc("axes", prop_cycle=(cycler("color", ["k", "r", "tab:orange", "y", "g", "c", "b", "m"])))
    plt.plot(x_axis, pca_lst, marker=".", linestyle=":", label="PCA")
        
    ipca_dic = {}
    
    for i in ipca_init_list:
        ipca_dic["{0}".format(i)]= np.loadtxt(input_path + "ipca_" + array_prefix + "_" + str(i) + ".txt").tolist()
        if len(ipca_dic["{0}".format(i)]) != len(x_axis):
            ipca_dic["{0}".format(i)] = [None] * (len(x_axis)-len(ipca_dic["{0}".format(i)])) + ipca_dic["{0}".format(i)]
        plt.plot(x_axis, ipca_dic["{0}".format(i)], marker=".", linestyle=":", label= "IPCA " + str(i))
        
    plt.title(title)
    plt.xlabel("Rank")
    plt.ylabel("SNR")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), shadow=False, ncol=4)
    plt.grid(axis="y")
    plt.tight_layout()
    #plt.show()
    plt.savefig(output_path + name_out + ".png", dpi=1000)
    
    
def plot_snr_difference(ipca_init_list,
             pca_ipca_end,
             interval,
             input_path,
             output_path,
             array_prefix,
             name_out="out",
             title="IPCA/PCA SNR Difference"):
                 
    x_axis = []
    for counter, i in enumerate(range(ipca_init_list[0], pca_ipca_end+1)):       
        if counter == 0 or ((i % interval) == 0):
            x_axis.append(i)
    
    pca_lst = np.loadtxt(input_path + "pca_" + array_prefix + "_.txt")
    
    plt.rc("axes", prop_cycle=(cycler("color", ["k", "r", "tab:orange", "y", "g", "c", "b", "m"])))
    plt.plot(x_axis, np.asarray(pca_lst)-np.asarray(pca_lst), marker=".", linestyle=":", label="PCA")
        
    ipca_dic = {}
    my = False    
    
    for i in ipca_init_list:
        ipca_dic["{0}".format(i)]= np.loadtxt(input_path + "ipca_" + array_prefix + "_" + str(i) + ".txt").tolist()
        if len(ipca_dic["{0}".format(i)]) != len(x_axis):
            my = True
            diff = len(x_axis)-len(ipca_dic["{0}".format(i)])
            ipca_dic["{0}".format(i)] = np.asarray([0] * diff + ipca_dic["{0}".format(i)])
        ipca_dic["{0}".format(i)] = ipca_dic["{0}".format(i)]-np.asarray(pca_lst)
        if my:
            for j in range(diff):
                ipca_dic["{0}".format(i)][j] = None
        plt.plot(x_axis, np.asarray(ipca_dic["{0}".format(i)]), marker=".", linestyle=":", label= "IPCA " + str(i))
        
    plt.title(title)
    plt.xlabel("Rank")
    plt.ylabel("SNR Difference")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), shadow=False, ncol=4)
    plt.grid(axis="y")
    plt.tight_layout()
    #plt.ylim(-5, 10)
    #plt.show()
    plt.savefig(output_path + name_out + ".png", dpi=1000)