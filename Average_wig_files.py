###############################################
##Dmitry Sutormin, 2022##
##TopoIV Topo-Seq analysis##

####
#Script computes by-position average of a set of wig files.
#Script prepares a correlation matrix between the tracks.
####

###############################################

#######
#Packages to be imported.
#######

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm as cm
import scipy
import scipy.cluster.hierarchy as sch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


#Path to the working directory.
PWD="Data_analysis\WIG_MatP_Nolivos_no_lim\\"

#Dictionary of replicas 
#'Replica name' : 'Path to wig file'
Dict_of_replicas_1={'TopoIV_Cfx_1' : PWD + "Replicate_1_+Cfx+IP_div_+Cfx-IP_FE.wig",
                    'TopoIV_Cfx_2' : PWD + "Replicate_2_+Cfx+IP_div_+Cfx-IP_FE.wig",
                    'TopoIV_Cfx_3' : PWD + "Replicate_3_+Cfx+IP_div_+Cfx-IP_FE.wig",
                    }

Dict_of_replicas_2={'TopoIV_GyrA_S83L_Cfx_1' : PWD + "Replicate_1_GyrA_S83L_+Cfx+IP_div_+Cfx-IP_FE.wig",
                    'TopoIV_GyrA_S83L_Cfx_2' : PWD + "Replicate_2_GyrA_S83L_+Cfx+IP_div_+Cfx-IP_FE.wig",
                    'TopoIV_GyrA_S83L_Cfx_3' : PWD + "Replicate_3_GyrA_S83L_+Cfx+IP_div_+Cfx-IP_FE.wig",
                    }

Dict_of_replicas_3={'MatP_1' : PWD + "MatP_1_SRR1946729.wig",
                    'MatP_2' : PWD + "MatP_2_SRR1946730.wig",
                    'MatP_3' : PWD + "MatP_3_SRR1946731.wig",
                    }

#ID or short description of the track (will be the name of a track in IGV).
#Output_name='TopoIV_Cfx_av'
Output_name='MatP_Nolivos_av'
#ID of chromosome (for w3110_Mu_SGS: NC_007779.1_w3110_Mu)
Chromosome_name='NC_007779.1_w3110_Mu'
#Output path for the corralation matrix.
Outpath=PWD


#######
#Parses WIG file.
#######

def wig_parsing(wigfile):
    print('Now is processing: ' + str(wigfile))
    wigin=open(wigfile, 'r')
    NE_values=[]
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            NE_values.append(float(line[0]))
    wigin.close()
    return NE_values


#########
#Compute correlation matrix and draw heatmaps.
#########

#Plot diagonal correlation matrix.
def correlation_matrix(df, cor_method, title, outpath_folder, file_name):
    #Create heatmap.
    fig=plt.figure(figsize=(5,5), dpi=100)
    ax1=fig.add_subplot(111)
    cmap=cm.get_cmap('rainbow', 30)
    
    #Create correlation matrix and draw heatmap.
    df_cor_matrix=df.corr(method=cor_method)   
    df_cor_matrix.to_csv(outpath_folder+file_name+'.csv', sep='\t', header=True, index=True)
    color_ax=ax1.imshow(df_cor_matrix, interpolation="nearest", cmap=cmap, norm=None, vmin=-1, vmax=1)
    ax1.grid(True, which='minor', linestyle="--", linewidth=0.5, color="black")
    plt.title(title)
    labels=list(df)
    ax1.set_xticks(np.arange(len(labels)))
    ax1.set_yticks(np.arange(len(labels)))    
    ax1.set_xticklabels(labels, fontsize=12, rotation=90)
    ax1.set_yticklabels(labels, fontsize=12)
    ax1.set_ylim(sorted(ax1.get_xlim(), reverse=True)) #Solves a bug in matplotlib 3.1.1 discussed here: https://stackoverflow.com/questions/56942670/matplotlib-seaborn-first-and-last-row-cut-in-half-of-heatmap-plot
    
    #Create text annotation for heatmap pixels.
    for i in range(len(labels)):
        for j in range(len(labels)):
            text = ax1.text(i, j, round(df_cor_matrix[labels[i]][labels[j]], 2), ha="center", va="center", color="black")    
    
    #Add colorbar, make sure to specify tick locations to match desired ticklabels.
    #Full scale:[-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55, -0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])
    #fig.colorbar(cax, ticks=[-1.00, -0.90, -0.80, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00], shrink=0.7)
    #axins1=inset_axes(ax1, width="5%",  height="50%",  loc='upper right', bbox_to_anchor=(1.05, 0., 1, 1), bbox_transform=ax1.transAxes, borderpad=0)    #From here: https://matplotlib.org/3.1.1/gallery/axes_grid1/demo_colorbar_with_inset_locator.html 
    fig.colorbar(color_ax, ticks=[-1.00, -0.80, -0.60, -0.40, -0.20, 0.00, 0.20, 0.40, 0.60, 0.80, 1.00], shrink=0.5, panchor=(1.0, 0.5))
    plt.tight_layout()
    plt.savefig(outpath_folder+file_name+'.png', dpi=400, figsize=(5, 5))
    plt.show()
    plt.close()
    return df_cor_matrix


#########
#Plot correlation matrix.
#########

def correlation_matrix_plot(df_cor_matrix, cor_method, title, outpath_folder, file_name):
    fig=plt.figure(figsize=(5,5), dpi=100)
    ax1=fig.add_subplot(111)
    cmap=cm.get_cmap('rainbow', 30)
    #Create correlation matrix and heatmap.
    color_ax=ax1.imshow(df_cor_matrix, interpolation="nearest", cmap=cmap, norm=None, vmin=-1, vmax=1, aspect="equal")
    ax1.grid(True, which='minor', linestyle="--", linewidth=0.5, color="black")
    plt.title(title)
    #Label ticks.
    labels=list(df_cor_matrix)
    ax1.set_xticks(np.arange(len(labels)))
    ax1.set_yticks(np.arange(len(labels)))    
    ax1.set_xticklabels(labels, fontsize=12, rotation=90)
    ax1.set_yticklabels(labels, fontsize=12)
    ax1.set_ylim(sorted(ax1.get_xlim(), reverse=True)) #Solves a bug in matplotlib 3.1.1 discussed here: https://stackoverflow.com/questions/56942670/matplotlib-seaborn-first-and-last-row-cut-in-half-of-heatmap-plot
    #Create text annotation for heatmap pixels.
    for i in range(len(labels)):
        for j in range(len(labels)):
            text = ax1.text(i, j, round(df_cor_matrix[labels[i]][labels[j]], 2), ha="center", va="center", color="black")    
    #Add colorbar.
    #Full scale:[-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55, -0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])
    #axins1=inset_axes(ax1, width="5%",  height="50%",  loc='upper right', bbox_to_anchor=(1.05, 0., 1, 1), bbox_transform=ax1.transAxes, borderpad=0)    #From here: https://matplotlib.org/3.1.1/gallery/axes_grid1/demo_colorbar_with_inset_locator.html 
    fig.colorbar(color_ax, ticks=[-1.00, -0.80, -0.60, -0.40, -0.20, 0.00, 0.20, 0.40, 0.60, 0.80, 1.00], shrink=0.5, panchor=(1.0, 0.5))
    plt.tight_layout()
    plt.savefig(outpath_folder+file_name+'.png', dpi=400, figsize=(5, 5))
    plt.show()
    plt.close()
    return


#######
#Identify clusters in a corralation matrix (hierarchy clustering).
#Code stolen from https://github.com/TheLoneNut/CorrelationMatrixClustering/blob/master/CorrelationMatrixClustering.ipynb
#######

def Clustering(clust_matrix, outpath_folder, file_name):
    X = clust_matrix.values
    d = sch.distance.pdist(X)   # vector of pairwise distances
    L = sch.linkage(d, method='complete')
    ind = sch.fcluster(L, 0.5*d.max(), 'distance')
    columns = [clust_matrix.columns.tolist()[i] for i in list((np.argsort(ind)))]
    clust_matrix = clust_matrix.reindex(columns, axis=1)
    clust_matrix.to_csv(outpath_folder+file_name+'.csv', sep='\t', header=True, index=True)
    return clust_matrix


#######
#Average wig files and write WIG.
#######

def average_write_wig(dict_of_replicas, average_file_path, name, Chromosome_name):

    #Write file with avaraged data.
    average_out=open(average_file_path+name+'.wig', 'w')
    average_out.write('track type=wiggle_0 name="'+name+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name+' start=1 step=1\n')
    
    for i in range(len(dict_of_replicas[list(dict_of_replicas.keys())[0]])):
        av_data_position=[]
        for replica_name, replica_data in dict_of_replicas.items():
            av_data_position.append(replica_data[i])
        average_out.write(str(np.mean(av_data_position))+'\n')
    average_out.close()
    
    return


#######
#Wrapper function.
#######

def clustering_wrapper(Dict_of_replicas, corr_type, Outpath, Output_name, Chromosome_name):

    #Read input wig files.
    dict_of_replicas={}
    for replica_name, replica_path in Dict_of_replicas.items():
        dict_of_replicas[replica_name]=wig_parsing(replica_path)
    
    #Calculate correlation matrix, plot, clusterize, plot again.
    Correlation_matrix=correlation_matrix(pd.DataFrame(dict_of_replicas), corr_type, 'Correlation of samples', Outpath, Output_name+'_correlation_matrix')
    Correlation_matrix_clusterized=Clustering(Correlation_matrix, Outpath, Output_name+'_correlation_matrix_clusterized')
    correlation_matrix_plot(Correlation_matrix_clusterized, corr_type, 'Correlation of samples clusterized', Outpath, Output_name+'_correlation_matrix_clusterized') 
    
    #Write position-averaged WIG-file.
    average_write_wig(dict_of_replicas, Outpath, Output_name, Chromosome_name)
    return

clustering_wrapper(Dict_of_replicas_3, 'pearson', Outpath, Output_name, Chromosome_name)

