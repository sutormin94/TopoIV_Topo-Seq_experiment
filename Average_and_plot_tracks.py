###############################################
##Dmitry Sutormin, 2022##
##TopoIV Topo-Seq analysis##

#The script takes tetrades of WIG files contain N3E or N5E values: A+IP+, A+IP-, A-IP+, A-IP-.
#It smooths A+IP- and A-IP- tracks and divides A+IP+ and A-IP+ by them.
#Once obtains A+IP+_div and A-IP+_div the script performs Audic-Clavery
#statistic test and returns regions of A+IP+_div where i and i+5 positions are
#significantly higher than corresponding in A-IP+. These regions till now are called GCSs.
#GCSs are stored in the output TXT file. 
#Also two plots are generated: 1) signal coverage over the genome for treated and untreated samples;
#2) Motif expected to be under the GCSs.

#Requirements: TAB file with deletions.
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm as cm
import matplotlib.patheffects as PathEffects
import numpy as np
import scipy
from scipy import stats
import scipy.cluster.hierarchy as sch
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import GC as GC_count
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import locale

locale.getdefaultlocale()

#######
#Variables to be defined.
#######

#Path to the working directory
pwd="Data_analysis\\"
#Path to the file with regions to be omitted (e.g. deletions).
Deletions="TopoIV_Topo-Seq_experiment\Additional_genome_features\Deletions_w3110_G_Mu_SGS.broadPeak"
#Correlation method.
Correlation_type='pearson'

#######
#Opens and reads BED file with deletions coordinates.
#Example:
#GenomeID\tStart\tEnd
#NC_007779.1_w3110_Mu\t274500\t372148
#######

def deletions_info(del_path):
    del_ar=[]
    filein=open(del_path, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        del_ar.append([int(line[1]), int(line[2])])
    filein.close()
    return del_ar

#######
#Returns nearby NE value if current position falls into deleted region of the genome.
#######

def get_value(i, ends, deletions):
    if i<0: #coordinate is out of the left genome border (start)
        j=len(ends)+i
    elif i>=len(ends): #coordinate is out of the right genome border (end)
        j=i-len(ends)
    else: #coordinate is within the genome borders
        check_in_del=0
        for dl in deletions: #check if coordinate falls into deletion
            if dl[1]>=i>=dl[0]:
                j=dl[1]-dl[0]+i+1
                check_in_del=1
        if check_in_del==0:
            j=i
    return ends[j]

#######
#Parses WIG file with N3/5E values.
#Computes a total number of Ends.
#######

def wig_parsing(wigfile):
    print('Now is processing: ' + str(wigfile))
    wigin=open(wigfile, 'r')
    NE_values=[]
    Total_NE=0
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            NE_values.append(int(line[0]))
            Total_NE+=int(line[0])
    print('Total number of ends: ' + str(Total_NE))
    wigin.close()
    return NE_values, Total_NE

#######
#Returns smoothed N3/5E tracks.
#Smoothing using sliding window (default - 100000 nt).
#######

def Smoothing(ends, deletions):
    smoothed=[]
    #Calculating the value for the first genome position
    mean=0.0
    window=10000
    window_float=float(window)
    for i in range(-window, window):
        mean=mean + get_value(i, ends, deletions)
    mean=mean/(2*window_float)
    smoothed.append(mean)
    #Calculating values for the part of the genome remains
    for i in range(1, len(ends)):
        mean=mean + (get_value(i+window, ends, deletions) - get_value(i-window, ends, deletions))/(2*window_float)
        smoothed.append(mean)
    return smoothed

#######
#Returns A+IP+/smoothed(A+IP-) and A-IP+/smoothed(A-IP-) tracks ready for GCSs calling.
#######

def norm_smooth_divide(ex_file_path, cont_file_path, un_ex_file_path, un_cont_file_path, deletions):
    #WIG parsing, total NE counting (for further normalization on reads number)
    treated_experiment=wig_parsing(ex_file_path) #+A+IP
    treated_control=wig_parsing(cont_file_path) #+A-IP
    untreated_experiment=wig_parsing(un_ex_file_path) #-A+IP
    untreated_control=wig_parsing(un_cont_file_path) #-A-IP
    #Normalization on the coverage of neutral region - NR (does not influenced by replication). 
    #NR is set to be 0:2500000 bp.
    te=np.mean(treated_experiment[0][0:2500000])
    tc=np.mean(treated_control[0][0:2500000]) 
    ue=np.mean(untreated_experiment[0][0:2500000])
    uc=np.mean(untreated_control[0][0:2500000])
    norm_array=[te, tc, ue, uc]
    print(norm_array)
    Min_total_NE=min(norm_array)
    print('Min_total_NE: ' + str(Min_total_NE))
    treated_experiment_norm=[1.0 * x * Min_total_NE/norm_array[0] for x in treated_experiment[0]] #+A+IP norm
    print('Passed!')
    treated_control_norm=[1.0 * x * Min_total_NE/norm_array[1] for x in treated_control[0]] #+A-IP norm
    print('Passed!')
    untreated_experiment_norm=[1.0 * x * Min_total_NE/norm_array[2] for x in untreated_experiment[0]] #-A+IP norm
    print('Passed!')
    untreated_control_norm=[1.0 * x * Min_total_NE/norm_array[3] for x in untreated_control[0]] #-A-IP norm
    print('Passed!')
    #Control samples smoothing: A+IP- and A-IP-
    un_experiment_norm_sm=Smoothing(untreated_experiment_norm, deletions) #-A+IP norm sm 
    un_control_norm_sm=Smoothing(untreated_control_norm, deletions) #-A-IP norm sm
    #Pairwise division: +A+IP/-A+IP and +A-IP/-A-IP
    ends_divide_IP=[] #+A+IP/-A+IP
    ends_divide_mock=[] #+A-IP/-A-IP
    for i in range (len(treated_experiment_norm)):
        if un_experiment_norm_sm[i]!=0:
            ends_divide_IP.append(treated_experiment_norm[i]/un_experiment_norm_sm[i])
        else:
            ends_divide_IP.append(0)
            
        if un_control_norm_sm[i]!=0:
            ends_divide_mock.append(treated_control_norm[i]/un_control_norm_sm[i])
        else:
            ends_divide_mock.append(0)
            
    return ends_divide_IP, ends_divide_mock

#######
#Plots the enrichment signal over the genome: +A+IP/smoothed(-A+IP) and +A-IP/smoothed(-A-IP)
#######

def plot_IP_signal(fname, IP_nd_ends, deletions, path_out):
    #Some hack to avoid some bug in matplotlib (OverflowError: In draw_path: Exceeded cell block limit)
    #See: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit
    mpl.rcParams['agg.path.chunksize']=10000
    #Regions to be masked (e.g. deletions).  
    mask_array=[]
    for k in range(len(IP_nd_ends)):
        check_in_del=0
        for dl in deletions:
            if dl[1]>=k>=dl[0]:
                mask_array.append(True)
                check_in_del=1
        if check_in_del==0:
            mask_array.append(False)
    IPed=np.ma.masked_array(IP_nd_ends, mask=mask_array)
    #Plotting the distribution of the signal around the genome for IPed and mock samples.
    xcoord=np.arange(0,4647999)
    rRNA_operons=[3424400, 3690771, 2725847, 4212776, 3466047, 223771, 3597167]
    Macrodomains=[0, 400000, 1200000, 2000000, 2600000, 3500000, 4400000]
    MatS_sites=[1139217, 1158978, 1264672, 1325557, 1345163, 1404637, 1454347, 1462801, 1502925, 1508331, 1517770, 1539623, 
                1566098, 1603288, 1622782, 1633140, 1667991, 1706818, 1747967, 1768029, 1841642, 1854015, 1919680]
    
    x_labels=list(np.around(np.array(np.arange(0,4.7,0.2)),2)) + ['dif'] + ['OriC']
    x_labels_coords=list(np.arange(0,4647999,200000)) + [1593653] + [3711950]  
    plt.figure(figsize=(16, 8), dpi=100)
    plt.suptitle(fname, fontsize=20)
    ax=plt.subplot() 
    ax.plot(xcoord, IPed, '-', label='+A+IP/smoothed(-A+IP)', color='blue', linewidth=1)    
    ax.set_xlabel('Genome position, nt', size=17)
    ax.set_ylabel('Cleavage activity enrichment, N3E', size=17)
    ax.set_xticks(x_labels_coords)
    ax.set_xticklabels(x_labels)
    ## Label dif, matS sites. Make ter inset.
    ## Label Ori, rRNA operons. Make ori inset.
    ax.legend(loc='upper right')
    plt.show()
    plt.savefig(path_out+fname+'_IP_cleavage_signal_N3E.png', dpi=300, figsize=(16, 8))
    plt.savefig(path_out+fname+'_IP_cleavage_signal_N3E.svg', dpi=300, figsize=(16, 8))
    plt.close()
    return

#######
#Plots the enrichment signal over the genome: +A+IP/smoothed(-A+IP) and +A-IP/smoothed(-A-IP)
#######

def plot_IP_and_mock_signals(fname, IP_nd_ends, mock_nd_ends, deletions, path_out):
    #Some hack to avoid some bug in matplotlib (OverflowError: In draw_path: Exceeded cell block limit)
    #See: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit
    mpl.rcParams['agg.path.chunksize']=10000
    #Regions to be masked (e.g. deletions).  
    mask_array=[]
    for k in range(len(IP_nd_ends)):
        check_in_del=0
        for dl in deletions:
            if dl[1]>=k>=dl[0]:
                mask_array.append(True)
                check_in_del=1
        if check_in_del==0:
            mask_array.append(False)
    IPed=np.ma.masked_array(IP_nd_ends, mask=mask_array)
    mock=np.ma.masked_array(mock_nd_ends, mask=mask_array)
    #Plotting the distribution of the signal around the genome for IPed and mock samples.
    xcoord=np.arange(0,4647999)
    x_labels=list(np.around(np.array(np.arange(0,4.7,0.2)),2)) + ['dif'] + ['OriC']
    x_labels_coords=list(np.arange(0,4647999,200000)) + [1593653] + [3711950]  
    plt.figure(figsize=(16, 8), dpi=100)
    plt.suptitle(fname, fontsize=20)
    ax=plt.subplot()
    ax.plot(xcoord, IPed, '-', label='+A+IP/smoothed(-A+IP)', color='blue', linewidth=1)
    ax.plot(xcoord, mock, '-', label='+A-IP/smoothed(-A-IP)', color='orange', linewidth=1)   
    ax.set_xlabel('Genome position, nt', size=17)
    ax.set_ylabel('Cleavage activity enrichment, N3E', size=17)
    ax.set_xticks(x_labels_coords)
    ax.set_xticklabels(x_labels)    
    ax.legend(loc='upper right')
    plt.show()
    plt.savefig(path_out+fname+'_IP_and_Mock_cleavage_signal_N3E.png', dpi=300, figsize=(16, 8))
    plt.savefig(path_out+fname+'_IP_and_Mock_cleavage_signal_N3E.svg', dpi=300, figsize=(16, 8))
    plt.close()
    return

#########
#Mask regions (deletions and multiplicated genes).
#########

def mask_dataframe(Replicas_dataframe, regions_to_mask):
    print('Dataframe size before masking: ' + str(Replicas_dataframe.shape))
      
    #Remove items falling into delelted or masked regions.
    maska=(Replicas_dataframe.index>=0)
    len_of_masked=0
    for j in range(len(regions_to_mask)):
        len_of_masked+=(regions_to_mask[j][1]-regions_to_mask[j][0])
        maska=maska & (~((Replicas_dataframe.index<regions_to_mask[j][1]) & (Replicas_dataframe.index>regions_to_mask[j][0])))
    Replicas_dataframe_masked=Replicas_dataframe[maska]  
    print('Length of masked regions: ' + str(len_of_masked))
    print('Dataframe size after masking: ' + str(Replicas_dataframe_masked.shape))
    return Replicas_dataframe_masked

#########
#Compute correlation matrix and draw heatmaps.
#########

def make_correlation_matrix_plot(df, cor_method, title, outpath_folder, file_name):
    fig=plt.figure(figsize=(5,5), dpi=100)
    ax1=fig.add_subplot(111)
    cmap=cm.get_cmap('rainbow', 30)
    #Create correlation matrix and heatmap.
    df_cor_matrix=df.corr(method=cor_method)
    df_cor_matrix.to_csv(outpath_folder+file_name+'.csv', sep='\t', header=True, index=True)
    color_ax=ax1.imshow(df_cor_matrix, interpolation="nearest", cmap=cmap, norm=None, vmin=-1, vmax=1, aspect="equal")
    ax1.grid(True, which='minor', linestyle="--", linewidth=0.5, color="black")
    plt.title(title)
    #Label ticks.
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
    #Add colorbar.
    #Full scale:[-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55, -0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])
    #axins1=inset_axes(ax1, width="5%",  height="50%",  loc='upper right', bbox_to_anchor=(1.05, 0., 1, 1), bbox_transform=ax1.transAxes, borderpad=0)    #From here: https://matplotlib.org/3.1.1/gallery/axes_grid1/demo_colorbar_with_inset_locator.html 
    fig.colorbar(color_ax, ticks=[-1.00, -0.80, -0.60, -0.40, -0.20, 0.00, 0.20, 0.40, 0.60, 0.80, 1.00], shrink=0.5, panchor=(1.0, 0.5))
    plt.tight_layout()
    plt.savefig(outpath_folder+file_name+'.png', dpi=300, figsize=(5, 5))
    plt.savefig(outpath_folder+file_name+'.svg', dpi=300, figsize=(5, 5))
    plt.show()
    plt.close()
    return df_cor_matrix

#######
#Write WIG.
#######

def write_wig(Averaged_track, average_file_path, name, Chromosome_name):

    #Write file with avaraged data.
    average_out=open(average_file_path+name+'.wig', 'w')
    average_out.write('track type=wiggle_0 name="'+name+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name+' start=1 step=1\n')
    
    for i in range(len(Averaged_track)):
        average_out.write(str(Averaged_track[i])+'\n')
        
    average_out.close()
    
    return

#######
#Average wig files.
#######

def average_wig(dict_of_replicas, average_file_path, name, Chromosome_name):
    #Average tracks.
    Averaged_track=[]
    for i in range(len(dict_of_replicas[list(dict_of_replicas.keys())[0]])):
        av_data_position=[]
        for replica_name, replica_data in dict_of_replicas.items():
            av_data_position.append(replica_data[i])
        Averaged_track.append(np.mean(av_data_position))
    #Write averaged track.    
    write_wig(Averaged_track, average_file_path, name, Chromosome_name)
    
    return Averaged_track

#######
#Wraps all the functions: data normalization, GCSs calling, 
#plotting, writing the data.
#######

def Wrapper(data_dictionary, corr_type, experiment_type, chromosome_name, deletions_inpath, path_out):
    
    #Parsing deletions
    deletions=deletions_info(deletions_inpath)
    
    dict_of_replicas_IP={}
    dict_of_replicas_Mock={}
    
    for Tet_ID, tetrade_dictionary in data_dictionary.items():
    
        #Define samples within the tetrade.
        print('Now we are working with: ' + str(Tet_ID))
        treated_control=tetrade_dictionary['A+IP-']
        treated_experiment=tetrade_dictionary['A+IP+']
        untreated_control=tetrade_dictionary['A-IP-']
        untreated_experiment=tetrade_dictionary['A-IP+']
    
        #Obtain pairwise divided tracks: A+IP+/A+IP- and A-IP+/A-IP-.
        ends_fitting=norm_smooth_divide(treated_experiment, treated_control, untreated_experiment, untreated_control, deletions)
        IP_norm_div=ends_fitting[0]
        mock_norm_div=ends_fitting[1]    
        
        #Plotting the distribution of the signal around the genome for treated and untreated samples.
        plot_IP_signal(Tet_ID, IP_norm_div, deletions, path_out)  
        plot_IP_and_mock_signals(Tet_ID, IP_norm_div, mock_norm_div, deletions, path_out) 
        
        dict_of_replicas_IP[Tet_ID]=IP_norm_div
        dict_of_replicas_Mock[Tet_ID]=mock_norm_div
    
    #Prepare dataframe for tracks correlation.
    Replicas_dataframe=pd.DataFrame(dict_of_replicas_IP, columns=list(dict_of_replicas_IP.keys()))
    Replicas_dataframe_masked=mask_dataframe(Replicas_dataframe, deletions)
    
    #Correlate tracks and plot.
    Correlation_matrix=make_correlation_matrix_plot(Replicas_dataframe_masked, corr_type, 'Correlation of samples', path_out, f'TopoIV_{experiment_type}_{corr_type}_tracks_correlation_matrix_masked')
    
    #Average biological replicates.
    IP_norm_div_av=average_wig(dict_of_replicas_IP, path_out, f'TopoIV_{experiment_type}_N3E_FE_av', chromosome_name)
    Mock_norm_div_av=average_wig(dict_of_replicas_Mock, path_out, f'TopoIV_{experiment_type}_N3E_FE_av', chromosome_name)
    
    #Plot averaged tracks. Plotting the distribution of the signal around the genome for treated and untreated samples.
    plot_IP_signal(experiment_type, IP_norm_div_av, deletions, path_out)  
    plot_IP_and_mock_signals(experiment_type, IP_norm_div_av, Mock_norm_div_av, deletions, path_out)   

    return



Data_dictionary_gyrA_S83L={'Replic_1_gyrA_S83L' : {'A+IP+': pwd + "WIG\gyrA-S83L_replic_1_Cfx_IP_N3E.wig", 
                                                   'A+IP-': pwd + "WIG\gyrA-S83L_replic_1_Cfx_no_IP_N3E.wig",
                                                   'A-IP+': pwd + "WIG\gyrA-S83L_replic_1_no_Cfx_IP_N3E.wig", 
                                                   'A-IP-': pwd + "WIG\gyrA-S83L_replic_1_no_Cfx_no_IP_N3E.wig",},
                           'Replic_2_gyrA_S83L' : {'A+IP+': pwd + "WIG\gyrA-S83L_replic_2_Cfx_IP_N3E.wig",
                                                   'A+IP-': pwd + "WIG\gyrA-S83L_replic_2_Cfx_no_IP_N3E.wig",
                                                   'A-IP+': pwd + "WIG\gyrA-S83L_replic_2_no_Cfx_IP_N3E.wig", 
                                                   'A-IP-': pwd + "WIG\gyrA-S83L_replic_2_no_Cfx_no_IP_N3E.wig",},
                           'Replic_3_gyrA_S83L' : {'A+IP+': pwd + "WIG\gyrA-S83L_replic_3_Cfx_IP_N3E.wig",
                                                   'A+IP-': pwd + "WIG\gyrA-S83L_replic_3_Cfx_no_IP_N3E.wig",
                                                   'A-IP+': pwd + "WIG\gyrA-S83L_replic_3_no_Cfx_IP_N3E.wig", 
                                                   'A-IP-': pwd + "WIG\gyrA-S83L_replic_3_no_Cfx_no_IP_N3E.wig",}
                           }

Data_dictionary_wt={'Replic_1_wt' :  {'A+IP+': pwd + "WIG\wt_replic_1_Cfx_IP_N3E.wig",
                                      'A+IP-': pwd + "WIG\wt_replic_1_Cfx_no_IP_N3E.wig",
                                      'A-IP+': pwd + "WIG\wt_replic_1_no_Cfx_IP_N3E.wig", 
                                      'A-IP-': pwd + "WIG\wt_replic_1_no_Cfx_no_IP_N3E.wig",},
                    'Replic_2_wt' :  {'A+IP+': pwd + "WIG\wt_replic_2_Cfx_IP_N3E.wig", 
                                      'A+IP-': pwd + "WIG\wt_replic_2_Cfx_no_IP_N3E.wig",
                                      'A-IP+': pwd + "WIG\wt_replic_2_no_Cfx_IP_N3E.wig",
                                      'A-IP-': pwd + "WIG\wt_replic_2_no_Cfx_no_IP_N3E.wig",},
                    'Replic_3_wt' :  {'A+IP+': pwd + "WIG\wt_replic_3_Cfx_IP_N3E.wig",
                                      'A+IP-': pwd + "WIG\wt_replic_3_Cfx_no_IP_N3E.wig",
                                      'A-IP+': pwd + "WIG\wt_replic_3_no_Cfx_IP_N3E.wig", 
                                      'A-IP-': pwd + "WIG\wt_replic_3_no_Cfx_no_IP_N3E.wig",}
                    }


#Type of the experiment.
experiment_type="Cfx_S83L"
#ID of chromosome (for w3110_Mu_SGS: NC_007779.1_w3110_Mu)
Chromosome_name='NC_007779.1_w3110_Mu'
#Output folder
Path_for_output=f'{pwd}Cleavage_tracks_final\{experiment_type}\\'
if not os.path.exists(Path_for_output):
    os.makedirs(Path_for_output)
    
Wrapper(Data_dictionary_gyrA_S83L, Correlation_type, experiment_type, Chromosome_name, Deletions, Path_for_output)
Wrapper(Data_dictionary_wt, Correlation_type, experiment_type, Chromosome_name, Deletions, Path_for_output)


print('Script ended its work succesfully!')
