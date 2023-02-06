###############################################
##Dmitry Sutormin, 2022##
##Topo IV Topo-Seq analysis##

#Takes WIG files with information about the cumulative distribution of some 
#protein over transcription units (TUs). Plots this information.
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.cluster.hierarchy as sch
from scipy import stats
from Bio import SeqIO
from Bio.SeqUtils import GC as GC_count
import matplotlib.patheffects as PathEffects
from matplotlib_venn import venn2, venn3, venn3_circles
from matplotlib import cm as cm
import collections
from collections import OrderedDict
import pandas as pd
from pandas import DataFrame


#Path to the working directory.
PWD='Data_analysis\Metagene_plots_repr_transcripts\\'

#Half-window width will be used to smooth signal.
Sm_window=100

#Name of the signal to plotted (protein or smth.).
Signal_name_TopoIV_wt='TopoIV_wt'
#Dictionary of pathes to input WIG data.
Wig_data_in_dict_transcripts_TopoIV_wt={'TopoIV wt All TUs'  :  PWD + 'Signal_of_TUs_wig\All_TUs_1672\Signal_TopoIV_Cfx_over_All_TUs_1672_width_15000bp_gb_5000bp.wig',
                                        'TopoIV wt HETU'     :  PWD + 'Signal_of_TUs_wig\HETU_no_rRNA_200\Signal_TopoIV_Cfx_over_HETU_no_rRNA_200_width_15000bp_gb_5000bp.wig',
                                        'TopoIV wt LETU'     :  PWD + 'Signal_of_TUs_wig\LETU_no_rRNA_200\Signal_TopoIV_Cfx_over_LETU_no_rRNA_200_width_15000bp_gb_5000bp.wig',  
                                        }

#Name of the signal to plotted (protein or smth.).
Signal_name_TopoIV_mut='TopoIV_mut'
Wig_data_in_dict_transcripts_TopoIV_mut={'TopoIV mut All TUs' :  PWD + 'Signal_of_TUs_wig\All_TUs_1672\Signal_TopoIV_GyrA_S83L_Cfx_over_All_TUs_1672_width_15000bp_gb_5000bp.wig',
                                         'TopoIV mut HETU'    :  PWD + 'Signal_of_TUs_wig\HETU_no_rRNA_200\Signal_TopoIV_GyrA_S83L_Cfx_over_HETU_no_rRNA_200_width_15000bp_gb_5000bp.wig',
                                         'TopoIV mut LETU'    :  PWD + 'Signal_of_TUs_wig\LETU_no_rRNA_200\Signal_TopoIV_GyrA_S83L_Cfx_over_LETU_no_rRNA_200_width_15000bp_gb_5000bp.wig',  
                                         }

#Name of the signal to plotted (protein or smth.).
Signal_name_MukB='MukB'
Wig_data_in_dict_transcripts_MukB={'MukB All TUs'       :  PWD + 'Signal_of_TUs_wig\All_TUs_1672\Signal_MukB_over_All_TUs_1672_width_15000bp_gb_5000bp.wig',
                                   'MukB HETU'          :  PWD + 'Signal_of_TUs_wig\HETU_no_rRNA_200\Signal_MukB_over_HETU_no_rRNA_200_width_15000bp_gb_5000bp.wig',
                                   'MukB LETU'          :  PWD + 'Signal_of_TUs_wig\LETU_no_rRNA_200\Signal_MukB_over_LETU_no_rRNA_200_width_15000bp_gb_5000bp.wig', 
                                   }

#Name of the signal to plotted (protein or smth.).
Signal_name_gyrase='DNA_gyrase'
Wig_data_in_dict_transcripts_Gyrase={'Gyrase All TUs'     :  PWD + 'Signal_of_TUs_wig\All_TUs_1672\Signal_Gyrase_over_All_TUs_1672_width_15000bp_gb_5000bp.wig',
                                     'Gyrase HETU'        :  PWD + 'Signal_of_TUs_wig\HETU_no_rRNA_200\Signal_Gyrase_over_HETU_no_rRNA_200_width_15000bp_gb_5000bp.wig',
                                     'Gyrase LETU'        :  PWD + 'Signal_of_TUs_wig\LETU_no_rRNA_200\Signal_Gyrase_over_LETU_no_rRNA_200_width_15000bp_gb_5000bp.wig',                                                                           
                                      }




#Set type to choose plotting parameters: genes, operons or transcripts.
Set_type="transcripts"

#######
#Checks if directory exists and if not creates.
#######

def Dir_check_create(some_path):
    if not os.path.exists(some_path):
        os.makedirs(some_path)    
    return

#Output path.
Out_path=f'{PWD}\Figures\Plot_combinations\\'
Dir_check_create(Out_path)


#######
#Parses WIG file with FE over TUs.
#######

def wig_FE_over_genes_parsing(name, wigfile):
    print('Now is processing: ' + str(name) + ' ' + str(wigfile))
    wigin=open(wigfile, 'r')
    NE_values=[]
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] in ['track']:
            ww_l=line[2].split('=')[1].rstrip('"').lstrip('"').split('_')
            win_width=int(ww_l[0])
            length=int(ww_l[1])
        if line[0] not in ['track', 'fixedStep']:
            NE_values.append(float(line[0]))
    wigin.close()
    print(f'Window width: {win_width}, length of TU: {length}')
    return NE_values, win_width, length

#######
#Returns smoothed tracks.
#######

def Smoothing(ends, window):
    smoothed=[]
    #Calculating the value for the first position
    sm=0.0
    window_float=float(window)
    sm+=np.mean(ends[:2*window])
    smoothed.append(sm)
    #Calculating values for the part of the array remains
    for i in range(len(ends)-2*window):
        sm+=(ends[i+2*window]-ends[i])/(2*window_float)
        smoothed.append(sm)
    return smoothed



#######
#Plot the signal for different groups of genes together.
#######


def plot_FE_TUs_groups(wig_in_dict, sm_window, output_path, set_name, set_type):
    #Number of genes in sets.
    TU_sets_v={'TopoIV wt All TUs' : 1672, 'TopoIV wt HETU' : 200, 'TopoIV wt LETU' : 200,
               'TopoIV mut All TUs' : 1672, 'TopoIV mut HETU' : 200, 'TopoIV mut LETU' : 200,
               'MukB All TUs' : 1672, 'MukB HETU' : 200, 'MukB LETU' : 200,
               'Gyrase All TUs' : 1672, 'Gyrase HETU' : 200, 'Gyrase LETU' : 200,
               }
    
    #FE averaged WIG parsing
    dict_of_wigs={}
    win_width=15000
    length=5000
    for name, file in wig_in_dict.items():
        data=wig_FE_over_genes_parsing(name, file)
        win_width=data[1]
        length=data[2]
        dict_of_wigs[name]=data[0]        
    positions=np.arange(-win_width, win_width+length, 1)    
    print(positions[0], positions[-1])
       
    
    #Plot FE over genes.
    plt.figure(figsize=(10, 6), dpi=100) #Def size - 10, 6; Closer look - 3, 6
    plot1=plt.subplot(111)
    ##Transcription units below.
    if set_name=='TopoIV_wt':
        #TopoIV wt 
        plot1.plot(positions, np.array(dict_of_wigs['TopoIV wt LETU']),    linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'LETU ({TU_sets_v["TopoIV wt LETU"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['TopoIV wt All TUs']), linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=f'All TUs ({TU_sets_v["TopoIV wt All TUs"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17       
        plot1.plot(positions, np.array(dict_of_wigs['TopoIV wt HETU']),    linestyle='-', color='#b08642', linewidth=2.5, alpha=1, label=f'HETU ({TU_sets_v["TopoIV wt HETU"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25         
    elif set_name=='TopoIV_mut':
        #TopoIV mut 
        plot1.plot(positions, np.array(dict_of_wigs['TopoIV mut LETU']),    linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'LETU ({TU_sets_v["TopoIV mut LETU"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['TopoIV mut All TUs']), linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=f'All TUs ({TU_sets_v["TopoIV mut All TUs"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17       
        plot1.plot(positions, np.array(dict_of_wigs['TopoIV mut HETU']),    linestyle='-', color='#b08642', linewidth=2.5, alpha=1, label=f'HETU ({TU_sets_v["TopoIV mut HETU"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25                 
    elif set_name=='MukB':
        #MukB
        plot1.plot(positions, np.array(dict_of_wigs['MukB LETU']),    linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'LETU ({TU_sets_v["MukB LETU"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['MukB All TUs']), linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=f'All TUs ({TU_sets_v["MukB All TUs"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17       
        plot1.plot(positions, np.array(dict_of_wigs['MukB HETU']),    linestyle='-', color='#b08642', linewidth=2.5, alpha=1, label=f'HETU ({TU_sets_v["MukB HETU"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25         
    elif set_name=='DNA_gyrase':
        #Gyrase
        plot1.plot(positions, np.array(dict_of_wigs['Gyrase LETU']),    linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'LETU ({TU_sets_v["Gyrase LETU"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Gyrase All TUs']), linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=f'All TUs ({TU_sets_v["Gyrase All TUs"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17       
        plot1.plot(positions, np.array(dict_of_wigs['Gyrase HETU']),    linestyle='-', color='#b08642', linewidth=2.5, alpha=1, label=f'HETU ({TU_sets_v["Gyrase HETU"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25         
         
                    
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()   
    plot1.set_xticks(ticks, minor='False')
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    ticks_lables[ticks.index(length)]='TE'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1) 
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axvline(length, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.set_yticks([1], minor='True')
    plot1.grid(axis='y', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axhline(1, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.legend(fontsize=15, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'{set_name} fold enrichment', size=20)    
    plt.savefig(f'{output_path}\\{set_name}_FE_over_{set_type}_{win_width}bp_nd_with_body_{length}bp.png', dpi=400, figsize=(10, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6
    plt.show()
    plt.close()     
    
        
    dict_of_wigs_sm={}
    for name, wig in dict_of_wigs.items():
        dict_of_wigs_sm[name]=Smoothing(wig, sm_window)
    positions_sm=np.arange(-win_width+sm_window, win_width+length-(sm_window-1), 1)  
    print(positions_sm[0], positions_sm[-1])
        
    #Plot smoothed FE over genes.
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)
    ##Transcription units below.
    if set_name=='TopoIV_wt':
        #TopoIV wt 
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['TopoIV wt LETU']),    linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'LETU ({TU_sets_v["TopoIV wt LETU"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['TopoIV wt All TUs']), linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=f'All TUs ({TU_sets_v["TopoIV wt All TUs"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17       
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['TopoIV wt HETU']),    linestyle='-', color='#b08642', linewidth=2.5, alpha=1, label=f'HETU ({TU_sets_v["TopoIV wt HETU"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25         
    elif set_name=='TopoIV_mut':
        #TopoIV mut 
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['TopoIV mut LETU']),    linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'LETU ({TU_sets_v["TopoIV mut LETU"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['TopoIV mut All TUs']), linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=f'All TUs ({TU_sets_v["TopoIV mut All TUs"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17       
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['TopoIV mut HETU']),    linestyle='-', color='#b08642', linewidth=2.5, alpha=1, label=f'HETU ({TU_sets_v["TopoIV mut HETU"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25                 
    elif set_name=='MukB':
        #MukB
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['MukB LETU']),    linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'LETU ({TU_sets_v["MukB LETU"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['MukB All TUs']), linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=f'All TUs ({TU_sets_v["MukB All TUs"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17       
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['MukB HETU']),    linestyle='-', color='#b08642', linewidth=2.5, alpha=1, label=f'HETU ({TU_sets_v["MukB HETU"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25         
    elif set_name=='DNA_gyrase':
        #Gyrase
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Gyrase LETU']),    linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'LETU ({TU_sets_v["Gyrase LETU"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Gyrase All TUs']), linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=f'All TUs ({TU_sets_v["Gyrase All TUs"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17       
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Gyrase HETU']),    linestyle='-', color='#b08642', linewidth=2.5, alpha=1, label=f'HETU ({TU_sets_v["Gyrase HETU"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25         
      
        
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    ticks_lables[ticks.index(length)]='TE'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1)
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)   
    plot1.axvline(0, color='black', linestyle=':', alpha=0.7, linewidth=1.5)
    plot1.axvline(length, color='black', linestyle=':', alpha=0.7, linewidth=1.5)    
    plot1.set_yticks([1], minor='True')
    plot1.grid(axis='y', which='minor', color='grey', linestyle='--', alpha=0.7, linewidth=1)   
    plot1.axhline(1, color='black', linestyle=':', alpha=0.7, linewidth=1.5)
    plot1.legend(fontsize=14.5, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'{set_name} fold enrichment', size=20)     
    plt.savefig(f'{output_path}\\{set_name}_FE_over_{set_type}_smoothed_{win_width}bp_nd_with_body_{length}bp_smoothed_{2*sm_window}bp.png', dpi=400, figsize=(10, 6), transparent=True)       
    plt.savefig(f'{output_path}\\{set_name}_FE_over_{set_type}_smoothed_{win_width}bp_nd_with_body_{length}bp_smoothed_{2*sm_window}bp.svg', dpi=400, figsize=(10, 6), transparent=True)   
    plt.show()
    plt.close()    
    return

#TopoIV wt
plot_FE_TUs_groups(Wig_data_in_dict_transcripts_TopoIV_wt, Sm_window, Out_path, Signal_name_TopoIV_wt, Set_type)
#TopoIV mut
plot_FE_TUs_groups(Wig_data_in_dict_transcripts_TopoIV_mut, Sm_window, Out_path, Signal_name_TopoIV_mut, Set_type)
#MukB
plot_FE_TUs_groups(Wig_data_in_dict_transcripts_MukB, Sm_window, Out_path, Signal_name_MukB, Set_type)
#DNA-gyrase
plot_FE_TUs_groups(Wig_data_in_dict_transcripts_Gyrase, Sm_window, Out_path, Signal_name_gyrase, Set_type)
