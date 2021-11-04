###############################################
##Dmitry Sutormin, 2019##
##TopoA ChIP-Seq analysis##

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
PWD='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Signal_over_TUs\Representative_transcripts\\'

#Name of the signal to plotted (protein or smth.).
Signal_name='TopoIV_MukB'
#Half-window width will be used to smooth signal.
Sm_window=100
#Dictionary of pathes to input WIG data.

Wig_data_in_dict_transcripts_TopoIV={'TopoIV All TUs' :      PWD + 'Signal_of_TUs_wig\All_TUs_no_dps_1660\Signal_TopoIV_Cfx_over_All_TUs_no_dps_1660_width_15000bp_gb_5000bp.wig',
                                     'TopoIV HETU'    :      PWD + 'Signal_of_TUs_wig\HETU_no_dps_rfa_200\Signal_TopoIV_Cfx_over_HETU_no_dps_rfa_200_width_15000bp_gb_5000bp.wig',
                                     'TopoIV LETU'    :      PWD + 'Signal_of_TUs_wig\LETU_no_dps_200\Signal_TopoIV_Cfx_over_LETU_no_dps_200_width_15000bp_gb_5000bp.wig',                                                                           
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
Out_path=f'{PWD}\Figures\Plot_combinations\{Signal_name}'
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
    TU_sets_v={'All genes' : 4119, 'HEG 270' : 269, 'LEG 270' : 270, 'HEG 370' : 369, 'LEG 370' : 370,
               'All operons' : 2327, 'HEO 144' : 143, 'LEO 144' : 144, 'HEO 186' : 185, 'LEO 186' : 186, 'LAO 27' : 27, 'SAO 27' : 27, 'rRNA_operons' : 7, 
               'All TUs no dps 1660 CTD-Rif-' : 1660, 'All TUs no dps 1660 CTD-Rif+' : 1660, 'All TUs no dps 1660 CTD+Rif-' : 1660, 'All TUs no dps 1660 CTD+Rif+' : 1660,
               'HETU no dps rfa 200 CTD-Rif-' : 200, 'HETU no dps rfa 200 CTD-Rif+' : 200, 'HETU no dps rfa 200 CTD+Rif-' : 200, 'HETU no dps rfa 200 CTD+Rif+' : 200, 
               'LETU no dps 200 CTD-Rif-' : 200, 'LETU no dps 200 CTD-Rif+' : 200, 'LETU no dps 200 CTD+Rif-' : 200, 'LETU no dps 200 CTD+Rif+' : 200, 
               'rRNA 7 CTD-Rif-' : 7, 'rRNA 7 CTD-Rif+' : 7, 'rRNA 7 CTD+Rif-' : 7, 'rRNA 7 CTD+Rif+' : 7, 
               'tRNA 49 CTD-Rif-' : 39, 'tRNA 49 CTD-Rif+' : 39, 'tRNA 49 CTD+Rif-' : 39,'tRNA 49 CTD+Rif+' : 39,
               'All TUs no dps 1660 RpoC'     : 1660, 'HETU no dps rfa 200 RpoC'     : 200, 'LETU no dps 200 RpoC'     : 200, 
               'All TUs no dps 1660 RpoC Rif' : 1660, 'HETU no dps rfa 200 RpoC Rif' : 200, 'LETU no dps 200 RpoC Rif' : 200, 
               'All TUs no dps 1660 RpoB'     : 1660, 'HETU no dps rfa 200 RpoB'     : 200, 'LETU no dps 200 RpoB'     : 200, 
               'All TUs no dps 1660 Rif-'     : 1660, 'HETU no dps rfa 200 Rif-'     : 200, 'LETU no dps 200 Rif-'     : 200,
               'All TUs no dps 1660 Rif+'     : 1660, 'HETU no dps rfa 200 Rif+'     : 200, 'LETU no dps 200 Rif+'     : 200,
               'TopoIV All TUs' : 1660, 'TopoIV HETU' : 1660, 'TopoIV LETU' : 1660,
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
    ##Genes below
    if set_type=='genes':
        #LEG_270
        plot1.plot(positions, dict_of_wigs['LEG 270'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LEG ({TU_sets_v["LEG 270"]})', zorder=6)
        plot1.plot(positions, dict_of_wigs['LEG 270 Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LEG Rif ({TU_sets_v["LEG 270"]})', zorder=5)   
        #LEG_370
        #plot1.plot(positions, dict_of_wigs['LEG_370'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'LEG ({TU_sets_v["LEG_370"]})', zorder=6)
        #plot1.plot(positions, dict_of_wigs['tRNA Rif'], linestyle='--', color='#E692A9', linewidth=1, alpha=0.8, label='tRNA Rif (86)', zorder=5)        
        #All_genes
        plot1.plot(positions, dict_of_wigs['All genes'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All genes ({TU_sets_v["All genes"]})', zorder=10)
        plot1.plot(positions, dict_of_wigs['All genes Rif'], linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All genes Rif ({TU_sets_v["All genes"]})', zorder=9)
        #HEG_270
        plot1.plot(positions, dict_of_wigs['HEG 270'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HEG ({TU_sets_v["HEG 270"]})', zorder=8)
        plot1.plot(positions, dict_of_wigs['HEG 270 Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HEG ({TU_sets_v["HEG 270"]})', zorder=7)
        #HEG_370
        #plot1.plot(positions, dict_of_wigs['HEG_370'], linestyle='-', color='#e4d1b4', linewidth=1, alpha=1.5, label=f'HEG ({TU_sets_v["HEG_370"]})', zorder=4)
        #plot1.plot(positions, dict_of_wigs['ncRNA Rif'], linestyle='--', color='#FFC000', linewidth=0.8, alpha=0.8, label='ncRNA Rif (18)', zorder=3)
    ##Operons below.
    elif set_type=="operons":
        #LEO_144
        #plot1.plot(positions, dict_of_wigs['LEO 144'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LEO ({TU_sets_v["LEO 144"]})', zorder=6)
        #plot1.plot(positions, dict_of_wigs['LEO_144_Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LEO Rif ({TU_sets_v["LEO_144_Rif"]})', zorder=5)   
        #LEO_186 or SAO_27
        #plot1.plot(positions, dict_of_wigs['SAO 27'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'SAO ({TU_sets_v["SAO 27"]})', zorder=6)
        plot1.plot(positions, dict_of_wigs['LEO 186'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1.5, label=f'LEO ({TU_sets_v["LEO 186"]})', zorder=6)
        plot1.plot(positions, dict_of_wigs['LEO 186 Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=0.8, label=f'LEO Rif ({TU_sets_v["LEO 186"]})', zorder=5)        
        #All_operons
        plot1.plot(positions, dict_of_wigs['All operons'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All operons ({TU_sets_v["All operons"]})', zorder=10)
        plot1.plot(positions, dict_of_wigs['All operons Rif'], linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All operons Rif ({TU_sets_v["All operons"]})', zorder=9)
        #HEO_144
        #plot1.plot(positions, dict_of_wigs['HEO 144'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HEO no rRNA, tRNA ({TU_sets_v["HEO 144"]})', zorder=8)
        #plot1.plot(positions, dict_of_wigs['HEO_no_tRNA_rRNA_144_Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HEO no rRNA, tRNA Rif ({TU_sets_v["HEO_no_tRNA_rRNA_144_Rif"]})', zorder=7)
        #HEO_186 or rRNA operons or LAO_27
        #plot1.plot(positions, dict_of_wigs['LAO 27'], linestyle='-', color='#e4d1b4', linewidth=1.5, alpha=1, label=f'LAO ({TU_sets_v["LAO 27"]})', zorder=4)
        #plot1.plot(positions, dict_of_wigs['rRNA_operons'], linestyle='-', color='#e4d1b4', linewidth=1.5, alpha=1, label=f'rRNA operons ({TU_sets_v["rRNA_operons"]})', zorder=4)
        plot1.plot(positions, dict_of_wigs['HEO 186'], linestyle='-', color='#b08642', linewidth=1.5, alpha=1.5, label=f'HEO ({TU_sets_v["HEO 186"]})', zorder=4)
        plot1.plot(positions, dict_of_wigs['HEO 186 Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=0.8, label=f'HEO Rif ({TU_sets_v["HEO 186"]})', zorder=3)  
        ##Transcription units below.
    elif set_type=="transcripts":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        #LETU, no dps
        #plot1.plot(positions,  np.array(dict_of_wigs['LETU no dps 200 CTD-Rif-'])-0.18, linestyle='-', color='#B6B8BD', linewidth=1, alpha=1, label=f'LETU Rif- ({TU_sets_v["LETU no dps 200 CTD-Rif-"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions,  np.array(dict_of_wigs['LETU no dps 200 CTD-Rif+'])-0.1, linestyle='-', color='#B6B8BD', linewidth=1, alpha=1, label=f'LETU CTD-/Rif+ ({TU_sets_v["LETU no dps 200 CTD-Rif+"]})', zorder=5) #Def linewidth=0.8; #R123 -0.1; R12 +0; R3 -0.45    
        #plot1.plot(positions,  np.array(dict_of_wigs['LETU no dps 200 CTD+Rif-'])-0.11, linestyle='-', color='#B6B8BD', linewidth=1, alpha=1, label=f'LETU CTD+ ({TU_sets_v["LETU no dps 200 CTD+Rif-"]})', zorder=6) #Def linewidth=1.5; #R123 -0.1; R123 +0.15; R123 -0.15
        #plot1.plot(positions,  np.array(dict_of_wigs['LETU no dps 200 CTD+Rif+'])-0.2, linestyle='-', color='#B6B8BD', linewidth=1, alpha=1, label=f'LETU CTD/Rif ({TU_sets_v["LETU no dps 200 CTD+Rif+"]})', zorder=5) #Def linewidth=1; #R23 -0.2; R23 +0.25; R23 -0.25         
        #All_TUs no dps
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no dps 1660 CTD-Rif-'])-0.18, linestyle='-', color='#333738', linewidth=1, alpha=1, label=f'All TUs Rif- ({TU_sets_v["All TUs no dps 1660 CTD-Rif-"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no dps 1660 CTD-Rif+'])-0.1, linestyle='--', color='#333738', linewidth=0.8, alpha=0.8, label=f'All TUs CTD-/Rif+ ({TU_sets_v["All TUs no dps 1660 CTD-Rif+"]})', zorder=9) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no dps 1660 CTD+Rif-'])-0.11, linestyle='-', color='#333738', linewidth=1, alpha=1, label=f'All TUs CTD+ ({TU_sets_v["All TUs no dps 1660 CTD+Rif-"]})', zorder=10) #Def linewidth=2; #R123 -0.25; R123 -0.25; R123 -0.27
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no dps 1660 CTD+Rif+'])-0.2, linestyle='-', color='#333738', linewidth=1, alpha=1, label=f'All TUs CTD/Rif ({TU_sets_v["All TUs no dps 1660 CTD+Rif+"]})', zorder=9) #Def linewidth=1; #R23 -0.25; R23 -0.22; R23 -0.24
        #HETU, no dps rfa
        #plot1.plot(positions, np.array(dict_of_wigs['HETU no dps rfa 200 CTD-Rif-'])-0.18, linestyle='-', color='#b08642', linewidth=1, alpha=1, label=f'HETU Rif- ({TU_sets_v["HETU no dps rfa 200 CTD-Rif-"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25
        #plot1.plot(positions, np.array(dict_of_wigs['HETU no dps rfa 200 CTD-Rif+'])-0.1, linestyle='-', color='#b08642', linewidth=1, alpha=1, label=f'HETU CTD-/Rif+ ({TU_sets_v["HETU no dps rfa 200 CTD-Rif+"]})', zorder=7) #Def linewidth=0.8 #R123 -0.05; R12 -0.03; R3 -0.32
        #plot1.plot(positions, np.array(dict_of_wigs['HETU no dps rfa 200 CTD+Rif-'])-0.11, linestyle='-', color='#b08642', linewidth=1, alpha=1, label=f'HETU CTD+ ({TU_sets_v["HETU no dps rfa 200 CTD+Rif-"]})', zorder=8) #Def linewidth=1.5 #R123 -0.20; R123 -0.25; R123 -0.28
        #plot1.plot(positions, np.array(dict_of_wigs['HETU no dps rfa 200 CTD+Rif+'])-0.2, linestyle='-', color='#b08642', linewidth=1, alpha=1, label=f'HETU CTD/Rif ({TU_sets_v["HETU no dps rfa 200 CTD+Rif+"]})', zorder=7) #Def linewidth=1 #R23 -0.15; R23 -0.20; R23 -0.24 
        #rRNA
        #plot1.plot(positions, dict_of_wigs['rRNA 7 CTD-Rif-'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=0.8, label=f'rRNA CTD-Rif- ({TU_sets_v["rRNA 7 CTD-Rif-"]})', zorder=8) #Def linewidth=1.5 #R123 -0; R12 -0; R3 -0
        #plot1.plot(positions, dict_of_wigs['rRNA 7 CTD-Rif+'], linestyle='--', color='#333738', linewidth=1, alpha=1, label=f'rRNA CTD-Rif+ ({TU_sets_v["rRNA 7 CTD-Rif-"]})', zorder=7) #Def linewidth=1 #R123 -0; R12 -0; R3 -0
        #plot1.plot(positions, dict_of_wigs['rRNA 7 CTD+Rif-'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'rRNA CTD+Rif- ({TU_sets_v["rRNA 7 CTD-Rif-"]})', zorder=8) #Def linewidth=1.5 #R123 -0; R123 -0; R123 -0
        #plot1.plot(positions, dict_of_wigs['rRNA 7 CTD+Rif+'], linestyle='--', color='#e4d1b4', linewidth=1, alpha=1, label=f'rRNA CTD+Rif+ ({TU_sets_v["rRNA 7 CTD-Rif-"]})', zorder=7) #Def linewidth=1 #R23 -0; R23 -0; R23 -0 
        #tRNA
        #plot1.plot(positions, np.array(dict_of_wigs['tRNA 49 CTD-Rif-'])-0.1, linestyle='-', color='#757d8b', linewidth=1.5, alpha=0.8, label=f'tRNA CTD-Rif- ({TU_sets_v["tRNA 49 CTD-Rif-"]})', zorder=8) #Def linewidth=1.5 #R123 -0.1; R12 -0.1; R3 -0.1
        #plot1.plot(positions, np.array(dict_of_wigs['tRNA 49 CTD-Rif+'])-0.05, linestyle='--', color='#333738', linewidth=1, alpha=1, label=f'tRNA CTD-Rif+ ({TU_sets_v["tRNA 49 CTD-Rif-"]})', zorder=7) #Def linewidth=1 #R123 -0.1; R12 -0.05; R3 -0.05
        #plot1.plot(positions, np.array(dict_of_wigs['tRNA 49 CTD+Rif-'])-0.25, linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'tRNA CTD+Rif- ({TU_sets_v["tRNA 49 CTD-Rif-"]})', zorder=8) #Def linewidth=1.5 #R123 -0.2; R123 -0.25; R123 -0.25
        #plot1.plot(positions, np.array(dict_of_wigs['tRNA 49 CTD+Rif+'])-0.15, linestyle='--', color='#e4d1b4', linewidth=1, alpha=1, label=f'tRNA CTD+Rif+ ({TU_sets_v["tRNA 49 CTD-Rif-"]})', zorder=7) #Def linewidth=1 #R23 -0.1; R23 -0.15; R23 -0.15   
        #DNA-gyrase.
        #plot1.plot(positions, np.array(dict_of_wigs['LETU no dps 200 Rif-'])-0.1, linestyle='-', color='#B6B8BD', linewidth=1, alpha=1, label=f'LETU ({TU_sets_v["LETU no dps 200 Rif-"]})', zorder=6) #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions, np.array(dict_of_wigs['LETU no dps 200 Rif+'])+0.05, linestyle='--', color='#B6B8BD', linewidth=1, alpha=0.8, label=f'LETU Rif ({TU_sets_v["LETU no dps 200 Rif+"]})', zorder=5) #R123 -0.1; R12 -0; R3 -0.45             
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no dps 1660 Rif-'])-0.1, linestyle='-', color='#333738', linewidth=1, alpha=1, label=f'All TUs ({TU_sets_v["All TUs no dps 1660 Rif-"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no dps 1660 Rif+'])+0.05, linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All TUs Rif ({TU_sets_v["All TUs no dps 1660 Rif+"]})', zorder=9) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42        
        #plot1.plot(positions, np.array(dict_of_wigs['HETU no dps rfa 200 Rif-'])-0.1, linestyle='-', color='#b08642', linewidth=1, alpha=1, label=f'HETU ({TU_sets_v["HETU no dps rfa 200 Rif-"]})', zorder=8) #R123 -0.05; R12 -0.05; R3 -0.25
        #plot1.plot(positions, np.array(dict_of_wigs['HETU no dps rfa 200 Rif+'])+0.05, linestyle='--', color='#b08642', linewidth=1, alpha=0.8, label=f'HETU Rif ({TU_sets_v["HETU no dps rfa 200 Rif+"]})', zorder=7) #R123 -0.05; R12 -0.03; R3 -0.32        
        #RpoC RNAP
        #plot1.plot(positions, np.array(dict_of_wigs['LETU no dps 200 RpoC'])-0.22, linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'LETU ({TU_sets_v["LETU no dps 200 RpoC"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no dps 1660 RpoC'])-0.22, linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=f'All TUs ({TU_sets_v["All TUs no dps 1660 RpoC"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17       
        #plot1.plot(positions, np.array(dict_of_wigs['HETU no dps rfa 200 RpoC'])-0.22, linestyle='-', color='#b08642', linewidth=2.5, alpha=1, label=f'HETU ({TU_sets_v["HETU no dps rfa 200 RpoC"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25         
        #plot1.plot(positions, np.array(dict_of_wigs['LETU no dps 200 RpoB']), linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'LETU ({TU_sets_v["LETU no dps 200 RpoB"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no dps 1660 RpoB']), linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=f'All TUs ({TU_sets_v["All TUs no dps 1660 RpoB"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17       
        #plot1.plot(positions, np.array(dict_of_wigs['HETU no dps rfa 200 RpoB']), linestyle='-', color='#b08642', linewidth=2.5, alpha=1, label=f'HETU ({TU_sets_v["HETU no dps rfa 200 RpoB"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25         
        #plot1.plot(positions, np.array(dict_of_wigs['LETU no dps 200 RpoC Rif'])+0.28,     linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'LETU ({TU_sets_v["LETU no dps 200 RpoC Rif"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no dps 1660 RpoC Rif'])+0.28, linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=f'All TUs ({TU_sets_v["All TUs no dps 1660 RpoC Rif"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17       
        #plot1.plot(positions, np.array(dict_of_wigs['HETU no dps rfa 200 RpoC Rif'])+0.28, linestyle='-', color='#b08642', linewidth=2.5, alpha=1, label=f'HETU ({TU_sets_v["HETU no dps rfa 200 RpoC Rif"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25         
        #TopoIV
        plot1.plot(positions, np.array(dict_of_wigs['TopoIV LETU']),     linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'LETU ({TU_sets_v["LETU no dps 200 RpoC Rif"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['TopoIV All TUs']), linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=f'All TUs ({TU_sets_v["All TUs no dps 1660 RpoC Rif"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17       
        plot1.plot(positions, np.array(dict_of_wigs['TopoIV HETU']), linestyle='-', color='#b08642', linewidth=2.5, alpha=1, label=f'HETU ({TU_sets_v["HETU no dps rfa 200 RpoC Rif"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25         
        
                    
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
    plot1.set_ylabel(f'TopoIV fold enrichment', size=20)
    #plot1.set_title(f'{set_name} signal over {set_type}', size=20)     
    plt.savefig(f'{output_path}\\{set_name}_FE_over_{set_type}_TopoIV_{win_width}bp_nd_with_body_{length}bp.png', dpi=400, figsize=(10, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6
    plt.show()
    plt.close()     
    
    #Smoothing.
    dict_of_wigs_sm={}
    for name, wig in dict_of_wigs.items():
        dict_of_wigs_sm[name]=Smoothing(wig, sm_window)
    positions_sm=np.arange(-win_width+sm_window, win_width+length-(sm_window-1), 1)  
    print(positions_sm[0], positions_sm[-1])
        
    #Plot smoothed FE over genes.
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)
    ##Genes below
    if set_type=='genes':
        #LEG_270
        plot1.plot(positions_sm, dict_of_wigs_sm['LEG 270'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LEG ({TU_sets_v["LEG 270"]})', zorder=6)
        plot1.plot(positions_sm, dict_of_wigs_sm['LEG 270 Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LEG Rif ({TU_sets_v["LEG 270"]})', zorder=5)   
        #LEG_370
        #plot1.plot(positions_sm, dict_of_wigs_sm['LEG_370'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'LEG ({TU_sets_v["LEG_370"]})', zorder=6)
        #plot1.plot(positions_sm, dict_of_wigs_sm['tRNA Rif'], linestyle='--', color='#E692A9', linewidth=1, alpha=0.8, label='tRNA Rif (86)', zorder=5)        
        #All_genes
        plot1.plot(positions_sm, dict_of_wigs_sm['All genes'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All genes ({TU_sets_v["All genes"]})', zorder=10)
        plot1.plot(positions_sm, dict_of_wigs_sm['All genes Rif'], linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All genes Rif ({TU_sets_v["All genes"]})', zorder=9)
        #HEG_270
        plot1.plot(positions_sm, dict_of_wigs_sm['HEG 270'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HEG ({TU_sets_v["HEG 270"]})', zorder=8)
        plot1.plot(positions_sm, dict_of_wigs_sm['HEG 270 Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HEG Rif ({TU_sets_v["HEG 270"]})', zorder=7)
        #HEG_370
        #plot1.plot(positions_sm, dict_of_wigs_sm['HEG_370'], linestyle='-', color='#e4d1b4', linewidth=1, alpha=1.5, label=f'HEG ({TU_sets_v["HEG_370"]})', zorder=4)
        #plot1.plot(positions_sm, dict_of_wigs_sm['ncRNA Rif'], linestyle='--', color='#FFC000', linewidth=0.8, alpha=0.8, label='ncRNA Rif (18)', zorder=3)  
    ##Operons below.
    elif set_type=="operons":
        #LEO_144
        #plot1.plot(positions_sm, dict_of_wigs_sm['LEO 144'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LEO ({TU_sets_v["LEO 144"]})', zorder=6)
        #plot1.plot(positions_sm, dict_of_wigs_sm['LEO_144_Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LEO Rif ({TU_sets_v["LEO_144_Rif"]})', zorder=5)   
        #LEO_186 or SAO 27
        #plot1.plot(positions_sm, dict_of_wigs_sm['SAO 27'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'SAO ({TU_sets_v["SAO 27"]})', zorder=6)
        plot1.plot(positions_sm, dict_of_wigs_sm['LEO 186'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1.5, label=f'LEO ({TU_sets_v["LEO 186"]})', zorder=6)
        plot1.plot(positions_sm, dict_of_wigs_sm['LEO 186 Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=0.8, label=f'LEO Rif ({TU_sets_v["LEO 186"]})', zorder=5)        
        #All_operons
        plot1.plot(positions_sm, dict_of_wigs_sm['All operons'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All operons ({TU_sets_v["All operons"]})', zorder=10)
        plot1.plot(positions_sm, dict_of_wigs_sm['All operons Rif'], linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All operons Rif ({TU_sets_v["All operons"]})', zorder=9)
        #HEO_144
        #plot1.plot(positions_sm, dict_of_wigs_sm['HEO 144'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HEO no rRNA, tRNA ({TU_sets_v["HEO 144"]})', zorder=8)
        #plot1.plot(positions_sm, dict_of_wigs_sm['HEO_no_tRNA_rRNA_144_Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HEO no rRNA, tRNA Rif ({TU_sets_v["HEO_no_tRNA_rRNA_144_Rif"]})', zorder=7)
        #HEO_186 or rRNA operons or LAO_27
        #plot1.plot(positions_sm, dict_of_wigs_sm['LAO 27'], linestyle='-', color='#e4d1b4', linewidth=1.5, alpha=1, label=f'LAO ({TU_sets_v["LAO 27"]})', zorder=4)
        #plot1.plot(positions_sm, dict_of_wigs_sm['rRNA_operons'], linestyle='-', color='#e4d1b4', linewidth=1.5, alpha=1, label=f'rRNA operons ({TU_sets_v["rRNA_operons"]})', zorder=4)
        plot1.plot(positions_sm, dict_of_wigs_sm['HEO 186'], linestyle='-', color='#b08642', linewidth=1.5, alpha=1.5, label=f'HEO ({TU_sets_v["HEO 186"]})', zorder=4)
        plot1.plot(positions_sm, dict_of_wigs_sm['HEO 186 Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=0.8, label=f'HEO Rif ({TU_sets_v["HEO 186"]})', zorder=3)  
    ##Transcription units below.
    elif set_type=="transcripts":
        #LETU, no dps
        #plot1.plot(positions_sm,  np.array(dict_of_wigs_sm['LETU no dps 200 CTD-Rif-'])-0.18, linestyle='-', color='#B6B8BD', linewidth=2, alpha=1, label=f'LETU Rif- ({TU_sets_v["LETU no dps 200 CTD-Rif-"]})', zorder=20)
        #plot1.plot(positions_sm,  np.array(dict_of_wigs_sm['LETU no dps 200 CTD-Rif+'])-0.1, linestyle='-', color='#B6B8BD', linewidth=2, alpha=1, label=f'LETU Rif+ ({TU_sets_v["LETU no dps 200 CTD-Rif+"]})', zorder=5)
        #plot1.plot(positions_sm,  np.array(dict_of_wigs_sm['LETU no dps 200 CTD+Rif-'])-0.11, linestyle='-', color='#B6B8BD', linewidth=2, alpha=1, label=f'LETU CTD+ ({TU_sets_v["LETU no dps 200 CTD+Rif-"]})', zorder=6)
        #plot1.plot(positions_sm,  np.array(dict_of_wigs_sm['LETU no dps 200 CTD+Rif+'])-0.2, linestyle='-', color='#B6B8BD', linewidth=2, alpha=1, label=f'LETU CTD/Rif ({TU_sets_v["LETU no dps 200 CTD+Rif+"]})', zorder=5)          
        #All_TUs no dps
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs no dps 1660 CTD-Rif-'])-0.18, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs Rif- ({TU_sets_v["All TUs no dps 1660 CTD-Rif-"]})', zorder=19)
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs no dps 1660 CTD-Rif+'])-0.1, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs Rif+ ({TU_sets_v["All TUs no dps 1660 CTD-Rif+"]})', zorder=9)
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs no dps 1660 CTD+Rif-'])-0.11, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs CTD+ ({TU_sets_v["All TUs no dps 1660 CTD+Rif-"]})', zorder=10)
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs no dps 1660 CTD+Rif+'])-0.2, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs CTD/Rif ({TU_sets_v["All TUs no dps 1660 CTD+Rif+"]})', zorder=9)        
        #HETU, no dps rfa
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU no dps rfa 200 CTD-Rif-'])-0.18, linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'HETU Rif- ({TU_sets_v["HETU no dps rfa 200 CTD-Rif-"]})', zorder=18)
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU no dps rfa 200 CTD-Rif+'])-0.1, linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'HETU Rif+ ({TU_sets_v["HETU no dps rfa 200 CTD-Rif+"]})', zorder=7)
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU no dps rfa 200 CTD+Rif-'])-0.11, linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'HETU CTD+ ({TU_sets_v["HETU no dps rfa 200 CTD+Rif-"]})', zorder=8)
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU no dps rfa 200 CTD+Rif+'])-0.2, linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'HETU CTD/Rif ({TU_sets_v["HETU no dps rfa 200 CTD+Rif+"]})', zorder=7)        
        #rRNA
        #plot1.plot(positions_sm, dict_of_wigs_sm['rRNA 7 CTD-Rif-'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=0.8, label=f'rRNA CTD-Rif- ({TU_sets_v["rRNA 7 CTD-Rif-"]})', zorder=8) #Def linewidth=1.5
        #plot1.plot(positions_sm, dict_of_wigs_sm['rRNA 7 CTD-Rif+'], linestyle='--', color='#333738', linewidth=1, alpha=1, label=f'rRNA CTD-Rif+ ({TU_sets_v["rRNA 7 CTD-Rif-"]})', zorder=7) #Def linewidth=1
        #plot1.plot(positions_sm, dict_of_wigs_sm['rRNA 7 CTD+Rif-'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'rRNA CTD+Rif- ({TU_sets_v["rRNA 7 CTD-Rif-"]})', zorder=8) #Def linewidth=1.5
        #plot1.plot(positions_sm, dict_of_wigs_sm['rRNA 7 CTD+Rif+'], linestyle='--', color='#e4d1b4', linewidth=1, alpha=1, label=f'rRNA CTD+Rif+ ({TU_sets_v["rRNA 7 CTD-Rif-"]})', zorder=7) #Def linewidth=1        
        #tRNA
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['tRNA 49 CTD-Rif-'])-0.1, linestyle='-', color='#757d8b', linewidth=1.5, alpha=0.8, label=f'tRNA CTD-Rif- ({TU_sets_v["tRNA 49 CTD-Rif-"]})', zorder=8) #Def linewidth=1.5
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['tRNA 49 CTD-Rif+'])-0.05, linestyle='--', color='#333738', linewidth=1, alpha=1, label=f'tRNA CTD-Rif+ ({TU_sets_v["tRNA 49 CTD-Rif-"]})', zorder=7) #Def linewidth=1  
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['tRNA 49 CTD+Rif-'])-0.25, linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'tRNA CTD+Rif- ({TU_sets_v["tRNA 49 CTD-Rif-"]})', zorder=8) #Def linewidth=1.5
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['tRNA 49 CTD+Rif+'])-0.15, linestyle='--', color='#e4d1b4', linewidth=1, alpha=1, label=f'tRNA CTD+Rif+ ({TU_sets_v["tRNA 49 CTD-Rif-"]})', zorder=7) #Def linewidth=1           
        #DNA-gyrase.
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['LETU no dps 200 Rif-'])-0.1, linestyle='-', color='#B6B8BD', linewidth=2, alpha=1, label=f'LETU ({TU_sets_v["LETU no dps 200 Rif-"]})', zorder=6) #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['LETU no dps 200 Rif+'])+0.05, linestyle='--', color='#B6B8BD', linewidth=1.2, alpha=1, label=f'LETU Rif ({TU_sets_v["LETU no dps 200 Rif+"]})', zorder=5) #R123 -0.1; R12 -0; R3 -0.45             
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs no dps 1660 Rif-'])-0.1, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs ({TU_sets_v["All TUs no dps 1660 Rif-"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs no dps 1660 Rif+'])+0.05, linestyle='--', color='#333738', linewidth=1.2, alpha=1, label=f'All TUs Rif ({TU_sets_v["All TUs no dps 1660 Rif+"]})', zorder=9) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42        
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU no dps rfa 200 Rif-'])-0.1, linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'HETU ({TU_sets_v["HETU no dps rfa 200 Rif-"]})', zorder=8) #R123 -0.05; R12 -0.05; R3 -0.25
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU no dps rfa 200 Rif+'])+0.05, linestyle='--', color='#b08642', linewidth=1.2, alpha=1, label=f'HETU Rif ({TU_sets_v["HETU no dps rfa 200 Rif+"]})', zorder=7) #R123 -0.05; R12 -0.03; R3 -0.32        
        #RpoC RNAP
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['LETU no dps 200 RpoC'])-0.22, linestyle='-', color='#B6B8BD', linewidth=2, alpha=1, label=f'LETU ({TU_sets_v["LETU no dps 200 RpoC"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs no dps 1660 RpoC'])-0.22, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs ({TU_sets_v["All TUs no dps 1660 RpoC"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17       
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU no dps rfa 200 RpoC'])-0.22, linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'HETU ({TU_sets_v["HETU no dps rfa 200 RpoC"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25         
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['LETU no dps 200 RpoB']), linestyle='-', color='#B6B8BD', linewidth=2, alpha=1, label=f'LETU ({TU_sets_v["LETU no dps 200 RpoB"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs no dps 1660 RpoB']), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs ({TU_sets_v["All TUs no dps 1660 RpoB"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17       
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU no dps rfa 200 RpoB']), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'HETU ({TU_sets_v["HETU no dps rfa 200 RpoB"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25         
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['LETU no dps 200 RpoC Rif'])+0.28,     linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'LETU ({TU_sets_v["LETU no dps 200 RpoC Rif"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs no dps 1660 RpoC Rif'])+0.28, linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=f'All TUs ({TU_sets_v["All TUs no dps 1660 RpoC Rif"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17       
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU no dps rfa 200 RpoC Rif'])+0.28, linestyle='-', color='#b08642', linewidth=2.5, alpha=1, label=f'HETU ({TU_sets_v["HETU no dps rfa 200 RpoC Rif"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25         
        #TopoIV
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['TopoIV LETU']),     linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'LETU ({TU_sets_v["LETU no dps 200 RpoC Rif"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['TopoIV All TUs']), linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=f'All TUs ({TU_sets_v["All TUs no dps 1660 RpoC Rif"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17       
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['TopoIV HETU']), linestyle='-', color='#b08642', linewidth=2.5, alpha=1, label=f'HETU ({TU_sets_v["HETU no dps rfa 200 RpoC Rif"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25         
        
        
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
    plot1.set_ylabel(f'TopoIV fold enrichment', size=20)
    #plot1.set_title(f'{set_name} signal over {set_type}', size=20)     
    plt.savefig(f'{output_path}\\{set_name}_FE_over_{set_type}_TopoIV_smoothed_{win_width}bp_nd_with_body_{length}bp_smoothed_{2*sm_window}bp.png', dpi=400, figsize=(10, 6), transparent=True)       
    plt.savefig(f'{output_path}\\{set_name}_FE_over_{set_type}_TopoIV_smoothed_{win_width}bp_nd_with_body_{length}bp_smoothed_{2*sm_window}bp.svg', dpi=400, figsize=(10, 6), transparent=True)   
    plt.show()
    plt.close()    
    return


#EcTopoI
#plot_FE_TUs_groups(Wig_data_in_dict_transcripts, Sm_window, Out_path, Signal_name, Set_type)
#RNA-pol
#plot_FE_TUs_groups(Wig_data_in_dict_transcripts_RNApol, Sm_window, Out_path, Signal_name, Set_type)
#DNA-gyrase
#plot_FE_TUs_groups(Wig_data_in_dict_transcripts_gyrase, Sm_window, Out_path, Signal_name, Set_type)
#TopoIV
plot_FE_TUs_groups(Wig_data_in_dict_transcripts_TopoIV, Sm_window, Out_path, Signal_name, Set_type)


#######
#Closer look into TUs, colored area under TSS and TU body.
#######


def plot_FE_TUs_groups_local(wig_in_dict, sm_window, output_path, set_name, set_type):
    #Number of genes in sets.
    TU_sets_v={'All TUs no dps 1660 CTD-Rif-' : 1660, 'All TUs no dps 1660 CTD-Rif+' : 1660, 'All TUs no dps 1660 CTD+Rif-' : 1660, 'All TUs no dps 1660 CTD+Rif+' : 1660,
               'HETU no dps rfa 200 CTD-Rif-' : 200, 'HETU no dps rfa 200 CTD-Rif+' : 200, 'HETU no dps rfa 200 CTD+Rif-' : 200, 'HETU no dps rfa 200 CTD+Rif+' : 200, 
               'LETU no dps 200 CTD-Rif-' : 200, 'LETU no dps 200 CTD-Rif+' : 200, 'LETU no dps 200 CTD+Rif-' : 200, 'LETU no dps 200 CTD+Rif+' : 200}
    
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
    plt.figure(figsize=(6, 6), dpi=100) #Def size - 10, 6; Closer look - 3, 6
    plot1=plt.subplot(111)
    
    ##Transcription units below.
    if set_type=="transcripts":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        #LETU, no dps
        #plot1.plot(positions,  np.array(dict_of_wigs['LETU no dps 200 CTD-Rif-'])-0.02, linestyle='-', color='#B6B8BD', linewidth=1, alpha=1, label=f'LETU Rif- ({TU_sets_v["LETU no dps 200 CTD-Rif-"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions,  np.array(dict_of_wigs['LETU no dps 200 CTD-Rif+'])-0.1, linestyle='-', color='#B6B8BD', linewidth=1, alpha=1, label=f'LETU CTD-/Rif+ ({TU_sets_v["LETU no dps 200 CTD-Rif+"]})', zorder=5) #Def linewidth=0.8; #R123 -0.1; R12 +0; R3 -0.45    
        #plot1.plot(positions,  np.array(dict_of_wigs['LETU no dps 200 CTD+Rif-'])-0.25, linestyle='--', color='#B6B8BD', linewidth=0.8, alpha=0.8, label=f'LETU CTD+Rif- ({TU_sets_v["LETU no dps 200 CTD+Rif-"]})', zorder=6) #Def linewidth=1.5; #R123 -0.1; R123 +0.15; R123 -0.15
        #plot1.plot(positions,  np.array(dict_of_wigs['LETU no dps 200 CTD+Rif+'])-0.25, linestyle='--', color='#B6B8BD', linewidth=0.8, alpha=0.8, label=f'LETU CTD+Rif+ ({TU_sets_v["LETU no dps 200 CTD+Rif+"]})', zorder=5) #Def linewidth=1; #R23 -0.2; R23 +0.25; R23 -0.25         
        #All_TUs no dps
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no dps 1660 CTD-Rif-'])-0.18, linestyle='-', color='#333738', linewidth=1, alpha=1, label=f'All TUs ({TU_sets_v["All TUs no dps 1660 CTD-Rif-"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        #plot1.fill_between(positions, 1, np.array(dict_of_wigs['All TUs no dps 1660 CTD-Rif-'])-0.18, where=((positions>-1000) & (positions<1000)), facecolor='#43c287', alpha=1, interpolate=True)        
        #plot1.fill_between(positions, 1, np.array(dict_of_wigs['All TUs no dps 1660 CTD-Rif-'])-0.18, where=((positions>1000) & (positions<5000)), facecolor='#7ce0ff', alpha=1, interpolate=True)         
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no dps 1660 CTD-Rif+'])-0.1, linestyle='--', color='#333738', linewidth=0.8, alpha=0.8, label=f'All TUs Rif ({TU_sets_v["All TUs no dps 1660 CTD-Rif+"]})', zorder=9) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42
        #plot1.fill_between(positions, 1, np.array(dict_of_wigs['All TUs no dps 1660 CTD-Rif+'])-0.1, where=((positions>-1000) & (positions<1000)), facecolor='#43c287', alpha=1, interpolate=True)        
        #plot1.fill_between(positions, 1, np.array(dict_of_wigs['All TUs no dps 1660 CTD-Rif+'])-0.1, where=((positions>1000) & (positions<5000)), facecolor='#7ce0ff', alpha=1, interpolate=True)          
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no dps 1660 CTD+Rif-'])-0.11, linestyle='-', color='#333738', linewidth=1, alpha=1, label=f'All TUs CTD ({TU_sets_v["All TUs no dps 1660 CTD-Rif+"]})', zorder=10) #Def linewidth=2; #R123 -0.25; R123 -0.25; R123 -0.27
        #plot1.fill_between(positions, 1, np.array(dict_of_wigs['All TUs no dps 1660 CTD+Rif-'])-0.11, where=((positions>-1000) & (positions<1000)), facecolor='#43c287', alpha=1, interpolate=True)        
        #plot1.fill_between(positions, 1, np.array(dict_of_wigs['All TUs no dps 1660 CTD+Rif-'])-0.11, where=((positions>1000) & (positions<5000)), facecolor='#7ce0ff', alpha=1, interpolate=True)           
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no dps 1660 CTD+Rif+'])-0.2, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs CTD/Rif ({TU_sets_v["All TUs no dps 1660 CTD+Rif+"]})', zorder=9) #Def linewidth=1; #R23 -0.25; R23 -0.22; R23 -0.24
        #plot1.fill_between(positions, 1, np.array(dict_of_wigs['All TUs no dps 1660 CTD+Rif+'])-0.2, where=((positions>-1000) & (positions<1000)), facecolor='#43c287', alpha=1, interpolate=True)        
        #plot1.fill_between(positions, 1, np.array(dict_of_wigs['All TUs no dps 1660 CTD+Rif+'])-0.2, where=((positions>1000) & (positions<5000)), facecolor='#7ce0ff', alpha=1, interpolate=True)          
        #HETU, no dps rfa
        plot1.plot(positions, np.array(dict_of_wigs['HETU no dps rfa 200 CTD-Rif-'])-0.18, linestyle='-', color='#b08642', linewidth=1, alpha=1, label=f'HETU ({TU_sets_v["HETU no dps rfa 200 CTD-Rif-"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25
        plot1.fill_between(positions, 1, np.array(dict_of_wigs['HETU no dps rfa 200 CTD-Rif-'])-0.18, where=((positions>-1000) & (positions<1000)), facecolor='#43c287', alpha=1, interpolate=True)        
        plot1.fill_between(positions, 1, np.array(dict_of_wigs['HETU no dps rfa 200 CTD-Rif-'])-0.18, where=((positions>1000) & (positions<5000)), facecolor='#7ce0ff', alpha=1, interpolate=True)          
        #plot1.plot(positions, np.array(dict_of_wigs['HETU no dps rfa 200 CTD-Rif+'])-0.1, linestyle='-', color='#b08642', linewidth=1, alpha=1, label=f'HETU Rif ({TU_sets_v["HETU no dps rfa 200 CTD-Rif+"]})', zorder=7) #Def linewidth=0.8 #R123 -0.05; R12 -0.03; R3 -0.32
        #plot1.fill_between(positions, 1, np.array(dict_of_wigs['HETU no dps rfa 200 CTD-Rif+'])-0.1, where=((positions>-1000) & (positions<1000)), facecolor='#43c287', alpha=1, interpolate=True)        
        #plot1.fill_between(positions, 1, np.array(dict_of_wigs['HETU no dps rfa 200 CTD-Rif+'])-0.1, where=((positions>1000) & (positions<5000)), facecolor='#7ce0ff', alpha=1, interpolate=True)         
        #plot1.plot(positions, np.array(dict_of_wigs['HETU no dps rfa 200 CTD+Rif-'])-0.11, linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'HETU CTD+Rif- ({TU_sets_v["HETU no dps rfa 200 CTD+Rif-"]})', zorder=8) #Def linewidth=1.5 #R123 -0.20; R123 -0.25; R123 -0.28
        #plot1.fill_between(positions, 1, np.array(dict_of_wigs['HETU no dps rfa 200 CTD+Rif-'])-0.11, where=((positions>-1000) & (positions<1000)), facecolor='#43c287', alpha=1, interpolate=True)        
        #plot1.fill_between(positions, 1, np.array(dict_of_wigs['HETU no dps rfa 200 CTD+Rif-'])-0.11, where=((positions>1000) & (positions<5000)), facecolor='#7ce0ff', alpha=1, interpolate=True)             
        #plot1.plot(positions, np.array(dict_of_wigs['HETU no dps rfa 200 CTD+Rif+'])-0.2, linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'HETU CTD/Rif ({TU_sets_v["HETU no dps rfa 200 CTD+Rif+"]})', zorder=7) #Def linewidth=1 #R23 -0.15; R23 -0.20; R23 -0.24  
        #plot1.fill_between(positions, 1, np.array(dict_of_wigs['HETU no dps rfa 200 CTD+Rif+'])-0.2, where=((positions>-1000) & (positions<1000)), facecolor='#43c287', alpha=1, interpolate=True)        
        #plot1.fill_between(positions, 1, np.array(dict_of_wigs['HETU no dps rfa 200 CTD+Rif+'])-0.2, where=((positions>1000) & (positions<5000)), facecolor='#7ce0ff', alpha=1, interpolate=True)         
        #DNA-gyrase.
        #plot1.plot(positions, np.array(dict_of_wigs['LETU 244 no ybiI'])-0.12, linestyle='-', color='#B6B8BD', linewidth=1, alpha=1, label=f'LETU ({TU_sets_v["LETU 244 no ybiI"]})', zorder=6) #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions, np.array(dict_of_wigs['LETU 244 no ybiI Rif'])+0.05, linestyle='--', color='#B6B8BD', linewidth=0.8, alpha=0.8, label=f'LETU Rif ({TU_sets_v["LETU 244 no ybiI Rif"]})', zorder=5) #R123 -0.1; R12 -0; R3 -0.45             
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs'])-0.12, linestyle='-', color='#333738', linewidth=1, alpha=1, label=f'All TUs ({TU_sets_v["All TUs"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs Rif'])+0.05, linestyle='--', color='#333738', linewidth=0.8, alpha=0.8, label=f'All TUs Rif ({TU_sets_v["All TUs Rif"]})', zorder=9) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42        
        #plot1.plot(positions, np.array(dict_of_wigs['HETU 321 no ompX'])-0.12, linestyle='-', color='#b08642', linewidth=1, alpha=1, label=f'HETU ({TU_sets_v["HETU 321 no ompX"]})', zorder=8) #R123 -0.05; R12 -0.05; R3 -0.25
        #plot1.plot(positions, np.array(dict_of_wigs['HETU 321 no ompX Rif'])+0.05, linestyle='--', color='#b08642', linewidth=0.8, alpha=0.8, label=f'HETU Rif ({TU_sets_v["HETU 321 no ompX Rif"]})', zorder=7) #R123 -0.05; R12 -0.03; R3 -0.32        
        
                    
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    plot1.set_xticks(ticks+list(range(-1000, 0, 500))+list(range(5000, 6000, 500)))
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
    plot1.set_xlim(-1000, 6000) 
    plot1.grid(axis='y', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axhline(1, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.legend(fontsize=15, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'EcTopoI fold enrichment', size=20)
    #plot1.set_title(f'{set_name} signal over {set_type}', size=20)     
    plt.savefig(f'{output_path}\\{set_name}_FE_over_{set_type}_HETU_EcTopoI_noCTD_noRif_346_TU_closer_{win_width}bp_nd_with_body_{length}_bp.png', dpi=400, figsize=(6, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6
    plt.show()
    plt.close()     
    
    #Smoothing.
    dict_of_wigs_sm={}
    for name, wig in dict_of_wigs.items():
        dict_of_wigs_sm[name]=Smoothing(wig, sm_window)
    positions_sm=np.arange(-win_width+sm_window, win_width+length-(sm_window-1), 1)  
    print(positions_sm[0], positions_sm[-1])
        
    #Plot smoothed FE over genes.
    plt.figure(figsize=(6, 6), dpi=100)
    plot1=plt.subplot(111)
 
    ##Transcription units below.
    if set_type=="transcripts":
        #LETU, no dps
        #plot1.plot(positions_sm,  np.array(dict_of_wigs_sm['LETU no dps 200 CTD-Rif-'])-0.02, linestyle='-', color='#B6B8BD', linewidth=2, alpha=1, label=f'LETU Rif- ({TU_sets_v["LETU no dps 200 CTD-Rif-"]})', zorder=20)
        #plot1.plot(positions_sm,  np.array(dict_of_wigs_sm['LETU no dps 200 CTD-Rif+'])-0.1, linestyle='-', color='#B6B8BD', linewidth=2, alpha=1, label=f'LETU Rif+ ({TU_sets_v["LETU no dps 200 CTD-Rif+"]})', zorder=5)
        #plot1.plot(positions_sm,  np.array(dict_of_wigs_sm['LETU no dps 200 CTD+Rif-'])-0.15, linestyle='--', color='#B6B8BD', linewidth=0.8, alpha=0.8, label=f'LETU CTD+Rif- ({TU_sets_v["LETU no dps 200 CTD+Rif-"]})', zorder=6)
        #plot1.plot(positions_sm,  np.array(dict_of_wigs_sm['LETU no dps 200 CTD+Rif+'])-0.25, linestyle='-.', color='#B6B8BD', linewidth=0.8, alpha=0.8, label=f'LETU CTD+Rif+ ({TU_sets_v["LETU no dps 200 CTD+Rif+"]})', zorder=5)          
        #All_TUs no dps
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs no dps 1660 CTD-Rif-'])-0.18, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs ({TU_sets_v["All TUs no dps 1660 CTD-Rif-"]})', zorder=19)
        #plot1.fill_between(positions_sm, 1, np.array(dict_of_wigs_sm['All TUs no dps 1660 CTD-Rif-'])-0.18, where=((positions_sm>-1000) & (positions_sm<1000)), facecolor='#43c287', alpha=1, interpolate=True)        
        #plot1.fill_between(positions_sm, 1, np.array(dict_of_wigs_sm['All TUs no dps 1660 CTD-Rif-'])-0.18, where=((positions_sm>1000) & (positions_sm<5000)), facecolor='#7ce0ff', alpha=1, interpolate=True)          
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs no dps 1660 CTD-Rif+'])-0.1, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs Rif ({TU_sets_v["All TUs no dps 1660 CTD-Rif+"]})', zorder=9)
        #plot1.fill_between(positions_sm, 1, np.array(dict_of_wigs_sm['All TUs no dps 1660 CTD-Rif+'])-0.1, where=((positions_sm>-1000) & (positions_sm<1000)), facecolor='#43c287', alpha=1, interpolate=True)        
        #plot1.fill_between(positions_sm, 1, np.array(dict_of_wigs_sm['All TUs no dps 1660 CTD-Rif+'])-0.1, where=((positions_sm>1000) & (positions_sm<5000)), facecolor='#7ce0ff', alpha=1, interpolate=True)  
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs no dps 1660 CTD+Rif-'])-0.11, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs CTD ({TU_sets_v["All TUs no dps 1660 CTD+Rif-"]})', zorder=10)
        #plot1.fill_between(positions_sm, 1, np.array(dict_of_wigs_sm['All TUs no dps 1660 CTD+Rif-'])-0.11, where=((positions_sm>-1000) & (positions_sm<1000)), facecolor='#43c287', alpha=1, interpolate=True)        
        #plot1.fill_between(positions_sm, 1, np.array(dict_of_wigs_sm['All TUs no dps 1660 CTD+Rif-'])-0.11, where=((positions_sm>1000) & (positions_sm<5000)), facecolor='#7ce0ff', alpha=1, interpolate=True)        
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs no dps 1660 CTD+Rif+'])-0.2, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs CTD/Rif ({TU_sets_v["All TUs no dps 1660 CTD+Rif+"]})', zorder=9)   
        #plot1.fill_between(positions_sm, 1, np.array(dict_of_wigs_sm['All TUs no dps 1660 CTD+Rif+'])-0.2, where=((positions_sm>-1000) & (positions_sm<1000)), facecolor='#43c287', alpha=1, interpolate=True)        
        #plot1.fill_between(positions_sm, 1, np.array(dict_of_wigs_sm['All TUs no dps 1660 CTD+Rif+'])-0.2, where=((positions_sm>1000) & (positions_sm<5000)), facecolor='#7ce0ff', alpha=1, interpolate=True)        
        #HETU, no dps rfa
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU no dps rfa 200 CTD-Rif-'])-0.18, linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'HETU Rif- ({TU_sets_v["HETU no dps rfa 200 CTD-Rif-"]})', zorder=18)
        plot1.fill_between(positions_sm, 1, np.array(dict_of_wigs_sm['HETU no dps rfa 200 CTD-Rif-'])-0.18, where=((positions_sm>-1000) & (positions_sm<1000)), facecolor='#43c287', alpha=1, interpolate=True)        
        plot1.fill_between(positions_sm, 1, np.array(dict_of_wigs_sm['HETU no dps rfa 200 CTD-Rif-'])-0.18, where=((positions_sm>1000) & (positions_sm<5000)), facecolor='#7ce0ff', alpha=1, interpolate=True)         
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU no dps rfa 200 CTD-Rif+'])-0.1, linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'HETU Rif ({TU_sets_v["HETU no dps rfa 200 CTD-Rif+"]})', zorder=7)
        #plot1.fill_between(positions_sm, 1, np.array(dict_of_wigs_sm['HETU no dps rfa 200 CTD-Rif+'])-0.1, where=((positions_sm>-1000) & (positions_sm<1000)), facecolor='#43c287', alpha=1, interpolate=True)        
        #plot1.fill_between(positions_sm, 1, np.array(dict_of_wigs_sm['HETU no dps rfa 200 CTD-Rif+'])-0.1, where=((positions_sm>1000) & (positions_sm<5000)), facecolor='#7ce0ff', alpha=1, interpolate=True)          
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU no dps rfa 200 CTD+Rif-'])-0.11, linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'HETU CTD ({TU_sets_v["HETU no dps rfa 200 CTD+Rif-"]})', zorder=8)
        #plot1.fill_between(positions_sm, 1, np.array(dict_of_wigs_sm['HETU no dps rfa 200 CTD+Rif-'])-0.11, where=((positions_sm>-1000) & (positions_sm<1000)), facecolor='#43c287', alpha=1, interpolate=True)        
        #plot1.fill_between(positions_sm, 1, np.array(dict_of_wigs_sm['HETU no dps rfa 200 CTD+Rif-'])-0.11, where=((positions_sm>1000) & (positions_sm<5000)), facecolor='#7ce0ff', alpha=1, interpolate=True)            
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU no dps rfa 200 CTD+Rif+'])-0.2, linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'HETU CTD/Rif  ({TU_sets_v["HETU no dps rfa 200 CTD+Rif+"]})', zorder=7)    
        #plot1.fill_between(positions_sm, 1, np.array(dict_of_wigs_sm['HETU no dps rfa 200 CTD+Rif+'])-0.2, where=((positions_sm>-1000) & (positions_sm<1000)), facecolor='#43c287', alpha=1, interpolate=True)        
        #plot1.fill_between(positions_sm, 1, np.array(dict_of_wigs_sm['HETU no dps rfa 200 CTD+Rif+'])-0.2, where=((positions_sm>1000) & (positions_sm<5000)), facecolor='#7ce0ff', alpha=1, interpolate=True)         
        #DNA-gyrase.
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['LETU 244 no ybiI'])-0.12, linestyle='-', color='#B6B8BD', linewidth=1.5, alpha=1, label=f'LETU ({TU_sets_v["LETU 244 no ybiI"]})', zorder=6) #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['LETU 244 no ybiI Rif'])+0.05, linestyle='--', color='#B6B8BD', linewidth=0.8, alpha=1, label=f'LETU Rif ({TU_sets_v["LETU 244 no ybiI Rif"]})', zorder=5) #R123 -0.1; R12 -0; R3 -0.45             
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs'])-0.12, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs ({TU_sets_v["All TUs"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs Rif'])+0.05, linestyle='--', color='#333738', linewidth=1.2, alpha=1, label=f'All TUs Rif ({TU_sets_v["All TUs Rif"]})', zorder=9) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42        
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU 321 no ompX'])-0.12, linestyle='-', color='#b08642', linewidth=1.5, alpha=1, label=f'HETU ({TU_sets_v["HETU 321 no ompX"]})', zorder=8) #R123 -0.05; R12 -0.05; R3 -0.25
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU 321 no ompX Rif'])+0.05, linestyle='--', color='#b08642', linewidth=0.8, alpha=1, label=f'HETU Rif ({TU_sets_v["HETU 321 no ompX Rif"]})', zorder=7) #R123 -0.05; R12 -0.03; R3 -0.32        
         
        
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    plot1.set_xticks(ticks+list(range(-1000, 0, 500))+list(range(5000, 6000, 500)))
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
    plot1.set_xlim(-1000, 6000) 
    plot1.grid(axis='y', which='minor', color='grey', linestyle='--', alpha=0.7, linewidth=1)   
    plot1.axhline(1, color='black', linestyle=':', alpha=0.7, linewidth=1.5)
    plot1.legend(fontsize=14.5, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'EcTopoI fold enrichment', size=20)
    #plot1.set_title(f'{set_name} signal over {set_type}', size=20)     
    plt.savefig(f'{output_path}\\{set_name}_FE_over_{set_type}_HETU_EcTopoI_noCTD_noRif_346_TU_closer_smoothed_{win_width}bp_nd_with_body_{length}bp_smoothed_{2*sm_window}bp.png', dpi=400, figsize=(6, 6), transparent=True)   
    plt.show()
    plt.close()    
    return


#EcTopoI
#plot_FE_TUs_groups_local(Wig_data_in_dict_transcripts, Sm_window, Out_path, Signal_name, Set_type)


#######
#Closer look into TUstart and TUend.
#######


def plot_FE_TUs_start_end(wig_in_dict, sm_window, output_path, set_name, set_type):
    #Number of genes in sets.
    TU_sets_v={'All TUs no dps 1660 CTD-Rif-' : 1660, 'All TUs no dps 1660 CTD-Rif+' : 1660, 'All TUs no dps 1660 CTD+Rif-' : 1660, 'All TUs no dps 1660 CTD+Rif+' : 1660,
               'HETU no dps rfa 200 CTD-Rif-' : 200, 'HETU no dps rfa 200 CTD-Rif+' : 200, 'HETU no dps rfa 200 CTD+Rif-' : 200, 'HETU no dps rfa 200 CTD+Rif+' : 200, 
               'LETU no dps 200 CTD-Rif-' : 200, 'LETU no dps 200 CTD-Rif+' : 200, 'LETU no dps 200 CTD+Rif-' : 200, 'LETU no dps 200 CTD+Rif+' : 200,
               'All TUs no dps 1660 RpoC' : 1660, 'HETU no dps rfa 200 RpoC' : 200, 'LETU no dps 200 RpoC' : 200,
               'All TUs no dps 1660 RpoC Rif' : 1660, 'HETU no dps rfa 200 RpoC Rif' : 200, 'LETU no dps 200 RpoC Rif' : 200, 
               'All TUs no dps 1660 RpoB' : 1660, 'HETU no dps rfa 200 RpoB' : 200, 'LETU no dps 200 RpoB' : 200}
    
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
    plt.figure(figsize=(3, 6), dpi=100) #Def size - 10, 6; Closer look - 3, 6
    plot1=plt.subplot(111)
    
    ##Transcription units below.
    if set_type=="transcripts":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        #LETU, no dps
        #plot1.plot(positions,  np.array(dict_of_wigs['LETU no dps 200 CTD-Rif-'])-0.18, linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'LETU ({TU_sets_v["LETU no dps 200 CTD-Rif-"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions,  np.array(dict_of_wigs['LETU no dps 200 CTD-Rif+'])-0.1, linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'LETU CTD-/Rif+ ({TU_sets_v["LETU no dps 200 CTD-Rif+"]})', zorder=5) #Def linewidth=0.8; #R123 -0.1; R12 +0; R3 -0.45    
        #plot1.plot(positions,  np.array(dict_of_wigs['LETU no dps 200 CTD+Rif-'])-0.11, linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'LETU CTD ({TU_sets_v["LETU no dps 200 CTD+Rif-"]})', zorder=6) #Def linewidth=1.5; #R123 -0.1; R123 +0.15; R123 -0.15
        #plot1.plot(positions,  np.array(dict_of_wigs['LETU no dps 200 CTD+Rif+'])-0.2, linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'LETU CTD/Rif ({TU_sets_v["LETU no dps 200 CTD+Rif+"]})', zorder=5) #Def linewidth=1; #R23 -0.2; R23 +0.25; R23 -0.25         
        #All_TUs no dps
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no dps 1660 CTD-Rif-'])-0.18, linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=f'All TUs ({TU_sets_v["All TUs no dps 1660 CTD-Rif-"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17       
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no dps 1660 CTD-Rif+'])-0.1, linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=f'All TUs Rif ({TU_sets_v["All TUs no dps 1660 CTD-Rif+"]})', zorder=9) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42         
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no dps 1660 CTD+Rif-'])-0.11, linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=f'All TUs CTD ({TU_sets_v["All TUs no dps 1660 CTD+Rif-"]})', zorder=10) #Def linewidth=2; #R123 -0.25; R123 -0.25; R123 -0.27
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no dps 1660 CTD+Rif+'])-0.2, linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=f'All TUs CTD/Rif ({TU_sets_v["All TUs no dps 1660 CTD+Rif+"]})', zorder=9) #Def linewidth=1; #R23 -0.25; R23 -0.22; R23 -0.24
        #HETU, no dps rfa
        #plot1.plot(positions, np.array(dict_of_wigs['HETU no dps rfa 200 CTD-Rif-'])-0.18, linestyle='-', color='#b08642', linewidth=2.5, alpha=1, label=f'HETU ({TU_sets_v["HETU no dps rfa 200 CTD-Rif-"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25         
        #plot1.plot(positions, np.array(dict_of_wigs['HETU no dps rfa 200 CTD-Rif+'])-0.1, linestyle='-', color='#b08642', linewidth=2.5, alpha=1, label=f'HETU Rif ({TU_sets_v["HETU no dps rfa 200 CTD-Rif+"]})', zorder=7) #Def linewidth=0.8 #R123 -0.05; R12 -0.03; R3 -0.32       
        #plot1.plot(positions, np.array(dict_of_wigs['HETU no dps rfa 200 CTD+Rif-'])-0.11, linestyle='-', color='#b08642', linewidth=2.5, alpha=1, label=f'HETU CTD ({TU_sets_v["HETU no dps rfa 200 CTD+Rif-"]})', zorder=8) #Def linewidth=1.5 #R123 -0.20; R123 -0.25; R123 -0.28
        #plot1.plot(positions, np.array(dict_of_wigs['HETU no dps rfa 200 CTD+Rif+'])-0.2, linestyle='-', color='#b08642', linewidth=2.5, alpha=1, label=f'HETU CTD/Rif ({TU_sets_v["HETU no dps rfa 200 CTD+Rif+"]})', zorder=7) #Def linewidth=1 #R23 -0.15; R23 -0.20; R23 -0.24  
        #RpoC RNAP
        #plot1.plot(positions,  np.array(dict_of_wigs['LETU no dps 200 RpoC'])-0.22, linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'LETU ({TU_sets_v["LETU no dps 200 RpoC"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no dps 1660 RpoC'])-0.22, linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=f'All TUs ({TU_sets_v["All TUs no dps 1660 RpoC"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17       
        #plot1.plot(positions, np.array(dict_of_wigs['HETU no dps rfa 200 RpoC'])-0.22, linestyle='-', color='#b08642', linewidth=2.5, alpha=1, label=f'HETU ({TU_sets_v["HETU no dps rfa 200 RpoC"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25         
        #plot1.plot(positions, np.array(dict_of_wigs['LETU no dps 200 RpoB']), linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'LETU ({TU_sets_v["LETU no dps 200 RpoB"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no dps 1660 RpoB']), linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=f'All TUs ({TU_sets_v["All TUs no dps 1660 RpoB"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17       
        #plot1.plot(positions, np.array(dict_of_wigs['HETU no dps rfa 200 RpoB']), linestyle='-', color='#b08642', linewidth=2.5, alpha=1, label=f'HETU ({TU_sets_v["HETU no dps rfa 200 RpoB"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25         
        plot1.plot(positions, np.array(dict_of_wigs['LETU no dps 200 RpoC Rif'])+0.28,     linestyle='-', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'LETU ({TU_sets_v["LETU no dps 200 RpoC Rif"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['All TUs no dps 1660 RpoC Rif'])+0.28, linestyle='-', color='#333738', linewidth=2.5, alpha=1, label=f'All TUs ({TU_sets_v["All TUs no dps 1660 RpoC Rif"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17       
        plot1.plot(positions, np.array(dict_of_wigs['HETU no dps rfa 200 RpoC Rif'])+0.28, linestyle='-', color='#b08642', linewidth=2.5, alpha=1, label=f'HETU ({TU_sets_v["HETU no dps rfa 200 RpoC Rif"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25         
        
                    
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    #plot1.set_xticks(ticks+list(range(-300, 200, 20))) #Start
    plot1.set_xticks(ticks+list(range(length-200, length+300, 20))) #End    
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
    #plot1.set_xlim(-300, 200) #Start
    plot1.set_xlim(length-200, length+300) #End
    #plot1.legend(fontsize=15, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'EcTopoI fold enrichment', size=20)  
    plt.savefig(f'{output_path}\\{set_name}_FE_over_{set_type}_RpoC_Rif_TU_end_{win_width}bp_nd_with_body_{length}_bp.png', dpi=400, figsize=(3, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6
    plt.show()
    plt.close()     
    
    return


#EcTopoI
#plot_FE_TUs_start_end(Wig_data_in_dict_transcripts, Sm_window, Out_path, Signal_name, Set_type)
#RNAP
#plot_FE_TUs_start_end(Wig_data_in_dict_transcripts_RNApol, Sm_window, Out_path, Signal_name, Set_type)


#######
#Plot the signal for different groups of genes together.
#######


def plot_FE_ChIP_groups(wig_in_dict, sm_window, output_path, set_name, set_type):
    #Number of genes within sets.
    TU_sets_v={'RNAP_Borukhov' : 1660, 'RNAP_Borukhov_HETU' : 200, 'RpoD_Myers' : 1660, 'RpoS_Seo' : 1660, 
               'TopA_Sutormin_All_TUs' : 1660, 'TopA_Sutormin_HETU' : 200, 'TopA_Sutormin_LETU' : 200,
               'Gyrase_Sutormin_All_TUs' : 1660, 'Gyrase_Sutormin_HETU' : 200, 'Gyrase_Rif_Sutormin_All_TUs' : 1660, 
               'TopA_delta11_All_TUs' : 1660, 'TopA_delta14_All_TUs' : 1660, 'TopA_delta30_All_TUs' : 1660, 
               'TopA_delta11_HETU' : 200, 'TopA_delta14_HETU' : 200, 'TopA_delta30_HETU' : 200, 
               'TopA_delta11_LETU' : 200, 'TopA_delta14_LETU' : 200, 'TopA_delta30_LETU' : 200,
               'TopA_Y319F_plus_All_TUs' : 1660, 'TopA_Y319F_plus_HETU' : 200, 'TopA_Y319F_plus_LETU' : 200}
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
    if set_type=="transcripts":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        #Normalizing coefficients: 1.20 (Mean), 2.11 (STD)
        plot1.plot(positions, (np.array(dict_of_wigs['RNAP_Borukhov'])-1.20)/2.11/2, linestyle='-', color='#B6B8BD', linewidth=2, alpha=1, label=f'All TUs RNAP (1660)', zorder=3) #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions, ((np.array(dict_of_wigs['RNAP_Borukhov_HETU'])-1.20)/2.11)/2, linestyle='-', color='#B6B8BD', linewidth=2, alpha=1, label=f'RNAP HETU (200)', zorder=3) #R123 -0; R12 +0.03; R3 -0.07        
        #Normalizing coefficients: 1.04 (Mean), 3.17 (STD)
        #plot1.plot(positions, (np.array(dict_of_wigs['RpoD_Myers'])-1.04)/3.17, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'RpoD', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        #Normalizing coefficients: 1.18 (Mean), 9.51 (STD)
        #plot1.plot(positions, (np.array(dict_of_wigs['RpoS_Seo'])-1.18)/9.51, linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'RpoS', zorder=1) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        #Normalizing coefficients: 1.26 (Mean), 2.34 (STD)
        plot1.plot(positions, ((np.array(dict_of_wigs['TopA_Sutormin_All_TUs'])-1.26)/2.34)+0.04, linestyle='-', color='#B08642', linewidth=2, alpha=1, label=f'All TUs EcTopoI (1660)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions, ((np.array(dict_of_wigs['TopA_Sutormin_HETU'])-1.26)/2.34)+0.07, linestyle='-', color='#B08642', linewidth=2, alpha=1, label=f'EcTopoI HETU (200)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions, (np.array(dict_of_wigs['TopA_Sutormin_LETU'])-1.26)/2.34, linestyle='-', color='#B08642', linewidth=2, alpha=1, label=f'EcTopoI LETU (200)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42          
        #plot1.plot(positions, np.array(dict_of_wigs['Gyrase_Rif_Sutormin_All_TUs'])+0.05, linestyle='-', color='#B6B8BD', linewidth=2, alpha=1, label=f'Gyrase +Rif (1660)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #Normalizing coefficients: 1.14 (Mean), 0.58 (STD)
        plot1.plot(positions, ((np.array(dict_of_wigs['Gyrase_Sutormin_All_TUs'])-1.14)/0.58)+0.05, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs Gyrase (1660)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions, ((np.array(dict_of_wigs['Gyrase_Sutormin_HETU'])-1.14)/0.58)-0.12, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Gyrase HETU (200)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42 
        #plot1.plot(positions, (np.array(dict_of_wigs['TopA_delta11_All_TUs'])-1.043)/0.711, linestyle='-', color='#B6B8BD', linewidth=2, alpha=1, label=f'All TUs EcTopoI delta11 (1660)', zorder=3) #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions, (np.array(dict_of_wigs['TopA_delta14_All_TUs'])-1.164)/1.080, linestyle='-', color='#B08642', linewidth=2, alpha=1, label=f'All TUs EcTopoI delta14 (1660)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions, (np.array(dict_of_wigs['TopA_delta30_All_TUs'])-1.140)/0.800, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs EcTopoI delta30 (1660)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions, (np.array(dict_of_wigs['TopA_delta11_HETU'])-1.043)/0.711, linestyle='-', color='#B6B8BD', linewidth=2, alpha=1, label=f'HETU EcTopoI delta11 (200)', zorder=3) #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions, (np.array(dict_of_wigs['TopA_delta14_HETU'])-1.164)/1.080, linestyle='-', color='#B08642', linewidth=2, alpha=1, label=f'HETU EcTopoI delta14 (200)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions, (np.array(dict_of_wigs['TopA_delta30_HETU'])-1.140)/0.800, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'HETU EcTopoI delta30 (200)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions, (np.array(dict_of_wigs['TopA_delta11_LETU'])-1.043)/0.711, linestyle='-', color='#B6B8BD', linewidth=2, alpha=1, label=f'LETU EcTopoI delta11 (200)', zorder=3) #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions, (np.array(dict_of_wigs['TopA_delta14_LETU'])-1.164)/1.080, linestyle='-', color='#B08642', linewidth=2, alpha=1, label=f'LETU EcTopoI delta14 (200)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions, (np.array(dict_of_wigs['TopA_delta30_LETU'])-1.140)/0.800, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'LETU EcTopoI delta30 (200)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42    
        #plot1.plot(positions, (np.array(dict_of_wigs['TopA_Y319F_plus_All_TUs'])-1.21)/1.80, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs EcTopoI Y319Fplus (1660)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions, (np.array(dict_of_wigs['TopA_Y319F_plus_HETU'])-1.21)/1.80,    linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs EcTopoI Y319Fplus (200)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions, (np.array(dict_of_wigs['TopA_Y319F_plus_LETU'])-1.21)/1.80,    linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs EcTopoI Y319Fplus (200)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
              
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    #plot1.set_xticks(ticks+list(range(-300, 200, 20))) #Start
    #plot1.set_xticks(ticks+list(range(4800, 5300, 20))) #End
    plot1.set_xticks(ticks, minor='False')
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    ticks_lables[ticks.index(length)]='TE'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1) 
    plot1.set_xticks([0, length], minor='True')
    #plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axvline(length, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.set_yticks([1], minor='True')
    #plot1.set_xlim(-300, 200) #Start
    #plot1.set_xlim(4800, 5300) #End
    #plot1.grid(axis='y', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axhline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.legend(fontsize=15, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'Fold enrichment, norm', size=20)
    #plot1.set_yscale('log')
    #plot1.set_title(f'{set_name} signal over {set_type}', size=20)     
    plt.savefig(f'{output_path}\\{set_name}_All_TUs_Norm_FE_over_{set_type}_All_TUs_EcTopoI_vs_Gyrase_vs_RNAP_{win_width}bp_nd_with_body_{length}_bp_end.png', dpi=400, figsize=(10, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6
    plt.show()
    plt.close()     
    
    #Smoothing.
    dict_of_wigs_sm={}
    for name, wig in dict_of_wigs.items():
        dict_of_wigs_sm[name]=Smoothing(wig, sm_window)
    positions_sm=np.arange(-win_width+sm_window, win_width+length-(sm_window-1), 1)  
    print(positions_sm[0], positions_sm[-1])
        
    #Plot smoothed FE over genes.
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)
    
    
    ##Transcription units below.
    if set_type=="transcripts":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        #Normalizing coefficients: 1.20 (Mean), 2.11 (STD)
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['RNAP_Borukhov'])-1.20)/2.11/2, linestyle='-', color='#B6B8BD', linewidth=2, alpha=1, label=f'All TUs RNAP (1660)', zorder=3) #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions_sm, ((np.array(dict_of_wigs_sm['RNAP_Borukhov_HETU'])-1.20)/2.11)/2, linestyle='-', color='#B6B8BD', linewidth=2, alpha=1, label=f'RNAP HETU (200)', zorder=3) #R123 -0; R12 +0.03; R3 -0.07
        #Normalizing coefficients: 1.04 (Mean), 3.17 (STD)
        #plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['RpoD_Myers'])-1.04)/3.17, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'RpoD', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        #Normalizing coefficients: 1.18 (Mean), 9.51 (STD)
        #plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['RpoS_Seo'])-1.18)/9.51, linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'RpoS', zorder=1) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        #Scaling coefficients: 1.26 (Mean), 2.34 (STD)
        plot1.plot(positions_sm, ((np.array(dict_of_wigs_sm['TopA_Sutormin_All_TUs'])-1.26)/2.34)+0.04, linestyle='-', color='#B08642', linewidth=2, alpha=1, label=f'All TUs EcTopoI (1660)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions_sm, ((np.array(dict_of_wigs_sm['TopA_Sutormin_HETU'])-1.26)/2.34)+0.07, linestyle='-', color='#B08642', linewidth=2, alpha=1, label=f'EcTopoI HETU (200)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['TopA_Sutormin_LETU'])-1.26)/2.34, linestyle='-', color='#B08642', linewidth=2, alpha=1, label=f'EcTopoI LETU (200)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Gyrase_Rif_Sutormin_All_TUs'])+0.05, linestyle='-', color='#B6B8BD', linewidth=2, alpha=1, label=f'Gyrase +Rif (1660)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #Scaling coefficients: 1.14 (Mean), 0.58 (STD)
        plot1.plot(positions_sm, ((np.array(dict_of_wigs_sm['Gyrase_Sutormin_All_TUs'])-1.14)/0.58)+0.05, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs Gyrase (1660)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions_sm, ((np.array(dict_of_wigs_sm['Gyrase_Sutormin_HETU'])-1.14)/0.58)-0.12, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Gyrase HETU (200)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['TopA_delta11_All_TUs'])-1.043)/0.711, linestyle='-', color='#B6B8BD', linewidth=2, alpha=1, label=f'All TUs EcTopoI delta11 (1660)', zorder=3) #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['TopA_delta14_All_TUs'])-1.164)/1.080, linestyle='-', color='#B08642', linewidth=2, alpha=1, label=f'All TUs EcTopoI delta14 (1660)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['TopA_delta30_All_TUs'])-1.140)/0.800, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs EcTopoI delta30 (1660)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['TopA_delta11_HETU'])-1.043)/0.711, linestyle='-', color='#B6B8BD', linewidth=2, alpha=1, label=f'HETU EcTopoI delta11 (200)', zorder=3) #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['TopA_delta14_HETU'])-1.164)/1.080, linestyle='-', color='#B08642', linewidth=2, alpha=1, label=f'HETU EcTopoI delta14 (200)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['TopA_delta30_HETU'])-1.140)/0.800, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'HETU EcTopoI delta30 (200)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['TopA_delta11_LETU'])-1.043)/0.711, linestyle='-', color='#B6B8BD', linewidth=2, alpha=1, label=f'LETU EcTopoI delta11 (200)', zorder=3) #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['TopA_delta14_LETU'])-1.164)/1.080, linestyle='-', color='#B08642', linewidth=2, alpha=1, label=f'LETU EcTopoI delta14 (200)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['TopA_delta30_LETU'])-1.140)/0.800, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'LETU EcTopoI delta30 (200)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['TopA_Y319F_plus_All_TUs'])-1.21)/1.80, linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs EcTopoI Y319Fplus (1660)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['TopA_Y319F_plus_HETU'])-1.21)/1.80,    linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs EcTopoI Y319Fplus (200)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        #plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['TopA_Y319F_plus_LETU'])-1.21)/1.80,    linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'All TUs EcTopoI Y319Fplus (200)', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42  
        
        
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    ticks_lables[ticks.index(length)]='TE'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1)
    plot1.set_xticks([0, length], minor='True')
    #plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)   
    plot1.axvline(0, color='black', linestyle=':', alpha=0.7, linewidth=1.5)
    plot1.axvline(length, color='black', linestyle=':', alpha=0.7, linewidth=1.5)    
    plot1.set_yticks([1], minor='True')
    #plot1.grid(axis='y', which='minor', color='grey', linestyle='--', alpha=0.7, linewidth=1)   
    plot1.axhline(0, color='black', linestyle=':', alpha=0.7, linewidth=1.5)
    plot1.legend(fontsize=14.5, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'Fold enrichment, norm', size=20)
    #plot1.set_yscale('log')
    #plot1.set_title(f'{set_name} signal over {set_type}', size=20)     
    plt.savefig(f'{output_path}\\{set_name}_All_TUs_Norm_FE_over_{set_type}_All_TUs_EcTopoI_vs_Gyrase_vs_RNAP_{win_width}bp_nd_with_body_{length}bp_smoothed_{2*sm_window}bp.png', dpi=400, figsize=(10, 6), transparent=True)
    plt.savefig(f'{output_path}\\{set_name}_All_TUs_Norm_FE_over_{set_type}_All_TUs_EcTopoI_vs_Gyrase_vs_RNAP_{win_width}bp_nd_with_body_{length}bp_smoothed_{2*sm_window}bp.svg', dpi=400, figsize=(10, 6), transparent=True) 
    plt.close()    
    return

#plot_FE_ChIP_groups(Wig_data_in_dict_transcripts_tr_assoc, Sm_window, Out_path, Signal_name, Set_type)
