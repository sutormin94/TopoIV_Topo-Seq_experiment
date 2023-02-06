###############################################
##Dmitry Sutormin, 2022##
##Topo-Seq analysis##

#The script takes results of scanning procedure (WIG file), returns scores for
#TCSs (writes TAB files contain coordinate\tN3E\tScore info), 
#computes Pearson correlation between N3E and score, plots (Score, N3E) scatter plots.

###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import numpy as np
import scipy
from scipy.stats import pearsonr

#######
#Variables to be defined.
#######

#Path to the input WIG score file.
path_to_res_score_files="Data_analysis\TCSs_analysis\Cfx_and_S83L_Cfx\\Unique_GCSs_E_coli_w3110_G_Mu_genome_score.wig"

#Input data - TCSs, TAB.
path_to_TCSs_files={'Cfx wt'  : "Data_analysis\TCSs_analysis\TCSs_sets\\Cfx_trusted_GCSs.txt",
                    'Cfx mut' : "Data_analysis\TCSs_analysis\TCSs_sets\\Cfx_S83L_trusted_GCSs.txt"}

#Output data - TCSs, TAB (with score info added).
path_to_TCSs_sets_with_score="Data_analysis\TCSs_analysis\TCSs_sets\\"
if not os.path.exists(path_to_TCSs_sets_with_score):
    os.makedirs(path_to_TCSs_sets_with_score)
Outpath_to_TCSs_files={'Cfx wt'  : path_to_TCSs_sets_with_score + "TopoIV_TCSs_wt_h_s.txt",
                       'Cfx mut' : path_to_TCSs_sets_with_score + "TopoIV_TCSs_gyrA_S83L_h_s.txt"}

#Output directory for N3E-score plots.
N3E_score_plot_path=path_to_TCSs_sets_with_score

#######
#Parsing and preparing score files.
#######

def Parsing_score(path_s):
    #Parse input WIG score file.
    filein=open(path_s, "r")
    score_ar=[]
    for line in filein:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            score_ar.append(float(line[0]))
    filein.close()
    return score_ar

#######
#Trusted GCSs data parsing.
#Calculate and return score for GCSs.
#Write GCSs files with scores info obtained.
#######

def GCSs_parsing_score_returning_info_writing(input_dict, score_ar, output_dict):
    GCSs_sets_dict={}
    for k, v in input_dict.items():
        GCSs_dict={}
        filein=open(v, 'r')
        fileout=open(output_dict[k], 'w')
        fileout.write('GCSs_coordinate\tN3E\tScore\n')
        for line in filein:
            line=line.rstrip().split('\t')
            if line[0] not in ['GCSs_coordinate']:
                GCSs_dict[int(line[0])]=[float(line[1]), (score_ar[int(line[0])-1+1]+score_ar[int(line[0])-1+4])/2]
                fileout.write(str(line[0]) + '\t' + str(line[1]) + '\t' + str((score_ar[int(line[0])-1+1]+score_ar[int(line[0])-1+4])/2) + '\n')
            else:
                continue
        filein.close()
        fileout.close()
        GCSs_sets_dict[k]=GCSs_dict
        print('Number of trusted GCSs for ' + str(k) + ' : ' + str(len(GCSs_dict)))
    return GCSs_sets_dict

#######
#Calculate N3E-score Pearson correlation.
#######

def correlation_h_s(GCSs_sets_dict):
    HS_dict={}
    for ab, GCSs_set in GCSs_sets_dict.items():
        #[[N3E data], [Score data]]
        HS_dict[ab]=[] 
        N3E=[]
        Score=[]
        for k, v in GCSs_set.items():
            N3E.append(v[0])
            Score.append(v[1])
        #Sorting.
        Score, N3E = (list(t) for t in zip(*sorted(zip(Score, N3E))))
        pearson_cor=scipy.stats.pearsonr(N3E, Score)
        print('Paerson correlation (N3E, score) for ' + str(ab) + ' : ' + str(pearson_cor))
        HS_dict[ab].append(N3E)
        HS_dict[ab].append(Score)
    return HS_dict, pearson_cor

#######
#Visualize (Score, N3E) dependencies.
#######

def Plot_N3E_score(HS_dict, plot_path):
    for dataset_name, dataset in HS_dict.items():
        fig=plt.figure(figsize=(15,15), dpi=100)

        plot0=plt.subplot2grid((1,1),(0,0))    
        fit=np.polyfit(HS_dict[dataset_name][1], HS_dict[dataset_name][0], 1)
        print(fit)
        fit_fn=np.poly1d(fit)  
        plot0.plot(HS_dict[dataset_name][1], HS_dict[dataset_name][0], 'o', fillstyle='none', color='#7FCE79', markeredgecolor='#7FCE79', markersize=2, alpha=0.8)
        plot0.plot(HS_dict[dataset_name][1], fit_fn(HS_dict[dataset_name][1]), '--k', label='y='+str(round(fit[0], 3))+'x+'+str(round(fit[1], 3)))
        plot0.set_xlabel('TCSs score', size=17) 
        plot0.set_ylabel('TCSs N3E', size=17)
        plot0.legend(loc='upper right', fontsize=17)
        plot0.set_title(dataset_name, size=18)
        
        plt.tight_layout()
        plt.savefig(plot_path + dataset_name + "_TCSs_N3E_score_correlations.png", dpi=400, figsize=(15, 15)) 
        plt.close()
        
    return

#######
#Visualize Score, N3E distributions.
#######

def Plot_N3E_score_distribution(HS_dict, score_ar, plot_path):
    for dataset_name, dataset in HS_dict.items():
        #matplotlib.rc('text', usetex = True)
        fig=plt.figure(figsize=(15,15), dpi=100)
        #Dataset N3E
        plot0=plt.subplot2grid((2,3),(0,0), rowspan=1, colspan=1)     
        plot0.hist(HS_dict[dataset_name][0], color='#7FCE79', edgecolor='black', alpha=0.8)
        plot0.annotate('Mean N3E='+str(round(np.mean(HS_dict[dataset_name][0]),2)), xy=(0.45, 0.9), xycoords='axes fraction', size=15)
        plot0.set_xlabel('TCSs N3E', size=17)
        plot0.set_ylabel('Number of TCSs', size=17)
        plot0.set_title(dataset_name, size=18)
        #Dataset sequence score
        plot1=plt.subplot2grid((2,3),(1,0), rowspan=1, colspan=1)     
        plot1.hist(HS_dict[dataset_name][1], color='#7FCE79', edgecolor='black', alpha=0.8)
        plot1.annotate('Mean score=\n'+str(round(np.mean(HS_dict[dataset_name][1]),2)), xy=(0.60, 0.8), xycoords='axes fraction', size=15)
        plot1.set_xlabel('TCSs score', size=17)
        plot1.set_ylabel('Number of TCSs', size=17)
        #General sequence score
        plot8=plt.subplot2grid((2,3),(0,1), rowspan=2, colspan=2)  
        plot8.hist(score_ar, color='#5762ff', edgecolor='black', alpha=0.8)
        plot8.annotate('Mean score='+str(round(np.mean(score_ar),2)), xy=(0.65, 0.9), xycoords='axes fraction',  weight="bold", size=18)
        plot8.set_xlabel('Score', size=17) 
        plot8.set_ylabel('Number of genome positions', size=17)
        plot8.set_title('Scores distribution for $\it{E. coli}$ DY330 MuSGS genome', size=18)
        plt.tight_layout()
        plt.savefig(plot_path + dataset_name + "_TCSs_N3E_score_distributions.png", dpi=400, figsize=(15, 15)) 
        plt.close()    
    return

#######
#Wrapps all the functions together.
#######

def Wrapper(path_s, input_dict, output_dict, plot_path):
    score_track=Parsing_score(path_s)
    GCSs_info_dict=GCSs_parsing_score_returning_info_writing(input_dict, score_track, output_dict)
    Sorted_for_plot_dict, pearson_cor=correlation_h_s(GCSs_info_dict)
    Plot_N3E_score(Sorted_for_plot_dict, plot_path)
    Plot_N3E_score_distribution(Sorted_for_plot_dict, score_track, plot_path)
    return
    
Wrapper(path_to_res_score_files, path_to_TCSs_files, Outpath_to_TCSs_files, N3E_score_plot_path)   

print('Script ended its work succesfully!')