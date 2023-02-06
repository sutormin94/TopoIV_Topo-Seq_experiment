###############################################
##Dmitry Sutormin, 2022##
##TopoIV Topo-Seq analysis##

#The script takes tetrades of WIG files contain N3E or N5E values: A+IP+, A+IP-, A-IP+, A-IP-.
#It smooths A+IP- and A-IP- tracks and divides A+IP+ and A-IP+ by them.
#Alongside N3E genomic tracks, density of TCSs or GCSs is plotted, and transcription levels.
#A plot is generated: signal coverage over the genome for treated and untreated samples.

#Requirements: TAB file with deletions.
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import stats
from Bio import SeqIO
from Bio.SeqUtils import GC as GC_count
import matplotlib.patheffects as PathEffects
import locale

locale.getdefaultlocale()

#######
#Variables to be defined.
#######

#Path to the working directory
pwd="Data_analysis\\"
#Path to the file with regions to be omitted (e.g. deletions).
Deletions="TopoIV_Topo-Seq_experiment\Additional_genome_features\Deletions_w3110_G_Mu_SGS.broadPeak"
#Path to the reference genome
Genome="TopoIV_Topo-Seq_experiment\Additional_genome_features\E_coli_w3110_G_Mu.fasta"


#######
#Opens and reads FASTA file with reference genome
#######

def read_genome(genome_path):
    #Opens FASTA file with the reference genome.
    genome=open(genome_path, 'r')
    for record in SeqIO.parse(genome, "fasta"):
        genomefa=str(record.seq)
    return genomefa

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
#Read filewith TCSs data.
#######

def read_TCSs(TCSs_path):
    TCSs_ar=[]
    filein=open(TCSs_path, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        TCSs_ar.append([int(line[1]), int(line[2]), float(line[6])])
    filein.close()
    return TCSs_ar

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
            Pos_value=float(line[0])
            NE_values.append(Pos_value)
            Total_NE+=Pos_value
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
        #if i%100000==0:
        #    print(i)
        mean=mean + (get_value(i+window, ends, deletions) - get_value(i-window, ends, deletions))/(2*window_float)
        smoothed.append(mean)
    return smoothed

#######
#Returns A+IP+/smoothed(A+IP-) and A-IP+/smoothed(A-IP-) tracks ready for GCSs calling.
#######

def norm_smooth_devide(ex_file_path, cont_file_path, un_ex_file_path, un_cont_file_path, deletions):
    #WIG parsing, total NE counting (for further normalization on reads number)
    treated_experiment=wig_parsing(ex_file_path) #+A+IP
    treated_control=wig_parsing(cont_file_path) #+A-IP
    untreated_experiment=wig_parsing(un_ex_file_path) #-A+IP
    untreated_control=wig_parsing(un_cont_file_path) #-A-IP
    #Normalization on the coverage of neutral region - NR (does not influenced by replication). 
    #NR is set to be 0:2500000 bp.
    #Adds pseudocounts to avoid zero values
    te=np.mean(treated_experiment[0][0:2500000])
    tc=np.mean(treated_control[0][0:2500000]) 
    ue=np.mean(untreated_experiment[0][0:2500000])
    uc=np.mean(untreated_control[0][0:2500000])
    norm_array=[te, tc, ue, uc]
    print(norm_array)
    Min_total_NE=min(norm_array)
    print('Min_total_NE: ' + str(Min_total_NE))
    treated_experiment_norm=[1.0 * (x + 1) * Min_total_NE/norm_array[0] for x in treated_experiment[0]] #+A+IP norm
    print('Passed!')
    treated_control_norm=[1.0 * (x + 1) * Min_total_NE/norm_array[1] for x in treated_control[0]] #+A-IP norm
    print('Passed!')
    untreated_experiment_norm=[1.0 * (x + 1) * Min_total_NE/norm_array[2] for x in untreated_experiment[0]] #-A+IP norm
    print('Passed!')
    untreated_control_norm=[1.0 * (x + 1) * Min_total_NE/norm_array[3] for x in untreated_control[0]] #-A-IP norm
    print('Passed!')
    #Control samples smoothing: A+IP- and A-IP-
    un_experiment_norm_sm=Smoothing(untreated_experiment_norm, deletions) #-A+IP norm sm 
    un_control_norm_sm=Smoothing(untreated_control_norm, deletions) #-A-IP norm sm
    #Pairwise division: +A+IP/-A+IP and +A-IP/-A-IP
    ends_divide_IP=[] #+A+IP/-A+IP
    ends_divide_mock=[] #+A-IP/-A-IP
    for i in range (len(treated_experiment_norm)):
        if treated_experiment_norm[i]!=0 and un_experiment_norm_sm[i]!=0:
            ends_divide_IP.append(treated_experiment_norm[i]/un_experiment_norm_sm[i])
        else:
            ends_divide_IP.append(0)
        if treated_control_norm[i]!=0 and un_control_norm_sm[i]!=0:
            ends_divide_mock.append(treated_control_norm[i]/un_control_norm_sm[i])
        else:
            ends_divide_mock.append(0) 
    return ends_divide_IP, ends_divide_mock, un_experiment_norm_sm, un_control_norm_sm


#######
#Genomic track binning.
#######

def bin_track(bin_edges, track_to_bin):
    Bin_start=0
    N3E_binned=[]
    for i in range(len(bin_edges)):
        Bin_end=int(bin_edges[i])
        if Bin_start<=len(track_to_bin)-1 and Bin_end<=len(track_to_bin)-1:
            N3E_mean_bin=np.mean(track_to_bin[Bin_start:Bin_end])
        elif Bin_start<=len(track_to_bin)-1 and Bin_end>len(track_to_bin)-1:
            N3E_mean_bin=np.mean(track_to_bin[Bin_start:])
        elif Bin_start>len(track_to_bin)-1 and Bin_end>len(track_to_bin)-1:
            N3E_mean_bin=0
        N3E_binned.append(N3E_mean_bin)
        Bin_start=Bin_end    
    
    return N3E_binned


#######
#Plots the enrichment signal over the genome: +A+IP/smoothed(-A+IP) and +A-IP/smoothed(-A-IP)
#######

def plot_enrichment_signal(fname, IP_nd_ends, mock_nd_ends, un_IP_sm, un_mock_sm, deletions, TCSs_coords, RNA_seq_data, path_out):
    #Some hack to avoid some bug in matplotlib (OverflowError: In draw_path: Exceeded cell block limit)
    #See: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit
    mpl.rcParams['agg.path.chunksize']=10000
    #Scaling smoothed tracks to make them visible on the plot.
    max_element=max(IP_nd_ends+mock_nd_ends) #Max N3E value of experimental tracks
    max_element_IP_sm=max(un_IP_sm)
    max_element_mock_sm=max(un_mock_sm)    
    if max_element<1000:
        un_IP_sm=[(max_element/2)*x/max_element_IP_sm for x in un_IP_sm]
        un_mock_sm=[(max_element/2)*x/max_element_mock_sm for x in un_mock_sm]
    else:
        un_IP_sm=[(max_element/4)*x/max_element_IP_sm for x in un_IP_sm]
        un_mock_sm=[(max_element/4)*x/max_element_mock_sm for x in un_mock_sm]        
    
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
    un_IPed=np.ma.masked_array(un_IP_sm, mask=mask_array)
    un_mock=np.ma.masked_array(un_mock_sm, mask=mask_array)
    RNA_seq_data_msk=np.ma.masked_array(RNA_seq_data, mask=mask_array)
    
    
    #Get left-most TCSs positions.
    TSCs_ar=[]
    for TCS in TCSs_coords:
        TCS_left=TCS[0]
        TSCs_ar.append(TCS_left)
        
    #Unweighted TCSs binning - count the number of TCSs in genomic bins.
    bin_width=40000
    chrom_len=4647999
    num_bin=chrom_len//bin_width
    hist, bin_edges=np.histogram(TSCs_ar, bins=num_bin+1) #Taken from this discussion: https://stackoverflow.com/questions/38772640/rebinning-a-list-of-numbers-in-python
    hist=list(hist)+[0]
    bin_edges=np.array(list(bin_edges))+bin_width
    print(bin_edges)
    
    #Weighted TCSs binning - calculated the average cleavage signal in genomic bins (average N3E enrichment).
    TCSs_coords_sorted=sorted(TCSs_coords, key=lambda x: x[0])
    print(TCSs_coords_sorted)
    Bin_start=0
    Start_TCSs_index=0
    Hist_uw=[]
    Hist_w_av=[]
    Hist_w=[]
    for i in range(len(bin_edges)):
        Bin_end=bin_edges[i]
        TCSs_count=0
        TCSs_N3E=0
        for j in range(len(TCSs_coords_sorted)-Start_TCSs_index):
            Start_index=Start_TCSs_index+j
            TCSs_coord=TCSs_coords_sorted[Start_index][0]
            if Bin_start<=TCSs_coord<Bin_end:
                TCSs_count+=1
                TCSs_N3E+=TCSs_coords_sorted[Start_index][2]
            elif TCSs_coord>Bin_end:
                break
        
        if TCSs_count>0:
            TCSs_av_N3E=TCSs_N3E/TCSs_count
        else:
            TCSs_av_N3E=0
        print(Bin_start, Bin_end, TCSs_count, TCSs_av_N3E)    
        Hist_uw.append(TCSs_count)
        Hist_w_av.append(TCSs_av_N3E)
        Hist_w.append(TCSs_N3E)
        Bin_start=Bin_end
        
    #Bin N3E data for IP sample.
    N3E_binned=bin_track(bin_edges, IP_nd_ends)
        
    #Bin RNA-Seq data.
    RNA_seq_binned=bin_track(bin_edges, RNA_seq_data_msk)
    
    #Plotting the distribution of the signal around the genome for IPed and mock samples.
    xcoord=np.arange(0,chrom_len)
    fig, plots=plt.subplots(6, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1, 1, 1, 1, 1]}, figsize=(10, 8))
        
    plots[0].plot(xcoord, IPed, '-', label='+A+IP/smoothed(-A+IP)', color='black', linewidth=0.5)
    plots[0].plot(xcoord, mock, '-', label='+A-IP/smoothed(-A-IP)', color='orange', linewidth=0.5)
    plots[0].plot(xcoord, un_IPed, '-', label='smoothed(-A+IP)', color='#5bbdff', linewidth=3)
    plots[0].plot(xcoord, un_mock, '-', label='smoothed(-A-IP)', color='#ed781f', linewidth=3)  
    plots[0].set_ylim([-20,500])
    plots[0].set_ylabel('Topo IV N3E FE', size=17)
    plots[0].legend(loc='upper right', frameon=False)
    plots[1].step(bin_edges, hist, color='black')    
    plots[2].step(bin_edges, Hist_w_av, color='black')  
    plots[3].step(bin_edges, Hist_w, color='black')  
    plots[4].step(bin_edges, N3E_binned, color='black') 
    plots[5].step(bin_edges, RNA_seq_binned, color='black') 
    plots[5].set_yscale('log')
    plots[5].set_xlabel('Genome position, Mb', size=17)     
    plt.show()
    plt.savefig(path_out+fname+'_signal_enrichment_binned.png', dpi=300, figsize=(16, 8))
    plt.savefig(path_out+fname+'_signal_enrichment_binned.svg', dpi=300, figsize=(16, 8))
    plt.close()
    return

#######
#Wraps all the functions: data normalization, GCSs calling, 
#plotting, writing the data.
#######

def Wrapper_func(tetrade_dictionary, deletions_inpath, genome_path, tcss_coord_path, rna_seq_data, path_out):
    #Parsing deletions
    deletions=deletions_info(deletions_inpath)
    
    #Parsing TCSs coords.
    TCSs_coords=read_TCSs(tcss_coord_path)    
    
    #Read RNA-Seq data.
    RNA_seq_data=wig_parsing(rna_seq_data)[0]
    
    #Define samples within the tetrade.
    Tet_ID=tetrade_dictionary['Tetrade name']
    print('Now we are working with: ' + str(Tet_ID))
    treated_control=tetrade_dictionary['A+IP-']
    treated_experiment=tetrade_dictionary['A+IP+']
    untreated_control=tetrade_dictionary['A-IP-']
    untreated_experiment=tetrade_dictionary['A-IP+']
    
    #Obtain pairwise divided tracks: A+IP+/A+IP- and A-IP+/A-IP-.
    ends_fitting=norm_smooth_devide(treated_experiment, treated_control, untreated_experiment, untreated_control, deletions)
    IP_norm_div=ends_fitting[0]
    mock_norm_div=ends_fitting[1]
    un_IP_sm=ends_fitting[2]
    un_mock_sm=ends_fitting[3]    
    
    #Plotting the distribution of the signal around the genome for treated and untreated.
    plot_enrichment_signal(Tet_ID, IP_norm_div, mock_norm_div, un_IP_sm, un_mock_sm, deletions, TCSs_coords, RNA_seq_data, path_out)    

    return



Tetrade_replic_1_10mkM_gyrS83L={'A+IP+': pwd + "WIG\gyrA-S83L_replic_1_Cfx_IP_N3E.wig",
                                'A+IP-': pwd + "WIG\gyrA-S83L_replic_1_Cfx_no_IP_N3E.wig",
                                'A-IP+': pwd + "WIG\gyrA-S83L_replic_1_no_Cfx_IP_N3E.wig",   
                                'A-IP-': pwd + "WIG\gyrA-S83L_replic_1_no_Cfx_no_IP_N3E.wig",
                                'Tetrade name': 'TopoIV_Replic_1_Cfx_10_mkM_gyrA_S83L'
                                }

Tetrade_Cfx_Gyrase={'A+IP+': "TopoIV_Topo-Seq_experiment\Additional_genome_features\\Gyrase_Cfx_IP_Mu_10mkM_1_edt_N3E.wig", 
                    'A+IP-': "TopoIV_Topo-Seq_experiment\Additional_genome_features\\Gyrase_Cfx_no_IP_Mu_10mkM_1_edt_N3E.wig",
                    'A-IP+': "TopoIV_Topo-Seq_experiment\Additional_genome_features\\Gyrase_no_Cfx_IP_Mu_1_edt_N3E.wig", 
                    'A-IP-': "TopoIV_Topo-Seq_experiment\Additional_genome_features\\Gyrase_no_Cfx_no_IP_Mu_1_edt_N3E.wig",
                    'Tetrade name': 'Gyrase_Replic_1_Cfx_10_mkM'
                    }

Tetrade_RifCfx_Gyrase={'A+IP+': "TopoIV_Topo-Seq_experiment\Additional_genome_features\\Gyrase_RifCfx_IP_Mu_122mkM_10mkM_1_edt_for_rev_depth.wig", 
                       'A+IP-': "TopoIV_Topo-Seq_experiment\Additional_genome_features\\Gyrase_RifCfx_no_IP_Mu_122mkM_10mkM_1_edt_N3E.wig",
                       'A-IP+': "TopoIV_Topo-Seq_experiment\Additional_genome_features\\Gyrase_no_Cfx_IP_Mu_1_edt_N3E.wig", 
                       'A-IP-': "TopoIV_Topo-Seq_experiment\Additional_genome_features\\Gyrase_no_Cfx_no_IP_Mu_1_edt_N3E.wig",
                       'Tetrade name': 'Gyrase_Replic_1_RifCfx_122mkM_10_mkM'
                       }

#Output folder
Path_for_output=f'{pwd}Genomic_tracks_plots\\'
if not os.path.exists(Path_for_output):
    os.makedirs(Path_for_output)  
    
#Path to TCSs data.
TCSs_coord_path="TopoIV_Topo-Seq_experiment\Additional_genome_features\TopoIV_TCSs_final_set_Cfx_wt_and_gyrA-S83L_N3E.BroadPeak"
GCSs_Cfx_coord_path="TopoIV_Topo-Seq_experiment\Additional_genome_features\Gyrase_Cfx_10mkM_trusted_GCSs_N3E_score.BroadPeak"
GCSs_RifCfx_coord_path="TopoIV_Topo-Seq_experiment\Additional_genome_features\Gyrase_RifCfx_trusted_GCSs_N3E_score.BroadPeak"

#Path to RNA-Seq data.
RNA_seq_data="TopoIV_Topo-Seq_experiment\Additional_genome_features\Sutormin_RNA_Seq_Exponential_av.wig"

Wrapper_func(Tetrade_replic_1_10mkM_gyrS83L, Deletions, Genome, TCSs_coord_path, RNA_seq_data, Path_for_output)
Wrapper_func(Tetrade_Cfx_Gyrase, Deletions, Genome, GCSs_Cfx_coord_path, RNA_seq_data, Path_for_output)
Wrapper_func(Tetrade_RifCfx_Gyrase, Deletions, Genome, GCSs_RifCfx_coord_path, RNA_seq_data, Path_for_output)


print('Script ended its work succesfully!')
