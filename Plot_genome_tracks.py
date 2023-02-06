###############################################
##Dmitry Sutormin, 2022##
##TopoIV Topo-Seq analysis##

#The script takes diades of WIG files contain N3E or N5E values: A+IP+, A-IP+.
#It smoothes A-IP+ track and divides the A+IP+ track by it.
#Once obtains A+IP+_div, the script plots signal over the genomic positions.
#Alongside N3E genomic tracks, density of TCSs and GCSs is plotted, matS sites locations, transcription levels, and enrichment of MatP.
#A plot is generated: signal intensity over the genome for treated sample alongside with TCSs and GCSs density, 
#MatP enrichment, matS-sites locations, and transcription levels.


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

def norm_smooth_devide(ex_file_path, un_ex_file_path, deletions):
    #WIG parsing, total NE counting (for further normalization on reads number)
    treated_experiment=wig_parsing(ex_file_path) #+A+IP
    untreated_experiment=wig_parsing(un_ex_file_path) #-A+IP
    #Normalization on the coverage of neutral region - NR (does not influenced by replication). 
    #NR is set to be 0:2500000 bp.
    #Adds pseudocounts to avoid zero values
    te=np.mean(treated_experiment[0][0:2500000])
    ue=np.mean(untreated_experiment[0][0:2500000])
    norm_array=[te, ue]
    print(norm_array)
    Min_total_NE=min(norm_array)
    print('Min_total_NE: ' + str(Min_total_NE))
    treated_experiment_norm=[1.0 * (x + 1) * Min_total_NE/norm_array[0] for x in treated_experiment[0]] #+A+IP norm
    print('Passed!')
    untreated_experiment_norm=[1.0 * (x + 1) * Min_total_NE/norm_array[1] for x in untreated_experiment[0]] #-A+IP norm
    print('Passed!')
    #Control samples smoothing: A+IP- and A-IP-
    un_experiment_norm_sm=Smoothing(untreated_experiment_norm, deletions) #-A+IP norm sm 
    #Pairwise division: +A+IP/-A+IP and +A-IP/-A-IP
    ends_divide_IP=[] #+A+IP/-A+IP
    for i in range (len(treated_experiment_norm)):
        if treated_experiment_norm[i]!=0 and un_experiment_norm_sm[i]!=0:
            ends_divide_IP.append(treated_experiment_norm[i]/un_experiment_norm_sm[i])
        else:
            ends_divide_IP.append(0)

    return ends_divide_IP


#######
#Mask zero position in genomic tracks.
#######

def mask_track(track_to_mask, deletions):
    #Regions to be masked (e.g. deletions).  
    mask_array=[]
    for k in range(len(track_to_mask)):
        check_in_del=0
        for dl in deletions:
            if dl[1]>=k>=dl[0]:
                mask_array.append(True)
                check_in_del=1
        if check_in_del==0:
            mask_array.append(False)
    track_masked=np.ma.masked_array(track_to_mask, mask=mask_array)
    
    return track_masked


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

def plot_enrichment_signal(fname, IP_nd_ends, MatP_track, RNA_Seq_track, deletions, TCSs_coords, matS_coords, path_out):
    #Some hack to avoid some bug in matplotlib (OverflowError: In draw_path: Exceeded cell block limit)
    #See: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit
    mpl.rcParams['agg.path.chunksize']=10000
    #Regions to be masked (e.g. deletions). 
    IPed=mask_track(IP_nd_ends, deletions)
    MatP_track_masked=mask_track(MatP_track, deletions)
    RNA_Seq_track_masked=mask_track(RNA_Seq_track, deletions)
    
    #Get left-most TCSs positions.
    TSCs_ar=[]
    for TCS in TCSs_coords:
        TCS_left=TCS[0]
        TSCs_ar.append(TCS_left)
        
    # Bin TCSs.
    bin_width=50000
    chrom_len=4647999
    num_bin=chrom_len//bin_width
    hist, bin_edges=np.histogram(TSCs_ar, bins=num_bin+1) #Taken from this discussion: https://stackoverflow.com/questions/38772640/rebinning-a-list-of-numbers-in-python
    hist=list(hist)+[0]
    bin_edges=np.array(list(bin_edges))+bin_width
    print(bin_edges)    
    
    #Bin RNA-Seq data.
    RNA_seq_binned=bin_track(bin_edges, RNA_Seq_track_masked)
        
    #Get left-most matS positions.
    matS_ar=[]
    for matS in matS_coords:
        matS_left=matS[0]
        matS_ar.append(matS_left)

    #Plotting the distribution of the signal around the genome for IPed sample.
    xcoord=np.arange(0,chrom_len)
    fig, plots=plt.subplots(4, sharex=True, figsize=(5, 7)) 
    
    plots[0].plot(xcoord, IPed, '-', label='Topo IV N3E FE', color='black', linewidth=0.25) 
    plots[0].scatter(TSCs_ar, [-12.5]*len(TSCs_ar), s=2, marker='o', c='black', alpha=0.1, label='TCSs') 
    plots[0].set_xlim([800000,2300000])
    plots[0].set_ylim([-25,200])
    plots[0].set_ylabel('Signal enrichment', size=6)
    plots[0].legend(loc='upper right', fontsize=6, frameon=False)
    plots[1].step(bin_edges, hist, color='black')
    plots[1].set_xlim([800000,2300000])
    plots[1].set_ylim([0,100])
    plots[2].step(bin_edges, RNA_seq_binned, color='black') 
    plots[2].set_xlim([800000,2300000])
    plots[2].set_yscale('log')
    plots[3].plot(xcoord, MatP_track_masked, '-', label='MatP FE', color='#00A0E3', linewidth=0.25) 
    plots[3].scatter(matS_ar, [-35]*len(matS_ar), s=2, marker='o', c='#00A0E3', alpha=1, label='matS') 
    plots[3].set_xlim([800000,2300000])
    plots[3].set_ylim([-70,650])
    plots[3].set_xlabel('Genome position, nt', size=6)
    plots[3].set_ylabel('Signal enrichment', size=6)
    plots[3].legend(loc='upper right', fontsize=6, frameon=False)
    plt.show()
    plt.savefig(path_out+fname+'_signal_enrichment.png', dpi=300, figsize=(16, 8))
    plt.savefig(path_out+fname+'_signal_enrichment.svg', dpi=300, figsize=(16, 8))
    plt.savefig(path_out+fname+'_signal_enrichment.pdf', dpi=300, figsize=(16, 8))
    plt.close()
    return

#######
#Wraps all the functions: data normalization, GCSs calling, 
#plotting, writing the data.
#######

def Wrapper_func(tetrade_dictionary, deletions_inpath, genome_path, matP_track_path, rna_seq_data, tcss_coord_path, matS_coord_path, path_out):
    #Parsing deletions
    deletions=deletions_info(deletions_inpath)
    
    #Parsing TCSs coords.
    TCSs_coords=deletions_info(tcss_coord_path)
    
    #Parsing matS coords.
    matS_coords=deletions_info(matS_coord_path)    
    
    #Read MatP data.
    MatP_track=wig_parsing(matP_track_path)[0]  
    
    #Read RNA-Seq data.
    RNA_Seq_track=wig_parsing(rna_seq_data)[0]     
    
    #Define samples within the tetrade.
    Tet_ID=tetrade_dictionary['Tetrade name']
    print('Now we are working with: ' + str(Tet_ID))
    treated_experiment=tetrade_dictionary['A+IP+']
    untreated_experiment=tetrade_dictionary['A-IP+']
    
    #Obtain pairwise divided tracks: A+IP+/A+IP- and A-IP+/A-IP-.
    IP_norm_div=norm_smooth_devide(treated_experiment, untreated_experiment, deletions)
    
    #Plotting the distribution of the signal around the genome for treated and untreated.
    plot_enrichment_signal(Tet_ID, IP_norm_div, MatP_track, RNA_Seq_track, deletions, TCSs_coords, matS_coords, path_out)    

    return



#Paths to the WIG files contain N3E or N5E that forms a tetrade: A+IP+, A+IP-, A-IP+, A-IP-.
#Name of the set (e.g. Cfx, RifCfx, Micro, Oxo and so on).
Tetrade_replic_1_10mkM_gyrS83L={'A+IP+': pwd + "WIG\gyrA-S83L_replic_1_Cfx_IP_N3E.wig", 
                                'A-IP+': pwd + "WIG\gyrA-S83L_replic_1_no_Cfx_IP_N3E.wig", 
                                'Tetrade name': 'Replic_1_Cfx_10_mkM_gyrA_S83L'
                                }

Tetrade_replic_1_Gyrase={'A+IP+': "TopoIV_Topo-Seq_experiment\Additional_genome_features\Gyrase_Cfx_IP_Mu_10mkM_1_edt_N3E.wig", 
                         'A-IP+': "TopoIV_Topo-Seq_experiment\Additional_genome_features\Gyrase_no_Cfx_IP_Mu_1_edt_N3E.wig", 
                         'Tetrade name': 'Gyrase_Replic_1_Cfx_10_mkM'
                         }

#Path to MatP data.
MatP_track_path="TopoIV_Topo-Seq_experiment\Additional_genome_features\MatP_FE_av.wig"

#Path to RNA-Seq data.
RNA_seq_data="TopoIV_Topo-Seq_experiment\Additional_genome_features\Sutormin_RNA_Seq_Exponential_av.wig"

#Path to TCSs data.
TCSs_coord_path="TopoIV_Topo-Seq_experiment\Additional_genome_features\TopoIV_TCSs_final_set_Cfx_wt_and_gyrA-S83L_N3E.BroadPeak"
GCSs_coord_path="TopoIV_Topo-Seq_experiment\Additional_genome_features\Gyrase_Cfx_10mkM_trusted_GCSs_N3E_score.BroadPeak"


#Path to matS data.
matS_coord_path="TopoIV_Topo-Seq_experiment\Additional_genome_features\MatS_sites_Mercier_W3110_Mu_SGS.BroadPeak"

#Output folder
Path_for_output=f'{pwd}Genomic_tracks_plots\\'
if not os.path.exists(Path_for_output):
    os.makedirs(Path_for_output)  
    
Wrapper_func(Tetrade_replic_1_10mkM_gyrS83L, Deletions, Genome, MatP_track_path, RNA_seq_data, TCSs_coord_path, matS_coord_path, Path_for_output)
Wrapper_func(Tetrade_replic_1_Gyrase, Deletions, Genome, MatP_track_path, RNA_seq_data, GCSs_coord_path, matS_coord_path, Path_for_output)


print('Script ended its work succesfully!')
