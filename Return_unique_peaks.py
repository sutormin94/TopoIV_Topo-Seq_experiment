###############################################
##Dmitry Sutormin, 2022##
##TopoIV Topo-Seq analysis##

#Takes narrowPeaks files with ChIP-Seq peaks identified in different biological replicas,
#return narrowPeak containing only unique peaks (non-overlaping with any other peak).
###############################################

#######
#Packages to be imported.
#######

import numpy as np
import Bio
from Bio import SeqIO
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles


#Path to the working directory with NarrowPeak files.
PWD_peaks="TopoIV_Topo-Seq_experiment\Additional_genome_features\\"
Set_name="TopoIV Topo-Seq, ChIP-Seq, and NorfIP"
Peaks_data={'TopoIV_ChIP_Seq' : PWD_peaks + "Sayyed_ChIPSeq_TCSs_trusted_w3110_G_Mu.BroadPeak",
            'TopoIV_Topo_Seq' : PWD_peaks + "TopoIV_TCSs_final_set_Cfx_wt_and_gyrA-S83L_N3E.BroadPeak",             
            'TopoIV_NorfIP'   : PWD_peaks + "Sayyed_NorfIP_TCSs_trusted_w3110_G_Mu.BroadPeak",
            }

#Path to the reference genome (e.g. E_coli_w3110_G_Mu.fasta).
Genome="TopoIV_Topo-Seq_experiment\Additional_genome_features\E_coli_w3110_G_Mu.fasta"
#Outpath.
Path_out="Data_analysis\\NorfIP_vs_ChIP_Seq_vs_Topo_Seq\\"
    
    
#######
#Opens and reads FASTA file with reference genome
#######

def read_genome(genome_path):
    #Opens FASTA file with the reference genome.
    genome=open(genome_path, 'r')
    for record in SeqIO.parse(genome, "fasta"):
        genomefa=str(record.seq)
        genome_id=record.name
    return len(genomefa), genomefa, genome_id

#######
#Opens and reads BED or narrowPeak files.
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
#Indicate where peaks occures by addition of 1 to these positions to genome-length array.
#######

def Indicate_where_peaks(genome_ar, peaks_ar, indicator):
    for peak in peaks_ar:
        for i in range (peak[1]-peak[0]):
            genome_ar[peak[0]+i]+=indicator
    return genome_ar


#######
#Find reproducible regions in genome-length array. Use for identification of peaks shared between all replicas.
#######  

def Find_rep_peaks(genome_ar, thr):
    peak=0
    rep_peaks_ar=[]
    for i in range(len(genome_ar)):
        if genome_ar[i]<thr and peak==0: #We are not in peak.
            continue
        elif genome_ar[i]>=thr and peak==0: #We are at left peak border.
            peak=1
            current_peak=[i]
            continue
        elif genome_ar[i]>=thr and peak==1: #We are within a peak.
            continue
        elif genome_ar[i]<thr and peak==1: #We are at the right peak border.
            peak=0
            current_peak.append(i)
            rep_peaks_ar.append(current_peak)
            continue
    return rep_peaks_ar


#######
#Find unique regions in genome-length array.
####### 

def Find_unique_peaks(genome_ar, thr):
    peak=0
    unique_peaks_ar=[]
    for i in range(len(genome_ar)):
        #Rules how to enter a peak.
        if genome_ar[i]==0 and peak==0: #We are not in peak.
            continue       
        elif genome_ar[i]==thr and peak==0: #We are at left peak border.
            peak=1
            current_peak=[i]
            continue
        elif (genome_ar[i] in [3, 3-thr]) and (peak==0): #We are at left peak border. But the peak is wrong.
            peak=2
            continue 
        
        #Rules how go deeper into the right peak (not wrong one).
        elif genome_ar[i]==thr and peak==1: #We are within a right peak.
            continue
        elif genome_ar[i]==0 and peak==1: #We are at the right border of a right peak.
            peak=0
            current_peak.append(i)
            unique_peaks_ar.append(current_peak)
            continue
        elif genome_ar[i]==3 and peak==1: #We are at the left border of a shared peak.
            peak=2
            continue
        elif genome_ar[i]==(3-thr) and peak==1: #We are at the precise junction of peaks of two different types.
            peak=0
            current_peak.append(i)
            unique_peaks_ar.append(current_peak)
            continue
        
        #Rules how to go deeper into the wrong peak.
        elif (genome_ar[i] in [3, 3-thr]) and (peak==2): #We are within the shared peak or a peak of a wrong type.
            continue
        elif genome_ar[i]==thr and peak==2: #We are within a unique portion of a peak included into the shared peak.
            continue
        elif genome_ar[i]==0 and peak==2: #We are at the shared peak right border - data should be discarded.
            peak=0
            continue                   
    
    return unique_peaks_ar

        
#######
#Write reproducible peaks in broadPeak format.
#######   

def write_bed(rep_peaks_ar, chrom_name, outpath):
    fileout=open(outpath, 'w')
    for i in range(len(rep_peaks_ar)):
        fileout.write(chrom_name+'\t'+str(rep_peaks_ar[i][0])+'\t'+str(rep_peaks_ar[i][1])+'\tPeak_'+str(i)+'\t10\t.\t-1\t-1\t-1\n')
    fileout.close()
    return

#######
#Identifies reproducible peaks with threshold (number of samples in which a peak should be present) given.
#######  

def overlap_call(rep_data, genome_length):
    #Create template genome-long array.
    genome_ar=[0]*genome_length
    #Indicates peaks.
    i=0
    for name, peaks_ar in rep_data.items():
        indicator=2**i
        genome_ar=Indicate_where_peaks(genome_ar, peaks_ar, indicator)
        i+=1
    #Identify reproducible peaks. Threshold: 2**0 + 2**1=3
    Rep_peaks_array=Find_rep_peaks(genome_ar, 3)
    
    #Identify unique peaks for set1. Threshold: 1
    Unique_set1=Find_unique_peaks(genome_ar, 1)
    
    #Identify unique peaks for set2. Threshold: 2
    Unique_set2=Find_unique_peaks(genome_ar, 2)
    
    return Rep_peaks_array, Unique_set1, Unique_set2

    
#######
#Wrapper: takes peaks from different biological replicas,
#Identifies reproducible regions, writes broadPeak file with reproducible peaks.
#######    

def Wrapper(reps_dict, set_name, genome_path, outpath):
    
    #Reads genome fasta.
    genome_length, genome_seq, chrom_name=read_genome(genome_path)
    
    #Reads replicas data.
    rep_data={}
    for name, rep_path in reps_dict.items():
        rep_data[name]=deletions_info(rep_path)
    
    #Create template genome-long array, Indicate peaks, identify reproducible and unique peaks by pairwise comparision of regions.
    keys_list=list(rep_data.keys())
    Rep_peaks_array12, Unique_peaks_1, Unique_peaks_2=overlap_call({keys_list[0] : rep_data[keys_list[0]], keys_list[1] : rep_data[keys_list[1]]}, genome_length)
    Unuque_peaks_ar=[Unique_peaks_1, Unique_peaks_2]
    print(f'Number of {keys_list[0]} peaks: {len(rep_data[keys_list[0]])}')
    print(f'Number of {keys_list[1]} peaks: {len(rep_data[keys_list[1]])}')
    print(f'Number of shared peaks: {len(Rep_peaks_array12)}')
    print(f'Number of {keys_list[0]} unique peaks: {len(Unique_peaks_1)}')
    print(f'Number of {keys_list[1]} unique peaks: {len(Unique_peaks_2)}')
    
    if len(keys_list)>=3:
        Rep_peaks_array13, Unique_peaks_1, Unique_peaks_3=overlap_call({keys_list[0] : rep_data[keys_list[0]], keys_list[2] : rep_data[keys_list[2]]}, genome_length)
        Rep_peaks_array23, Unique_peaks_2, Unique_peaks_3=overlap_call({keys_list[1] : rep_data[keys_list[1]], keys_list[2] : rep_data[keys_list[2]]}, genome_length)
        Unuque_peaks_ar.append(Unique_peaks_3)
        if len(keys_list)>=4:
            Rep_peaks_array14, Unique_peaks_1, Unique_peaks_4=overlap_call({keys_list[0] : rep_data[keys_list[0]], keys_list[3] : rep_data[keys_list[3]]}, genome_length)
            Rep_peaks_array24, Unique_peaks_2, Unique_peaks_4=overlap_call({keys_list[1] : rep_data[keys_list[1]], keys_list[3] : rep_data[keys_list[3]]}, genome_length)
            Rep_peaks_array34, Unique_peaks_3, Unique_peaks_4=overlap_call({keys_list[2] : rep_data[keys_list[2]], keys_list[3] : rep_data[keys_list[3]]}, genome_length)
            Unuque_peaks_ar.append(Unique_peaks_4)
    
          
    #Write unique peaks.
    for i in range(len(keys_list)):
        write_bed(Unuque_peaks_ar[i], chrom_name, f'{outpath}{keys_list[i]}_not_overlapping_peaks.narrowPeak')

    return
            
Wrapper(Peaks_data, Set_name, Genome, Path_out)       
