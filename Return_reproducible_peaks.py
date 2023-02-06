###############################################
##Dmitry Sutormin, 2022##
##TopoIV ChIP-Seq analysis##

#Takes narrowPeaks files with ChIP-Seq peaks identified in different biological replicas 
#(or for different proteins if one wants to investigate the colocalization of signals),
#return narrowPeak only with reproducible (overlapping) peaks.
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
#Threshold for reproducible peaks calling (must not exceed number of replicas).
Threshold=int(3)
#Outpath.
Path_out="Data_analysis\\NorfIP_vs_ChIP_Seq_vs_Topo_Seq\\"
#Output name.
Output_name="TopoIV_Topo_ChIP_Norf"
    
    
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

def Indicate_where_peaks(genome_length, peaks_ar):
    #Create template genome-long array.
    genome_ar=[0]*genome_length    
    for peak in peaks_ar:
        for i in range (peak[1]-peak[0]):
            genome_ar[peak[0]+i]+=1
    #Remove overlappings.
    for i in range(len(genome_ar)):
        if genome_ar[i]>0:
            genome_ar[i]=1
    return genome_ar

#######
#Find reproducible regions in genome-length array. Use for identification of peaks shared between all replicas.
#######  

def Find_rep_peaks(genome_ar, thr):
    peak=0
    rep_peaks_ar=[]
    for i in range(len(genome_ar)):
        
        #Rules how to enter a peak.
        if genome_ar[i]<thr and peak==0: #We are not in peak.
            continue
        elif genome_ar[i]>=thr and peak==0: #We are at left peak border.
            peak=1
            current_peak=[i]
            continue
        
        #Rules how go deeper into the peak.
        elif genome_ar[i]>=thr and peak==1: #We are within a peak.
            continue
        elif genome_ar[i]<thr and peak==1: #We are at the right peak border.
            peak=0
            current_peak.append(i)
            rep_peaks_ar.append(current_peak)
            continue
    return rep_peaks_ar

#######
#Find reproducible and unique regions by pairwise comparison of coordinates.
####### 

def regs_overlap(name_1, reg_ar_1, name_2, reg_ar_2):
    
    Not_overlapped_peak=0
    Overlapped_peaks=0
    for peak_1 in reg_ar_1:
        peak_1_start=peak_1[0]
        peak_1_end=peak_1[1]
        overlap=0
        for peak_2 in reg_ar_2:
            peak_2_start=peak_2[0]
            peak_2_end=peak_2[1] 
            if (peak_2_end>=peak_1_start>=peak_2_start) or (peak_2_end>=peak_1_end>=peak_2_start) or (peak_1_end>=peak_2_start>=peak_1_start) or (peak_1_end>=peak_2_end>=peak_1_start):
                overlap=1
        if overlap==0:
            Not_overlapped_peak+=1
        else:
            Overlapped_peaks+=1
        
    print(f'Total number of {name_1} peaks: {len(reg_ar_1)}') 
    print(f'Total number of {name_2} peaks: {len(reg_ar_2)}')
    print(f'Number of {name_1} peaks overlapping with {name_2} peaks: {Overlapped_peaks}')
    print(f'Number of {name_1} peaks NOT overlapping with {name_2} peaks: {Not_overlapped_peak}')
    
    return

#######
#Count peaks shared between 2 samples. Non-symmetric because different possible overlappings between peaks of two samples.
#######  

def Find_rep_peaks_pairwise(rep_data, keys_list1, keys_list2, peaks_shared_ar):
    peaks_ar1=rep_data[keys_list1] 
    peaks_ar2=rep_data[keys_list2]
    #Number of peaks shared.
    num_peak12_shared=0
    num_peak21_shared=0
    
    #Count the number of peaks of sample 1 shared with sample 2.
    for peak1 in peaks_ar1:
        shared=0
        for shared_peak in peaks_shared_ar:
            shared_peak_med=int((shared_peak[0]+shared_peak[1])/2)
            if peak1[1]>=shared_peak_med>=peak1[0]:
                shared=1
        if shared==1:
            num_peak12_shared+=1
                
    #Count the number of peaks of sample 2 shared with sample 1.
    for peak2 in peaks_ar2:
        shared=0
        for shared_peak in peaks_shared_ar:
            shared_peak_med=int((shared_peak[0]+shared_peak[1])/2)
            if peak2[1]>=shared_peak_med>=peak2[0]:
                shared=1
        if shared==1:
            num_peak21_shared+=1            
    print("Number of peaks of sample ", keys_list1, " shared with peaks of sample ", keys_list2, ": ", num_peak12_shared)
    print("Number of peaks of sample ", keys_list2, " shared with peaks of sample ", keys_list1, ": ", num_peak21_shared)
    return num_peak12_shared, num_peak21_shared 
        
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

def overlap_call(rep_data, thr, genome_length):
    genome_ar_init=np.array([0]*genome_length)
    #Indicates peaks.
    for name, peaks_ar in rep_data.items():
        genome_ar=np.array(Indicate_where_peaks(genome_length, peaks_ar))
        genome_ar_init+=genome_ar
    #Identify reproducible peaks.
    Rep_peaks_array=Find_rep_peaks(genome_ar_init, thr)
    return Rep_peaks_array

#######
#Calculate matrix of Jaccardian distancies from matrix comtains number of shared peaks.
#######  

def num_to_jaccard_dist(matrix_num):
    matrix_jac=np.zeros((len(matrix_num), len(matrix_num)))
    #Square matrix assumed.
    for i in range(len(matrix_num)):
        for j in range(len(matrix_num)):
            Jaccard_dist=matrix_num[i,j]/(matrix_num[i,i] + matrix_num[j,j] - matrix_num[i,j])
            matrix_jac[i,j]=Jaccard_dist     
    return matrix_jac

#######
#Draw heatmap.
#######  

def draw_heatmap(Matrix_of_shared, keys_list, Title, File_name, pics_outpath):
    #Visualize data with heatmap.
    fig=plt.figure(figsize=(4,4), dpi=100)
    ax=fig.add_subplot(111)
    #Based on: https://matplotlib.org/gallery/images_contours_and_fields/image_annotated_heatmap.html#sphx-glr-gallery-images-contours-and-fields-image-annotated-heatmap-py
    im = ax.imshow(Matrix_of_shared)
    ax.set_xticks(np.arange(len(keys_list)))
    ax.set_yticks(np.arange(len(keys_list)))
    ax.set_xticklabels(keys_list)
    ax.set_yticklabels(keys_list)  
    ax.set_ylim(sorted(ax.get_xlim(), reverse=True)) #Solves a bug in matplotlib 3.1.1 discussed here: https://stackoverflow.com/questions/56942670/matplotlib-seaborn-first-and-last-row-cut-in-half-of-heatmap-plot
    for i in range(len(keys_list)):
        for j in range(len(keys_list)):
            if isinstance(Matrix_of_shared[i, j], int):
                text = ax.text(j, i, Matrix_of_shared[i, j], ha="center", va="center", color="w")
            else:
                text = ax.text(j, i, round(Matrix_of_shared[i, j],2), ha="center", va="center", color="w")
    ax.set_title(Title)
    fig.tight_layout()
    plt.show()
    plt.savefig(pics_outpath+File_name, dpi=400, figsize=(4, 4))
    return
    
#######
#Wrapper: takes peaks from different biological replicas,
#Identifies reproducible regions, writes broadPeak file with reproducible peaks.
#######    

def Wrapper(reps_dict, set_name, thr, genome_path, outpath, pathout_name):
    
    #Reads genome fasta.
    genome_length, genome_seq, chrom_name=read_genome(genome_path)
    
    #Reads replicas data.
    rep_data={}
    for name, rep_path in reps_dict.items():
        rep_data[name]=deletions_info(rep_path)
        
    #Count reproduced and unique peaks by euristic pairwise comparison of coordinates.
    for name_1, reg_ar_1 in rep_data.items():
        for name_2, reg_ar_2 in rep_data.items():
            regs_overlap(name_1, reg_ar_1, name_2, reg_ar_2)
    
    #Create template genome-long array, Indicates peaks, Identify reproducible peaks.
    Rep_peaks_array=overlap_call(rep_data, thr, genome_length)
    keys_list=list(rep_data.keys())
    Rep_peaks_array12=overlap_call({keys_list[0] : rep_data[keys_list[0]], keys_list[1] : rep_data[keys_list[1]]}, 2, genome_length)
    if len(keys_list)>=3:
        Rep_peaks_array13=overlap_call({keys_list[0] : rep_data[keys_list[0]], keys_list[2] : rep_data[keys_list[2]]}, 2, genome_length)
        Rep_peaks_array23=overlap_call({keys_list[1] : rep_data[keys_list[1]], keys_list[2] : rep_data[keys_list[2]]}, 2, genome_length)
        if len(keys_list)>=4:
            Rep_peaks_array14=overlap_call({keys_list[0] : rep_data[keys_list[0]], keys_list[3] : rep_data[keys_list[3]]}, 2, genome_length)
            Rep_peaks_array24=overlap_call({keys_list[1] : rep_data[keys_list[1]], keys_list[3] : rep_data[keys_list[3]]}, 2, genome_length)
            Rep_peaks_array34=overlap_call({keys_list[2] : rep_data[keys_list[2]], keys_list[3] : rep_data[keys_list[3]]}, 2, genome_length) 
    
    #Count peaks shared between pairs of samples.
    Num_peaks12_shared, Num_peaks21_shared=Find_rep_peaks_pairwise(rep_data, keys_list[0], keys_list[1], Rep_peaks_array12)
    if len(keys_list)>=3:
        Num_peaks13_shared, Num_peaks31_shared=Find_rep_peaks_pairwise(rep_data, keys_list[0], keys_list[2], Rep_peaks_array13)
        Num_peaks23_shared, Num_peaks32_shared=Find_rep_peaks_pairwise(rep_data, keys_list[1], keys_list[2], Rep_peaks_array23)
        if len(keys_list)>=4:
            Num_peaks14_shared, Num_peaks41_shared=Find_rep_peaks_pairwise(rep_data, keys_list[0], keys_list[3], Rep_peaks_array14)
            Num_peaks24_shared, Num_peaks42_shared=Find_rep_peaks_pairwise(rep_data, keys_list[1], keys_list[3], Rep_peaks_array24) 
            Num_peaks34_shared, Num_peaks43_shared=Find_rep_peaks_pairwise(rep_data, keys_list[2], keys_list[3], Rep_peaks_array34)
    
    #Create 2d array for visualization.
    Matrix_of_shared=np.array([[len(rep_data[keys_list[0]]), Num_peaks12_shared],
                               [Num_peaks21_shared, len(rep_data[keys_list[1]])]])    
    if len(keys_list)==3:
        Matrix_of_shared=np.array([[len(rep_data[keys_list[0]]), Num_peaks12_shared, Num_peaks13_shared],
                                   [Num_peaks21_shared, len(rep_data[keys_list[1]]), Num_peaks23_shared],
                                   [Num_peaks31_shared, Num_peaks32_shared, len(rep_data[keys_list[2]])]])
        
    if len(keys_list)==4:
        Matrix_of_shared=np.array([[len(rep_data[keys_list[0]]), Num_peaks12_shared, Num_peaks13_shared, Num_peaks14_shared],
                                   [Num_peaks21_shared, len(rep_data[keys_list[1]]), Num_peaks23_shared, Num_peaks24_shared],
                                   [Num_peaks31_shared, Num_peaks32_shared, len(rep_data[keys_list[2]]), Num_peaks34_shared],
                                   [Num_peaks41_shared, Num_peaks42_shared, Num_peaks43_shared, len(rep_data[keys_list[3]])]])
     
    #Visualize data with heatmap.
    draw_heatmap(Matrix_of_shared, keys_list, "Number of peaks shared between\n" + set_name + " samples", pathout_name+'_heatmap_num.png', outpath) 
    
    #Compute Jaccardian distances between samples, visualize data with heatmap.
    matrix_jac=num_to_jaccard_dist(Matrix_of_shared)
    draw_heatmap(matrix_jac, keys_list, "Jaccardian distances between\n" + set_name + " samples", pathout_name+'_heatmap_jac_dist.png', outpath)   
    
    #Create Venn diagram represents replicas overlapping.
    plt.figure(figsize=(4,4))   
    keys_list=list(rep_data.keys())
    if len(keys_list)==2:
        venn2(subsets=(len(rep_data[keys_list[0]])-len(Rep_peaks_array), len(rep_data[keys_list[1]])-len(Rep_peaks_array), len(Rep_peaks_array)), set_labels=(keys_list[0], keys_list[1]),
              set_colors=('#8ec7f2', '#ff8692'), alpha = 0.7)
        venn2_circles(subsets=(len(rep_data[keys_list[0]])-len(Rep_peaks_array), len(rep_data[keys_list[1]])-len(Rep_peaks_array), len(Rep_peaks_array)), linestyle='solid')    
        plt.show()
    elif len(keys_list)==3:
        #123-overlap.
        Rep_peaks_array_123=overlap_call(rep_data, 3, genome_length)
        #12-overlap.
        Rep_peaks_array_12=overlap_call({keys_list[0]: rep_data[keys_list[0]], keys_list[1]: rep_data[keys_list[1]]}, 2, genome_length)
        #13-overlap.
        Rep_peaks_array_13=overlap_call({keys_list[0]: rep_data[keys_list[0]], keys_list[2]: rep_data[keys_list[2]]}, 2, genome_length)
        #23-overlap.
        Rep_peaks_array_23=overlap_call({keys_list[1]: rep_data[keys_list[1]], keys_list[2]: rep_data[keys_list[2]]}, 2, genome_length)        
        #Plot.
        Only_1=len(rep_data[keys_list[0]])-len(Rep_peaks_array_12)-len(Rep_peaks_array_13)+len(Rep_peaks_array_123)
        Only_2=len(rep_data[keys_list[1]])-len(Rep_peaks_array_12)-len(Rep_peaks_array_23)+len(Rep_peaks_array_123)
        Only_3=len(rep_data[keys_list[2]])-len(Rep_peaks_array_13)-len(Rep_peaks_array_23)+len(Rep_peaks_array_123)
        Only_12=len(Rep_peaks_array_12)-len(Rep_peaks_array_123)
        Only_13=len(Rep_peaks_array_13)-len(Rep_peaks_array_123)
        Only_23=len(Rep_peaks_array_23)-len(Rep_peaks_array_123)
        Only_123=len(Rep_peaks_array_123)
        print(Only_1, Only_2, Only_12, Only_3, Only_13, Only_23, Only_123)
        venn3(subsets=(Only_1, Only_2, Only_12, Only_3, Only_13, Only_23, Only_123), set_labels=(keys_list[0], keys_list[1], keys_list[2]))
        venn3_circles(subsets=(Only_1, Only_2, Only_12, Only_3, Only_13, Only_23, Only_123), linestyle='solid')    
        plt.show()   
    else:
        print('Impossible to draw Venn diagram. Currently only 2 and 3 sets can be processed.')
    plt.savefig(outpath+pathout_name+'_venn.png', dpi=400, figsize=(4, 4))
    plt.savefig(outpath+pathout_name+'_venn.svg', dpi=400, figsize=(4, 4))
    
    #Write reproducible peaks.
    write_bed(Rep_peaks_array, chrom_name, outpath+pathout_name+'.narrowPeak')
    return
            
Wrapper(Peaks_data, Set_name, Threshold, Genome, Path_out, Output_name)       
    