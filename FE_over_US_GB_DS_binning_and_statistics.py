###############################################
##Dmitry Sutormin, 2022##
##TopoIV Topo-Seq analysis##

#Script computes Fold Enrichment (FE) over upstream (US)
#and downstream (DS) regions of transcription units (TUs)
#and over TUs bodies.
#Script is dedicated for noisy data. To handle the data it performes data binning for smoothing.
#Also it keeps signal data for all TUs to contstruct confidential interval.
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import random
from Bio import SeqIO
from scipy import stats


#Path to the directory with input files.
PWD='Data_analysis\\'
#Path to TUs groups file.
TUs_groups_path="TopoIV_Topo-Seq_experiment\Additional_genome_features\\"
#Path to the input annotation, type of annotation and name of TUs set.       
##1##
Path_to_annotation_1=TUs_groups_path + 'DY330_RNA-Seq_transcripts_representative_EP_del_cor.txt'
Type_of_annot_1='broadPeak'             
Genes_set_name_1='All_TUs_1672' 
##2##                                   
Path_to_annotation_2=TUs_groups_path + 'DY330_RNA-Seq_transcripts_representative_no_rRNA_EP_del_cor_HETU_200.txt'
Type_of_annot_2='broadPeak'             
Genes_set_name_2='HETU_no_rRNA_200'         
##3##                                   
Path_to_annotation_3=TUs_groups_path + 'DY330_RNA-Seq_transcripts_representative_no_rRNA_EP_del_cor_LETU_200.txt'
Type_of_annot_3='broadPeak'             
Genes_set_name_3='LETU_no_rRNA_200' 
##4##                                   
Path_to_annotation_4=TUs_groups_path + 'DY330_RNA-Seq_transcripts_EP_del_cor_rRNA_7.txt'
Type_of_annot_4='broadPeak'             
Genes_set_name_4='rRNA_7'               
##5##                                   
Path_to_annotation_5=TUs_groups_path + 'DY330_RNA-Seq_transcripts_EP_del_cor_tRNA_49.txt'
Type_of_annot_5='broadPeak'             
Genes_set_name_5='tRNA_49'    

#Path to the file with regions to be omitted (e.g. deletions).
Deletions_inpath='TopoIV_Topo-Seq_experiment\Additional_genome_features\\Deletions_w3110_G_Mu_SGS.broadPeak'
#Width of US, DS regions.
Win_width=15000
#Length of GB.
Length=5000
#Bin width.
Bin_width=200

#Dictionary of pathes to input data.
Set_name='TopoIV_MukB_binned'
Path_to_input_files='Data_analysis\FE\\'
Dict_of_wigs_path={'TopoIV_Cfx'                    : Path_to_input_files + 'TopoIV_Cfx_av.wig',
                   'TopoIV_GyrA_S83L_Cfx'          : Path_to_input_files + 'TopoIV_GyrA_S83L_Cfx_av.wig',
                   'MukB'                          : Path_to_input_files + 'Nolivos_MukB_IP_av_FE.wig', 
                   'Gyrase'                        : Path_to_input_files + 'Sutormin_Gyrase_Cfx_10mkM_FE_av.wig'
                   }



#######
#Checks if directory exists and if not creates.
#######

def Dir_check_create(some_path):
    if not os.path.exists(some_path):
        os.makedirs(some_path)    
    return

#Path to the output directory.
Out_path='Data_analysis\Metagene_plots_repr_transcripts\\'

#Output path.
def create_out_dirs(out_path, genes_set_name, track_set_name):
    Dir_check_create(out_path)
    Dir_check_create(f'{out_path}\Figures\Plot_combinations\{track_set_name}\{genes_set_name}')
    Dir_check_create(f'{out_path}\Signal_of_TUs_tab\{track_set_name}\{genes_set_name}')
    Dir_check_create(f'{out_path}\Signal_of_TUs_wig\{track_set_name}\{genes_set_name}')    
    return

create_out_dirs(Out_path, Genes_set_name_1,  Set_name)
create_out_dirs(Out_path, Genes_set_name_2,  Set_name)
create_out_dirs(Out_path, Genes_set_name_3,  Set_name)
create_out_dirs(Out_path, Genes_set_name_4,  Set_name)
create_out_dirs(Out_path, Genes_set_name_5,  Set_name)


#######
#Parses WIG file with N3/5E values.
#Computes a total number of Ends.
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

#######
#Parsing gff file and preparing gene annotation.
#######

def parse_gff_annotation(gff_inpath, deletions_inpath):
    #Parsing deletions
    deletions=deletions_info(deletions_inpath)
    
    filein=open(gff_inpath, 'r')
    genes_annotation={'Gene': {},
                      'rRNA': {},
                      'tRNA': {},
                      'ncRNA': {}
                      }
    data_source={}
    for line in filein:
        line=line.rstrip()
        if line[0]!='#':
            line=line.split('\t')
            #What occurs in the annotation:
            if line[1] not in data_source:
                data_source[line[1]]={line[2]: 1}
            else:
                if line[2] not in data_source[line[1]]:
                    data_source[line[1]][line[2]]=1
                else:
                    data_source[line[1]][line[2]]+=1
            #Classify genes:
            #Protein coding genes.
            if line[1]=='ena' and line[2]=='gene':
                gene_start=int(line[3])
                gene_end=int(line[4])
                gene_strand=line[6]
                annot=line[8].split(';')
                if annot[1][0]=="N":
                    if annot[2].split('=')[1]=='protein_coding':
                        gene_name=annot[1].split('=')[1]
                        genes_annotation['Gene'][gene_name]=[gene_start, gene_end, gene_strand]
            #rRNA genes.
            elif line[1]=='ena' and line[2]=='rRNA_gene':
                gene_start=int(line[3])
                gene_end=int(line[4])
                gene_strand=line[6]
                annot=line[8].split(';')
                if annot[1][0]=="N":
                    if annot[2].split('=')[1]=='rRNA':
                        gene_name=annot[1].split('=')[1]
                        genes_annotation['rRNA'][gene_name]=[gene_start, gene_end, gene_strand] 
            #tRNA genes.
            elif line[1]=='ena' and line[2]=='tRNA_gene':
                gene_start=int(line[3])
                gene_end=int(line[4])
                gene_strand=line[6]
                annot=line[8].split(';')
                if annot[1][0]=="N":
                    if annot[2].split('=')[1]=='tRNA':
                        gene_name=annot[1].split('=')[1]
                        genes_annotation['tRNA'][gene_name]=[gene_start, gene_end, gene_strand]
            #Other non-coding RNAs.
            elif line[1]=='Rfam' and line[2]=='ncRNA_gene':
                gene_start=int(line[3])
                gene_end=int(line[4])
                gene_strand=line[6]
                annot=line[8].split(';')
                if annot[1][0]=="N":
                    if annot[2].split('=')[1]=='ncRNA':
                        gene_name=annot[1].split('=')[1]
                        genes_annotation['ncRNA'][gene_name]=[gene_start, gene_end, gene_strand]
    filein.close()            
    return genes_annotation, data_source


#######
#Reads annotation of particular set of genes .tab BroadPeak-like (determined on a basis of expression level).
#######

def parse_expression_annotation(annot_inpath):
    genes_annotation={}
    filein=open(annot_inpath, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        if line[0] not in ['GeneID', 'OperonID', 'TU_ID']:
            TU_name=line[1].lstrip('"').rstrip(';"')
            TU_start=int(line[2])
            TU_end=int(line[3])
            TU_strand=line[4]
            TU_expression=float(line[5].replace(',','.'))
            genes_annotation[TU_name]=[TU_start, TU_end, TU_strand, TU_expression]
    filein.close()            
    return genes_annotation


#######
#Write .wig file.
#######

def write_wig(ar, fileout_path, name):
    fileout=open(fileout_path, 'w')
    fileout.write(f'track type=wiggle_0 name="{name}" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom=NC_007779.1_w3110_Mu start=1 step=1\n')
    for point in ar:
        fileout.write(f'{point}\n')
    fileout.close()
    return

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
#Mask deleted regions in wig files, return mean and STD for normalization.
#######

def wig_signal_mask_mean_std(wig_data):
    
    #Mask deletions.
    deletions=[[274500, 372148], [793800, 807500], [1199000, 1214000]]
    mask=[0]*len(wig_data)
    for delition in deletions:
        mask[delition[0]:delition[1]]=[1]*(delition[1]-delition[0])  
    
    wig_data_masked=[]
    for i in range(len(mask)):
        if mask[i]==0:
            wig_data_masked.append(wig_data[i])
    
    Signal_mean_genome=np.mean(wig_data_masked)
    Signal_std_genome=np.std(wig_data_masked)  
    
    print(f'Global mean: {Signal_mean_genome}; Global std: {Signal_std_genome}')     
    
    return Signal_mean_genome, Signal_std_genome

#######
#Returns binned tracks.
#######

def Binning(data_ar, bin_width):
    Binned=[]
    #Calculating number of bins.
    N_bins=int(len(data_ar)/bin_width)
    for i in range(N_bins):
        bin_value=np.mean(data_ar[i*bin_width:(i+1)*bin_width])
        Binned.append(bin_value)
    if N_bins!=0:
        if len(data_ar)>((i+1)*bin_width):
            Binned.append(np.mean(data_ar[(i+1)*bin_width:]))
    elif N_bins==0:
        Binned.append(np.mean(data_ar))
    return Binned

#######
#Scale regions (gene bodies) to equal length: make long shorter and short longer.
#######

def scale_gene_body(ar, length):
    scaled=[]
    if len(ar)>length: #array should be shrinked
        #Determines positions to be taken (other positions will be discarded).
        positions_to_take=[]
        while len(positions_to_take)!=length:
            position=random.randint(0,len(ar)-1)
            if position not in positions_to_take:
                positions_to_take.append(position)
            else:
                continue
        positions_to_take.sort()
        for pos in positions_to_take:
            scaled.append(ar[pos])
    elif len(ar)<length:
        #Determine positions to be duplicated (other positions will be discarded).
        scaled=ar
        for i in range(length-len(ar)):
            position=random.randint(0,len(scaled))
            if position==0:
                scaled=scaled[:position+1]+scaled[position:position+1]+scaled[position+1:]
            else:
                scaled=scaled[:position]+scaled[position-1:position]+scaled[position:]        
    elif len(ar)==length:
        scaled=ar

    return scaled

#######
#Write .tab file with FE info for genes US, GB, and DS.
#######

def write_genes_FE(dict1, dict2, dict3, FE_track_name, path_out):
    fileout=open(path_out, 'w')
    fileout.write(f'Gene_name\tStart\tEnd\tStrand\t{FE_track_name}_FE_US\t{FE_track_name}_FE_GB\t{FE_track_name}_FE_DS\n')
    for k, v in dict1.items():
        fileout.write(f'{k}\t{v[1]}\t{v[2]}\t{v[3]}\t{v[0]}\t{dict2[k][0]}\t{dict3[k][0]}\n')
    fileout.close()
    return

#######
#Convert dictionary to array, discard keys.
#######

def dict_to_ar(dictionary):
    ar=[]
    for k,v in dictionary.items():
        ar.append(v[0]) 
    return ar

#######
#Computes standard error of mean.
#######

def compute_standard_error(ar):
    std_err=np.std(ar)/np.sqrt(len(ar))
    return std_err

#######
#Returns FE or Ded FE over the set of genes (US, GB, DS) - for each gene separately. For non-strand-specific data.
#######

def genes_and_FE_nss(gene_annotation, genes_set_name, tracks_group_name, FE_track, FE_track_name, out_path, deletions_inpath, win_width, length, bin_width):
    
    #Parsing deletions
    deletions=deletions_info(deletions_inpath) 

    #Calculate FE over genes.
    gene_US=[]
    gene_DS=[]
    gene_B=[]
    gene_US_mean_dict={}
    gene_DS_mean_dict={}
    gene_B_mean_dict={} 
    
    #Number of genes.
    Num_genes=len(gene_annotation)
    
    for gene_name, gene_info in gene_annotation.items():
        delited=0
        for deletion in deletions:
            if deletion[1]>=gene_info[0]>=deletion[0] or deletion[1]>=gene_info[1]>=deletion[0]:
                delited=1
        if delited==0:
            start=gene_info[0]
            end=gene_info[1]
            glen=len(FE_track)
            if gene_info[2]=='+':
                if start<win_width:
                    gene_US.append(FE_track[glen-(win_width-start):] + FE_track[:start])
                    gene_US_mean_dict[gene_name]=[np.mean(FE_track[glen-(win_width-start):] + FE_track[:start]), start, end, gene_info[2]]                   
                else:
                    gene_US.append(FE_track[start-win_width:start])
                    gene_US_mean_dict[gene_name]=[np.mean(FE_track[start-win_width:start]), start, end, gene_info[2]]                    
                if end+win_width>glen:
                    gene_DS.append(FE_track[end:] + FE_track[:end+win_width-glen])
                    gene_DS_mean_dict[gene_name]=[np.mean(FE_track[end:] + FE_track[:end+win_width-glen]), start, end, gene_info[2]]                    
                else:
                    gene_DS.append(FE_track[end:end+win_width])
                    gene_DS_mean_dict[gene_name]=[np.mean(FE_track[end:end+win_width]), start, end, gene_info[2]]
                
                gene_B.append(FE_track[start:end])
                gene_B_mean_dict[gene_name]=[np.mean(FE_track[start:end]), start, end, gene_info[2]]               
                
            elif gene_info[2]=='-':
                if start<win_width:
                    gene_DS.append(FE_track[:start][::-1] + FE_track[glen-(win_width-start):][::-1])                  
                    gene_DS_mean_dict[gene_name]=[np.mean(FE_track[:start][::-1] + FE_track[glen-(win_width-start):][::-1]), start, end, gene_info[2]]
                else:
                    gene_DS.append(FE_track[start-win_width:start][::-1])                                               
                    gene_DS_mean_dict[gene_name]=[np.mean(FE_track[start-win_width:start][::-1]), start, end, gene_info[2]]                  
                if end+win_width>glen:
                    gene_US.append(FE_track[:end+win_width-glen][::-1] + FE_track[end:][::-1])                         
                    gene_US_mean_dict[gene_name]=[np.mean(FE_track[:end+win_width-glen][::-1] + FE_track[end:][::-1]), start, end, gene_info[2]]               
                else:
                    gene_US.append(FE_track[end:end+win_width][::-1])                                                   
                    gene_US_mean_dict[gene_name]=[np.mean(FE_track[end:end+win_width][::-1]), start, end, gene_info[2]]                                                                                   
                
                gene_B.append(FE_track[start:end][::-1])
                gene_B_mean_dict[gene_name]=[np.mean(FE_track[start:end][::-1]), start, end, gene_info[2]]                                                         
                
    #Data binning.
    print(len(gene_US), len(gene_B), len(gene_DS))
    
    gene_US_binned=[]
    gene_DS_binned=[]
    gene_B_binned=[]   
    
    for i in range(len(gene_US)):
        gene_US_binned.append(Binning(gene_US[i], bin_width))
        gene_DS_binned.append(Binning(gene_DS[i], bin_width))
        gene_B_binned.append(Binning(gene_B[i], bin_width))
    
    #Scale GB length.
    print(f'GB F scaling in progress, it takes some time...')
    length_binned=int(length/bin_width)
    gene_B_binned_sc=[]
    for gene in gene_B_binned:
        gene_B_binned_sc.append(scale_gene_body(gene, length_binned))
       
    #Calculate mean, std.
    win_width_binned=int(win_width/bin_width)
    
    gene_US_binned_mean=[]
    gene_DS_binned_mean=[]
    gene_US_binned_std=[]
    gene_DS_binned_std=[]
  
    for i in range(win_width_binned):
        gene_US_bin_ar=[]
        gene_DS_bin_ar=[]
       
        for j in range(len(gene_US_binned)):
            gene_US_bin_ar.append(gene_US_binned[j][i])
            gene_DS_bin_ar.append(gene_DS_binned[j][i])
        
        gene_US_binned_mean.append(np.mean(gene_US_bin_ar))
        gene_US_binned_std.append(np.std(gene_US_bin_ar))
        gene_DS_binned_mean.append(np.mean(gene_DS_bin_ar))
        gene_DS_binned_std.append(np.std(gene_DS_bin_ar))    
        
    gene_B_binned_mean=[]
    gene_B_binned_std=[]   
    
    for i in range(length_binned):
        gene_B_bin_ar=[]
        
        for j in range(len(gene_B_binned_sc)):
            gene_B_bin_ar.append(gene_B_binned_sc[j][i])
        
        gene_B_binned_mean.append(np.mean(gene_B_bin_ar))
        gene_B_binned_std.append(np.std(gene_B_bin_ar))  

    #Write wig-like file with FE over US, GB, DS.
    print(f'Writing FE over TU, GB, DS...')
    gene_binned_mean=np.concatenate((gene_US_binned_mean, gene_B_binned_mean, gene_DS_binned_mean), axis=None)
    gene_binned_std=np.concatenate((gene_US_binned_std,  gene_B_binned_std,  gene_DS_binned_std),  axis=None)
    write_wig(gene_binned_mean, f'{out_path}\Signal_of_TUs_wig\{tracks_group_name}\{genes_set_name}\Mean_signal_{FE_track_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_bin_width_{bin_width}bp.wig', f'{win_width}_{length}_{bin_width}')
    write_wig(gene_binned_std,  f'{out_path}\Signal_of_TUs_wig\{tracks_group_name}\{genes_set_name}\STD_signal_{FE_track_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_bin_width_{bin_width}bp.wig',  f'{win_width}_{length}_{bin_width}')


    #Plot FE over US, GB, DS. 
    Upper_conf_interval=np.array(gene_binned_mean)+(np.array(gene_binned_std)/np.sqrt(Num_genes))
    Lower_conf_interval=np.array(gene_binned_mean)-(np.array(gene_binned_std)/np.sqrt(Num_genes))   
    
    print(f'Plotting FE over TU, GB, DS...')
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)  
    
    positions_bn=np.arange(-win_width+int(bin_width/2), win_width+length-int(bin_width/2)+1, bin_width)
    
    plot1.plot(positions_bn,  gene_binned_mean, linestyle='-', color='#c44733', linewidth=2.5, alpha=1, label='Coding strand') 
    plot1.plot(positions_bn,  Upper_conf_interval,  linestyle='-', color='#c44733', linewidth=0.8, alpha=0.6)
    plot1.plot(positions_bn,  Lower_conf_interval,  linestyle='-', color='#c44733', linewidth=0.8, alpha=0.6)
    plot1.fill_between(positions_bn, Lower_conf_interval, Upper_conf_interval, facecolor='#dca0ff', alpha=0.4, interpolate=True)       
    
    #plot1.set_ylim(-0.6, 0.5) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    ticks_lables[ticks.index(length)]='TE'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1)
    plot1.set_xticks([0, length], minor='True')
    plot1.set_yticks([1], minor='True') 
    plot1.axhline(1, color='black', linestyle='--', alpha=0.5, linewidth=0.5, zorder=20) 
    plot1.axvline(0, color='black', linestyle='--', alpha=0.5, linewidth=0.5, zorder=21)
    plot1.axvline(length, color='black', linestyle='--', alpha=0.5, linewidth=0.5, zorder=22)       
    plot1.legend(fontsize=12, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel('EcTopoI FE', size=20)    
    plt.savefig(f'{out_path}\Figures\Plot_combinations\{tracks_group_name}\{genes_set_name}\\{FE_track_name}_over_{genes_set_name}_{win_width}bp_bin_width_{bin_width}bp.png', dpi=400, figsize=(10, 6))   
    plt.close()      
  
    return gene_US_binned, gene_B_binned_sc, gene_DS_binned


#######
#Normalizes tracks and divides strand-specific tracks by a non-strand-specific track.
#######

def genes_FE_comb_plot(gene_annotation, genes_set_name, tracks_group_name, tracks_dict, norm_coefs_dict, out_path, win_width, length, bin_width):
    
    #Number of genes.
    Num_genes=len(gene_annotation)
    
    #Normalize and divide tracks, compute mean and std.
    #For US and DS regions.
    win_width_binned=int(win_width/bin_width)
    #For gene body regions.     
    length_binned=int(length/bin_width) 
    
       
    #Normalize tracks, calculate mean, std.
    #For US and DS regions.
    Mean_STD_dict={}
    
    for condition_name, condition_data in tracks_dict.items():
        
        gene_US_binned=condition_data[0]
        gene_DS_binned=condition_data[2]
         
        gene_US_binned_mean=[]
        gene_DS_binned_mean=[]
        gene_US_binned_std=[]
        gene_DS_binned_std=[]        
  
        for i in range(win_width_binned):
            gene_US_bin_ar=[]
            gene_DS_bin_ar=[]        
           
            for j in range(len(gene_US_binned)):
                gene_US_bin_ar.append(((gene_US_binned[j][i]-norm_coefs_dict[condition_name]['Mean'])/norm_coefs_dict[condition_name]['STD']))
                gene_DS_bin_ar.append(((gene_DS_binned[j][i]-norm_coefs_dict[condition_name]['Mean'])/norm_coefs_dict[condition_name]['STD']))            
            
            gene_US_binned_mean.append(np.mean(gene_US_bin_ar))
            gene_US_binned_std.append(np.std(gene_US_bin_ar))    
            gene_DS_binned_mean.append(np.mean(gene_DS_bin_ar))
            gene_DS_binned_std.append(np.std(gene_DS_bin_ar))    
        
        Mean_STD_dict[condition_name]={'US_mean' : gene_US_binned_mean, 'DS_mean' : gene_DS_binned_mean,
                                       'US_STD'  : gene_US_binned_std,  'DS_STD'  : gene_DS_binned_std}
    
    #For gene body regions.    
    for condition_name, condition_data in tracks_dict.items():
        
        gene_B_binned=condition_data[1]

        gene_B_binned_mean=[]
        gene_B_binned_std=[]         
    
        for i in range(length_binned):
            gene_B_bin_ar=[]
            
            for j in range(len(gene_B_binned)):
                gene_B_bin_ar.append(((gene_B_binned[j][i]-norm_coefs_dict[condition_name]['Mean'])/norm_coefs_dict[condition_name]['STD']))
            
            gene_B_binned_mean.append(np.mean(gene_B_bin_ar))
            gene_B_binned_std.append(np.std(gene_B_bin_ar))  
        
        Mean_STD_dict[condition_name]['GB_mean']=gene_B_binned_mean
        Mean_STD_dict[condition_name]['GB_STD']=gene_B_binned_std
  
     
    #Write wig-like file with FE over US, GB, DS for normalized tracks.
    print(f'Writing FE over TU, GB, DS...(normalized)')  
    Mean_STD_dict_concat={}
    for condition_name, condition_data in Mean_STD_dict.items():
        #Concatenate US, GB, DS.
        gene_binned_mean=np.concatenate((condition_data['US_mean'], condition_data['GB_mean'], condition_data['DS_mean']), axis=None)
        gene_binned_std=np.concatenate((condition_data['US_STD'], condition_data['GB_STD'], condition_data['DS_STD']),  axis=None)
        write_wig(gene_binned_mean,   f'{out_path}\Signal_of_TUs_wig\{tracks_group_name}\{genes_set_name}\Mean_normalized_signal_{tracks_group_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_bin_width_{bin_width}bp.wig',   f'{win_width}_{length}_{bin_width}')    
        write_wig(gene_binned_std,    f'{out_path}\Signal_of_TUs_wig\{tracks_group_name}\{genes_set_name}\STD_normalized_signal_{tracks_group_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp_bin_width_{bin_width}bp.wig',    f'{win_width}_{length}_{bin_width}')
        
        #Calculate upper and lower confidential intervals (+/- STE).
        Upper_conf_interval=np.array(gene_binned_mean)+(np.array(gene_binned_std)/np.sqrt(Num_genes))
        Lower_conf_interval=np.array(gene_binned_mean)-(np.array(gene_binned_std)/np.sqrt(Num_genes))         
        Mean_STD_dict_concat[condition_name]={'Mean' : gene_binned_mean, 'STD' : gene_binned_std, 
                                              'Upper_conf' : Upper_conf_interval, 'Lower_conf' : Lower_conf_interval}

    #Plot FE over US, GB, DS. Make metagene plots.
    print(f'Plotting FE over TU, GB, DS...')
    plt.figure(figsize=(7.5, 4.5), dpi=100)
    plot1=plt.subplot(111)  
    
    positions_bn=np.arange(-win_width+int(bin_width/2), win_width+length-int(bin_width/2)+1, bin_width)
    color_ar=['#6a65c7', '#c96458', '#333738', '#B6B8BD', '#d9ac3a', '#31ed57']
    
    i=0
    for condition_name, condition_data in Mean_STD_dict_concat.items():
        plot1.plot(positions_bn,  condition_data['Mean'], linestyle='-', color=color_ar[i], linewidth=1.5, alpha=1, label=f'EcTopoI FE {condition_name} {genes_set_name} ({Num_genes})') 
        plot1.plot(positions_bn,  condition_data['Upper_conf'],  linestyle='-', color=color_ar[i], linewidth=0.8, alpha=0.2)
        plot1.plot(positions_bn,  condition_data['Lower_conf'],  linestyle='-', color=color_ar[i], linewidth=0.8, alpha=0.2)
        plot1.fill_between(positions_bn, condition_data['Lower_conf'], condition_data['Upper_conf'], facecolor=color_ar[i], alpha=0.1, interpolate=True)        
        i+=1
    
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='$TU_{start}$'
    ticks_lables[ticks.index(length)]='$TU_{end}$'
    ticks_lables1=ticks_lables[:ticks_lables.index('$TU_{end}$')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1)
    plot1.tick_params(axis='x', which='major', labelsize=15, pad=1)
    plot1.set_xticks([0, length], minor='True')
    plot1.set_yticks([1], minor='True') 
    plot1.axhline(0, color='black', linestyle='--', alpha=0.5, linewidth=0.5, zorder=20) 
    plot1.axvline(0, color='black', linestyle='--', alpha=0.5, linewidth=0.5, zorder=21)
    plot1.axvline(length, color='black', linestyle='--', alpha=0.5, linewidth=0.5, zorder=22) 
    plot1.tick_params(axis='y', which='major', labelsize=15, pad=2)
    plot1.legend(fontsize=10, frameon=False, markerscale=5, handlelength=0.5)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel('EcTopoI FE, norm', size=20)    
    plt.tight_layout()
    plt.savefig(f'{out_path}\Figures\Plot_combinations\{tracks_group_name}\{genes_set_name}\\{tracks_group_name}_over_{genes_set_name}_{win_width}bp_bin_width_{bin_width}bp.png', dpi=400, figsize=(7.5, 4.5), transparent=True)   
    plt.savefig(f'{out_path}\Figures\Plot_combinations\{tracks_group_name}\{genes_set_name}\\{tracks_group_name}_over_{genes_set_name}_{win_width}bp_bin_width_{bin_width}bp.svg', dpi=400, figsize=(7.5, 4.5), transparent=True)       
    plt.close()    
    
    
    #Compare signal in US, GB, DS. Do statistics. Make barplots. 
    #Get US, DS, TSS, GB signal. Normalize data.
    TSS_halfwidth=int(200/bin_width)
    USDS_length=int(5000/bin_width)
    
    Barplot_data_dict={}
    for condition_name, condition_data in tracks_dict.items():
        gene_US_binned=condition_data[0]
        gene_GB_binned=condition_data[1]
        gene_DS_binned=condition_data[2]        
        
        US_means_ar=[]
        TSS_means_ar=[]
        GB_means_ar=[]
        DS_means_ar=[]   
        
        for i in range(len(gene_US_binned)):
            US_means_ar.append(np.mean(((np.array(gene_US_binned[i][-USDS_length-TSS_halfwidth:-TSS_halfwidth])-norm_coefs_dict[condition_name]['Mean'])/norm_coefs_dict[condition_name]['STD'])))
            TSS_means_ar.append(np.mean(((np.array(gene_US_binned[i][-TSS_halfwidth:]+gene_GB_binned[i][:TSS_halfwidth])-norm_coefs_dict[condition_name]['Mean'])/norm_coefs_dict[condition_name]['STD'])))
            GB_means_ar.append(np.mean(((np.array(gene_GB_binned[i][TSS_halfwidth:])-norm_coefs_dict[condition_name]['Mean'])/norm_coefs_dict[condition_name]['STD']))) 
            DS_means_ar.append(np.mean(((np.array(gene_DS_binned[i][TSS_halfwidth:TSS_halfwidth+USDS_length])-norm_coefs_dict[condition_name]['Mean'])/norm_coefs_dict[condition_name]['STD'])))  
            
        Barplot_data_dict[condition_name]=[US_means_ar, TSS_means_ar, GB_means_ar, DS_means_ar]
    
    #Prepare bar coordinates.    
    X_coords=[]
    X_coords_m=[]
    for i in range(4):
        group_coords=np.array(np.arange(1,len(tracks_dict)+1,1).tolist())
        X_coords+=list(group_coords+(i*(len(tracks_dict)+1)))
        X_coords_m.append(np.mean(group_coords)+(i*(len(tracks_dict)+1)))
    print(X_coords)
    print(X_coords_m)
    
    #Prepare data: mean, STE.
    Data_ar=[]
    for i in range(4):
        for condition_name, condition_data in Barplot_data_dict.items():
            Data_ar.append(condition_data[i])
            
    Mean_data_ar=[np.mean(x) for x in Data_ar]
    StEr_data_ar=[compute_standard_error(x) for x in Data_ar]
    
    #Prepare colors.
    color_list_base=['#6a65c7', '#c96458', '#333738', '#B6B8BD', '#d9ac3a']
    color_list=[]
    for i in range(4):
        color_list+=list(color_list_base[:len(tracks_dict)])
        
    #Prepare xticks labels.
    Regions=['US', 'TSS', 'TU body', 'DS']
        
    #Compare US, TSS, TUB.
    fig, plot_av=plt.subplots(1,1,figsize=(3,3), dpi=100)
      
    Bars=plot_av.bar(X_coords, Mean_data_ar, yerr=StEr_data_ar, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.9, color=color_list, edgecolor='k', linewidth=1)
    plot_av.set_ylabel('EcTopoI FE, norm', size=16)
    plot_av.set_xticks(X_coords_m)
    plot_av.set_xticklabels(Regions, rotation=0, size=11)     
    plt.legend(Bars[:len(tracks_dict)], list(tracks_dict.keys()), frameon=False, loc='best', markerscale=1, handlelength=0.7, handletextpad=0.3)  
    plt.tight_layout()
    plt.savefig(f'{out_path}\Figures\Plot_combinations\{tracks_group_name}\{genes_set_name}\\US_GB_DS_mean_{tracks_group_name}_over_{genes_set_name}_US_DS_width_{int(USDS_length*bin_width)}bp_TSS_width_{int(TSS_halfwidth*bin_width*2)}bp_bin_width_{bin_width}bp.svg', dpi=300, size=(3,3), transparent=True) 
    plt.savefig(f'{out_path}\Figures\Plot_combinations\{tracks_group_name}\{genes_set_name}\\US_GB_DS_mean_{tracks_group_name}_over_{genes_set_name}_US_DS_width_{int(USDS_length*bin_width)}bp_TSS_width_{int(TSS_halfwidth*bin_width*2)}bp_bin_width_{bin_width}bp.png', dpi=300, size=(3,3), transparent=True)   
    
    plt.close() 
    
    #Welch t-test. Compare signal for different regions.
    for i in range(4):
        for j in range(len(tracks_dict)):
            for k in range(len(tracks_dict)):
                if k>j:
                    Intervals_stat=stats.ttest_ind(Data_ar[(i*len(tracks_dict))+j], Data_ar[(i*len(tracks_dict))+k], equal_var=False)
                    print(f'Test difference between: {list(tracks_dict.keys())[j]} {Regions[i]} and {list(tracks_dict.keys())[k]} {Regions[i]}')
                    print(f'Sample size: {len(Data_ar[(i*len(tracks_dict))+j])}, Sample size: {len(Data_ar[(i*len(tracks_dict))+k])}')
                    print(f'\nT-test FE Mean1={round(np.mean(Data_ar[(i*len(tracks_dict))+j]),3)}; Mean2={round(np.mean(Data_ar[(i*len(tracks_dict))+k]),3)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n')    

    return

#######
#Wrapper: reads data and gene annotation, computes signal over US, GB, DS, plots and writes the output.
#######


def Wrapper_signal_over_TUs(dict_of_wigs_path, path_to_annotation, type_of_annot, genes_set_name, set_name, deletions_inpath, win_width, length, bin_width, out_path):
    
    #Reads annotation.
    print(f'Now working with {path_to_annotation}')
    if type_of_annot=='gff':
        gene_annotation=parse_gff_annotation(path_to_annotation, deletions_inpath)[0]['Gene']
    elif type_of_annot=='broadPeak':
        gene_annotation=parse_expression_annotation(path_to_annotation)    
              
    #Reads input data in wig files (not strand-specific).
    #Calculate normalization coefficients.
    dict_of_wigs={}
    norm_coefs_dict={}
    for wig_name, wig_data in dict_of_wigs_path.items():
        wig_ar=wig_parsing(wig_data)
        dict_of_wigs[wig_name]=wig_ar
        Global_mean, Global_STD=wig_signal_mask_mean_std(wig_ar)
        norm_coefs_dict[wig_name]={'Mean' : Global_mean, 'STD' : Global_STD}
        
    #Calculate and plot signal over TUs (not strand-specific).
    binned_data_dict={}
    for FE_track_name, FE_track in dict_of_wigs.items():
        binned_data_dict[FE_track_name]=genes_and_FE_nss(gene_annotation, genes_set_name, set_name, FE_track, FE_track_name, out_path, deletions_inpath, win_width, length, bin_width)
    
    #Normalize tracks, divide strand-specific tracks by not strand-specific.
    genes_FE_comb_plot(gene_annotation, genes_set_name, set_name, binned_data_dict, norm_coefs_dict, out_path, win_width, length, bin_width)
    
    return

Wrapper_signal_over_TUs(Dict_of_wigs_path, Path_to_annotation_1,  Type_of_annot_1,  Genes_set_name_1,  Set_name, Deletions_inpath, Win_width, Length, Bin_width, Out_path)
Wrapper_signal_over_TUs(Dict_of_wigs_path, Path_to_annotation_2,  Type_of_annot_2,  Genes_set_name_2,  Set_name, Deletions_inpath, Win_width, Length, Bin_width, Out_path)
Wrapper_signal_over_TUs(Dict_of_wigs_path, Path_to_annotation_3,  Type_of_annot_3,  Genes_set_name_3,  Set_name, Deletions_inpath, Win_width, Length, Bin_width, Out_path)
Wrapper_signal_over_TUs(Dict_of_wigs_path, Path_to_annotation_4,  Type_of_annot_4,  Genes_set_name_4,  Set_name, Deletions_inpath, Win_width, Length, Bin_width, Out_path)
Wrapper_signal_over_TUs(Dict_of_wigs_path, Path_to_annotation_5,  Type_of_annot_5,  Genes_set_name_5,  Set_name, Deletions_inpath, Win_width, Length, Bin_width, Out_path)
