###############################################
##Dmitry Sutormin, 2022##
##TopoIV ChIP-Seq analysis##

#Converts coordinates from one reference genome to coordinates for
#another reference genome for a set of regions in the BroadPeak format.
###############################################

#######
#Packages to be imported.
#######

import os
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

#######
#Variables to be defined.
#######

#Path to initial reference genome to extract sequences of regions to be transfered.
Init_ref_genome_path="TopoIV_Topo-Seq_experiment\Additional_genome_features\E_coli_K-12_MG1655_NC_000913.2.fasta"

#Path to file with coordinates to be transfered in a BroadPeak format.
Regions_to_transfer_BroadPeak_path="TopoIV_Topo-Seq_experiment\Additional_genome_features\Sayyed_NorfIP_TCSs_trusted_MG1655.BroadPeak"

#Regions set name.
Regs_name="Sayyed_NorfIP_TCSs_trusted"

#Path to the second reference genome.
Scnd_ref_genome_path="TopoIV_Topo-Seq_experiment\Additional_genome_features\E_coli_w3110_G_Mu.fasta"

#Path to store temp and final files.
Outpath="Data_analysis\Sayyed_TopoIV\\"

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
#BroadPeak-formalized intervals parsing (NAPs sites or BIMEs, etc) and filtering intervals that are not deleted.
#######

def broadpeak_pars(intervals_path):
    filein=open(intervals_path, 'r')
    reg_dict={}       
    for line in filein:
        line=line.rstrip().split('\t')
        int_start=int(line[1])
        int_end=int(line[2])
        reg_name=line[3]   
        reg_dict[reg_name]=[int_start, int_end]        
    print("Number of regions: " + str(len(reg_dict)))
    filein.close()
    return reg_dict

#######
#Extract sequences from a reference sequence with an array of coordinates.
#######

def extract_seq_write(regs_dict, genomefa, outpath, regs_name):
    fileout=open(outpath+regs_name+'.fasta', 'w')
    for reg_name, coord_pair in regs_dict.items():
        int_start=coord_pair[0]
        int_end=coord_pair[1]
        int_seq=genomefa[int_start:int_end]
        fileout.write(f'>{reg_name}\n{int_seq}\n')
    fileout.close()
    return

#######
#Blast sequences extracted from the initial reference sequence in the second reference sequence.
#######

def blast_seq(ref2_genome_path, outpath, regs_name):
    
    #Create BLAST database from the second reference genome.
    Make_ref2_genome_db=NcbimakeblastdbCommandline(dbtype="nucl", input_file=ref2_genome_path)    
    print('Making blast database: ' + str(Make_ref2_genome_db))
    Make_ref2_genome_db()    
    
    #Run BLAST for promoters.
    Regs_blast=NcbiblastnCommandline(query=f'{outpath}{regs_name}.fasta', db=ref2_genome_path, out=f'{outpath}{regs_name}_blast_result.txt', outfmt=6)  
    print('Run blast of the extracted sequences: ' + str(Regs_blast))
    Regs_blast()     

    return


#######
#Read BLAST output.
#######

def read_blast_output(regs_dict, regs_name, outpath):
    
    #Read BLAST results.
    Blast_output=open(f'{outpath}\{regs_name}_blast_result.txt', 'r')
    
    Regions_coords_dict={}
    
    Names_list=[]
    for line in Blast_output:
        line=line.rstrip().split('\t')
        
        if int(line[8]) < int(line[9]):
            start=int(line[8])-1 #Correction for 0-based enumeration of python (BLAST is 1-based)
            end=int(line[9])     #Correction for slice method in python (end of the interval is not included)
            strand="+"
        else:
            start=int(line[9])-1 #Correction for 0-based enumeration of python (BLAST is 1-based)
            end=int(line[8])     #Correction for slice method in python (end of the interval is not included)
            strand="-"
            
        Site_name=line[0]
        Chromosome_ID=line[1]
        Reg_start, Reg_end=regs_dict[Site_name]
        Reg_len=Reg_end-Reg_start+1
        if Site_name in Names_list:
            i+=1
            if float(line[2])>0.95 and int(line[3])>Reg_len*0.95:
                Regions_coords_dict[f'{Site_name}_{i}']=[start, end, strand, Chromosome_ID]
        else:
            i=0
            Names_list.append(Site_name)
            if float(line[2])>0.95 and int(line[3])>Reg_len*0.95:
                Regions_coords_dict[f'{Site_name}_{i}']=[start, end, strand, Chromosome_ID]                    
                                   
    return Regions_coords_dict


#######
#Write BroadPeak format.
#######

def write_BroadPeak(Regions_coords_dict, outpath, regs_name):
    
    fileout=open(f'{outpath}{regs_name}_w3110_G_Mu.BroadPeak', 'w')
    for site_name, site_data in Regions_coords_dict.items():
        fileout.write(f'{site_data[3]}\t{site_data[0]}\t{site_data[1]}\t{site_name}\t10\t.\t1\t-1\t-1\n')
    fileout.close()
    
    return


def Wrapper_function(ref1_genome_path, ref2_genome_path, intervals_path, regs_name, outpath):
    
    genomefa1=read_genome(ref1_genome_path)
    Regs_dict=broadpeak_pars(intervals_path)
    extract_seq_write(Regs_dict, genomefa1, outpath, regs_name)
    blast_seq(ref2_genome_path, outpath, regs_name)
    Regions_coords_dict=read_blast_output(Regs_dict, regs_name, outpath)
    write_BroadPeak(Regions_coords_dict, outpath, regs_name)
    
    return

Wrapper_function(Init_ref_genome_path, Scnd_ref_genome_path, Regions_to_transfer_BroadPeak_path, Regs_name, Outpath)