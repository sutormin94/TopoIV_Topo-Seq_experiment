###############################################
##Dmitry Sutormin, 2022##
##TopoIV Topo-Seq analysis##

#The script takes sets of trusted TCSs as input and plots motifs using the sequences under the TCSs.
#Also it writes sequences and motif to the files.
###############################################

#######
#Packages to be imported.
#######

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from Bio import SeqIO
from Bio import SeqUtils
from Bio.SeqUtils import GC as GC_count
from scipy.stats import binom_test

#######
#Variables to be defined.
#######

print('Variables to be defined:')

#Input data - GCSs, TAB.
path_to_TCSs_files={'Cfx_S83L': "Data_analysis\GCSs_analysis\GCSs_sets\Cfx_S83L_trusted_GCSs.txt",
                    'Cfx'     : "Data_analysis\GCSs_analysis\GCSs_sets\Cfx_trusted_GCSs.txt",
                    }

#Path to the genome FASTA.
Genome_path="TopoIV_Topo-Seq_experiment\Additional_genome_features\E_coli_w3110_G_Mu.fasta"

#Path for the output.
Output_path="Data_analysis\GCSs_analysis\Motif\\"
if not os.path.exists(Output_path):
        os.makedirs(Output_path)

#######
#Trusted GCSs data parsing.
#######

def trusted_GCSs_parsing(input_dict):
        GCSs_sets_dict={}
        for k, v in input_dict.items():
                ar=[]
                filein=open(v, 'r')
                for line in filein:
                        line=line.rstrip().split('\t')
                        if line[0] not in ['GCSs_coordinate']:
                                ar.append(int(line[0]))
                        else:
                                continue
                filein.close()
                GCSs_sets_dict[k]=ar
                print('Number of trusted TCSs for ' + str(k) + ' : ' + str(len(ar)))
        return GCSs_sets_dict

#######
#Genome sequence parsing.
#######

def genome_seq(genome_path):
        genome=open(genome_path, 'r')
        for record in SeqIO.parse(genome, "fasta"):
                genome_sequence=str(record.seq)
        genome.close()
        print('Whole genome average GC: ' + str(SeqUtils.GC(genome_sequence)))
        print('Whole genome length: ' + str(len(genome_sequence)))        
        return genome_sequence

#######
#Returns list of DNA seqs under the GCSs. Seqs have len=win_width.
#Writes sequences under the GCSs to file.
#######

def return_seqs(GCS_coords, win_range, genomefa, filepath_full_len, filepath_6bp_LOGO): 
        fileout=open(filepath_full_len, 'w')
        fileout_6bp_LOGO=open(filepath_6bp_LOGO, 'w')
        seqs=[]
        for i in range(len(GCS_coords)):
                seq=genomefa[int(GCS_coords[i]- win_range[0] - 1):int(GCS_coords[i]+ win_range[1] - 1)]
                seq_6bp_LOGO=genomefa[int(GCS_coords[i]-1):int(GCS_coords[i]-1+6)]
                seqs.append(seq)
                fileout.write('>'+str(GCS_coords[i])+'\n'+str(seq)+'\n')
                fileout_6bp_LOGO.write('>'+str(GCS_coords[i])+'\n'+str(seq_6bp_LOGO)+'\n')
        fileout.close()
        fileout_6bp_LOGO.close()
        print('Number of sequences (TCSs) analysing: ' + str(len(seqs)))
        return seqs

#######
#PFM construction.
#Scans sequences stack by columns, counts the number of particular letters.
#Returns a range of PFMs - "positional frequencies matrixes" .
#######

def make_PFM(seqs_list):
        matrix=[]
        template=seqs_list[0]
        for i in range(len(template)):
                column=[0, 0, 0, 0]
                for j in range(len(seqs_list)):
                        if seqs_list[j][i] == str('A'):
                                column[0] = column[0] + 1
                        elif seqs_list[j][i] == str('T'):
                                column[1] = column[1] + 1
                        elif seqs_list[j][i] == str('G'):
                                column[2] = column[2] + 1
                        elif seqs_list[j][i] == str('C'):
                                column[3] = column[3] + 1
                matrix.append(column)
        #Returns a range of PFMs.
        GC_percent = []
        GT_percent = []
        CT_percent = []
        A_percent = []
        T_percent = []
        G_percent = []
        C_percent = []
        for i in range(len(matrix)):
                GC = float((int(matrix[i][2]) + int(matrix[i][3]))) / (
                        int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
                GT = float((int(matrix[i][1]) + int(matrix[i][2]))) / (
                        int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
                CT = float((int(matrix[i][1]) + int(matrix[i][3]))) / (
                        int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
                A = float((int(matrix[i][0]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
                T = float((int(matrix[i][1]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
                G = float((int(matrix[i][2]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
                C = float((int(matrix[i][3]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
                GC_percent.append(GC)
                GT_percent.append(GT)
                CT_percent.append(CT)
                A_percent.append(A)
                T_percent.append(T)
                G_percent.append(G)
                C_percent.append(C)
        return {'Num_seqs': len(seqs_list), 'A': A_percent, 'T': T_percent, 'G': G_percent, 'C': C_percent, 'CT': CT_percent, 'GT': GT_percent, 'GC': GC_percent}

#######
#Writes PFM data to file.
#######

def write_motif(ar, filepath, coord_shift):
        fileout=open(filepath, 'w')
        fileout.write("#X\tY\n")
        for i in range(len(ar)):
                fileout.write(str((-coord_shift/2)+1+i) + '\t' + str(ar[i])+'\n')
        fileout.close()
        return

#######
#Plotting the motif.
#Matrix type - type of the PFM to plot.
#######

def Plotting(PFMs_set, title, matrix_type, write_out, win_width):
        x_axis=[]
        for i in range(len(PFMs_set['Cfx'])):
                x_axis.append(-(win_width/2)+1+i)      
        print('len(x_axis)=' + str(len(x_axis)))
        ax_range = [-win_width/2, win_width/2, 0.3, 1.0]
        plt.figure(dpi=100, figsize=(16, 6))
        plt.suptitle(str(title), fontsize=20)
        plot1 = plt.subplot()
        plot1.set_xticks([0], minor=True)
        plot1.xaxis.grid(True, which='minor', linewidth=0.5, linestyle='--', alpha=1)            
        #Cfx_S83L
        plot1.plot(x_axis, PFMs_set['Cfx_S83L'], color='#7FCE79', linewidth=4, alpha=0.6)
        plot1.plot(x_axis, PFMs_set['Cfx_S83L'], color='#454F24', linewidth=1, alpha=0.6)
        plot1.plot(x_axis, PFMs_set['Cfx_S83L'], 'o', fillstyle='none', color='#7FCE79', markeredgecolor='#454F24', markersize=2, alpha=0.6)        
        #Cfx
        plot1.plot(x_axis, PFMs_set['Cfx'], color='#ff878b', linewidth=4, alpha=0.6)
        plot1.plot(x_axis, PFMs_set['Cfx'], color='#7D212B', linewidth=1, alpha=0.6)
        plot1.plot(x_axis, PFMs_set['Cfx'], 'o', fillstyle='none', color='#ff878b', markeredgecolor='#4D590D', markersize=2, alpha=0.6)       
        #Tracks annotation
        plot1.annotate('Cfx S83L', xytext=(-75, 0.8), xy=(40, 0.85), color='#7FCE79', weight="bold", size=15)
        plot1.annotate('Cfx', xytext=(-75, 0.75), xy=(40, 0.85), color='#ff878b', weight="bold", size=15)        
        #Set axis parameters
        plot1.tick_params(axis='both', direction='in', bottom='on', top='on', left='on', right='on')
        plot1.axis(ax_range)
        plot1.set_xlim(-win_width/2, win_width/2)
        plot1.set_xticks(np.concatenate((np.arange(-(win_width/2)+5, (win_width/2)+2, 10), [0, 3, -63, -17, 20, 66])), minor=False)
        plot1.tick_params(axis='both', which='major', labelsize=14)
        plot1.set_xlabel('Position, nt', size=17)
        plot1.set_ylabel(str(matrix_type)+'%', size=17)
        #plt.show()
        plt.savefig(write_out, dpi=400, figsize=(16, 6)) 
        plt.close()
        return

#######
#Plotting the motif with statistic.
#Matrix type - type of the PFM to plot.
#######

def Plotting_stat(GC_PFM, num_seq, title, matrix_type, genome_sequence, write_out, win_width):
        #GC statistics module
        #Counts average GC% over the whole genome
        GC_genome=GC_count(genome_sequence)/100
        print('GC% of the reference genome: ' + str(GC_genome))
        
        #Counts GC% p-value in the particular pwm column.
        #Returns p-value array and auxiliary Zero array for plotting.
        alignment_thick=num_seq
        pValue=[]
        Zero=[]
        for i in range(len(GC_PFM)):
                pValue.append(binom_test(float(GC_PFM[i]) * alignment_thick, n=alignment_thick, p=GC_genome))
                Zero.append(1)
        #Plotting   
        x_axis=[]
        for i in range(len(GC_PFM)):
                x_axis.append(-(win_width/2)+1+i)      
        print('len(x_axis)=' + str(len(x_axis)))
        ax_range = [-win_width/2, win_width/2, 0.3, 1.0]
        plt.figure(dpi=100, figsize=(16, 6))
        plt.suptitle(str(title), fontsize=20)
        plot1 = plt.subplot()
        plot1.set_xticks([0], minor=True)
        plot1.xaxis.grid(True, which='minor', linewidth=0.5, linestyle='--', alpha=1)            
        #GC% pwm plotting
        plot1.plot(x_axis, GC_PFM, color='green', linewidth=1)
        plot1.plot(x_axis, GC_PFM, 'go', markersize=3)
        plot1.axis(ax_range) 
        plot1.annotate(matrix_type+'%', xytext=(65, 0.65), xy=(40, 0.85), color='green', weight="bold", size=15)
        txt=plot1.annotate('p-value', xytext=(80, 0.60), xy=(-105, 0.64), color='cyan', weight="bold", size=15)
        txt.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='black')])   
        plot1.set_xlabel('Position, nt', size=17)
        plot1.set_ylabel(matrix_type+'%', size=17)                 
        #Set axis parameters
        plot1.tick_params(axis='both', direction='in', bottom='on', top='on', left='on', right='on')
        plot1.set_xlim(-win_width/2, win_width/2)
        plot1.set_xticks(np.concatenate((np.arange(-(win_width/2)+5, (win_width/2)+2, 10), [0, 3, -63, -17, 20, 66])))
        #p-value plotting
        plot2=plot1.twinx()
        plot2.plot(x_axis, pValue, 'k', linewidth=0.5, alpha=0.6)
        plot2.fill_between(x_axis, pValue, Zero, color='cyan', alpha=0.2)
        plot2.set_yticks(np.arange(0, 1.01, 0.01), minor=False)
        plot2.set_yscale('log')
        plot2.set_yticks([0.005], minor=True)
        plot2.yaxis.grid(True, which='minor', linewidth=1, linestyle='--', alpha=0.8)
        plot2.annotate('Confidence level = 0.005', xytext=(45, 0.0025), xy=(40, 0.8), color='black', size=15)
        plot2.set_ylim(0.0000001, 1.0)
        plot2.set_xlim(-win_width/2, win_width/2)
        plot2.set_ylabel('p-value, logarithmic scale', size=17)        
        #plt.show()
        plt.savefig(write_out, dpi=400, figsize=(16, 6)) 
        plt.close()
        return

#######
#Wraps all the functions together.
#######

def wrap_function(TCSs_input, genome_input_path, output_path):
        win_width=170
        win_range=[(win_width/2)-2, (win_width/2)+2]
        PFM_type='GC'
        plot_title='TopoIV motifs obtained for different Topo-Seq conditions'
        TCSs_dict=trusted_GCSs_parsing(TCSs_input)
        genome_sequence=genome_seq(genome_input_path)
        dict_of_PFMs={}
        for k, v in TCSs_dict.items():
                sequences_list=return_seqs(v, win_range, genome_sequence, output_path+str(k)+'_sequences_under_TCSs_full.fasta', output_path+str(k)+'_sequences_under_TCSs_6bp_LOGO.fasta')
                PFMs=make_PFM(sequences_list)
                write_motif(PFMs[PFM_type], output_path+str(k)+'_GC_pfm.txt', win_width)
                dict_of_PFMs[k]=PFMs[PFM_type]
                Plotting_stat(PFMs[PFM_type], PFMs['Num_seqs'], 'TopoIV motif statistic for '+k, PFM_type, genome_sequence, output_path+PFM_type+'_topoIV_motif_trusted_TCSs_statistic_'+str(k)+'.png', win_width)
        Plotting(dict_of_PFMs, plot_title, PFM_type, output_path+'TopoIV_motif_trusted_TCSs_Cfx'+str(PFM_type)+'.png', win_width)
        Plotting(dict_of_PFMs, plot_title, PFM_type, output_path+'TopoIV_motif_trusted_TCSs_Cfx'+str(PFM_type)+'.svg', win_width)
        return

wrap_function(path_to_TCSs_files, Genome_path, Output_path)

print('Script ended its work succesfully!') 