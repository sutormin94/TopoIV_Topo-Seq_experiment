###############################################
##Dmitry Sutormin, 2022##
##TopoIV Topo-Seq analysis##

#The script takes raw TCSs data, returns only trusted TCSs, 
#computes GCSs shared between different conditions, 
#draws Venn diagrams of the sets overlappings, 
#writes TCSs sets.
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib.pyplot as plt
import collections
from matplotlib_venn import venn2, venn3, venn3_circles, venn2_circles
import numpy as np

#######
#Variables to be defined.
#######

print('Variables to be defined:')

#Path to the working directory
pwd="Data_analysis\\"

#Input data
path_to_cfx_S83L_replicas={'Cfx_1_10_S83L': pwd + "TCSs_analysis\TCSs_calling_final_0.05_220\Replic_1_Cfx_10_mkM_gyrA_S83L\Replic_1_Cfx_10_mkM_gyrA_S83L_raw_GCSs_called.txt",
                           'Cfx_2_10_S83L': pwd + "TCSs_analysis\TCSs_calling_final_0.05_220\Replic_2_Cfx_10_mkM_gyrA_S83L\Replic_2_Cfx_10_mkM_gyrA_S83L_raw_GCSs_called.txt",
                           'Cfx_3_10_S83L': pwd + "TCSs_analysis\TCSs_calling_final_0.05_220\Replic_3_Cfx_10_mkM_gyrA_S83L\Replic_3_Cfx_10_mkM_gyrA_S83L_raw_GCSs_called.txt"
                           }

path_to_cfx_replicas={'Cfx_1_10': pwd + "TCSs_analysis\TCSs_calling_final_0.05_220\Replic_1_Cfx_10_mkM\Replic_1_Cfx_10_mkM_raw_GCSs_called.txt",
                      'Cfx_2_10': pwd + "TCSs_analysis\TCSs_calling_final_0.05_220\Replic_2_Cfx_10_mkM\Replic_2_Cfx_10_mkM_raw_GCSs_called.txt",
                      'Cfx_3_10': pwd + "TCSs_analysis\TCSs_calling_final_0.05_220\Replic_3_Cfx_10_mkM\Replic_3_Cfx_10_mkM_raw_GCSs_called.txt"
                      }

path_to_all_replicas={'Cfx_1_10_S83L': pwd + "TCSs_analysis\TCSs_calling_final_0.05_220\Replic_1_Cfx_10_mkM_gyrA_S83L\Replic_1_Cfx_10_mkM_gyrA_S83L_raw_GCSs_called.txt",
                      'Cfx_2_10_S83L': pwd + "TCSs_analysis\TCSs_calling_final_0.05_220\Replic_2_Cfx_10_mkM_gyrA_S83L\Replic_2_Cfx_10_mkM_gyrA_S83L_raw_GCSs_called.txt",
                      'Cfx_3_10_S83L': pwd + "TCSs_analysis\TCSs_calling_final_0.05_220\Replic_3_Cfx_10_mkM_gyrA_S83L\Replic_3_Cfx_10_mkM_gyrA_S83L_raw_GCSs_called.txt",      
                      'Cfx_1_10'     : pwd + "TCSs_analysis\TCSs_calling_final_0.05_220\Replic_1_Cfx_10_mkM\Replic_1_Cfx_10_mkM_raw_GCSs_called.txt",
                      'Cfx_2_10'     : pwd + "TCSs_analysis\TCSs_calling_final_0.05_220\Replic_2_Cfx_10_mkM\Replic_2_Cfx_10_mkM_raw_GCSs_called.txt",
                      'Cfx_3_10'     : pwd + "TCSs_analysis\TCSs_calling_final_0.05_220\Replic_3_Cfx_10_mkM\Replic_3_Cfx_10_mkM_raw_GCSs_called.txt"
                      }

path_to_NorfIP_replicas={'NorfIP_1' : pwd + "Sayyed_TopoIV\\NorfIP_1_sites.txt",
                         'NorfIP_2' : pwd + "Sayyed_TopoIV\\NorfIP_2_sites.txt",
                         'NorfIP_3' : pwd + "Sayyed_TopoIV\\NorfIP_3_sites.txt",}

#Configuration of the output for the GCSs data in replicas.
Replicas_path_out="Data_analysis\TCSs_analysis\TCSs_sets\\"
if not os.path.exists(Replicas_path_out):
    os.makedirs(Replicas_path_out)
Cfx_S83L_name="Cfx_S83L"
Cfx_name="Cfx"
All_conditions_name="All_conditions_GCSs"
#Configuration of the output for GCSs trusted.
Cfx_S83L_path=Replicas_path_out + "Cfx_S83L_trusted_GCSs.txt"
Cfx_path=Replicas_path_out + "Cfx_trusted_GCSs.txt"
Cfx_S83L_Cfx_shared_path=Replicas_path_out + "Cfx_S83L_Cfx_shared_trusted_GCSs.txt"
#Outpath for Venn diagrams.
plot_outpath=Replicas_path_out

#Configuration of the output for NorfIP data.
Replicas_path_out_NorfIP=pwd + "Sayyed_TopoIV\\"
if not os.path.exists(Replicas_path_out_NorfIP):
    os.makedirs(Replicas_path_out_NorfIP)
NorfIP_name="NorfIP"
NorfIP_path=Replicas_path_out_NorfIP + "NorfIP_trusted_TCSs.txt"
NorfIP_plot_outpath=Replicas_path_out_NorfIP

#######
#Parsing raw GCSs coordinates, returns dictionary - GCSs_coordinate:N3E.
#######

def read_GCSs_file(GCSs_file_path):
    GCSs_dict={}
    GCSs_in=open(GCSs_file_path, 'r')
    for line in GCSs_in:
        line=line.rstrip().split('\t')
        if line[0] not in ['GCSs_coordinate']:
            GCSs_dict[int(line[0])]=float(line[1])
    GCSs_in.close()
    return GCSs_dict

#######
#Combines replicates into one GCSs table.
#######

def combine_replicates(replicas_dict, path_out, name):
    #Merges a range of replicates
    GCSs_replicas_dict={}
    names_ar=[]
    for key, value in replicas_dict.items(): #Iterates replicas
        names_ar.append(key)
        #Read file with raw GCSs
        Raw_GCSs_dict=read_GCSs_file(value)
        for k, v in Raw_GCSs_dict.items(): #Iterates raw GCSs
            #Table filling process initiation
            if len(names_ar)==1:
                GCSs_replicas_dict[k]=[v]
            #Table filling process continuing (the table already contains at least one GCSs set)
            else:
                #If GCSs is already in the table
                if k in GCSs_replicas_dict:
                    GCSs_replicas_dict[k].append(v)
                #If this is the first occurrence of the element in a NON empty table.
                else:
                    add_el=[]
                    for j in range(len(names_ar)-1):
                        add_el.append(0)
                    add_el.append(v)
                    GCSs_replicas_dict[k]=add_el
        #If table body line contains less elements than header does, hence add zero.
        for k, v in GCSs_replicas_dict.items():
            if len(v)<len(names_ar):
                GCSs_replicas_dict[k].append(0)
    #Sorting the list of dictionary keys.
    GCSs_replicas_dict_sorted=collections.OrderedDict(sorted(GCSs_replicas_dict.items()))
    #Writes merged GCSs data
    fileout=open(path_out + name + '_GCSs_replicates.txt', 'w')
    #Header
    fileout.write('GCSs_coordinate\t')
    for i in names_ar:
        fileout.write(str(i) + '_N3E\t')
    fileout.write('\n')
    #Body of the table
    for k, v in GCSs_replicas_dict_sorted.items():
        fileout.write(str(k) + '\t')
        for i in GCSs_replicas_dict_sorted[k]:
            fileout.write(str(i) + '\t')
        fileout.write('\n')
    fileout.close()
    return GCSs_replicas_dict
 
#Prepares TCSs table for all Topo-Seq conditions
combine_replicates(path_to_all_replicas, Replicas_path_out, All_conditions_name)  
#Prepares TCSs table for all NorfIP conditions
combine_replicates(path_to_NorfIP_replicas, Replicas_path_out_NorfIP, NorfIP_name)  


#######
#Returns only trusted GCSs - observed at least 2 times within 3 biological replicates.
#Data organization: 1. coordinate of GCSs, 2.-4. N3E values for biological replicates 1-3
#######

def trusted(ar):
    av_height=0
    ind=0
    for i in range(len(ar)):
        if ar[i]>0:
            ind=ind+1
            av_height=av_height+ar[i]
    if ind>1:
        return av_height/ind
    else:
        return "No signal"

def trusted_GCSs_calling(GCSs_dictionary):
    ar=[]
    for k, v in GCSs_dictionary.items():
        if trusted(v)!="No signal":
            ar.append([k, trusted(v)])
    return ar

def replicas_comb_trust_wrapper(replicas_dict, path_out, name):
    print('Now working with: ' + str(name))
    cur_GCSs_dict=combine_replicates(replicas_dict, path_out, name)
    cur_GCSs_trusted=trusted_GCSs_calling(cur_GCSs_dict)
    print('Number of trusted GCSs for ' + str(name) + ' : ' + str(len(cur_GCSs_trusted)))
    return cur_GCSs_trusted

#FIdentify trusted TCSs for Topo-Seq replicates.
Cfx_S83L=replicas_comb_trust_wrapper(path_to_cfx_S83L_replicas, Replicas_path_out, Cfx_S83L_name)
Cfx=replicas_comb_trust_wrapper(path_to_cfx_replicas, Replicas_path_out, Cfx_name)
Antibs_GCSs_sets=[Cfx_S83L, Cfx]
#Identify trusted TCSs for NorfIP replicates.
NorfIP_TCSs=replicas_comb_trust_wrapper(path_to_NorfIP_replicas, Replicas_path_out_NorfIP, NorfIP_name)


#######
#GCSs shared between pairs of antibiotics - Cfx and Cfx_S83L.
#######

def pairs_construction(ar1, ar2):
    double=[]
    for i in range(len(ar1)):
        for j in range(len(ar2)):
            if ar1[i][0]==ar2[j][0]:
                double.append([ar1[i][0], ar1[i][1], ar2[j][1]]) #GCSs coordinate, N3E_1, N3E_2 
    return double

Cfx_S83L_Cfx_shared_GCSs=pairs_construction(Cfx_S83L, Cfx)
print('Number of GCSs shared between Cfx_S83L and Cfx: ' + str(len(Cfx_S83L_Cfx_shared_GCSs)) + '\n')


#######
#Parses replicas, overlaps lists of GCSs, output data for Venn diagram construction.
#######

def replicates_parsing_to_list_and_overlapping(replicas_dict):
    #Parsing
    GCSs_dict={}
    for k, v in replicas_dict.items(): #Iterate replicas.
        GCSs_dict[k]=[]
        for c, h in read_GCSs_file(v).items(): #Iterate GCSs.
            GCSs_dict[k].append(c)
    return GCSs_dict

#######
#Venn diagram represents GCSs sets overlapping.
#description2: one, two, one_two
#description3: one, two, one_two, three, one_three, two_three, one_two_three
#######

venn_data_2=[len(Cfx)-len(Cfx_S83L_Cfx_shared_GCSs), len(Cfx_S83L)-len(Cfx_S83L_Cfx_shared_GCSs), len(Cfx_S83L_Cfx_shared_GCSs)]
venn2(subsets = (venn_data_2), set_labels = ("Ciprofloxacin", "Ciprofloxacin S83L"))
venn2_circles(subsets = (venn_data_2))
plt.savefig(plot_outpath+'Cfx_vs_Cfx_S83L_venn.png', dpi=320)
plt.savefig(plot_outpath+'Cfx_vs_Cfx_S83L_venn.svg', dpi=320)
plt.close()

Dict_of_cfx_S83L_replicas=replicates_parsing_to_list_and_overlapping(path_to_cfx_S83L_replicas)
venn3([set(Dict_of_cfx_S83L_replicas['Cfx_1_10_S83L']), set(Dict_of_cfx_S83L_replicas['Cfx_2_10_S83L']), set(Dict_of_cfx_S83L_replicas['Cfx_3_10_S83L'])], set_labels = ('Cfx_S83L_10mkM_1', 'Cfx_S83L_10mkM_2', 'Cfx_S83L_10mkM_3'))
venn3_circles([set(Dict_of_cfx_S83L_replicas['Cfx_1_10_S83L']), set(Dict_of_cfx_S83L_replicas['Cfx_2_10_S83L']), set(Dict_of_cfx_S83L_replicas['Cfx_3_10_S83L'])])
plt.savefig(plot_outpath+'Cfx_S83L_replicas_venn.png', dpi=320)
plt.savefig(plot_outpath+'Cfx_S83L_replicas_venn.svg', dpi=320)
plt.close()

Dict_of_cfx_replicas=replicates_parsing_to_list_and_overlapping(path_to_cfx_replicas)
venn3([set(Dict_of_cfx_replicas['Cfx_1_10']), set(Dict_of_cfx_replicas['Cfx_2_10']), set(Dict_of_cfx_replicas['Cfx_3_10'])], set_labels = ('Cfx_10mkM_1', 'Cfx_10mkM_2', 'Cfx_10mkM_3'))
venn3_circles([set(Dict_of_cfx_replicas['Cfx_1_10']), set(Dict_of_cfx_replicas['Cfx_2_10']), set(Dict_of_cfx_replicas['Cfx_3_10'])])
plt.savefig(plot_outpath+'Cfx_replicas_venn.png', dpi=320)
plt.savefig(plot_outpath+'Cfx_replicas_venn.svg', dpi=320)
plt.close()

Dict_of_NorfIP_replicas=replicates_parsing_to_list_and_overlapping(path_to_NorfIP_replicas)
venn3([set(Dict_of_NorfIP_replicas['NorfIP_1']), set(Dict_of_NorfIP_replicas['NorfIP_2']), set(Dict_of_NorfIP_replicas['NorfIP_3'])], set_labels = ('NorfIP_1', 'NorfIP_2', 'NorfIP_3'))
venn3_circles([set(Dict_of_NorfIP_replicas['NorfIP_1']), set(Dict_of_NorfIP_replicas['NorfIP_2']), set(Dict_of_NorfIP_replicas['NorfIP_3'])])
plt.savefig(NorfIP_plot_outpath+'NorfIP_replicas_venn.png', dpi=320)
plt.savefig(NorfIP_plot_outpath+'NorfIP_replicas_venn.svg', dpi=320)
plt.close()


#######
#GCSs sets average N3E estimation.
#######

def average_height(ar):
    av_he=0
    for i in range(len(ar)):
        peak_he=np.mean(ar[i][1:])
        av_he=av_he+peak_he
    return av_he/len(ar)

print('Cfx average GCSs N3E: '              + str(average_height(Cfx)))
print('Cfx_S83L average GCSs N3E: '         + str(average_height(Cfx_S83L)))
print('Cfx and Cfx_S83L average GCSs N3E: ' + str(average_height(Cfx_S83L_Cfx_shared_GCSs)))


#######
#Write down files with GCSs lists - trusted or shared.
#######

All_GCSs_sets={Cfx_S83L_path: Antibs_GCSs_sets[0], Cfx_path: Antibs_GCSs_sets[1], NorfIP_path: NorfIP_TCSs}

def write_GCSs_file(dictionary):
    for k, v in dictionary.items(): #Iterates lists to be written
        v.sort(key=lambda tup: tup[0])  #Sorting lists by the zero elements of the sublists they consist of 
        fileout=open(k, 'w')
        fileout.write('GCSs_coordinate\tN3E\n')
        for i in range(len(v)):
            fileout.write(str(v[i][0]) + '\t' + str(np.mean(v[i][1:])) + '\n')
        fileout.close()
    return

write_GCSs_file(All_GCSs_sets)


def write_Cfx_S83L_Cfx_shared_GCSs(ar, path):
    fileout=open(path, 'w')
    fileout.write('GCSs_coordinate\tCfx_S83L_N3E\tCfx_N3E\n')
    ar.sort(key=lambda tup: tup[0])
    for i in range(len(ar)):
        fileout.write(str(ar[i][0]) + '\t' + str(ar[i][1]) + '\t' + str(ar[i][2]) + '\n')
    fileout.close()
    return
    
write_Cfx_S83L_Cfx_shared_GCSs(Cfx_S83L_Cfx_shared_GCSs, Cfx_S83L_Cfx_shared_path)
 
print('Script ended its work succesfully!') 
 
