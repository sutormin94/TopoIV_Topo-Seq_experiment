###############################################
##Dmitry Sutormin, 2022##
##Microscopy data analysis##

#Analysis of cell length and DNA content distributions.
###############################################

#######
#Packages to be imported.
#######

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns 
from scipy.stats import norm, shapiro, normaltest, mannwhitneyu


#######
#Import data.
#######

#Path to raw data table.
Microscopy_data_path="TopoIV_Topo-Seq_experiment\Additional_genome_features\Supplementary_Table_3_microscopy_raw_quantification_data.xlsx"
#Name of a worksheet with cell length data.
WS_name_cell_length="Cell_length"
#Name of a worksheet with DNA content data.
WS_name_DNA_content="DNA_content"


#######
#Fit data.
#######

def fitted_data(data):
    avg=np.mean(data)
    var=np.var(data)
    
    pdf_x=np.linspace(0,10,100)
    pdf_y=1.0/np.sqrt(2*np.pi*var)*np.exp(-0.5*(pdf_x-avg)**2/var)
    return pdf_x, pdf_y


##########
## Read cell length data and plot for E. coli parC-SPA and E. coli parC-SPA gyrA-S83L strains treated and not treated with Cfx.
#########

def plot_cell_length_dist(data_inpath, sheetname, out_plot_path):
    
    cl_data=pd.read_excel(data_inpath, sheet_name=sheetname, header=0, index_col=0)
    
    #Import data.
    wt_no_Cfx=cl_data.loc[:, 'wt-Cfx, mkm'].dropna()
    wt_Cfx=cl_data.loc[:, 'wt+Cfx, mkm'].dropna()
    mut_no_Cfx=cl_data.loc[:, 'gyrA-S83L-Cfx, mkm'].dropna()
    mut_Cfx=cl_data.loc[:, 'gyrA-S83L+Cfx, mkm'].dropna()
    
    #Prepare data.
    Data_dict={'wt -Cfx'        : [wt_no_Cfx,  "#f5d81a", r'$wt -Cfx$'], 
               'wt +Cfx'        : [wt_Cfx,     "#5bff7c", r'$wt +Cfx$'],
               'gyrA-S83L -Cfx' : [mut_no_Cfx, "#878787", r'$gyrA-S83L -Cfx$'], 
               'gyrA-S83L +Cfx' : [mut_Cfx,    "#212DA7", r'$gyrA-S83L +Cfx$']}
    
    #Compare the data statistically.
    #Test for normality.
    print('\nTest if data is normally distributed using normaltest')
    p_value_thr=1e-3
    for dataset_name, dataset_data in Data_dict.items():
        statistics, p_value=normaltest(dataset_data[0])
        if p_value<p_value_thr:
            print(f'{dataset_name}, number of cells={len(dataset_data[0])}, mean length={round(np.mean(dataset_data[0]),2)}um: {statistics}, {p_value} - not normal')
        else:
            print(f'{dataset_name}, number of cells={len(dataset_data[0])}, mean length={round(np.mean(dataset_data[0]),2)}um: {statistics}, {p_value} - normal')
            
    #Compare means. Use non-parametric Mann-Whitney U test. Null hypothesis: sub-populations came from the same general population - means are the same.
    print('\nCompare datasets means by mannwhitneyu')
    for dataset_name_1, dataset_data_1 in Data_dict.items():
        for dataset_name_2, dataset_data_2 in Data_dict.items():
            if dataset_name_1!=dataset_name_2:
                statistics, p_value=mannwhitneyu(dataset_data_1[0], dataset_data_2[0], alternative='two-sided')
                if p_value<p_value_thr:
                    print(f'{dataset_name_1} (mean={round(np.mean(dataset_data_1[0]),2)}um) vs {dataset_name_2} (mean={round(np.mean(dataset_data_2[0]),2)}um): {statistics}, {p_value} - means are different')
                else:
                    print(f'{dataset_name_1} (mean={round(np.mean(dataset_data_1[0]),2)}um) vs {dataset_name_2} (mean={round(np.mean(dataset_data_2[0]),2)}um): {statistics}, {p_value} - means are equal')
    
    #Calculate fraction of long cells.    
    print('\nDataset name\tNumber of long cells\tTotal number of cells\tFraction of long cells')
    wt_mean=np.mean(Data_dict['wt -Cfx'][0])
    for dataset_name, dataset_data in Data_dict.items():
        Number_of_long_cells=dataset_data[0][dataset_data[0]>2*wt_mean].count()
        Total_number_of_cells=len(dataset_data[0])
        Fraction_of_long_cells=float(Number_of_long_cells)/Total_number_of_cells
        print(dataset_name, Number_of_long_cells, Total_number_of_cells, Fraction_of_long_cells)
    
    #Plot data.
    xticks_ar=[1,2,4,6,8,10.1]
    xtickslabels_ar=[1,2,4,6,8,'>10']
    bins_ar=list(np.linspace(0,10.2,50))
    fig, plot_1=plt.subplots(1,1,figsize=(4,2.5), dpi=100)   
    
    for dataset_name, dataset_data in Data_dict.items():
        clipped_array=np.clip(dataset_data[0], bins_ar[0], bins_ar[-1])
        weights=np.ones_like(clipped_array)/(len(clipped_array)) #Taken from https://stackoverflow.com/questions/42481698/probability-density-histogram-with-matplotlib-doesnt-make-sense     
        plot_1.hist(clipped_array, bins=bins_ar, weights=weights, color=dataset_data[1], alpha=0.3, label=f'{dataset_data[2]} ({len(dataset_data[0])})')
        plot_1.hist(clipped_array, bins=bins_ar, weights=weights, histtype=u'step', color=dataset_data[1])
       
    plot_1.axvline(2*wt_mean, color='black', linestyle='--', linewidth=1)
    plot_1.set_xticks(xticks_ar)
    plot_1.set_xticklabels(xtickslabels_ar)
    plot_1.set_xlabel(r'Cell length, $\mu$m', size=12)
    plot_1.set_ylabel('Fraction', size=12)
    plot_1.set_xlim([1,10.5])
    
    plot_1.spines["top"].set_visible(False)
    plot_1.spines["right"].set_visible(False)
    
    plt.legend(fontsize=8, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3)
    plt.tight_layout()
    plt.show()
    
    plt.savefig(out_plot_path+'.png', dpi=300)
    plt.savefig(out_plot_path+'.svg', dpi=300)
    
    return

#Path to output plot.
Plot_outpath_cell_length="Data_analysis\Microscopy\Cells_length_measurement_distribution_wt_vs_gyAS83L_mut"

plot_cell_length_dist(Microscopy_data_path, WS_name_cell_length, Plot_outpath_cell_length)



##########
## Read DNA content data and plot for E. coli parC-SPA and E. coli parC-SPA gyrA-S83L strains treated and not treated with Cfx.
#########

def plot_DNA_content_dist(data_inpath, sheetname, out_plot_path):
    
    cl_data=pd.read_excel(data_inpath, sheet_name=sheetname, header=0, index_col=0)
    Scaling_factor=1000000
    
    #Import data.
    wt_no_Cfx=cl_data.loc[:, 'wt-Cfx'].dropna()/Scaling_factor
    wt_Cfx=cl_data.loc[:, 'wt+Cfx'].dropna()/Scaling_factor
    mut_no_Cfx=cl_data.loc[:, 'gyrA-S83L-Cfx'].dropna()/Scaling_factor
    mut_Cfx=cl_data.loc[:, 'gyrA-S83L+Cfx'].dropna()/Scaling_factor
    
    #Prepare data.
    Data_dict={'wt -Cfx'        : [wt_no_Cfx,  "#f5d81a", r'$wt -Cfx$'], 
               'wt +Cfx'        : [wt_Cfx,     "#5bff7c", r'$wt +Cfx$'],
               'gyrA-S83L -Cfx' : [mut_no_Cfx, "#878787", r'$gyrA-S83L -Cfx$'], 
               'gyrA-S83L +Cfx' : [mut_Cfx,    "#212DA7", r'$gyrA-S83L +Cfx$']}
    
    #Compare the data statistically.
    #Test for normality.
    print('\nTest if data is normally distributed using normaltest')
    p_value_thr=1e-3
    for dataset_name, dataset_data in Data_dict.items():
        statistics, p_value=normaltest(dataset_data[0])
        if p_value<p_value_thr:
            print(f'{dataset_name}, number of cells={len(dataset_data[0])}, mean DNA content={round(np.mean(dataset_data[0]),2)} r.u.: {statistics}, {p_value} - not normal')
        else:
            print(f'{dataset_name}, number of cells={len(dataset_data[0])}, mean DNA content={round(np.mean(dataset_data[0]),2)} r.u.: {statistics}, {p_value} - normal')
            
    #Compare means. Use non-parametric Mann-Whitney U test. Null hypothesis: sub-populations came from the same general population - means are the same.
    print('\nCompare datasets means by mannwhitneyu')
    for dataset_name_1, dataset_data_1 in Data_dict.items():
        for dataset_name_2, dataset_data_2 in Data_dict.items():
            if dataset_name_1!=dataset_name_2:
                statistics, p_value=mannwhitneyu(dataset_data_1[0], dataset_data_2[0], alternative='two-sided')
                if p_value<p_value_thr:
                    print(f'{dataset_name_1} (mean={round(np.mean(dataset_data_1[0]),2)}um) vs {dataset_name_2} (mean={round(np.mean(dataset_data_2[0]),2)} r.u.): {statistics}, {p_value} - means are different')
                else:
                    print(f'{dataset_name_1} (mean={round(np.mean(dataset_data_1[0]),2)}um) vs {dataset_name_2} (mean={round(np.mean(dataset_data_2[0]),2)} r.u.): {statistics}, {p_value} - means are equal')
    
    #Calculate fraction of cells with a high DNA content. 
    print('\nDataset name\tNumber of cells with a high DNA content\tTotal number of cells\tFraction of cells with a high DNA content')
    wt_mean=np.mean(Data_dict['wt -Cfx'][0])
    for dataset_name, dataset_data in Data_dict.items():
        Number_of_long_cells=dataset_data[0][dataset_data[0]>2*wt_mean].count()
        Total_number_of_cells=len(dataset_data[0])
        Fraction_of_long_cells=float(Number_of_long_cells)/Total_number_of_cells
        print(dataset_name, Number_of_long_cells, Total_number_of_cells, Fraction_of_long_cells)
    
    #Plot data.
    xticks_ar=[1,2,4,6,8.1]
    xtickslabels_ar=[1,2,4,6,'>8']
    bins_ar=list(np.linspace(0,8.2,40))
    fig, plot_1=plt.subplots(1,1,figsize=(3.6,2.5), dpi=100)   
    
    for dataset_name, dataset_data in Data_dict.items():
        clipped_array=np.clip(dataset_data[0], bins_ar[0], bins_ar[-1])
        weights=np.ones_like(clipped_array)/(len(clipped_array)) #Taken from https://stackoverflow.com/questions/42481698/probability-density-histogram-with-matplotlib-doesnt-make-sense     
        plot_1.hist(clipped_array, bins=bins_ar, weights=weights, color=dataset_data[1], alpha=0.3, label=f'{dataset_data[2]} ({len(dataset_data[0])})')
        plot_1.hist(clipped_array, bins=bins_ar, weights=weights, histtype=u'step', color=dataset_data[1])
       
    plot_1.axvline(2*wt_mean, color='black', linestyle='--', linewidth=1)
    plot_1.set_xticks(xticks_ar)
    plot_1.set_xticklabels(xtickslabels_ar)
    plot_1.set_xlabel(r'DNA content, relative units', size=12)
    plot_1.set_ylabel('Fraction', size=12)
    plot_1.set_xlim([0,8.5])
    
    plot_1.spines["top"].set_visible(False)
    plot_1.spines["right"].set_visible(False)
    
    plt.legend(fontsize=7.5, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, loc='upper right')
    plt.tight_layout()
    plt.show()
    
    plt.savefig(out_plot_path+'.png', dpi=300)
    plt.savefig(out_plot_path+'.svg', dpi=300)
    
    return

#Path to output plot.
Plot_outpath_DNA_content="Data_analysis\Microscopy\DNA_content_measurement_distribution_wt_vs_gyAS83L_mut"

plot_DNA_content_dist(Microscopy_data_path, WS_name_DNA_content, Plot_outpath_DNA_content)