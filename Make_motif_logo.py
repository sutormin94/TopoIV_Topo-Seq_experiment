###############################################
##Dmitry Sutormin, 2022##
##TopoIV Topo-Seq analysis##

####
#Script creates and plots a web-logo for a multiple alignment of DNA sequences.
####

###############################################

#######
#Packages to be imported.
#######

import weblogo
from weblogo import LogoOptions, LogoFormat, eps_formatter, read_seq_data, LogoData, png_formatter, pdf_formatter

Alig_inpath1="Data_analysis\TCSs_analysis\Combined_motif\Combined_motif_originated_sequences_30.fasta" 
Out_path1="Data_analysis\TCSs_analysis\Combined_motif\Combined_motif_originated_sequences_30.pdf"

Alig_inpath2="Data_analysis\TCSs_analysis\Cfx_and_S83L_Cfx\\Unique_GCSs_Combined_motif_originated_sequences_40.fasta" 
Out_path2="Data_analysis\TCSs_analysis\Cfx_and_S83L_Cfx\\Unique_GCSs_TopoIV_Cfx_and_S83L_Cfx_combined_motif_originated_sequences_40.pdf"


#######
#Creates motif logo.
#######

def Create_logo(alig_inpath, out_path):
    MFA_data=open(alig_inpath)
    MFA_seqs=read_seq_data(MFA_data)
    logodata=LogoData.from_seqs(MFA_seqs)
    logooptions=LogoOptions(yaxis_scale=1.8, pad_right=True, stacks_per_line=170)
    logooptions.show_errorbars=False
    logoformat=LogoFormat(logodata, logooptions)
    pdf=weblogo.logo_formatter.pdf_formatter(logodata, logoformat)
    logout=open(out_path, 'wb')
    logout.write(pdf)
    logout.close()
    return

Create_logo(Alig_inpath1, Out_path1)
Create_logo(Alig_inpath2, Out_path2)