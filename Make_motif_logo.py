import weblogo
from weblogo import LogoOptions, LogoFormat, eps_formatter, read_seq_data, LogoData, png_formatter, pdf_formatter

Alig_inpath="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoIV-Topo-Seq\Data_analysis\GCSs_analysis\Motif\Cfx_sequences_under_TCSs_full.fasta" 

Out_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\TopoIV-Topo-Seq\Data_analysis\GCSs_analysis\Motif\LOGO_Cfx_sequences_under_TCSs_full.pdf"


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

Create_logo(Alig_inpath, Out_path)