# TopoIV_Topo-Seq_experiment

Exploring Topoisomerase IV cleavage sites across *E. coli W3110* genome

Topoisomerase IV (Topo IV) is a type II topoisomerase that unlinks catenated DNA molecules in bacterial cells and also relaxes positive supercoiling.
It operates with two DNA segments - so called G and T. During catalysis the enzyme introduces temporal double-stranded break into G-segment, transfers T-segment through it and religates the gap. 
5'-ends of the DNA break are stabilized by formation of an intermediate covalent complex between DNA and Topo IV.
Topo-Seq is a ChIP-Seq-like approach that exploits formation of these intermediates to map the Topo IV cleavage sites (TCSs) with a single-base precision.

This repository contains a set of bash and python scripts which were used for Topo-Seq data analysis and visualization. 
Raw sequencing data and some processed files (coverage depth WIG, N3E WIG, GCSs lists) can be retrieved from GEO datasets with accession GSE224562.


## Raw_reads_processing.sh

Shell script that makes QC of the reads before and after the trimming procedure. 
Than script maps trimmed and paired reads to the reference genome, prepares sorted and 
indexed BAM-files suitable for visualization with IGV.

**Requirements:** factqc, trimmomatic, bwa mem, samtools, shell

**Input:** Raw reads files (FASTQ), Genome file (FASTA)

**Output:** FastQC reports, SAM files, sorted and indexed BAM files

######################

## SAM_to_coverage_and_N5E_N3E.py

Script takes SAM files as an input, performs QC filtering of reads relying on the alignment quality and a presence of the partner: 
only reads pairs that have a score<256 are stored. Than the script computes coverage depth for DNA chains separately and for both. 
Additionally it calculates N5E (number of DNA fragments starts) and N3E (number of DNA fragments ends) values for every genome position. 
Coverage depth, N3E and N5E info returns as WIG files.

**Requirements:** python 3

**Input:** SAM files, chromosome identificator (for output WIG)

**Output:** SAM files contain proper aligned reads only, TAB files with start coordinates of DNA fragments alignments and alignment lengths, a range of WIG files:
coverage depth for forward strand, reverse strand, and for both strands, N3E, N5E, and N3E+N5E data 

######################

## TCSs_calling_topoIV.py

The script takes WIG files tetrade that contain N3E values: A+IP+, A+IP-, A-IP+, A-IP-. It smooths A+IP- and A-IP- tracks and 
divides A+IP+ and A-IP+ by them. Once A+IP+_div and A-IP+_div are obtained the script performs Audic-Clavery statistic test (Audic & Claverie, 1997) 
and returns regions of A+IP+_div where i and i+5 positions are significantly higher than corresponding in A-IP+. 
These regions from now are called TCSs. GCSs are stored in the output TXT file. Also two plots are generated: 1) 
signal coverage over the genome for treated and untreated samples; 2) Motif expected to be under the TCSs.

**Requirements:** python 3

**Input:** TAB file with deletions coordinates, N3E or N5E values containing WIG files forms a tetrade (A+IP+, A+IP-, A-IP+, A-IP-), FASTA genome file

**Output:** Plot shows coverage depth normalization before TCSs calling, TAB file with raw TCSs coordinates and N3E values, Plot with raw motif

######################

## TCSs_filtering_and_overlapping_topoIV.py

The script takes raw TCSs data, returns only trusted TCSs (found in at least two biological replicates), finds TCSs shared between different conditions, 
draws Venn diagrams of the sets overlappings, writes TCSs sets.

**Requirements:** python 3

**Input:** TAB files with raw TCSs (triplicates)

**Output:** TAB files with trusted TCSs, TAB files with shared TCSs, 
Venn diagrams

######################

## Motifs_visualization_sequences_extraction.py

The script takes sets of trusted TCSs as input and plots motifs using the sequences under the TCSs.
Also it writes sequences and motif to files.

**Requirements:** python 3

**Input:** TAB files with trusted GCSs, FASTA genome file

**Output:** FASTA file with sequences under TCSs, TAB file with motif PFM (GC degenerate for further Fourier analysis), Topo IV motif plots constructed for trusted TCSs

######################

## Combined_motif_construction_scanning_plotting_writing.py

The script takes sets of trusted TCSs as input, filters TCSs with highest N3E, makes a combined set consists of these TCSs, returns sequences under them and constructs
PSSM matrix by the way getting rid of antibiotic-specific bias at positions forming the cleavage site. Than the script scans a sequence of interest with the PSSM, 
returns the results of scanning, plots combined motif and writes it in a GC% degenerate and in a non-degenerate forms.

**Requirements:** python 3

**Input:** TAB files with trusted TCSs, FASTA genome file, FASTA file with a sequence to scan using the motif constructed, 
dataset name and chromosome identificator (both for output WIG)

**Output:** WIG files with the scanning procedure results (score track), FASTA file with sequences under TCSs that were selected for combined motif construction,
Plot of the combined gyrase motif, TAB file with motif PFM (GC degenerate for further Fourier analysis) vertically oriented, 
TAB file with motif non-degenerate and corrected PFM horizontally oriented.

######################

## Return_GCSs_score_height_correlation.py

The script takes results of scanning procedure (WIG file) and 
returns score for GCSs (writes TAB files contain coordinate\tN3E\tScore info), 
computes Pearson correlation between N3E and score, plots (Score, N3E) scatter plots.

**Requirements:** python 3

**Input:** TAB files with trusted GCSs, WIG file with the scanning procedure results (score track).

**Output:** TAB files with trusted GCSs (score info added),
plot with distributions of N3E and score values for different GCSs sets and overall genome (the last is for score),
(score, N3E) scatter plot with 1d trend line, correlation (score, N3E)

######################

## TopoIV_combined_motif_vs_Gyrase_combined_motif.py

The script takes several PFMs (e.g. for Gyrase and Topo IV) and plots them together. 

**Requirements:** python 3

**Input:** TAB files with PFM data in the format: position, frequency of a GC pair.

**Output:** Plot with motifs.

######################

## Genome_intervals_analysis.py

The script analyzes sets of genome intervals (transcription units - TUs, BIME-1s, BIME-2s, IHF sites, Fis sites, H-NS sites, MatP sites, chromosomal macrodomains, etc.)
for the enrichment of TCSs (binomial test), compares their N3E and score with mean TCSs N3E and score (t-test or poisson test), 
compares intervals mean score with genome mean score (t-test).

**Note: ** Returns some warning messages due to the ommitting statistics (t-test) for too short sets of values.

**Requirements:** python 3

**Input:** TAB files with trusted TCSs and score info, WIG genome score file, TAB transcription units data files, TAB intervals data files

**Output:** TAB file with numbers of TCSs are associated with TUs compartments (USUS, USGB, GBDS, DSDS), 
TAB file with TCSs-TUs association analysis (TCSs number, TCSs N3E, TCSs score, TUs compartments score), TAB file with normalized numbers of TCSs are associated with TUs,
TAB file with the number of TCSs are associated with particular intervals set and statistics (TCSs number, TCSs N3E, TCSs score),
TAB file with the number of TCSs are associated with particular intervals (BIMEs-1, BIMEs-2), TAB file with intervals score statistics

######################

## Average_and_plot_tracks.py

Script takes all N3E tracks for all experimental and control samples of all replicates and create genome-wide plots of Topo IV cleavage activity.

**Requirements:** python 3

**Input:** WIG files with N3E data.

**Output:** WIG files with averaged N3E data. Genome-wide plots of Topo IV cleavage activity.

######################

## Average_wig_files.py

Takes a set of WIG files (organized as a dictionary) and computes by-position average WIG. Returns correlation matrix 
(file and heatmap) of pair-wize correlations between genomic tracks.

**Requirements:** Python 2 or 3

**Input:** WIG files to be averaged

**Output:** Averaged WIG file, correlation matrix (CSV), correlation matrix heatmap (PNG), clusterized correlation matrix (CSV), clusterized correlation matrix heatmap (PNG)

######################

## Bed_to_wig_converter.py

Script takes BED files with coverage depth and converts them to WIG format.

**Requirements:** Python 2 or 3

**Input:** BED files

**Output:** WIG files

######################

## Compute_fold_enrichment.py

Takes two WIG files (for IP and Mock control) and computes by-position fold enrichment. 

**Requirements:** Python 2 or 3

**Input:** WIG files (IP and Mock control)

**Output:** WIG file with Fold Enrichment

######################

## Convert_coordinates_broadpeak.py

Converts intervals coordinates in one reference genome to coordinates in
another reference genome for a set of regions in the BroadPeak format.

**Requirements:** Python 3, blast+

**Input:** Reference genome 1 (FASTA), Reference genome 2 (FASTA), BroadPeak file with intervals

**Output:** BroadPeak file with intervals

######################

## FE_over_US_GB_DS.py

Takes wig tracks of different continuous genomic features (Topo IV FE, gyrase FE, GC%, MukB ChiP-Seq, etc.). Computes signal over TUs upstream (US),
downstream (DS) and over TUs bodies (GB). 

**Requirements:** Python 3

**Input:** Files with signal data (WIG), genome annotation (GFF or BroadPeak), regions to be omitted (BroadPeak)

**Output:** WIG files with metagene signal over all TUs, TAB files with average signal for each of TUs, plot of average signal over all TUs, histogram of the signal over TUs

######################

## FE_over_US_GB_DS_binning_and_statistics.py

Takes wig tracks of different continuous genomic features (Topo IV FE, gyrase FE, GC%, MukB ChiP-Seq, etc.). Computes signal over TUs upstream (US),
downstream (DS) and over TUs bodies (GB). Script is dedicated for noisy data. 
To handle the data it performes data binning for smoothing. Also it keeps signal data for all TUs to contstruct confidential interval.
To plot different features together (e.g. ChIP-Seq data from different expriments or conditions) it normalizes the data.

**Requirements:** Python 3

**Input:** Files with signal data (WIG), genome annotation (GFF or BroadPeak), regions to be omitted (BroadPeak), size of a bin

**Output:** WIG files with metagene signal over all TUs, TAB files with average signal for each of TUs, 
metagen plot of average signal over all TUs, bar plot of average signal in US, DS regions and in TUs bodies

######################

## Make_motif_logo.py

Script creates and plots a web-logo for a multiple alignment of DNA sequences.

**Requirements:** Python 3

**Input:** Multiple alignment of DNA sequences in FASTA format.

**Output:** Web-logo plot.

######################

## Microscopy_data_analysis.py

Analysis of cell length and DNA content distributions.

**Requirements:** Python 3

**Input:** XLSX tables with raw data of cell measurments.

**Output:** Plot of cell length distribution. Plot of DNA content distribution.

######################

## Plot_genome_tracks.py

The script takes diades of WIG files contain N3E or N5E values: A+IP+, A-IP+.
It smoothes A-IP+ track and divides the A+IP+ track by it.
Once obtains A+IP+_div, the script plots signal over the genomic positions.
Alongside N3E genomic tracks, density of TCSs and GCSs is plotted, matS sites locations, transcription levels, and enrichment of MatP.

**Requirements:** Python 3

**Input:** BroadPeak files with coordinates of TCSs, GCSs, MatS sites.
WIG files with N3E tracks for Gyrase and Topo IV, fold enrichment track for MatP, transcription levels.

**Output:** Plot with genomic tracks.

######################

## Plot_genome_tracks_tetrade.py

The script takes tetrades of WIG files contain N3E or N5E values: A+IP+, A+IP-, A-IP+, A-IP-.
It smooths A+IP- and A-IP- tracks and divides A+IP+ and A-IP+ by them.
Alongside N3E genomic tracks, density of TCSs or GCSs is plotted, and transcription levels.
A plot is generated: signal coverage over the genome for treated and untreated samples.

**Requirements:** Python 3

**Input:** BroadPeak files with coordinates of TCSs or GCSs.
WIG files with N3E tracks for Gyrase or Topo IV, transcription levels.

**Output:** Plot with genomic tracks.

######################

## Plot_signal_over_transcription_units.py

Takes signal over transcriptions units in WIG format generated by FE_over_US_GB_DS.py, plots metagene signal over TUs upstream, downstream and TU bodies.

**Requirements:** Python 3

**Input:** WIG files with metagene signal over TUs sets

**Output:** Plots with metagene signal

######################

## Return_reproducible_peaks.py

Takes a dictionary of NarrowPeak or BroadPeak files with intervals.
Identifies overlapping regions and writes them as a NarrowPeak file.

**Requirements:** Python 2 or 3

**Input:** Peaks coordinates (NarrowPeaks), reference genome (FASTA)

**Output:** Reproducible peaks coordinates (NarrowPeak), heatmap (number of shared peaks), heatmap (Jaccardian distance between peak sets), Venn diagram (peak sets overlap)

######################

## Return_unique_peaks.py

Takes a dictionary of narrowPeak files with peak regions.
Identifies unique regions which do not overlap with other peaks, writes them as a narrowPeak file.

**Requirements:** Python 2 or 3

**Input:** Peaks coordinates (NarrowPeaks), reference genome (FASTA)

**Output:** Unique peaks coordinates (narrowPeak), heatmap (number of shared peaks), heatmap (Jaccardian distance between peak sets), Venn diagram (peak sets overlap)
