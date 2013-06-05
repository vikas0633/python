### load function scripts
A_hash_file.py - /Users/vikas0633/Desktop/script/python/ - a script that hash the second column using first column as key
B_hash_mRNA_IDs.py - /Users/vikas0633/Desktop/script/python/ - returns a uniq mRNA id hash
C_loadFasta.py - /Users/vikas0633/Desktop/script/python/ - script to load fasta sequences
D_longest_fasta_sequence_header.py - /Users/vikas0633/Desktop/script/python/ - script return headers of longest sequence
E_get_chr_size_gff3.py - script takes a gff3 file and returns max position for each chromosome - /Users/vikas0633/Desktop/script/python/

3.R Plant_week29 - calculating cluster distribution across the chromosome in R

### cluster scripts
4a.py plant_week29 - calculating clusters based on genotype (input IGV file)
4a1.r plant_week29 - for plotting results from previous step
4b.py plant_week30 - find clusters regulations and pattern based on size 
4c.py plant_week30 - from the clusters, make it to the inter-intra genic analysis
4d.py plant_week31 - calculating regulated sequences in the cluster

### taking out DNA and igv
5.py plant_week30 - make fasta files
5b.py plant_week30 - take out mapping positions from igv file using fasta file 				 containing sequences 

6.py plant/lotus_mRNA - a script for taking out cDNAs from transcripts file i.e. MG20 file

7.py plant/2011_week31 - a script for any fasta file which looks for a pattern and returns a count and possibility of being random


### miRNA Mapping
8a.py plant/2011_week31 - script for replacing 'U' to 'T'


### gap filling simulated
9.sh 2011_gap_filling/2011_week33/ script for gap filling project
9a.py 2011_gap_filling/2011_week33/ for taking out rep element
9b.py 2011_gap_filling/2011_week33/ for replacing genome by N
9c.py 2011_gap_filling/2011_week33/ for taking out sequence where only one read is mapped
9d.py 2011_gap_filling/2011_week33/ take out the rep elements and put these in the genomic region

### gap_filling real data
10a.py - take out N-region from the ref_genome
10b.py - remove any additional N-region which might be present in 10Kb flanking region
10c.py -  
10q.py - 2011_gap_filling/2011_week39/20110926 calculating insert_size/distances
10r.py - 2011_gap_filling/2011_week39/20110927 ### check if you have all elements
10s.py - 2011_gap_filling/2011_week39/20110927 ### make reverse complement multi-fasta
10t.py - 2011_gap_filling/2011_week39/20110928 ### make summary output for elements

### shell Script for smallRNA pipeline
12.sh - /plant/2011_week37/
12a.py - /plant/2011_week40/20111004/
12b.py - /plant/2011_week37/ for taking out sequences with particular pattern
12c.py - /plant/2011_week38/20110920/ make fasta file from profile
12d.py - /plant/2011_week38/20110920/ replace U by T in miRNA database
12e.py - /plant/2011_week38/20110920/ for counting mismatches
12f.py - /plant/2011_week38/20110920/ for miRNA mapping from counting
12g.py - /plant/2011_week38/20110923/ cluster predictions
12h.py - /plant/2011_week38/20110923/ genome region analysis cluster by cluster or position by position  
12i.py - /plant/2011_week38/20110924/ for counting regulated sequences in cluster
12j.py - /plant/2011_week40/20111004/ for counting unique size sequences and plotting the distribution from profile files
12k.py - /plant/2011_week40/20111004/ for counting unique size sequences and plotting the distribution from cluster files
12l.py - /plant/2011_week42/20111017 - for making profile from fastq files
12m.py - /plant/2011_week42/20111017 - remove reads which were mapped on repeats
12n.py - /plant/2011_week42/20111017 - for normalizing profiles
12o.py - /plant/2011_week42/20111019 - ### generate a text file with 0's of size total_no_of_libraries*total_no_of_libraries
12p.py - /plant/2011_week42/20111019 - for plotting expression data between two sets
12r.py - /plant/2011_week43/scripts -  for making sql batch script
mysql_batch - /plant/2011_week43/scripts - for saving file into the mysql database
12s.py - add anotation to the sequences
12t.py - script for making mySQl-add column batch file
12u.py - make a script which can make a fasta file with abundance in the header as of format Sequence_xAbundance
12v.py - script for making igv files for clustering and visualization
12w.py - script for filtering reads based on score
12x.py - script for making chromosome length file for miRNA predictions
12y.py - script for making file compatible for tasiRNA predictions
12z.py - script for taking out clusters of ta-siRNAs
12aa.py - script for combining sequences for profile with the mapped genomic sequences which contain genomic annotation
12ab.py - script for making header for the table with library wise abundances
12ac.py - script for making table the library wise abundances
12ad.py - script for plotting abundaces
12ae.py - python script for reaplcing libraries header
12af.py - miRNA to MySQL database 
12ag.py - tasi-RNA to MySQL database
12ah.py - add length coloumn to the mySQL database
12ai.py - make a file with unique abundances
12ai.py - 
12aj.py - script for spliting the fastq files by size
12ak.py - make artificial adapters
12al.py - parse the mirdeep 2 output to the mysql supported output
12am.py - script to find the miRNA* sequences in profiles - /Users/vgupta/Desktop/script/python
12an.py - script to add an identifier based non-redundancy - /Users/vgupta/Desktop/script/python



### gapfillRE
13_20110929_gapfillRE.sh - shell script processing other python scripts data
13_20110929_positive_control.sh - for running posistive control- a bit different as input comes from blast
13a.py -	for filtering of reads where both ends map to rep elements
13b.py - 	for making reference compatible, i.e. adding headers, removing small letters
13c.py - 	for taking out all gap positions
13d.py -	take out genomic sequence with flanking regions
13e.py -	remove additional N regions around targeted gap
13f.py -	take out hanging reads mapping on the flanking region
13g.py -	filter out pairs mapped to the flanking region
13h.py - 	filter out diretional reads i.e. for 5'&3', 5',3' 
13i.py - 	taking out top four condidate suitable for replacement of gap and make score table
13j.py -	reporting for gap regions which have no appropriate rep element for gap
13k.py - 	pick out best possible element from scores
13l.py -	print final list of elements with score
13m.py - 	count correctly inserted elements(only for positive controls)
13n.py - 	correct sequence name in fasta file (remove every thing after spaces), problem when mapping
13o.py - 	for taking out a particular fasta sequecnce
13p.py - 	script to remove pair mapped 
13q.py - 	take out all the contigs alraeady placed in the psuedomolecule
13r.py -  	add length of the contigs 
13s.py -    add distances from 5 prime and 3 prime ends

###bactrial project with Niels
14.py - for finding a gene in many genomes
14b.py - for finding a gene in many genomes using blast for unannotated genomes
14c.py - script for taking list of genes and concatanating these by species.

### genome-wide signatures
15.py - script to process the genome wide signature 
15a.py - script to add length and relavant columns

###svend's data
16.py

### common
spr.py open file and calculate spearman co-efficient between all columns

### counting correcting
correct_read.py - count reads that has been corrected by ECHO - /Users/vikas0633/Desktop/plant/2012_week7/package_0.2



### making patterns '\/_' for regulations
make_patterns.py - making patterns '\/_' for regulations



### using R from python
18_plot_sv.py - for plotting results obtinaed from the breakdancer-/u/vgupta/2011_genome_structure/2012_week7/20120215

### yasu's data
19_filter_markers.py - ~/Desktop/yasu/ - for filtering positions with the markers and storing these
19_merge_marker.py - ~/Desktop/yasu/ - script for merging different files based on some columns
19_remove_marker_positions.py - ~/Desktop/yasu/ - script for removing the existing markers and keeping only new SNPs 

### 28 accession data
20_compare_fq_mapping.py - script for comparing read-1 and read-2 mapped files to same reference
20_divide_on_adaptors.py - script for deviding fastq file based on different adapters (demultiplexing)
20_trim_reads.py - script for trimming the fastq reads and quality scores


### fastq script kit
20_compare_fq_mapping.py - script for comparing read-1 and read-2 mapped files to same reference
20_divide_on_adaptors.py - script for deviding fastq file based on different adapters (demultiplexing)
20_trim_reads.py - script for trimming the fastq reads and quality scores
20_compare_fq.py - script for counting common reads in two fastq files 
20d_count_mapped_fastq_inSam.py - script for counting common reads in two fastq files - script for counting the reads mapped


### function for making filtered fastq file using a infile(fasta_file) and fastq files
17_20120212_filter_fastq.py - /Users/vikas0633/Desktop/plant/2012_week7

### genomic toolkit
21a_remove_chacters.py - this script removes the any other character than ATGCN
21b_better_header.py - this script keeps only 4th field separated by '|'  - /Users/vikas0633/Desktop/plant/2012_week29/21b_better_header.py
21c_add_1_start.py - this script can add +1 to start position in a fasta file - /Users/vikas0633/Desktop/plant/2012_week29
21d_take_out_gene.py - this script takes out a sequence from fasta file given correct header name - /Users/vikas0633/Desktop/plant/2012_week29
21d_take_out_gene_list_headers.py - this script takes out a sequence from fasta file given correct header names in a file - /Users/vikas0633/Desktop/plant/2012_week29
21e_gff2gtf.py - this script converts gff, gff3 format to gtf format - /Users/vikas0633/Desktop/plant/2012_week30
gtf_to_gff.pl - this script converts gtf to gff3 format - /Users/vikas0633/Desktop/script/perl
81_parse.pl - script for calculating N50 value - /home/vgupta/short_projects/spider_gene_models/2012_week26/tranctula_palle_assembly_version1/
21f_merge_two_files.py - script for merging two files based on given columns - /Users/vikas0633/Desktop/plant/2012_week36
21g_para_gtf.py - script for calculating exon/intron/transcripts lengths - /Users/vikas0633/Desktop/script/python 
gff_convert.pl - script for inter-converting different gff formats - /Users/vikas0633/Desktop/script/perl
intersection of gene models - bedtools intersect
21h_calculate_seq_len.py - take a fasta file and print sequences in decreasing length -  /Users/vikas0633/Desktop/script/python
21h_plot_seq_len.py - take a fasta file and plot sequence length -  /Users/vikas0633/Desktop/script/python
21i_RMoutput2GTF.py - take a tab-formatted RMoutput file as parse it to make a gtf file - /Users/vikas0633/Desktop/script/python
21j_orf2fasta.py - script takes fasta file and output from orffinder and take out sequences with the orfs - /Users/vikas0633/Desktop/script/python
21k_make_input4_glimmerHMM.py - this scripts takes a gene structure file (gff/gtf) and makes a exon file parsable by glimmerHMM - /Users/vikas0633/Desktop/script/python
gff_to_genbank.py - Convert a GFF and associated FASTA file into GenBank format - /Users/vikas0633/Desktop/script/python
21l_pileup2GTF.py - script converts a pileup to a gtf file based on the coverage - /Users/vikas0633/Desktop/script/python
21m_gff2genestru.py - script creates input for gb format conversion script -  /Users/vikas0633/Desktop/script/python
21n_overlap_gff.py - takes two or more gff files merge the files where you see an overlap - /Users/vikas0633/Desktop/script/python
21n_intersect_gff.py - takes two or more gff files merge the files where you see an intersection - /Users/vikas0633/Desktop/script/python
21o_extract_seq_model.py - script takes out sequences/GTF models from given co-ordinate - /Users/vikas0633/Desktop/script/python
21p_filter_fasta.py - script to filter fasta file based on the length of the sequences - /Users/vikas0633/Desktop/script/python
21q_combine_GTF.py - This is the script for combining various annotations files - /Users/vikas0633/Desktop/script/python
21r_make_CDS.py- script to create CDS file from fasta (containing exon sequences generated by bedtools) and GTF/GFF3 file  - /Users/vikas0633/Desktop/script/python
21s_summary_eval.py - script for summarizing eval output - /Users/vikas0633/Desktop/script/python
21t_tau.py - script to add ORF to the gff file - /Users/vikas0633/Desktop/script/python
21u_make_gff2.py - script makes gff2 file for the TAU input, same as Stig's 26_parse.pl - /Users/vikas0633/Desktop/script/python
21v_format_gff3.py- script to format gff3 file in order to put in MySQL table - /Users/vikas0633/Desktop/script/python
21v2_format_gff3.py- script to format gff3 file in order to put in MySQL table sequal to 21v - /Users/vikas0633/Desktop/script/python
21w1_format_fasta.py- fasta file has duplicate entries - /Users/vikas0633/Desktop/script/python
21w1_format_orthoMCL.py - format OrthoMCL output - /Users/vikas0633/Desktop/script/python
21x_exon_repeat.py- find the exon Repeat over lap - /Users/vikas0633/Desktop/script/python
21y_strand_fasta.py - script takes a GFF3 file and correct fasta file if minus strand - /Users/vikas0633/Desktop/script/python
21z_foramt_IPR.py - script takes raw output from IPRScan and make non-redundant gene_ID\annotation  - /Users/vikas0633/Desktop/script/python
21aa_countMShit_in_GFF.py - script to count the uniq MS supported genes - /Users/vikas0633/Desktop/script/python
21ab_split_gff.py - script to split sorted GFF file based on contig/sequence/chro name -  /Users/vikas0633/Desktop/script/python
21ac_addType.py - script to add gene type - /Users/vikas0633/Desktop/script/python
21ad_makebed.py - script to make bed format file from the given column names - /Users/vikas0633/Desktop/script/python
21ae_correct_UTR.py - script to correct the UTR co-ordinates - /Users/vikas0633/Desktop/script/python
21af_format_protein_list_headers.py - script to get the corresponding headers between corrected and real fasta file - /Users/vikas0633/Desktop/script/python
21ag_cal_CSD_gene_overlap.py - script to calculate the CDS vs gene overlap - /Users/vikas0633/Desktop/script/python
21ag_cal_CSD_exon_overlap.py - script to calculate the CDS vs exon overlap - /Users/vikas0633/Desktop/script/python
21ah_find_longest_isoform.py - script was made for finding longest isoform in the spider protein set - /Users/vikas0633/Desktop/script/python
21ah_count_N_between_genes.py - script to count Ns between the genes - /Users/vikas0633/Desktop/script/python
21ai_modify_gene_names.py - script to modify gene names based on N counts - /Users/vikas0633/Desktop/script/python
21aj_add_mRNA.py - script to add dummy mRNAs if absent - /Users/vikas0633/Desktop/script/python
	
### transcripts
22a.py - /Users/vikas0633/Desktop/plant/2012_week29/22a.py - script for parsing tophat/cufflink generated GTF files against a target (-G cufflink) annotation file
22b.py - /Users/vikas0633/Desktop/plant/2012_week29/22b.py - normalize transcript profile table
22c.py - /Users/vikas0633/Desktop/plant/2012_week29/22c.py - script for making plots from profile tables generated using 22b
22d.py - /Users/vikas0633/Desktop/plant/2012_week29/22d.py - add profile tables to MySQL 
22e.py - /Users/vikas0633/Desktop/plant/2012_week29/22e.py - add annotations to profiles using fasta files
22f.py - /Users/vikas0633/Desktop/plant/2012_week29/22f.py - add annotations to profiles using two column formatted file
22g.py - /Users/vikas0633/Desktop/plant/2012_week29/22g.py - script to get pattern frequency from a profile table given a regulation, abundance and score cut-off 
22h.py - /Users/vikas0633/Desktop/plant/2012_week30/22h.py - script for finding complementary pattern between small RNAs and transcripts

### mysql
23a_mysql_header.py -/Users/vikas0633/Desktop/script/python/ - script for making headers for mysql tables

### Blast
24a_filter_blast.py - /Users/vikas0633/Desktop/script/python/ - script for filtering blast results 

### python plots
25a_plot_gene_freq.py - /Users/vikas0633/Desktop/script/python/ - script for plotting gene frequencies across each chromosome

### random scripts
26_summary_mirDeepP.py - /Users/vikas0633/Desktop/script/python/ - script for taking all the outputs from mirDeepP and putting it together


### Spider project
27_summary_MS_hit.py - /Users/vikas0633/Desktop/script/python/ - script for process MS hit text file
27_foramt_fasta_spider.py - /Users/vikas0633/Desktop/script/python/ - script to format the fasta headers according to the Thomas's explanations

### UNC RNA-seq project
28a_obo_parser.py - /Users/vikas0633/Desktop/script/python/ - script to obo file from the geneontology.org
28b_MSU_RAP_ids.py - /Users/vikas0633/Desktop/script/python/ - MSU id parser
28c_gff3_validator.py - /Users/vikas0633/Desktop/script/python/ - Script to validate a gff3 file


### 29. snpEff data analysis
29a_MakeGeneWideTable.py - /Users/vikas0633/Desktop/script/python/ - script to put the snpEff data togehter

### In general scripts
100_intersect_columns.py - /Users/vikas0633/Desktop/script/python/ - script to find non-overlapping entries between the two columns
21ab_split_gff.py - script to split sorted GFF file based on contig/sequence/chro name -  /Users/vikas0633/Desktop/script/python
101_filter_fastq_len.py - script to filter a fastq file based on read length - /Users/vikas0633/Desktop/script/python
102_flat2fasta_anno.py - script to make fasta file from the MySQL output - /Users/vikas0633/Desktop/script/python
103_sort_gff_blocks.py - /Users/vikas0633/Desktop/script/python/  - script to sort GFF3 file blocks
104_intersect_files_column.py - /Users/vikas0633/Desktop/script/python/ - script to print the desired columns given keys from the files
105_match_IDs_from_2gff3_files.py - /Users/vikas0633/Desktop/script/python/ - script will take two gff3 files and print out the corresponding mRNA IDs



