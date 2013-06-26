## Custom Made Modules
A_hash_file.py - a script that hash the second column using first column as key  <br />
B_hash_mRNA_IDs.py - returns a uniq mRNA id hash <br />
C_loadFasta.py - script to load fasta sequences <br />
D_longest_fasta_sequence_header.py - script return headers of longest sequence <br />
E_get_chr_size_gff3.py - script takes a gff3 file and returns max position for each chromosome <br />


## smallRNA Clustering Scripts
4a.py - calculating clusters based on genotype (input IGV file) <br />
4a1.r - for plotting results from previous step <br />
4b.py - find clusters regulations and pattern based on size <br /> 
4c.py - from the clusters, make it to the inter-intra genic analysis <br />
4d.py - calculating regulated sequences in the cluster <br />


## FASTA Handlers
5.py  - make fasta files <br />
5b.py - take out mapping positions from igv file using fasta file containing sequences <br />
6.py  - a script for taking out cDNAs from transcripts file i.e. MG20 file <br />
7.py  - a script for any fasta file which looks for a pattern and returns a count and possibility of being random <br />


## miRNA Mapping
8a.py - script for replacing 'U' to 'T' <br />


## Genome Gap Filling Simulation
9.sh - script for gap filling project <br />
9a.py - for taking out rep element <br />
9b.py - for replacing genome by N <br />
9c.py - for taking out sequence where only one read is mapped <br />
9d.py - take out the rep elements and put these in the genomic region <br />

## Gap Filling Real Data
10a.py - take out N-region from the ref_genome <br />
10b.py - remove any additional N-region which might be present in 10Kb flanking region <br />
10c.py -  <br />
10q.py - calculating insert_size/distances <br />
10r.py - check if you have all elements <br />
10s.py - make reverse complement multi-fasta <br />
10t.py - make summary output for elements <br />

### ShortRan Scripts
12.sh - /plant/2011_week37/ <br />
12a.py - /plant/2011_week40/20111004/ <br />
12b.py - for taking out sequences with particular pattern <br />
12c.py - make fasta file from profile <br />
12d.py - replace U by T in miRNA database <br />
12e.py - for counting mismatches <br />
12f.py - for miRNA mapping from counting <br />
12g.py - cluster predictions <br />
12h.py - genome region analysis cluster by cluster or position by position <br />
12i.py - for counting regulated sequences in cluster <br />
12j.py - for counting unique size sequences and plotting the distribution from profile files <br />
12k.py - for counting unique size sequences and plotting the distribution from cluster files <br />
12l.py - for making profile from fastq files <br />
12m.py - remove reads which were mapped on repeats <br />
12n.py - for normalizing profiles <br />
12o.py - generate a text file with 0's of size total_no_of_libraries*total_no_of_libraries <br />
12p.py - for plotting expression data between two sets <br />
12r.py - for making sql batch script <br />
12s.py - add anotation to the sequences <br />
12t.py - script for making mySQl-add column batch file <br />
12u.py - make a script which can make a fasta file with abundance in the header as of format Sequence-xAbundance <br />
12v.py - script for making igv files for clustering and visualization <br />
12w.py - script for filtering reads based on score <br />
12x.py - script for making chromosome length file for miRNA predictions <br />
12y.py - script for making file compatible for tasiRNA predictions <br />
12z.py - script for taking out clusters of ta-siRNAs <br />
12aa.py - script for combining sequences for profile with the mapped genomic sequences which contain genomic annotation <br />
12ab.py - script for making header for the table with library wise abundances <br />
12ac.py - script for making table the library wise abundances <br />
12ad.py - script for plotting abundances <br />
12ae.py - python script for reaplcing libraries header <br />
12af.py - miRNA to MySQL database <br />
12ag.py - tasi-RNA to MySQL database <br />
12ah.py - add length coloumn to the mySQL database <br />
12ai.py - make a file with unique abundances <br />
12ai.py - <br />
12aj.py - script for spliting the fastq files by size <br />
12ak.py - make artificial adapters <br />
12al.py - parse the mirdeep 2 output to the mysql supported output <br />
12am.py - script to find the miRNA* sequences in profiles - /Users/vgupta/Desktop/script/python <br />
12an.py - script to add an identifier based non-redundancy - /Users/vgupta/Desktop/script/python <br />
mysql_batch - for saving file into the mysql database <br />


### gapfillRE
13_20110929_gapfillRE.sh - shell script processing other python scripts data <br />
13_20110929_positive_control.sh - for running posistive control- a bit different as input comes from blast <br />
13a.py -	for filtering of reads where both ends map to rep elements <br />
13b.py - 	for making reference compatible, i.e. adding headers, removing small letters <br />
13c.py - 	for taking out all gap positions <br />
13d.py -	take out genomic sequence with flanking regions <br />
13e.py -	remove additional N regions around targeted gap <br />
13f.py -	take out hanging reads mapping on the flanking region <br />
13g.py -	filter out pairs mapped to the flanking region <br />
13h.py - 	filter out diretional reads i.e. for 5'&3', 5',3' <br /> 
13i.py - 	taking out top four condidate suitable for replacement of gap and make score table <br />
13j.py -	reporting for gap regions which have no appropriate rep element for gap <br />
13k.py - 	pick out best possible element from scores <br />
13l.py -	print final list of elements with score <br />
13m.py - 	count correctly inserted elements(only for positive controls) <br />
13n.py - 	correct sequence name in fasta file (remove every thing after spaces), problem when mapping <br />
13o.py - 	for taking out a particular fasta sequecnce <br />
13p.py - 	script to remove pair mapped <br />
13q.py - 	take out all the contigs alraeady placed in the psuedomolecule <br />
13r.py -  	add length of the contigs <br /> 
13s.py -    add distances from 5 prime and 3 prime ends <br />

### Bactrial Genome Project With Niels
14.py - for finding a gene in many genomes <br />
14b.py - for finding a gene in many genomes using blast for unannotated genomes <br />
14c.py - script for taking list of genes and concatanating these by species. <br />

### Genome-wide Signatures
15.py - script to process the genome wide signature <br />
15a.py - script to add length and relavant columns <br />

### Svend's Data 
16.py <br />

### SpearmanRank
spr.py open file and calculate spearman co-efficient between all columns <br />

### Counting Corrected Reads
correct_read.py - count reads that has been corrected by ECHO <br />

### Making Patterns '\/_' For Regulations
make_patterns.py - making patterns '\/_' for regulations <br />

### Using R From Python
18_plot_sv.py - for plotting results obtinaed from the breakdancer <br />

### Yasu's Data
19_filter_markers.py - for filtering positions with the markers and storing these <br />
19_merge_marker.py - script for merging different files based on some columns <br />
19_remove_marker_positions.py - script for removing the existing markers and keeping only new SNPs <br />

### 28 Accession Data
20_compare_fq_mapping.py - script for comparing read-1 and read-2 mapped files to same reference <br />
20_divide_on_adaptors.py - script for deviding fastq file based on different adapters (demultiplexing) <br />
20_trim_reads.py - script for trimming the fastq reads and quality scores <br />


### Fastq Script Kit
20_compare_fq_mapping.py - script for comparing read-1 and read-2 mapped files to same reference <br />
20_divide_on_adaptors.py - script for deviding fastq file based on different adapters (demultiplexing) <br />
20_trim_reads.py - script for trimming the fastq reads and quality scores <br />
20_compare_fq.py - script for counting common reads in two fastq files <br /> 
20d_count_mapped_fastq_inSam.py - script for counting common reads in two fastq files - script for counting the reads mapped <br />


### Function For Making Filtered Fastq File 
17_20120212_filter_fastq.py - /Users/vikas0633/Desktop/plant/2012_week7 <br />

### Genomic Toolkit
21a_remove_chacters.py - this script removes the any other character than ATGCN
21b_better_header.py - this script keeps only 4th field separated by '|' <br />
21c_add_1_start.py - this script can add +1 to start position in a fasta file  <br />
21d_take_out_gene.py - this script takes out a sequence from fasta file given correct header name <br />
21d_take_out_gene_list_headers.py - this script takes out a sequence from fasta file given correct header names in a file <br />
21e_gff2gtf.py - this script converts gff, gff3 format to gtf format <br />
gtf_to_gff.pl - this script converts gtf to gff3 format <br />
81_parse.pl - script for calculating N50 value <br />
21f_merge_two_files.py - script for merging two files based on given columns <br />
21g_para_gtf.py - script for calculating exon/intron/transcripts lengths <br />
gff_convert.pl - script for inter-converting different gff formats <br />
intersection of gene models - bedtools intersect <br />
21h_calculate_seq_len.py - take a fasta file and print sequences in decreasing length <br />
21h_plot_seq_len.py - take a fasta file and plot sequence length <br />
21i_RMoutput2GTF.py - take a tab-formatted RMoutput file as parse it to make a gtf file <br />
21j_orf2fasta.py - script takes fasta file and output from orffinder and take out sequences with the orfs <br />
21k_make_input4_glimmerHMM.py - this scripts takes a gene structure file (gff/gtf) and makes a exon file parsable by glimmerHMM <br />
gff_to_genbank.py - Convert a GFF and associated FASTA file into GenBank format <br />
21l_pileup2GTF.py - script converts a pileup to a gtf file based on the coverage <br />
21m_gff2genestru.py - script creates input for gb format conversion script <br />
21n_overlap_gff.py - takes two or more gff files merge the files where you see an overlap <br />
21n_intersect_gff.py - takes two or more gff files merge the files where you see an intersection <br />
21o_extract_seq_model.py - script takes out sequences/GTF models from given co-ordinate <br />
21p_filter_fasta.py - script to filter fasta file based on the length of the sequences <br />
21q_combine_GTF.py - This is the script for combining various annotations files <br />
21r_make_CDS.py- script to create CDS file from fasta (containing exon sequences generated by bedtools) and GTF/GFF3 file <br />
21s_summary_eval.py - script for summarizing eval output <br />
21t_tau.py - script to add ORF to the gff file <br />
21u_make_gff2.py - script makes gff2 file for the TAU input, same as Stig's 26_parse.pl <br />
21v_format_gff3.py- script to format gff3 file in order to put in MySQL table <br />
21v2_format_gff3.py- script to format gff3 file in order to put in MySQL table sequal to 21v <br />
21w1_format_fasta.py- fasta file has duplicate entries <br /> 
21w1_format_orthoMCL.py - format OrthoMCL output <br />
21x_exon_repeat.py- find the exon Repeat over lap <br />
21y_strand_fasta.py - script takes a GFF3 file and correct fasta file if minus strand <br />
21z_foramt_IPR.py - script takes raw output from IPRScan and make non-redundant gene_ID\annotation <br />
21aa_countMShit_in_GFF.py - script to count the uniq MS supported genes <br />
21ab_split_gff.py - script to split sorted GFF file based on contig/sequence/chro name <br />
21ac_addType.py - script to add gene type <br />
21ad_makebed.py - script to make bed format file from the given column names <br />
21ae_correct_UTR.py - script to correct the UTR co-ordinates <br />
21af_format_protein_list_headers.py - script to get the corresponding headers between corrected and real fasta file <br />
21ag_cal_CSD_gene_overlap.py - script to calculate the CDS vs gene overlap <br /> 
21ah_find_longest_isoform.py - script was made for finding longest isoform in the spider protein set <br />
21ah_count_N_between_genes.py - script to count Ns between the genes <br />
21ai_modify_gene_names.py - script to modify gene names based on N counts <br />
21aj_add_mRNA.py - script to add dummy mRNAs if absent <br />
21ak_remove_redundant.py - script to remove the redundant node gene models <br />
21al_correct_strand.py - this script takes strand from CDS and assigns the same to mRNA, exons and UTRs, GFF3 files <br />
21am_update_GFF3_fasta.py - this script updates GFF3 and fasta given a different file <br />

	
### Transcripts Handlers
22a.py - script for parsing tophat/cufflink generated GTF files against a target (-G cufflink) annotation file <br />
22b.py - normalize transcript profile table <br />
22c.py - script for making plots from profile tables generated using 22b <br />
22d.py - add profile tables to MySQL <br />
22e.py - add annotations to profiles using fasta files <br />
22f.py - add annotations to profiles using two column formatted file <br />
22g.py - script to get pattern frequency from a profile table given a regulation, abundance and score cut-off  <br />
22h.py - script for finding complementary pattern between small RNAs and transcripts <br />

### MYSQL
23a_mysql_header.py - script for making headers for mysql tables <br />

### Blast
24a_filter_blast.py - script for filtering blast results  <br />

### Python Plots
25a_plot_gene_freq.py - script for plotting gene frequencies across each chromosome <br />

### MirDeepP Summary
26_summary_mirDeepP.py - script for taking all the outputs from mirDeepP and putting it together <br />


### Spider Project
27_summary_MS_hit.py - script for process MS hit text file <br />
27_foramt_fasta_spider.py - script to format the fasta headers according to the Thomas's explanations <br />

### UNC RNA-seq project
28a_obo_parser.py - script to obo file from the geneontology.org <br />
28b_MSU_RAP_ids.py - MSU id parser <br />
28c_gff3_validator.py - Script to validate a gff3 file <br />


### 29. snpEff data analysis
29a_MakeGeneWideTable.py - script to put the snpEff data togehter <br />
29b_MakeGeneWideTableUnique.py - script to summarize snpEff data <br />

### General Scripts
100_intersect_columns.py - script to find non-overlapping entries between the two columns <br />
21ab_split_gff.py - script to split sorted GFF file based on contig/sequence/chro name <br />
101_filter_fastq_len.py - script to filter a fastq file based on read length <br />
102_flat2fasta_anno.py - script to make fasta file from the MySQL output <br />
103_sort_gff_blocks.py - script to sort GFF3 file blocks <br />
104_intersect_files_column.py - script to print the desired columns given keys from the files <br />
105_match_IDs_from_2gff3_files.py - script will take two gff3 files and print out the corresponding mRNA IDs <br />
106_filter_out_against_genelist.py - this script will filter out the genes which are in the list <br />
107_ParameterGFF3.py - Script to calculate Gene, mRNA, exon, CDS count, Total length and average length <br />
108_filterExactOverlapGFF3.py - script to filter overlapping start/end genemodels <br />
109_AddPhaseGFF3.py - script to add Phase <br />
