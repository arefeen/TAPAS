# TAPAS: Tool for Alternative Polyadenylation Site Analysis

Accoring to the central dogma of molecular biology, a pre-mRNA is synthesized from the coding sequence of a gene during the transcriptional process. This pre-mRNA is coverted into a (mature) mRNA by the post-transcriptional process. The post-transcriptional process consists of three major steps. One of them is the addition of polyadenylation (polyA) tail using the polyadenylation pocess, which in turn consists of two steps: the cleavage at the 3' end of a pre-mRNA and the addition of a polyA tail at the cleavage site. But, due to the effect of certain *cis*-acting elements and *trans*-acting factors alternative cleavage sites can be formed from in a pre-mRNA. More precisely, a single pre-mRNA can often produce more than one mRNA with 3' untranslated regions (3' UTRs) of different lengths. TAPAS is a tool for detecting such alternative polyadenylation cleavage sites (APSs) within a gene from RNA-Seq data. If two biological samples with multiple replicates are given, TAPAS can indentify differentially expressed APSs between the samples. Moreover, its differential analysis has been extended to discover the shortening/lengthening of 3' UTRs within a gene.


### Requirements
1. The tool runs on linux machine.
2. samtools 1.3
3. R 3.1.3
4. matrixStats, locfit and stats packages of R.

## APS detection
APS_detection of TAPAS detects novel APSs of genes.

USAGE

	./APS_detection {OPTIONS}

OPTIONS

	-ref <annotation_file_name>	An annotation file is given using this option. 
					e.g. Human annotation file is given in Finding_APS directory for reference.

	-cov <coverage_file_name>	A read coverage file is provided using this option.
					Samtools is used to have the read coverage.

	-l <int>			Read length

	-o <output_file_name>		Output file name is given using this option.

	-p <double>			A penalty value can be provided using this option.
					If nothing is given, the value is determined from the 
					read coverage of the 3' UTR frame.

EXAMPLE

	./APS_detection -ref refFlat_sf.txt -cov coverage_read_50.txt -l 76 -o expression_with_cp_read_50.txt


### Output of APS detection
The output file consists of six columns: gene name, chromosome name, strand of the gene, detected APSs, abundance of those APSs, read count of those APSs respectively. <br />
Note: The abundance of each detected APS = read count of the APS / length of the 3' UTR (that contains the APS)  

	
## Differential APS analysis
Diff_APS_Analysis of TAPAS does differential analyses between two biological samples.

USAGE

	./Diff_APS_Analysis {OPTIONS}

OPTIONS

	-C1 <cond1_f1,cond1_f2,cond1_f3,..>	Comma separated file names of condition 1 are given using this option. Each of these files
						is the APS detection file, outputted by the first part of TAPAS (outputted by APS_detection).
	
	-C2 <cond2_f1,cond2_f2,cond2_f3,..>	Comma separated file names of condition 2 are given using this option. Each of these files 
                                                is the APS detection file, outputted by the first part of TAPAS (outputted by APS_detection).

	-a <annotation_file_name>		An annotation file is given using this option. This file is similar to the annotation file
						of the APS detection analysis.

	-cutoff	<int>				Cutoff value is given using this option. This parameter is explained in TAPAS manuscript.
						Default value: 70

	-type	<d/s>				Type of differential analysis. d -> differential APS analysis, s -> shortening/lengthening
						event analysis.

	-o <output_file_name>			Ouput file name is given using this option. 
						Default: for differential APS analysis "diff_result_final.txt", 
							for shortening/lengthening event analysis "decision_output.txt"

EXAMPLE

	./Diff_APS_Analysis -C1 cond1_r1.txt,cond1_r2.txt,cond1_r3.txt,cond1_r4.txt,cond1_r5.txt,cond1_r6.txt -C2 cond2_r1.txt,cond2_r2.txt,cond2_r3.txt,cond2_r4.txt,cond2_r5.txt,cond2_r6.txt -a refFlat_sf.txt -cutoff 70 -type s -o deci_output.txt
					

### Output of differential analysis
For differential APS analysis, the output file contains eight columns: chromosome name, gene name, strand, APS, log2 fold change, p-value, adjusted p-value, decision. <br />
For shortening/lengthening event analysis, the output file contains six columns: chromosome name, gene name, strand, shorter APS, longer APS, log2 fold change. It also produces
differentially expressed APS file. This file contains eight columns: chromosome name, gene name, strand, APS, fold change, log2 fold change, p-value, adjusted p-value.
P-value is always adjusted using BH method.


## Helping command

### How to calculate coverage using samtools

Command to calculate coverage from bam file:

	samtools sort accepted_hits.bam -o sorted_accepted_hits.bam
	samtools index -b sorted_accepted_hits.bam
	samtools view -b sorted_accepted_hits.bam > accepted_reads.bam

	samtools depth accepted_reads.bam > read_coverage.txt
