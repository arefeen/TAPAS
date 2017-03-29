# TAPAS: Tool for Alternative Polyadenylation Site Analysis

Accoring to the central dogma of molecular biology, a pre-mRNA is synthesized from the coding sequence of a gene during the transcriptional process. This pre-mRNA is coverted into a (mature) mRNA by the post-transcriptional process. The post-transcriptional process consists of three major steps. One of them is the addition of polyadenylation (polyA) tail using the polyadenylation pocess, which in turn consists of two steps: the cleavage at the 3' end of a pre-mRNA and the addition of a polyA tail at the cleavage site. But, due to the effect of certain *cis*-acting elements and *trans*-acting factors alternative cleavage sites can be formed from in a pre-mRNA. More precisely, a single pre-mRNA can often produce more than one mRNA with 3' untranslated regions (3' UTRs) of different lengths. TAPAS is a tool for detecting such alternative polyadenylation cleavage sites (APSs) within a gene from RNA-Seq data. If two biological samples with multiple replicates are given, TAPAS can indentify differentially expressed APSs between the samples. Moreover, its differential analysis has been extended to discover the shortening/lengthening of 3' UTRs within a gene.


### Requirements
1. gcc version 4.8.3 20140911
2. samtools 1.3
3. R 3.1.3
4. matrixStats, locfit and stats packages of R

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


### Ouput of APS detection
The output file consists of six columns: gene name, chromosome name, strand of the gene, detected APSs, abundance of those APSs, read count of those APSs respectively.
Note: The abundance of each detected APS = read count of the APS / length of the 3' UTR (that contains the APS)  

	
## Differential APS analysis
