# TAPAS: Tool for Alternative Polyadenylation Site Analysis

Accoring to the central dogma of molecular biology, a pre-mRNA is synthesized from the coding sequence of a gene during the transcriptional process. This pre-mRNA is coverted to a mature mRNA by the post-transcriptional process. Post-transcriptional process contains three major steps. One of these steps is the addition of polyadenylation (polyA) tail using the polyadenylation pocess. Therefore, the polyadenylation process has two steps: cleavage the 3' end of a pre-mRNA and addition of a polyA tail at the cleavage site. But, due to certain *cis*-acting elements and *trans*-acting factors alternative cleavage sites can be formed from a pre-mRNA. More precisely, a single pre-mRNA can produce more than one mRNA with 3' untranslated regions (3' UTRs) of different lengths. TAPAS is a RNA-Seq based tool for detecting these alternative polyadenylation cleavage sites (APSs) with in a gene. If two biological samples with multiple replicates are given, TAPAS indentifies the differentially expressed APSs between these samples. It can also extend its differential analysis to identify shortening/lengthening event genes between the samples.

### Requirements
1. gcc version 4.8.3 20140911
2. samtools 1.3
3. R 3.1.3
4. matrixStats, locfit and stats packages of R

## APS detection
APS_detection of TAPAS detects novel APSs of genes.

USAGE	./APS_detection {OPTIONS}

OPTIONS

	-ref <annotation_file_name>	An annotation file is given using this option. 
					Human annotation file is given in Finding_APS directory.

	-cov <coverage_file_name>	A read coverage file is provided using this option.
					Samtools is used to have the read coverage.

	-l <int>			Read length

	-o <output_file_name>		Output file name is given using this option.

	-p <double>			A penalty value can be provided using this option.
					If nothing is given, this value is determined from the coverage
					of the 3' UTR frame.

## Differential APS analysis
