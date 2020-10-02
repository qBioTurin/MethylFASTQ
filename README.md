# Methyl-Seq

MethylFASTQ is a Python tool to simulate bisulfite sequencing data in a
highly customizable way. 

## Features

MethylFASTQ produces two kind of files.
* FASTQ file(s) that contains the sequenced reads. The id of each FASTQ record identifies the true mapping position of that read.
* methylation call file that contains the true methylation call of each sequenced cytosine

It allows to simulate both whole genome bisulfite sequencing (WGBS) and targeted bisulfite sequencing experiments.

## Options

* -i fasta-file : path of the FASTA file containing the genome to be sequenced
* -o output-file-path : directory where to save the artificial dataset
* --seq {single_end,paired_end} : selected sequencing mode
* --lib {directional,non_directional} : selected library mode
* --chr chromosome-id [chromosome-id ...] : list of FASTA id of chromosomes to be sequenced
* --regions target-regions : path of the .bed file containing chromosome regions to be sequenced
* --coverage coverage : selected depth of coverage
* --fragment fragment-size : size of produced fragments in fragmentation step
* --read read-length : length of produced reads
* --processes num-processes : maximum number of producer processes
* --buffer buffer-size :  buffer size of a single process
* --cg CG-methylation-probability : probability that a C in CG context is methylated
* --chg CHG-methylation-probability : probability that a C in CHG context is methylated (H = {A, C, T})
* --chh CHH-methylation-probability : probability that a C in CHH context is methylated
* --snp snp-probability : probability that a single nucleotide is a SNP (spontaneous mutation)
* --error sequencing-error-probability : probability that a single nucleotide is a sequencing error
* --maxq max-phred-score : maximum quality score in the produced reads
* --minq min-phred-score : minimum phred score in the produced reads (not implemented yet)

## Targeted sequencing 

If you want to use targeted sequencing mode, you have to provide a region file describing the chromosome regions to be sequenced. 
The region file must be a TAB-delimited file with 3+ columns describing
* the chromosome ID 
* the starting bp of the region to be sequenced
* the last  bp of the region to be sequenced 

Check the region file in the example folder. 