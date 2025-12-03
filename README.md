# rp1b

## Planning - Writing functions for Making mutations writing code for:
Cleaning FASTA files into base string
Making mutations (SNPs, indels as defined functions)
x30 depth with 100bp read sections on the mutated genomes - WIP 
Use minimap (mapped to original reference genome)
Apply bcftools
Find precision and recall of SNP caller

## Use code to make the pipeline (nextflow):
Map reads
Run 2 different variant callers (bcftools and snippy)
Combine results into 1 VCF file

## Validation:
Test code on Zam’s genomes and compare against my own validation set
Run pipeline on one real E coli FASTQ and reference
Assign each VCF record a score of how much you trust it (accuracy)

## Discussion
Discuss and evaluate performance of tool and what did/didn’t work
Comment on precision /recall of tool on simulated data
Is the combined VCF better than either?
Walk through how it performed on real data
Does tview support your high/low scores?

Submit in github

