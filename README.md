
# BASS (Bulk ATAC-Seq Snakemake pipeline) Introduction

Written by: Austin Montgomery

Date: 14 May 2020

## Pipeline Summary
1. Raw read QC: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. Trim reads (including adapters): [Trim Galore!](https://github.com/FelixKrueger/TrimGalore)
3. Map reads to chosen genome: [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
4. Remove PCR duplicates: [Picard](https://broadinstitute.github.io/picard/)
5. Remove poor quality reads and merge condition replicates: [Samtools](http://www.htslib.org/doc/samtools.html)
    1. Sort reads `samtools sort`
    2. Create alignment index `samtools index`
    3. Filter reads `samtools view -b -@ 4 -f 2 -h -q 30 chr2L chr2R chr3L chr3R chrX`
        - `-f 2` Only keep reads mapped in proper pair
        - `-h` Output header
        - `-q 30` Only keep reads with alignment score of >= 30
        - `chr2L chr2R chr3L chr3R chrX` Only keep reads from these chromosomes
6. Call narrow peaks on each replicate and merged alignment: [macs2](https://github.com/taoliu/MACS)
    1. Call reads `-f BAM -g dm --nomodel --shift -37 --extsize 73 --call-summits --keep-dup all`
        - `-f BAM` BAM format (This only takes reads from R1)
        - `-g dm` Drosophila genome
        - `--nomodel` Skip building the model
        - `--shift` Extend read back 37 bps
        - `--extsize 73` Extend read forward 73 bps
        - `--call-summits` Generate summits file
        - `--keep-dup all` Use any duplicates
7. Annotate peak files and find motifs: [Homer](http://homer.ucsd.edu/homer/)
    1. `findMotifsGenome.pl`
    2. `annotatePeaks.pl`
8. Create bigwig files: [deepTools](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html)
    1. `bamCoverage --extend-reads`
9. Concatenate and merge peak files to find consensus peakset: [bedtools](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)
    1. `cat {input.a} {input.b} | sort -k1,1 -k2,2n | bedtools merge -i - > {output}`
        - `cat {input.a} {input.b}` Add the different condition peaksets together into one peakfile
        - `sort -k1,1 -2,2n` Sort the peakfile by genomic position
        - `bedtools merge -i -` Merge peaks that overlap with each other
10. Create required .saf file to generate count matrix [awk](https://www.gnu.org/software/gawk/manual/gawk.html)
    1. `awk  '{{OFS="\t";print $1"."$2+1"."$3, $1, $2+1, $3, "."}}' {input} > {output}`
        - Take input file and print required fields to new file
11. Generate count matrix from consensus peakset: [featureCounts](https://www.rdocumentation.org/packages/Rsubread/versions/1.22.2/topics/featureCounts)
    1. [scripts/featureCounts.R](https://github.com/bigmonty12/BASS/blob/master/scripts/featureCounts.R)
12. Find differential peaks between conditions, PCA, plots: [deseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
    1. [scripts/deseq2-init.R](https://github.com/bigmonty12/BASS/blob/master/scripts/deseq2-init.R)
    2. [scripts/deseq2.R](https://github.com/bigmonty12/BASS/blob/master/scripts/deseq2.R)
    3. [scripts/plot-pca.R](https://github.com/bigmonty12/BASS/blob/master/scripts/plot-pca.R)
13. Annotate differential peaks: [Homer](http://homer.ucsd.edu/homer/)
    1. `annotatePeaks.pl`
14. Create QC file for raw read, peak calling, and other results [MultiQC](https://github.com/ewels/MultiQC)

## How to Run

i. Navigate to pipeline directory

ii. Make sure snakemake is available for analysis (Through conda or other download)

iii. Run:
```
snakemake --use-conda --cores 20 --directory /chosen/project/directory
```
iv. To generate a report of analysis:
```
snakemake --report chosen_project-report.html --directory /chosen/project/directory
```

## Current Issues
- Sometimes the featureCounts.Rds object isn't formatted correctly. Current fix is to delete the object and rerun the pipeline from featureCounts step
- I'm not sold on the current peak-calling method. I am evaluating other methods (HMMRATAC, Genrich, other macs2 parameters) to determine the best.
- There are additional checks and methods I plan to add (also might switch over to nextflow)
