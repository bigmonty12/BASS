# bulk-atac-snake

## Current Pipeline Steps
- Trim (TrimGalore!)
- Align (Bowtie2 | Samtools View)
- Filter (Samtools View)
- Remove Duplicates (Picard)
- Samtools Index
- Peak Calling (macs2)
- Create bigwig (deeptools)
- Merge across conditions to find consensus peakset
- QC (Fastqc, multiqc, samtools stat)
- Differential peak analysis (DEseq2)
- PCA

## Developmental Plans
- Footprinting analysis (rgt-hint?)
- Peak annotation (homer? ChIPSeeker?)
- Visualizations (Correlation/Heatmap, FRiP/PeakCounts, ReadLengthDistribution) 
