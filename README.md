# bulk-atac-snake

## Current Pipeline Steps
- Trim (TrimGalore!)
- Align (Bowtie2 | Samtools View)
- Filter (Samtools View)
- Remove Duplicates (Picard)
- Samtools Index
- Peak Calling (macs2)
- Create bigwig (deeptools)
- QC (Fastqc, multiqc, samtools stat)
## Developmental Plans
- Merge across conditions to find consensus peakset
- Convert consensus peakset to .saf
- Differential peak analysis (DEseq2)
- Footprinting analysis (rgt-hint?)
- Peak annotation (homer? ChIPSeeker?)
- Visualizations (PCA, FRiP, Peak Counts, other DEseq2 plots) 
