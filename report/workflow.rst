**Differential ATAC-seq peaks were found with the following workflow:**
    - Reads were trimmed w/ `TrimGalore`_
    - Reads were mapped onto {{ snakemake.config["ref"]["name"] }} w/ `Bowtie2`_
    - PCR duplicates were removed w/ Picard_
    - Replicates were merged w/ `Samtools`_
    - Peaks were called on each replicate and merged file w/ macs2_
    - Counts matrix was created w/ R package of `featureCounts`_
    - Differential analysis was done w/ `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_
    - Quality control was performed w/ FastQC_, Samtools_, and Picard_ and then aggregated via MultiQC_

.. _TrimGalore: https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
.. _Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
.. _Picard: https://broadinstitute.github.io/picard
.. _Samtools: http://samtools.sourceforge.net/
.. _macs2: https://github.com/taoliu/MACS
.. _featureCounts: http://bioinf.wehi.edu.au/featureCounts/
.. _FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _MultiQC: http://multiqc.info/
