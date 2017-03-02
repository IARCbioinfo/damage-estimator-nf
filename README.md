# damage-estimator-nf
Nextflow pipeline to run "Damage Estimator"
Cf. https://github.com/Ettwiller/Damage-estimator

# Tools needed to be installed prior to run
- samtools
- R with GGPLOT2 package

# Mandatory parameters
- bam_folder: folder containing .bam and .bam.bai files on which to run "Damage Estimator"
- de_path : location of folder containing damage estimator files (.pl and .r)
- genome_ref : the genome of reference

# Optional parameters
- cluster_options
- Q : Phred score quality threshold (Sanger encoding). Only keep the bases with a Q score above a given threshold (default 0)
- q : mapping quality. Only keep the reads that passes a given threshold (default 10)
- max_coverage_limit MAX (DEFAULT 100) : If a position has equal or more than MAX reads (R1 or R2), the position is not used to calculate the damage. This option is put in place in order to avoid high coverage regions of the genome being the main driver for the damage estimation program.
- min_coverage_limit MIN (DEFAULT 1) : If a position has equal or less than MIN reads (R1 or R2), the position is not used to calculate the damage. This option is put in place in order to calculate damage only in on-target regions (in cases of enrichment protocol such as exome ....)
- qualityscore MIN : Discard the match or mismatch if the base on a read has less than MIN base quality. Important parameters. The lower this limit is, the less the damage is apparent.

For exome bams, we recommend: -Q 20 -q 20 --max_coverage_limit 300 --min_coverage_limit 30

# Example:
nextflow run iarcbioinfo/damage-estimator-nf --bam_folder BAM/ --de_path /path/ --genome_ref ref.fasta 
