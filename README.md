# damage-estimator-nf
## Nextflow pipeline to run "Damage Estimator"

## Description
Cf. https://github.com/Ettwiller/Damage-estimator

## Dependencies
1. Nextflow: for installation procedures see the file config. (also for cluster-options for your cluster scheduler)
2. samtools
3. R with GGPLOT2 package
The tool writes in tmp folder so check that yours is specified in your .bash_profile (export TMPDIR=/data/tmp/, export TMP=/data/tmp)

## Input
  | Type      | Description     |
  |-----------|---------------|
  | bam folder    | Folder containing the bam files on which you want to run "Damage Estimator"|

## Parameters

* #### Mandatory
| Name               | Example value | Description     |
|--------------------|---------------|-----------------| 
| --bam_folder       |            PATH/FOLDER | folder containing .bam and .bam.bai files on which to run "Damage Estimator" (bams should preferably be generated by bwa mapping of Illumina paired-end sequencing) |
| --de_path          |            PATH/DE | location of folder containing damage estimator files (.pl and .r) |
| --ref              |            PATH/FILE | genome of reference (fasta file) |

 * #### Optional
| Name                 | Default value | Description     |
|----------------------|---------------|-----------------| 
| --Q                  |            0 | Phred score quality threshold (Sanger encoding). Only keep the bases with a Q score above a given threshold |
| --mq                 |            10 | mapping quality. Only keep the reads that passes a given threshold |
| --max_coverage_limit |            100 | If a position has equal or more than MAX reads (R1 or R2), the position is not used to calculate the damage. This option is put in place in order to avoid high coverage regions of the genome being the main driver for the damage estimation program. |
| --min_coverage_limit |            1 | If a position has equal or less than MIN reads (R1 or R2), the position is not used to calculate the damage. This option is put in place in order to calculate damage only in on-target regions (in cases of enrichment protocol such as exome ....) |
| --qualityscore       |            30 | Discard the match or mismatch if the base on a read has less than MIN base quality. Important parameters. The lower this limit is, the less the damage is apparent. |

For exome bams, we recommend: --Q 20 --mq 20 --max_coverage_limit 300 --min_coverage_limit 30

## Usage 

```
nextflow run iarcbioinfo/damage-estimator.nf --bam_folder BAM/ --de_path /path/ --genome_ref ref.fasta
```

## Output

  | Type      | Description     |
  |-----------|---------------|
  | Graph    | ...... |
  | SMR    | ...... |
  | Table    | ...... |
  
## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------| 
  | VOEGELE Catherine    |            voegelec@iarc.fr | Developer|
  
