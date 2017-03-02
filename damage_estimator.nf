#! /usr/bin/env nextflow

// IARC - C. VOEGELE - Last update 02 March 2017
// Cf. https://github.com/Ettwiller/Damage-estimator
// Usage: ./damage-estimator.nf --bam_folder BAM/ --de_path PATH/ --genome_ref ref.fasta

// Set parameters
params.help = null
params.cluster_options = ""
params.de="/appli57/damage_estimator/"
params.genome_ref ="/appli57/reference/GATK_Bundle/ucsc.hg19.fasta"
params.Q ="0"
params.mq = "10"
params.max_coverage_limit = "100"
params.min_coverage_limit = "1"
params.qualityscore ="30"

if (params.help) {
    log.info ''
    log.info '--------------------------------------------------'
    log.info 'NEXTFLOW: DAMAGE ESTIMATOR'
    log.info '--------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run iarcbioinfo/damage-estimator.nf --bam_folder BAM/ --de_path PATH/ --genome_ref ref.fasta'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --bam_folder   FOLDER                    Folder containing .bam and .bam.bai files'
    log.info '    --de_path	     FOLDER                    Folder containing Damage Estimator files (.pl, .r)'
    log.info '    --genome_ref   FILE                      Genome of reference file'
    log.info 'Options:'
    log.info '    --cluster_options	STRING		   Specific options for cluster scheduler'
    log.info '    --Q			NUMBER		   Phred score quality threshold'
    log.info '    --mq			NUMBER		   Mapping quality'
    log.info '    --max_coverage_limit	NUMBER		   Maximum coverage limit'
    log.info '    --min_coverage_limit	NUMBER		   Minimum coverage limit'
    log.info '    --qualityscore		NUMBER		   Quality score'
    log.info ''
    exit 1
}

// Check folder with bam files + bai files + pair them

bams = Channel.fromPath( params.bam_folder+'/*.bam' )
              .ifEmpty { error "Cannot find any bam file in: ${params.bam_folder}" }
    		.map {  path -> [ path.name.replace("bam",""), path ] }	
bais = Channel.fromPath( params.bam_folder+'/*.bam.bai' )
    		.map {  path -> [ path.name.replace("bam.bai",""), path ] }		

b_pair = bams
	.phase (bais)
	.map { pair1, pair2 -> [ pair1[1], pair2[1] ] }


////////// STEP 01 ################### Launch split mapped reads - takes between 3 and 4 hours
////////////////////////////////////// Parameter Q (OPTIONAL): Phred score quality threshold (Sanger encoding). Only keep the bases with a Q score above a given threshold (default 0).
////////////////////////////////////// Parameter q (OPTIONAL): mapping quality. Only keep the reads that passes a given threshold (default 10).

process split_mapped_reads {
	echo "split_map_reads"
	tag { bam_tag }
	memory '12GB'
	clusterOptions = params.cluster_options
	publishDir 'DE_output_smr'

input:
	file pair from b_pair
output:
	set val(bam_tag), file("${bam_tag}_file1") into smr_1
 	set val(bam_tag), file("${bam_tag}_file2") into smr_2
shell:
	bam_tag = pair[0].baseName
     	'''
	perl !{params.de}split_mapped_reads.pl --bam !{pair[0]} -genome !{params.genome_ref} -mpileup1 !{bam_tag}_file1 -mpileup2 !{bam_tag}_file2 -Q !{params.Q} -q !{params.mq}
       '''
}

////////// STEP 02 ################### Launch estimate damage
////////////////////////////////////// Parameter --max_coverage_limit (DEFAULT 100) : If a position has equal or more than MAX reads (R1 or R2), the position is not used to calculate the damage. This option is put in place in order to avoid high coverage regions of the genome being the main driver for the damage estimation program.
////////////////////////////////////// Parameter --min_coverage_limit (DEFAULT 1) : If a position has equal or less than MIN reads (R1 or R2), the position is not used to calculate the damage. This option is put in place in order to calculate damage only in on-target regions (in cases of enrichment protocol such as exome ....)
////////////////////////////////////// Parameter --qualityscore (DEFAULT 30): Discard the match or mismatch if the base on a read has less than MIN base quality. Important parameters. The lower this limit is, the less the damage is apparent.

process estimate_damage {
	echo "estimate_damage"
	tag { bam_tag }
	memory '12GB'
	clusterOptions = params.cluster_options
	publishDir 'DE_output_tables'

input:
	set val(bam_tag), file("${bam_tag}_file1") from smr_1
 	set val(bam_tag), file("${bam_tag}_file2") from smr_2
output:
	set val(bam_tag), file("${bam_tag}_estimate_damage.tab") into ed

shell:
     	'''
	perl !{params.de}estimate_damage.pl --mpileup1 !{bam_tag}_file1 --mpileup2 !{bam_tag}_file2 --id !{bam_tag} --qualityscore !{params.qualityscore} --max_coverage_limit !{params.max_coverage_limit} --min_coverage_limit !{params.min_coverage_limit} > !{bam_tag}_estimate_damage.tab
       '''
}

////////// STEP 03 ################### Generate plot

process generate_plot {
	echo "generate_plot"
	tag { bam_tag }
	clusterOptions = params.cluster_options
	publishDir 'DE_output_plots', mode: 'move'

input:
	set val(bam_tag), file("${bam_tag}_estimate_damage.tab") from ed
output:
	set val(bam_tag), file("${bam_tag}_estimate_damage.png") into ed_plot

shell:
     	'''
	Rscript !{params.de}plot_damage.R !{bam_tag}_estimate_damage.tab !{bam_tag}_estimate_damage.png
       '''
}
