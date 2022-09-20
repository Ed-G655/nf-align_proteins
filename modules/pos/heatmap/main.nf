#!/usr/bin/env nextflow
/*================================================================

								---- MODULE PIPELINE ---------

/*================================================================
The Aguilar Lab presents...

- A pipeline to realize a global aligment between FASTA query sequences and
    FASTA reference sequence file

==================================================================
Version: 0.1
Project repository:
==================================================================
Authors:

- Bioinformatics Design
 Jose Eduardo Garcia-Lopez (jeduardogl655@gmail.com)



- Bioinformatics Development
 Jose Eduardo Garcia-Lopez (jeduardogl655@gmail.com)


- Nextflow Port
 Jose Eduardo Garcia-Lopez (jeduardogl655@gmail.com)

///////////////////////////////////////////////////////////////

  Define pipeline Name
  This will be used as a name to include in the results and intermediates directory names
*/
pipeline_name = "nf-align_proteins.nf"

/*This directories will be automatically created by the pipeline to store files during the run
*/
results_dir = "${params.output_dir}/${pipeline_name}-results/"
intermediates_dir = "${params.output_dir}/${pipeline_name}-intermediate/"

/*================================================================/*

/* MODULE START */

/* POS-PLOT_HEATMAP_IDENTITY */

process HEATMAP {
	errorStrategy 'ignore'
	tag "$TSV"

	publishDir "${results_dir}/heatmap-plot/",mode:"copy"

	input:
	file TSV
	file R_script

	output:
	file "*"

	shell:
	"""
  Rscript --vanilla ${R_script}

	"""
	stub:
	"""
	      touch  ${TSV}.png
	"""
}
