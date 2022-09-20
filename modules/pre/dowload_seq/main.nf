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

/* Pre-DOWLOAD_SEQUENCES */

process DOWLOAD_SEQUENCES {
	tag "$Query_list"

	publishDir "${results_dir}/faa_files/",mode:"copy"

	input:
	each Query_list
	each Python_script

	output:
	file "*.faa"

	shell:
	"""
  python3 ${Python_script} ${Query_list}

	"""
	stub:
	"""
	      touch  ${Query_list}.faa
	"""
}
