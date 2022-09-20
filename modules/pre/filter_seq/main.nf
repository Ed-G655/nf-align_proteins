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

/* Pre-FILTER_SEQ */

process FILTER_SEQ {
	tag "$FAA_FILE"

	publishDir "${results_dir}/fim_fastas/",mode:"copy"

	input:
	each FAA_FILE

	output:
	file "*.fa"

	shell:
	"""
   grep "fimbrial" ${FAA_FILE} | tr -d ">" >  ${FAA_FILE.baseName}_fim.txt
	 seqtk subseq ${FAA_FILE} ${FAA_FILE.baseName}_fim.txt > ${FAA_FILE.baseName}.fa
	"""
	stub:
	"""
	      touch  ${FAA_FILE.baseName}.fa
	"""
}
