#!/usr/bin/env nextflow

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

Pre-processing:


Core-processing:
Core-align_seq

Pos-processing
Pos1-filter_hits
Pos2-plot_hits

Analysis:

================================================================*/

/* Define the help message as a function to call when needed *//////////////////////////////
def helpMessage() {
	log.info"""
  ==========================================
							${pipeline_name}

	The global aligment pipeline
  v${version}
  ==========================================

	Usage:

	nextflow run ${pipeline_name}.nf --ref_fa <path to input 1 > --query_fa <path to input 2>

		--ref_fa	<-  Directory whith reference FASTA files;

	  --query_fa	<- Directory whith with query FASTA files;

	  --output_dir     <- directory where results, intermediate and log files will be stored;
	      default: same dir where --query_fasta resides

	  -resume	   <- Use cached results if the executed project has been run before;
	      default: not activated
	      This native NF option checks if anything has changed from a previous pipeline execution.
	      Then, it resumes the run from the last successful stage.
	      i.e. If for some reason your previous run got interrupted,
	      running the -resume option will take it from the last successful pipeline stage
	      instead of starting over
	      Read more here: https://www.nextflow.io/docs/latest/getstarted.html#getstart-resume
	  --help           <- Shows Pipeline Information
	  --version        <- Show version
	""".stripIndent()
}

/*//////////////////////////////
  Define pipeline version
  If you bump the number, remember to bump it in the header description at the begining of this script too
*/
version = "0.1"

/*//////////////////////////////
  Define pipeline Name
  This will be used as a name to include in the results and intermediates directory names
*/
pipeline_name = "nf-align_proteins"

/*
  Initiate default values for parameters
  to avoid "WARN: Access to undefined parameter" messages
*/
params.ref_fa = false  //if no inputh path is provided, value is false to provoke the error during the parameter validation block
params.query_fa = false  //if no inputh path is provided, value is false to provoke the error during the parameter validation block
params.help = false //default is false to not trigger help message automatically at every run
params.version = false //default is false to not trigger version message automatically at every run

/*//////////////////////////////
  If the user inputs the --help flag
  print the help message and exit pipeline
*/
if (params.help){
	helpMessage()
	exit 0
}

/*//////////////////////////////
  If the user inputs the --version flag
  print the pipeline version
*/
if (params.version){
	println "${pipeline_name} v${version}"
	exit 0
}

/*//////////////////////////////
  Define the Nextflow version under which this pipeline was developed or successfuly tested
  Updated by iaguilar at MAY 2021
*/
nextflow_required_version = '20.01.0'
/*
  Try Catch to verify compatible Nextflow version
  If user Nextflow version is lower than the required version pipeline will continue
  but a message is printed to tell the user maybe it's a good idea to update her/his Nextflow
*/
try {
	if( ! nextflow.version.matches(">= $nextflow_required_version") ){
		throw GroovyException('Your Nextflow version is older than Pipeline required version')
	}
} catch (all) {
	log.error "-----\n" +
			"  This pipeline requires Nextflow version: $nextflow_required_version \n" +
      "  But you are running version: $workflow.nextflow.version \n" +
			"  The pipeline will continue but some things may not work as intended\n" +
			"  You may want to run `nextflow self-update` to update Nextflow\n" +
			"============================================================"
}

/*//////////////////////////////
  INPUT PARAMETER VALIDATION BLOCK
*/

/* Check if the input directory is provided
    if it was not provided, it keeps the 'false' value assigned in the parameter initiation block above
    and this test fails
*/
if ( !params.ref_fa | !params.query_fa) {
  log.error " Please provide the --ref_fa AND --query_fa \n\n" +
  " For more information, execute: nextflow run ${pipeline_name}.nf --help"
  exit 1
}

/*
Output directory definition
Default value to create directory is the parent dir of --input_dir
*/
params.output_dir = file(params.ref_fa).getParent() //!! maybe creates bug, should check

/*
  Results and Intermediate directory definition
  They are always relative to the base Output Directory
  and they always include the pipeline name in the variable pipeline_name defined by this Script

  This directories will be automatically created by the pipeline to store files during the run
*/
results_dir = "${params.output_dir}/${pipeline_name}-results/"
intermediates_dir = "${params.output_dir}/${pipeline_name}-intermediate/"

/*
Useful functions definition
*/

/*//////////////////////////////
  LOG RUN INFORMATION
*/
log.info"""
==========================================
The ${pipeline_name} pipeline
v${version}
==========================================
"""
log.info "--Nextflow metadata--"
/* define function to store nextflow metadata summary info */
def nfsummary = [:]
/* log parameter values beign used into summary */
/* For the following runtime metadata origins, see https://www.nextflow.io/docs/latest/metadata.html */
nfsummary['Resumed run?'] = workflow.resume
nfsummary['Run Name']			= workflow.runName
nfsummary['Current user']		= workflow.userName
/* string transform the time and date of run start; remove : chars and replace spaces by underscores */
nfsummary['Start time']			= workflow.start.toString().replace(":", "").replace(" ", "_")
nfsummary['Script dir']		 = workflow.projectDir
nfsummary['Working dir']		 = workflow.workDir
nfsummary['Current dir']		= workflow.launchDir
nfsummary['Launch command'] = workflow.commandLine
log.info nfsummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "\n\n--Pipeline Parameters--"
/* define function to store nextflow metadata summary info */
def pipelinesummary = [:]
/* log parameter values beign used into summary */
pipelinesummary['Input ref FASTA dir']			= params.ref_fa
pipelinesummary['Input query FASTA dir']			= params.query_fa
pipelinesummary['Results Dir']		= results_dir
pipelinesummary['Intermediate Dir']		= intermediates_dir
/* print stored summary info */
log.info pipelinesummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "==========================================\nPipeline Start"

/*//////////////////////////////
  PIPELINE START
*/

/* enable DSL2*/
nextflow.enable.dsl=2

/*
	READ GENERAL INPUTS
*/

	/* Load ref FASTAs files into channel */

// Define a function to get protein prefix
def get_fasta_prefix = { file -> file.baseName.replaceAll(/.fa/, "") }

/* Load ref FASTAs files into channel */
Channel
	.fromPath( "${params.ref_fa}*.ref" )
	.map{ file -> tuple(get_fasta_prefix(file) , file) }
	.groupTuple()
	.set{ ref_input }

	/* Load query FASTAs files into channel */
Channel
	.fromPath( "${params.query_fa}*.fa" )
	.map{ file -> tuple(file.baseName , file) }
	.groupTuple()
	.set{ query_input }

/*
 Load Python fileS
*/
/* Load vcf_split.py script */
	Channel
		.fromPath( "./modules/core/align-seq/globlal_alignment.py" )
		.set{ Python_script_1 }
/*
 Load R fileS
*/


	/*	  Import modules */
											/* PRE-processing */
include {	ALIGN_SEQ	} from './modules/core/align-seq/main.nf'


/*  main pipeline logic */
workflow  {
/* PRE-processing */
						// join FASTA channels by prefix
					//	ref_input.view()
						// query_input.view()
						GROUP_BY_PREFIX = ref_input.join(query_input)
						// PRE 1: ALIGN FASTA SEQUECES
						ALIGN_SEQ(GROUP_BY_PREFIX, Python_script_1)

}
