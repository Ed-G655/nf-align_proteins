#!/usr/bin/env python3

""" This scrip dowload a list of FASTA files from NCBI assembly using NCBI Edirect tools.

-The script takes a list with Organism IDs

Basic ideas:
    1. Read list with the IDs
    2. Parse list to dowload each FASTA

Autor:  Eduardo García López"""

## Import python libraries
import sys
import subprocess


## Read args from command line
    ## Uncomment For debugging only
    ## Comment for production mode only
#sys.argv = ("0", "test.txt")

##get IDs list
ID_list = sys.argv[1]

#Read file as list
ID_sequence = str(ID_list)

#content_list = my_file.readlines()
print(ID_list)


#Define output name
outName = ID_sequence.replace('\n', '') + ".faa"

print ("Dowload {} at {}".format(ID_sequence.replace('\n', ''), outName))
#Dedine command to convert BAM files
esearch = "esearch -db assembly -query {} | elink -target nuccore -name \
        assembly_nuccore_insdc | elink -target protein | efetch -format fasta \
        > {}".format(ID_sequence.replace('\n', ''), outName)
print ("The command used was: " + esearch)
#Pass command to shell
subprocess.call(esearch, shell=True)
