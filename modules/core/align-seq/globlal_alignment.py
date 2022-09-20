#!/usr/bin/env python3

""" This scrip realize a global aligment between FASTA query sequences and
    FASTA reference sequence file

The script takes two input files:
    1) A FASTA file (.fa.ref) with the reference sequences
    2) A FASTA file (.fa) with query sequences

Basic ideas:
    1. Grab FASTA files.
    2. Read through FASTA files, getthing one sequence at time.
    3. Run global aligment
    5. Calculate the score for each aligment.
    6- Write ouput as TSV file

Autor:  Eduardo García López"""

## Import python libraries
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.Align import substitution_matrices
import re
import sys

#Import blosum62 matrix
blosum62 = substitution_matrices.load("BLOSUM62")


## Read args from command line
    ## Uncomment For debugging only
    ## Comment for production mode only
#sys.argv = ("0", "test/data/test.fa.ref", "test/data/test_query.fa", "test/data/sample_out.tsv")



##get reference FASTA file from args
reference_fa = sys.argv[1]

##get the query FASTA file from args
query_fasta = sys.argv[2]

##get ouput file name from args
Output_file = sys.argv[3]

#Get query strain name
strain_ID = query_fasta.split('.fa')[0]


## Define FASTA sequences as lists
list_query = list(SeqIO.parse(str(query_fasta), "fasta"))
list_reference = list(SeqIO.parse(str(reference_fa), "fasta"))

#Get the lengt of FASTA lists
lengt_query = (len(list_query)) - 1 # Get the number of query sequences
lengt_reference = (len(list_reference)) - 1 # Get the number of reference faa sequences

## Write output files
results = open(str(Output_file),'w') # Ouput TSV dataframe

## Print output file header
results.write("referenceID\tstrain_ID\tquery_seqID\tquery_description\tidentity\tscore\tpos\tchanges\n")

#Print number of seqnces
print ("\n-------------------\n")
print ("Query sequences {}".format(lengt_query + 1))
print ("Reference sequences {}".format(lengt_reference + 1))
print ("\n--------------------\n")

#Define function to calculate identity between sequences
def calculate_identity(sequenceA, sequenceB):
        """
        Returns the percentage of identical characters between two sequences.
        Assumes the sequences are aligned.
        """
        # Get identical matches between two sequences
        identical_matches = pairwise2.align.globalxx(sequenceA, sequenceB, score_only=True)
        # Get global aligment with no gap penalties, identical characters have score of 1, otherwise 0.
        align_xx = pairwise2.align.globalxx(sequenceA, sequenceB)
        align = pairwise2.format_alignment(*align_xx[0]).split("\n")
        len_align = len(str(align[0]))
        indentity =  identical_matches /  len_align *100
        return (indentity)

reference = 0 # Define number of reference sequence to parse data

# Loop to align each query sequence with each reference seq in FAA file
while reference <= lengt_reference:
    query = 0 # Define number of query sequence to parse data
    reference_id = list_reference[reference].id # Get reference ID
    reference_seq = list_reference[reference].seq # Get reference sequence
# Loop through sequences
    while query <= lengt_query:
        query_id = list_query[query].id # Get query ID
        query_description = list_query[query].description
        query_seq = list_query[query].seq # Get query sequence
        # Print iterator message
        print (f"Aligning {reference_id} and {query_id}, query {query + 1} of {lengt_query + 1 }, ref {reference+1} of {lengt_reference+1}\n")
        # Get indentity
        indentity = calculate_identity(reference_seq, query_seq)
        #Align sequences with BLOSUM62 matrix
        alignment_b62 = pairwise2.align.globalds(reference_seq, query_seq, blosum62, -10, -0.5,  score_only=True)
        # Get global aligment with no gap penalties, identical characters have score of 1, otherwise 0.
        align_b62_f = pairwise2.align.globalds(reference_seq, query_seq, blosum62, -10, -0.5)
        align = pairwise2.format_alignment(*align_b62_f[0]).split("\n")

        # Save alignment_b62
        ref = str(align[0]) # Ref alignment sequence
        sequence = str(align[1])  # Notation alignment sequence ( | , . or " ")
        alt = str(align[2])     # Alt alignment sequence
        # Find all occurrences of mismatch in query sequence
        mismatch = '\\.'
        variants = [i.start() for i in re.finditer(mismatch, sequence)]
        # Find deletions on query sequence
        space = ' '
        deletions = [i.start() for i in re.finditer(space, sequence)]
        # Join list variants and deletions
        positions = variants + deletions
        if len(positions) > 0:
            #Annotate sequence changes
            changes = str([f"{ref[variant]}>{alt[variant]}" for variant in positions]).replace("[", "").replace("]", "")
            #Remove list characters
            res_pos = str(variants + deletions).replace("[", "").replace("]", "")
        elif len(positions) == 0:
            changes = "."
            res_pos = "."
        #Write results
        results.write(f"{reference_id}\t{strain_ID}\t{query_id}\t{query_description}\t{indentity}\t{alignment_b62}\t{res_pos}\t{changes}\n")

        query += 1
    reference += 1

## Close files
results.close()
