# SeqReadsProcessing (SRP)
Tools for Next Generation Sequencing reads processing suitable for both short Illumina and PacBio NanoPore long reads


# installation

`git clone https://github.com/XuewenWangUGA/SeqReadsProcessing/

cd SeqReadsProcessing `

# Dependency 
The Biopython is needed for some scripts. To install Biopython

`pip install biopython`
For more details, visit [Biopython](https://biopython.org/wiki/Download)

# Extracting reads from sequence IDs
seqExtract_fromID.py

FUNCTION:

This script will take a list of IDs from a file and then extract the ID associated reads from a file containing  a lots of reads/sequences in fasta format. The hash or dictionary algorithem is used so it is very fast. The minimum memory of computer should be bigger than the size of the read file.

USAGE:

`python3 seqExtract_fromID.py Options`

 Options: 
 
              -i input_idfile  The file name of id list, one sequence ID per line
 
              -s input_fasta_seqfile The sequence file containing all sequences in fasta format
 
              -o Outfile_extracted_seq The file name to save the extracted sequences in fasta format
 
               -h help
  
Testing example:

e.g. `python3 seqExtract_fromID.py -s testseq.fasta -i id.txt -o  testseq.fasta_extracted.fa`

 A statistical information will output to standout during run. e.g.
 
 Input id file is id.txt
 
Input sequence file:  testseq.fasta

Output seq file is testseq.fasta_extracted.fa

Warning: sequence is not available for IDs: >111 

Total # of id to bait:  3

Total # of extracted sequence:  2

Total # of not extracted sequence:      1

Result file with xxtracted sequence:    testseq.fasta_extracted.fa


# Long or short NGS reads length filtering
filterReads.py

FUNCTION:

This script will take the next generation sequencing reads, e.g. PacBio SMART sequencing reads, from the input file in fasta format
Then will filter the sequence based on the cutoff value of minimal sequence length in bp. The reads longer than this cutoff will be output to the result file.

Parameters: outseqf in_seqfile length_cutoff

USAGE: 

`python3 filterReads.py OUT_file IN_seq_file length_cutoff`

Example on testing data:

`python3 filterReads.py ./testseq.filtered.fasta ./testseq.fasta 1000`

Path information could be added before the input and output file name if files are not in current directory. A statistical summary will be provided for total reads, read length, N50 length before and after filtering.
e.g.
          
Total number of input sequences/reads:   7

Total length (bp) of input sequences/reads:      7313

Total number of sequences/reads after filtering:         3

Total length (bp) of filtered sequences/reads:   4558

N50 length (bp): 2201



# Filter out reads with unambiguous base N
filterFastqNReads.py

FUNCTION:
This script will take the next generation sequencing reads in fastq format, and then will filter out the sequence which has unknown base pair Ns. 
          The reads without unambiguous Ns will be kept in the result file while those with Ns will be discarded. This will retrieve the perfect reads without missing data.

USAGE: 

`filterFastqNReads.py [-h] -i INPUT [-o OUTPUT]`

optional arguments:

  -h, --help          show this help message and exit
  
  -i INPUT, --input INPUT
                      
                      The input file containing reads in fastq format
                        
  -o OUTPUT, --output OUTPUT
                      
                      The output file to save reads without base N in fastq format
  
   
Example on testing data:

`python3 filterFastqNReads.py -i testdata_illumina_1.fq -o testdata_illumina_1.filteredN.fq`

A statistical summary will be provided for total reads, read length, N50 length before and after filtering.



