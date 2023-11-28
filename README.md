# SeqReadsProcessing (SRP)
Tools for Next Generation Sequencing reads processing suitable for both short Illumina and PacBio NanoPore long reads

These tools are intitially developed by Xuewen Wang. Free to use for academic research and education. Any any other application, please require for license.


# Installation

`git clone https://github.com/XuewenWangUGA/SeqReadsProcessing/`

`cd SeqReadsProcessing `

# Dependency 
The Biopython is needed for some scripts. To install Biopython

`pip install biopython`

For more details, visit [Biopython](https://biopython.org/wiki/Download)

# Extracting reads from sequence IDs
seqExtract_fromID.py

FUNCTION:

This script will take a list of IDs from a file and then extract the ID associated reads from a file containing  a lots of reads/sequences in fasta format. The hash or dictionary algorithm is used so it is very fast. The minimum memory of computer should be bigger than the size of the read file.

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


# NGS reads length filtering and read statistical summary
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



# Get summary of bp for fastq.gz

filterReadsfastq_gz.py

FUNCTION:

This script will take the next generation sequencing reads in fastq.gz format, and then will filter out the sequence less than given length threshold of base pair.

`python3 scripts/filterReadsfastq_gz.py out.fastq merged.fastq.gz 1`


# get the summary for a genome assembly file

AssemblyStat.py

This script will genome assembly in .fasta or .fasta.gz format, and then will report summary for this assembly. This tool also can filter out the unwanted reads with length less than the cut off value.

`python3 AssemblyStat.py OUT_file_in_fasta IN_seq_file_in_fasta min_length_cutoff`

e.g.: python3 AssemblyStat.py ./testseq.filtered.fasta ./testseq.fasta 1000

for help:

`python3 AssemblyStat.py`

output: a filtered fasta file and statistical summary: 

    Total number of input sequences/reads:   334
    Total length (bp) of input sequences/reads:      493088611
    Total number of sequences/reads after filtering:         334
    Total length (bp) of filtered sequences/reads:   493088611
    N50 length (bp): 11300744
    L50 count of sequences : 15
    
    N80 length (bp): 4697535
    L80 count of sequences : 33
    
    N90 length (bp): 1254627
    L90 count of sequences : 56




# Clean reads with unambiguous base N
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


# Fastq to fasta convertion 
fq2fa.py

FUNCTION:
This script will take the next generation sequencing reads, e.g. PacBio SMART sequencing reads, from the input file in fastq format
Then will convert sequence into fasta format and write in the output file.

Parameters: outseqf in_seqfile

USAGE: 

`python3 fq2fa.py OUT_file IN_seq_file`

Example on testing data:

`python3 fq2fa.py  testseq.fasta testdata_illumina_1.fq`

Path information could be added before the input and output file name if files are not in current directroy.

Version 2021-June 18th

# Sequence length tool

calculate the lengths of fasta sequences from an input file and output to a result file

sequence_length.py is a tool to calculate the sequence(s) length(s) in a input file in fasta format.

USAGE:

 `python3 sequence_length.py  OUT_file IN_seq_file length_cutoff_in_bp`
 
e.g., `python3 sequence_length.py  Path/length_out_File_name Path/testseq.fasta 1000`

details in help:

python3 sequence_length.py -h
