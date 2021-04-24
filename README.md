# SeqReadsProcessing (SRP)
Tools for short and long next generation sequencing reads processing

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

Path information could be added before the input and output file name if files are not in current directory.


# Filter out reads with unambiguous base N
filterFastqNReads.py

FUNCTION:
This script will take the next generation sequencing reads in fastq format, and then will filter out the sequence which has unknown base pair Ns. 
          The reads without unambiguous Ns will be kept in the result file. This will retrieve the perfect reads without missing data.

USAGE: 

`filterFastqNReads.py [-h] -i INPUT [-o OUTPUT]`

optional arguments:

  -h, --help          show this help message and exit
  
  -i INPUT, --input INPUT
                      The input file containing reads in fastq format
                        
  -o OUTPUT, --output OUTPUT
                      The output file to save reads without base N in fastq format
  
   
Example on testing data:

python3 filterFastqNReads.py -i testdata_illumina_1.fq -o testdata_illumina_1.filteredN.fq
