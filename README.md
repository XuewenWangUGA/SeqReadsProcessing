# SeqReadsProcessing
Tools for short and long next generation sequencing reads processing

# Long or short reads filtering
filterReads.py

FUNCTION:
This script will take the next generation sequencing reads, e.g. from PacBio SMART sequencing, in the input file in fasta format
Then will filter the sequence based on the cutoff value of minimal sequence length in bp. The reads longer than this cutoff will be output to the result file.

Parameters: outseqf in_seqfile length_cutoff

USAGE: 

`python3 filterReads.py OUT_file IN_seq_file length_cutoff`

Example on testing data:

`python3 filterReads.py ./testseq.filtered.fasta ./testseq.fasta 1000`

Path information could be added before the input and output file name if files are not in current directroy.


Version 2020-Feb 10th
