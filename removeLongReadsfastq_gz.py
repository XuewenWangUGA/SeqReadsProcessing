#!/usr/bin/python
import sys
import gzip
import shutil

from Bio import SeqIO

'''
FUNCTION:
This script will take the next generation sequencing reads, e.g. PacBio SMART sequencing reads, from the input file in fasta format
Then will filter/remove the sequence longer than the cutoff value of minimal sequence length in bp. The reads shorter than this cutoff will be output to the result file.
Parameters: outseqf in_seqfile Max_length_cutoff
USAGE: 
`python3 RemoveLongReadsfastq_gz.py OUT_file IN_seq_file Max_length_cutoff`
Example on testing data:
`python3 RemoveLongReadsfastq_gz.py ./testseq.filtered.fasta ./testseq.fasta 5000`
Path information could be added before the input and output file name if files in compressed format .gz are not in current directroy.
Version 2024-May 26th
# use python version 3
'''
__author__ ="Xuewen Wang"


def gzip_file(source_file_path, dest_file_path):
    with open(source_file_path, 'rb') as f_in:
        with gzip.open(dest_file_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

usage=""" python3 removeLongReadsfastq_gz.py OUT_file_in_fastq.gz IN_seq_file_in_fastq.gz Max_length_cutoff
          e.g.: python3 removeLongReadsfastq_gz.py ./testseq.filtered.fastq.gz ./testseq.fastq.gz 5000
          """
print(usage)

lengths = []
sumlen = 0
sumraw=0
ct=0
ctgood=0



#outf = open(sys.argv[1], 'w')
with gzip.open(sys.argv[1], 'wt') as outf:
    with gzip.open(sys.argv[2], "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            ct +=1
            sumraw += len(record.seq)
            if len(record.seq) < int(sys.argv[3]):
                ctgood +=1
                lengths.append(len(record.seq))
                sumlen += len(record.seq)
                SeqIO.write(record, outf, "fastq")

#'Report stats'
n = sumlen/2
print("Total number of input sequences/reads:\t", ct)
print("Total length (bp) of input sequences/reads:\t", sumraw)
print("Total number of sequences/reads after filtering:\t", ctgood)
print("Total length (bp) of filtered sequences/reads:\t", sumlen)

'Count the N50 length'
total = 0
lengths.sort(reverse = True)
for x in range(len(lengths)):
     total += lengths[x]
     if total >= n:
          print("N50 length (bp): %i" % (lengths[x-1]))
          break 

