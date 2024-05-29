#!/usr/bin/python
import sys
import gzip

from Bio import SeqIO

'''
FUNCTION:
This script will take the next generation sequencing reads, e.g. PacBio SMART sequencing reads, from the input file in fasta format
Then will filter the sequence based on the cutoff value of minimal sequence length in bp. The reads longer than this cutoff will be output to the result file.
Parameters: outfileName in_seqfile.fastq.gz length_cutoff
USAGE: 
Example on testing data:
`python3 filterReadsReportLengthfastq_gz ./testseq.lengths.txt ./testseq.fastq 1000` 
Path information could be added before the input and output file name if files are not in current directroy.
Version 2024-April 1st
# use python version 3
'''
__author__ ="Xuewen Wang"

usage=""" python3 filterReadsReportLengthfasta_gz.py OUT_seqId_length.txt IN_seq_file_in_fastq.gz length_cutoff
          e.g.: python3 filterReadsReportLengthfastq_gz.py ./testseq.lengths.txt ./testseq.fastq.gz 1000 
          """
#print(usage)

seqIdLen = {}
sumlen = 0
sumraw=0
ct=0
ctgood=0

outf = open(sys.argv[1], 'w')
with gzip.open(sys.argv[2], "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        ct +=1
        sumraw += len(record.seq)
        seqlen=int(len(record.seq))
        if seqlen > int(sys.argv[3]):
            ctgood +=1
            seqid=record.id            
            #seqIdLen[seqid]=seqlen
            sumlen += seqlen
            outinfo=seqid+"\t"+str(seqlen)
            #print(outinfo)
            #SeqIO.write(record, outf, "fastq")
            outf.write(outinfo+"\n")

