#!/usr/bin/python
import sys
import gzip
from Bio import SeqIO

'''
FUNCTION:
This script will take the next generation sequencing reads, e.g. PacBio SMART sequencing reads, from the input file in fasta format
Then will filter the sequence based on the cutoff value of minimal sequence length in bp. The reads longer than this cutoff will be output to the result file.
Parameters: outseqf in_seqfile length_cutoff
USAGE: 
`python3 AssemblyStat.py OUT_file IN_seq_file length_cutoff`
Example on testing data:
`python3 AssemblyStat.py ./testseq.filtered.fasta ./testseq.fasta 1000`
Path information could be added before the input and output file name if files are not in current directroy.
Version 2020-Feb 10th
# use python version 3
'''
__author__ ="Xuewen Wang"

usage=""" python3 AssemblyStat.py OUT_file_in_fasta IN_seq_file_in_fasta.gz min_length_cutoff
          e.g.: python3 AssemblyStat.py ./testseq.filtered.fasta ./testseq.fasta 1000 
          """
print(usage)

lengths = []
sumlen = 0
sumraw=0
ct=0
ctgood=0

outf = open(sys.argv[1], 'w')
with gzip.open(sys.argv[2],"rt") as IN:
 for record in SeqIO.parse(IN, "fastq"):
 #for record in SeqIO.parse(sys.argv[2], "fasta"):
     ct +=1
     sumraw += len(record.seq)
     if int(len(record.seq)) > int(sys.argv[3]):
          ctgood +=1
          lengths.append(int(len(record.seq)))
          sumlen += int(len(record.seq))
          SeqIO.write(record, outf, "fasta")

'Report stats'
n = sumlen/2
print("Total number of input sequences/reads:\t", ct)
print("Total length (bp) of input sequences/reads:\t", sumraw)
print("Total number of sequences/reads after filtering:\t", ctgood)
print("Total length (bp) of filtered sequences/reads:\t", sumlen)

'Count the N50 length'
total = 0
seqcout=0
N50=True
N80=True
N90=True

lengths.sort(reverse = True)
for x in range(len(lengths)):
     total += lengths[x]
     seqcout += 1
     if (total >= n and N50):
          print("N50 length (bp): %i" % (lengths[x-1]))
          print("L50 count of sequences : %i\n" %seqcout)
          N50=False
     
     'N80 summary'
     if (total >= sumlen*0.8 and N80):
          print("N80 length (bp): %i" % (lengths[x-1]))
          print("L80 count of sequences : %i\n" %seqcout)
          N80=False
     
     'N90 summary'
     if (total >= sumlen*0.9 and N90):
          print("N90 length (bp): %i" % (lengths[x-1]))
          print("L90 count of sequences : %i\n" %seqcout)
          N90=False
          break





