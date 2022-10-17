#!/usr/bin/python
import sys

from Bio import SeqIO

'''
FUNCTION:
This script will take the sequencing in fasta format, e.g. genome_assembly.fa;
Then this tool will filter the sequence based on the cutoff value of minimal sequence length in bp. The reads longer than this cutoff will be used to calculate the length. The output file is a plain text file with results of the length of each sequence  with sequence ID , length in two collumns seprately by tab. The summary information of whole sequene input file will also be generated into terminal console.
Parameters: out_file_name in_seqfile length_cutoff
USAGE: 
`python3 sequence_length.py OUT_fileName IN_seq_file length_cutoff`
Example on testing data:
`python3 sequence_length.py ./out_lengthFile_name ./testseq.fasta 1000`
Path information could be added before the input and output file name if files are not in current directroy.
Version 2022-Oct 17th
# use python version 3
'''
__author__ ="Xuewen Wang"

usage=""" python3 sequence_length.py OUT_file IN_seq_file_in_fasta length_cutoff
          e.g.: python3 sequence_length.py ./out.length.txt ./testseq.fasta 1000 
          """
print(usage)

lengths = []
sumlen = 0
sumraw=0
ct=0
ctgood=0

outf = open(sys.argv[1], 'w')
outf.write("#ID\tLength(bp)\n")
for record in SeqIO.parse(sys.argv[2], "fasta"):
     ct +=1
     mylen=int(len(record.seq))
     sumraw += mylen
     if int(mylen) > int(sys.argv[3]):
          ctgood +=1
          lengths.append(mylen)
          sumlen += mylen
          out=record.id+"\t"+str(mylen)+"\n"
          outf.write(out)

'Report stats'
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
          
outf.close()
