#!/usr/bin/python
import sys

from Bio import SeqIO

'''
FUNCTION:
This script will take the next generation sequencing reads, e.g. PacBio SMART sequencing reads, from the input file in fasta format
Then will filter the sequence based on the cutoff value of minimal sequence length in bp. The reads longer than this cutoff will be output to the result file.
Parameters: outseqf in_seqfile length_cutoff
USAGE: 
`python3 filterReads.py OUT_file IN_seq_file length_cutoff`
Example on testing data:
`python3 filterReads.py ./testseq.filtered.fasta ./testseq.fasta 1000`
Path information could be added before the input and output file name if files are not in current directroy.
Version 2020-Feb 10th
# use python version 3
'''
__author__ ="Xuewen Wang"

usage=""" #python3 AssemblyStat.py OUT_file_in_fasta IN_seq_file_in_fasta min_length_cutoff
          e.g.: python3 AssemblyStat.py ./testseq.filtered.fasta ./testseq.fasta 1000 
          developed by PhD Xuewen Wang @ https://github.com/XuewenWangUGA, V1.1, 2025
          """
print(usage)

lengths = []
sumlen = 0
sumraw=0
ct=0
ctgood=0
maxlen=0
minlen=0
total_gc = 0
totalNcount=0
n_count=0
total_bases=0
overall_gc=0.00
length10k=10000
length10kCt=0
length1m=1000000
length1mCt=0

#outf = open(sys.argv[1], 'w')
for record in SeqIO.parse(sys.argv[2], "fasta"):
    ct += 1
    seq = record.seq.upper()
    length = len(seq)
    sumraw += length

    if ct==1:
         minlen=length
    elif length <minlen:
         minlen=length


    gc_count = seq.count("G") + seq.count("C")
    n_count=seq.count("N")
    gc_percent = 100 * gc_count / length if length > 0 else 0
    totalNcount += n_count
    total_gc += gc_count
    #print(f"{record.id}: {gc_percent:.2f}% GC")

    if length > length10k:
        length10kCt +=1
    if length >length1m:
        length1mCt +=1

    if length > maxlen:
        maxlen = length
    if length > int(sys.argv[3]):
        ctgood += 1
        lengths.append(length)
        sumlen += length
        #SeqIO.write(record, outf, "fasta")

#outf.close()

if sumraw > 0:
    overall_gc = 100 * total_gc / sumraw
    #print(f"\nOverall GC content: {overall_gc:.2f}%")


'Report stats'
n = sumlen/2
print("Total number of contigs/sequences:\t", ct)
print("Total length (bp) of input contigs/sequences:\t", sumraw)
print(f"Average GC content%:\t {overall_gc:.2f}")
#print("Total count of Ns:\t", totalNcount)
print("The longest sequence bp:\t", maxlen)
print("The shortest sequence bp:\t", minlen)
print(F"Sequence count longer than {length10k}:\t",length10kCt)
print(F"Sequence count longer than {length1m}:\t",length1mCt)


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
          print("N50 length (bp):\t%i" % (lengths[x-1]))
          print("L50 count of sequences:\t%i\n" %seqcout)
          N50=False
     
     'N80 summary'
     if (total >= sumlen*0.8 and N80):
          print("N80 length (bp):\t%i" % (lengths[x-1]))
          print("L80 count of sequences:\t%i\n" %seqcout)
          N80=False
     
     'N90 summary'
     if (total >= sumlen*0.9 and N90):
          print("N90 length (bp):\t%i" % (lengths[x-1]))
          print("L90 count of sequences:\t%i\n" %seqcout)
          N90=False
          break

print()
print("Total number of contigs/sequences after filtering:\t", ctgood)
print("Total length (bp) of after filtering:\t", sumlen)
print(f"Length filter cutoff bp:\t{sys.argv[3]}")
print()




