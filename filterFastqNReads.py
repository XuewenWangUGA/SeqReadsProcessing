#!/usr/bin/python
import sys
from Bio import SeqIO
import argparse

'''
filterFastqNReads.py
FUNCTION:
This script will take the next generation sequencing reads in fastq format. 
          Then it will filter out the sequence which has unknown base pair Ns. 
          The reads without N will be output to the result file. 
          This will retrieve the perfect reads without missing data.

HELP: python3 filterFastqNReads.py -h [--help]

USAGE: 

python3 filterFastqNReads.py -i input_fastq_file -o out_fastq_result_file

Example on testing data:

python3 filterFastqNReads.py -i testdata_illumina_1.fq -o testdata_illumina_1.filterN.fq

Path information could be added before the input and output file name if files are not in current directory.

Version 2021-April-10th

'''
__author__ ="Xuewen Wang"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",required=True, help="The input file containing reads in fastq format")
parser.add_argument("-o","--output", default='filtered.fq', help="The output file to save reads without base Ns in fastq format")
args=parser.parse_args()
ind=args.input
otd=args.output


lengths = []
sumlen = 0
sumraw=0
ct=0
ctgood=0
subword="N"

outf = open(otd, 'w')
for record in SeqIO.parse(ind, "fastq"):
     ct +=1
     sumraw += len(record.seq)
     if subword in record.seq:
          next
     else:     
          SeqIO.write(record, outf, "fastq")
          sumlen += len(record.seq)
          ctgood +=1

# Report stats
print("Total number of input sequences/reads:\t", ct)
print("Total length (bp) of input sequences/reads:\t", sumraw)
print("Total number of sequences/reads after filtering:\t", ctgood)
print("Total length (bp) of filtered sequences/reads:\t", sumlen)
print("Result:\t",otd)

