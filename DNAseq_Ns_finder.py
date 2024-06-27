#!/usr/bin/python
import sys
import gzip
from Bio import SeqIO

def find_continuous_ns(dna_sequence):
  """Finds all continuous Ns in a DNA sequence.

  Args:
    dna_sequence: A string representing the DNA sequence.

  Returns:
    A list of tuples, where each tuple represents a continuous N segment.
    Each tuple contains two values: the position of the first N in the
    segment (in 1-based indexing), and the length of the segment.
  """

  # Initialize the list of N segments.
  n_segments = []

  # Keep track of the start and end positions of the current N segment.
  start = None
  end = None

  # Iterate over the DNA sequence.
  for i, nucleotide in enumerate(dna_sequence):
    # start a new N segment.
    if nucleotide == "N" and start is None:
      start = i + 1

    # end the current N segment.
    elif nucleotide != "N" and start is not None:
      end = i
      n_segments.append((start, end - start + 1))
      start = None

  if start is not None:
    end = len(dna_sequence)
    n_segments.append((start, end - start + 1))

  # Return the list of N segments.
  return n_segments


  __author__ = "Xuewen Wang"

  usage = """ python3 DNAseq_Ns_finder.py DNAfastaFile >Ns_result.tsv
            e.g.: python3 DNAseq_Ns_finder.py DNA_Ns_test.fasta >Ns_result.tsv
            """
  print(usage)


if len(sys.argv) == 1:
  print("python3 DNAseq_Ns_finder.py DNAfastaFile >Ns_result.tsv")
  exit(1)

print(f"#ID\tStart_Pos\tLength_of_Ns(bp)")
#with gzip.open(sys.argv[1], "rt") as handle:
#fileinput="C:\snpADtestData\DNA_Ns_test.fasta"
with open(sys.argv[1], "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        id=record.id
        dna_sequence=record.seq
        n_segments = find_continuous_ns(dna_sequence)
        for n_segment in n_segments:
            print(f"{id}\t{n_segment[0]}\t{n_segment[1]}")
