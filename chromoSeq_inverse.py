#!/usr/bin/python
import sys
import getopt
import os
from Bio import SeqIO

xapp=os.path.basename(__file__)
__author__ = "Xuewen Wang"
def usage():
    version="1.0.0"
    print("Function: This tool programmed in python3 is designed to: \n"
          "\tinverse the sequences (-s) in a genome assembly;\n"
          "\treorder the sequences in the same given order in the idInfo file (-i)\n"
          "\tDrop off sequences if the id is not given"
    )
    print(f"Usage: python {xapp} [options]")
    print(f"e.g., python3 {xapp} -s dnaSeqtest.fasta -i idInvert.txt -o dnaInvertedSeq.fa")
    print("Options:")
    print("\t-h, --help: show this help message and exit")
    print("\t-s, --seqfile: string, required, input file name of the sequences in fasta format")
    print("\t-i, --idfile: string, required input file name containing inversion informaiton. format: id TAB + or - for keep and inverse the seq\n "
          "\t e.g. chr01\t+ "
          "\n\t chr02 -"
          "one id per line")
    print("\t-o, --outfile: the output file name to store the inversed sequences; optional.")
    print("\t-t, --threads: int, the number of parallelized computing threads, default 1. not implemented yet")
    print(f"Version: {version}, April,14th,2024")
    print(f"Support: {__author__},  xwang.kib@gmail.com")



# tested in python_3/3.9.5
# parse arguments
def xpar(argv):
    idfile = r'C:\macrohaptype_testData\testData\idInvert.txt'
    seqfile = r'C:\macrohaptype_testData\testData\dnaSeqtest.fasta'
    outfile = seqfile + "_inverted.fa"

    try:
        opts, args = getopt.getopt(argv, "hi:s:o:", ["idfile=", "seqfile=", "outfile="])
    except getopt.GetoptError  as err:
        print(err)
        usage()
        sys.exit(2)

    if not opts:
        # Print the usage message and exit
        print("Oops: No options provided")
        usage()
        sys.exit(0)

    for opt, arg in opts:
        if opt == '-h':
            print(' parameters: -s seq_file.fasta -i id_inverse_information.txt -o inversed_seq_outFile_prefix')
            usage()
            sys.exit(0)
        elif opt in ("-i", "--idfile"):
            idfile = arg
        elif opt in ("-s", "--seqfile"):
            seqfile = arg
        elif opt in ("-o", "--outfile"):
            outfile = arg
    print("Input id file is " + idfile)
    print("Input sequence file: ", seqfile)
    print('Output seq file is ' + outfile)
    print("")
    return idfile, seqfile, outfile


if __name__ == "__main__":
    idfile, seqfile, outfile = xpar(sys.argv[1:])
idf = idfile
seqf = seqfile
outf = outfile

idlist = []

# read in and hash seq
record_dict = SeqIO.index(seqf, "fasta")

# read in id file
OUTH = open(outf, "w")
ctok = 0
ctno = 0
with open(idf) as INID:
    for ids in INID:
        if ids.strip() == "":
            continue #skip empty lines
        ids, orientation= ids.strip("\n").split('\t')  #>1 - or +
        if (ids[0] == ">"):
            # if input id has >, then preprocess it first
            idlenx = len(ids)
            ids = ids[1:idlenx]
        idlist.append(ids)
        if (ids in record_dict):
            #extrct_seq = record_dict[ids].format("fasta")
            seq=record_dict[ids].seq
            if orientation == '-' :
                reverse_complement_seq = seq.reverse_complement()
                #update seq
                seq=reverse_complement_seq
            fastaSeq=">"+ids+"\n"+str(seq)
            OUTH.write(fastaSeq+"\n")
            ctok += 1
        else:
            print(f"Warning: sequence is not available for IDs: >{ids} \n")
            ctno += 1
total_ids = len(idlist)
print(f"Total # of ids in the id file:\t{total_ids}\n")
print(f"Total # of extracted sequence:\t{ctok}\n")
print(f"Total # of not extracted sequence:\t{ctno}\n")
print(f"Result file with inverted sequence:\t{outf}")
OUTH.close()