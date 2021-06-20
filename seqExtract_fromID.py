#!/usr/bin/python
import sys
import getopt
from Bio import SeqIO
''' this script programmed in python3 is designed to extract the seq with seqID from a fasta file
    The hash /dict method is used
    
    USAGE: 
    python3 seqExtract_fromID.py -s seq_file_fasta -i id_file_one_seqID_per_line -o out_file_extract_seq
    
    e.g.: python3 seqExtract_fromID.py -s testseq.fasta -i id.txt -o  testseq.fasta_extracted.fa
    the strategy used as previous perl script: seqExtract_fromID_rel.pl by Xuewen Wang
    Version:  Nov-9-2020
    Author: Xuewen Wang, Email: xwang.kib@gmail.com
'''
__author__="Xuewen Wang"

# module load python_3/3.6.6
#parse arguments
def xpar(argv):
    idfile = ''
    seqfile=''
    outfile = seqfile+"_extracted.fa"

    try:
        opts, args = getopt.getopt(argv, "hi:o:s:", ["idfile=", "seqfile=", "outfile="])
    except getopt.GetoptError:
        print("parameters: python3 seqExtract_fromID.py -i input_idfile -s fasta_seqfile -o Outfile_extracted_seq ")
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-h':
            print(' parameters: -i input_idfile -s input_fasta_seqfile -o Outfile_extracted_seq')
            sys.exit()
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
idf=idfile
seqf=seqfile
outf=outfile

idlist=[]

# read in and hash seq
record_dict=SeqIO.index(seqf, "fasta")

# read in id file
OUTH=open(outf,"w")
ctok=0
ctno=0
with open(idf) as INID:
    for ids in INID:
        ids=ids.strip("\n")
        if (ids[0] == ">"):
            idlenx = len(ids)
            ids = ids[1:idlenx]  
        idlist.append(ids) 
        if(ids in record_dict): 
            extrct_seq=record_dict[ids].format("fasta")
            OUTH.write(extrct_seq)
            ctok+=1
        else:
            print(f"Warning: sequence is not available for IDs: >{ids} \n")
            ctno += 1
total_ids=len(idlist)
print(f"Total # of id to bait:\t{total_ids}\n")
print(f"Total # of extracted sequence:\t{ctok}\n")
print(f"Total # of not extracted sequence:\t{ctno}\n")
print(f"Result file with xxtracted sequence:\t{outf}")
OUTH.close()



