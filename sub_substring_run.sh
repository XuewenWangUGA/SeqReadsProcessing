#!/usr/bin/bash
# to run the sub_substring

#step1: format refseq
refseq=".fasta"
perl formatfasta2oneline.pl -i $refseqrefseq  -o $refseq.formated.fasta
#step2: extract subseq based on range file
perl sub_substring.pl -o MHtartgetRegionSeq.fasta -s $refseq.formated.fasta -r range_MHtartgetRegionSeq.txt

#-r range file containing extraction >ID and positions, >ID name is same and sequences will be extracted
	#read range file: sequence range for extraction
	# format of rangefile.txt , tab delimited, have title exact title line
	# seqID	range_start	range_end	
	# >g123	19058478	19781595	
	# >g2	22173271	22591661