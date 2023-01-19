#!/usr/bin/perl -w
use strict;


# NOte: copy right is owned by Xuewen Wang, ";
# find the sequence in specified range (s)from a sequences file called eg.refseqfasta.txt and save in OUT file given by user;
#print "Copy right reserved 2010 Mar, script designed by Xuewen Wang.\n\n";

# perl sub_substring.pl -o outresult.txt -s inputsequence.fa -r range_positionfile.txt
#read arguments
my %commandin = @ARGV;
 if ((scalar @ARGV)%2 != 0){print "arguments must in pair";}
 my $outfile=$commandin{"-o"};
 my $refseq=$commandin{"-s"};
 my $rangefile=$commandin{"-r"};



# Output file
 	print  usageinfo();
        open (OUT, ">$outfile");  
       	print OUT usageinfo();
	print OUT formatedrealtime();
	print OUT "Source sequence file for extracting subsequences: $refseq\n";
	print OUT "Position data file for extracting subsequences: $rangefile\n";
 	print OUT "Result file of extracted subsequence : $outfile\n";


# read sequence file (refseq.fa), >ID and seq must be in one line only;
	# format of seqfile.fa , tab delimited
	# seqID	seq
	# >g123	attggatagatagatagtag
	# >g2	tttgagagtcgctgctgctgctgctg
	open (REFDNA, $refseq)|| die (" !!!!!!!!!!!filed open failed, pleasse check!!!!!!!!");
	
#read range file containing extraction >ID and positions, >ID name is same and sequences will be extracted
	#read range file: sequence range for extraction
	# format of rangefile.txt , tab delimited, have title exact title line
	# seqID	range_start	range_end	
	# >g123	19058478	19781595	
	# >g2	22173271	22591661	
        open (RANGE, $rangefile)|| die (" !!!!!!!!!!!filed open failed, pleasse check!!!!!!!!");
        my @rangedata=<RANGE>;
        close RANGE;
	

while (my $refseries=<REFDNA>){ 
	#1 while
#get reference seqeunce and ID from file refseq.txt in fasta
       (my $refID, my $Sequences) = split('\t', $refseries);
        chomp $refID;
        #print OUT $refID, "\n"; #ok
#retrieve sequence
         #$Sequences=join("",@refseries);
         $Sequences =~ tr/actg/ACTG/;
         $Sequences =~ s/\s//; #move white space
         $Sequences =~ s/[0-9]//g; #move number
	 chomp $Sequences;
	 #print OUT $Sequences, "\n"; #ok
	 
#extract sequence  
       	foreach my $series(@rangedata){
		chomp $series;		
		unless($series=~ /range_start/){ # position file must contain this word, here is to remove head line
			(my $seqID, my $startposition, my $endposition) = split('\t', $series);			
		      if ($seqID eq $refID){
			#$startposition=~ s/\s//g; #move any space
			#$endposition=~ s/\s//g;
			my $extrlength=$endposition-$startposition+1;
			print OUT "position error\n" if($extrlength<=0);
			#output ID and extracted sequence
			#print OUT "$seqID:$startposition", "_", "$endposition:$extrlength\t";  
			print OUT "$seqID\t", "$extrlength\t";
			  $startposition=$startposition-1; # the start position from 0 in perl instead of 1
 			print OUT uc substr($Sequences,$startposition,$extrlength),"\n";
		      } #end of if
		}#unless
      
	}# foreach /while	

} #END #1 while
close OUT;
close REFDNA;


###########################################
 sub usageinfo
 {# print program name and usage for help if no input available in command line
 my @usage=(); # show how to use the programme
 $usage[0]="Usage: \n";
 $usage[1]=" for    help: perl $0 ; \n";
 $usage[2]=" for running: perl $0 -o outresultfile.txt -s inputsequencefile.fa -r range_positionfile.txt\n  \n";
 $usage[3]="Author: Xuewen Wang\n";
 $usage[4]="year 2011\n\n";
 return @usage; 
 }
 
 
sub formatedrealtime {
my @weekday = ("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat");
my $local_time = gmtime();
my $runtime= "results yielded at $local_time\n";
return $runtime;
}

sub filename{      
	print "typing your output file name (^D to quit filename set):\n";	
	my $filename;				
   
                $filename=<stdin>; 
		chomp $filename;
		if ($filename eq "^D"){exit;}		
		if (-e $filename or $filename=~/^$/) {
			print "file name empty or already exist! Try again.\n";
			filename();					
		}else{		        	
			print "file name $filename accepted.\n";			
			return $filename;		
										
		}
} #end sub

