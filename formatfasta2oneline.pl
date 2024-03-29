#!/usr/bin/perl -w
use strict;
&usageinfo();
# usage example: perl formatfasta2oneline.pl -i input_seq.fasta  -o formated.input_seq.fasta
#programmed by Xuewen Wang, 2019 Nov; 
# Current version:  2019 Nov;
my %commandin = @ARGV;
if ((scalar @ARGV)%2 != 0){print "arguments must in pair";}
my $refin=$commandin{"-i"}; 
my $outputfile=$commandin{"-o"}; 
 unless(-e $refin){print "Sequence file does not exist in the working directoty\n";}
 open(OUT, ">$outputfile")|| die (" formated file writing failed, pleasse check\n");
 my $satfile=$outputfile.".sat1";
 open(OUTsat, ">$satfile")|| die (" sat1 file writing failed, pleasse check\n");
 my $filename=$refin;
open(FILE, $filename) || die("Couldn't read file $filename\n");  
local $/ = "\n>";  # read by FASTA tag >
my $seqformat="false";
my $seqno=0;
my $seqlen=0;
my $eachseqlen=0;
my $startpos =0;
my $overlap=0;
my $ii=0;
my $seq="";
my $testlength=0;

while (my $seqinfo = <FILE>) {
chomp $seqinfo;
     
    $seqinfo =~ s/^>*(.+)\n//;  
     my $seqID=">".$1; #   head line
	 my @seqIDcontn=split('\s',$seqID); #truncate seqID by first space
	 $seqID=$seqIDcontn[0];
	 $seqno ++ if($1);
	 $seqformat="TRUE" if($1);
	 $seqinfo=~ s/[0-9\n\s]//g;  
	 $eachseqlen=length $seqinfo;
	 $seqlen=$seqlen+$eachseqlen;
	 print OUT  $seqID,"\t";		    
	 print OUT  $seqinfo, "\n";
		
} #end while

&runtime(\*OUTsat); #for run log information
if($seqformat eq "TRUE"){
	print OUTsat "statistic summary of input sequence\(s\)",  "\n";
	print OUTsat "total input sequences#\t: ", $seqno, "\n"; 
	print OUTsat "total length (bp) after formatted:\t", $seqlen, "\n";
	print OUTsat "total chunking sites after formatted:\t", $ii, "\n";
}else{ 
	print OUTsat "input sequence format is not qualified. check the presence of sign '>' \n";
	exit;
}
close FILE;
close OUT;
close OUTsat;
sub usageinfo
 {
 my @usage=(); 
 $usage[0]="Function: format DNA sequences in a file to 1 line  with ID <tab> sequence for each sequence\n";
 $usage[1]=" for    help: perl $0 ; \n";
 $usage[2]=" for running: perl $0  -i [input DNA sequence file name]  -o [formated sequence file name]\n";
 $usage[3]="Author: Xuewen Wang\n";
 $usage[4]="year 2009-2019\n";
 unless(@ARGV){print @usage; exit;} 
 }
sub runtime() {
my $OUTfile=shift @_;
my $local_time = gmtime();
print {$OUTfile} "$0 was run and results were yielded at $local_time\n";
}
exit;

