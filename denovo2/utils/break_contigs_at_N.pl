#!/usr/bin/perl -w

#
# Program : break_contigs_at_N.pl
#
# Authors : Craig Cummings and Vrunda Sheth
#
# Purpose : Split sequences at runs of Ns
#
# Usage   : break_contigs_at_N.pl -i <input_file> -o <output_file> -n <max_allowed_consec_N>
#
# Output  : Multi-fasta file.  Contig IDs have a '_n' suffix indicating their order in the
#           original contig from which they were split.
#

use Getopt::Std;

getopt('ion', \%opts);
if (defined($opts{'i'}) && defined($opts{'o'}) && defined($opts{'n'}) ) {
    $infile               = $opts{'i'};
    $outfile              = $opts{'o'};
    $max_allowed_consec_N = $opts{'n'};    
} 
else {
    die "Usage: $0 -i <input_file> -o <output_file> -n <max_allowed_consec_N>\n\n";
}

$seq = '';

open(IN,$infile) or die "Can't open input file $infile: $!\n";
while(<IN>) {
    chomp;
    if ($_=~/^>/) {
	if ($seq ne '') {
	    $sequenceOf{$seqID} = $seq;
	}
	$seqID = $_;
	$seq   = '';
    }
    else {
	$seq .= $_;
    }
}
$sequenceOf{$seqID} = $seq;

open(OUT,">$outfile") or die;

foreach $id (keys %sequenceOf) {
    $count = 0;
    @subseqs = split(/N{$max_allowed_consec_N,}/, $sequenceOf{$id});
    foreach $subseq (@subseqs) {
	if ($subseq) {
	    print OUT "$id"."_"."$count\n";
	    for ($i=0; $i<=length($subseq); $i+=60) {
		print OUT substr($subseq, $i, 60), "\n";
	    }
	}
	$count++;
    }
}

		
		
			
				
