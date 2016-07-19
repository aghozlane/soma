#!/usr/bin/perl

#
# Program : assembly_stats.pl
#
# Authors : Vrunda Sheth and Craig Cummings
#
# Purpose : Calculate assembly statistics (e.g., # contigs, N50, base frequencies) for an
#           assembly in multi-Fasta format.
#
# Usage   : assembly_stats.pl multi-fasta_file  (output to STDOUT)
#

if (scalar(@ARGV) == 1) {
    $infile = shift(@ARGV);
} 
else {
    die "Usage: $0 <multi-fasta file>\n\n";
}

$count = 0;
$sum = 0;

open(IN,$infile) or die "Can't open input fasta file '$infile': $!\n";
while(<IN>) {
    chomp;
    if ($_ =~ /^>/) {
	$count++;
	if (defined($seq)) {
	    $combinedSeq .= $seq;
	    $sum += length($seq);
	    push(@lengths, length($seq));
	    $seq = '';
	}
    }
    else {
	$seq .= $_;
    }
}
$combinedSeq .= $seq;
$sum += length($seq);
push(@lengths, length($seq));

$average = int($sum / $count);

$add = 0;
$num = 0;

@lengths = sort {$b <=> $a} @lengths;
while ($add < ($sum/2)) {
    $add += $lengths[$num++];
}

$n50 = $lengths[$num-1];

@combinedSeq = split(//, $combinedSeq);
for ($i=0; $i<=$sum; $i++) {
    if (defined($combinedSeq[$i])) {
	$baseCount{$combinedSeq[$i]}++;
    }
}

printf("\nperc A               : %2.1d\n", ($baseCount{'A'}/$sum) * 100);
printf("perc C               : %2.1d\n", ($baseCount{'C'}/$sum) * 100);
printf("perc G               : %2.1d\n", ($baseCount{'G'}/$sum) * 100);
printf("perc T               : %2.1d\n", ($baseCount{'T'}/$sum) * 100);
if (defined($baseCount{'N'}) && ($baseCount{'N'} > 0)) {
    printf("perc N               : %2.1d\n", ($baseCount{'N'}/$sum) * 100);
} else {
    printf("perc N               : %2.1d\n", 0);
}

print "Sum contig length    : $sum\n";
print "Num contigs          : $count\n";
print "Mean contig length   : $average\n";
print "Median contig length : $lengths[int($#lengths/2)]\n";
print "N50 value            : $n50\n";
print "Max                  : $lengths[0]\n\n";
close(IN);


