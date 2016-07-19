#!/usr/bin/perl -w

#
# Program : parse_mummer_coord_file.pl
#
# Author  : Craig Cummings
#
# Purpose : Parse and summarize .coord files produced by the MUMmer show-coords program
#
# Usage   : parse_mummer_coord_file.pl <coord_file>
#           Coord file must have been produced using the -c (percent coverage columns)
#           and -l (sequence length columns) command line switches.
#

if (scalar(@ARGV) == 1) {
    $coordFile = shift(@ARGV);
} 
else {
    die "\nUsage: $0 <COORD_FILE>\n\n";
}

# Get reference sequence file name
open(COORDS, $coordFile) or die "Can't open coords file $coordFile: $!\n";
$firstLine = <COORDS>;		# sequence files are on first line of .coords file
@seqFiles = split(/\s+/, $firstLine);
$referenceFile = $seqFiles[0];
close(COORDS);


# Open reference sequence file and determine sequence length
$refSeq = '';
open(REF, $referenceFile) or die "Can't open reference genome file at $referenceFile: $!\n";
while (<REF>) {
    chomp;
    unless (($_ =~ /^>/) or ($_ !~ /\S+/)) {
	$refSeq .= $_;
    }
}
$genomeLength = length($refSeq);   


# Initialize coverage array
for ($i=0; $i<=$genomeLength; $i++) {
    $hitCoverage[$i] = 0;
}


# Parse hits in coords file
open(COORDS, $coordFile) or die;
while (<COORDS>) {
    chomp;
    if (/\[S1\]/) {		#  Header row
	unless ( /\[LEN 1\]/ && /\[LEN 2\]/ && /\[% IDY\]/ ) {
	    die "Coord file not in expected format.  Did you run show-coords with the ".
		"-c and -l options?\n\n";
	}
    }
    elsif (/^\s*\d+/) {
	s/\|//g;
	s/^\s+//;
	@aLine = split(/\s+/, $_);
	($refLeft, $refRight, $percID) = ($aLine[0], $aLine[1], $aLine[6]);
	for ($i=$refLeft; $i<=$refRight; $i++) {
	    $hitCoverage[$i]++;
	}
	$weightedPercId = ($refRight - $refLeft) * ($percID/100);
	$sumWeightedPercId += $weightedPercId;
	$numHits++;
	$numBases = $refRight - $refLeft - 1;
	$sumBases += $numBases;
    } 
    else {			# blank rows, other header rows, '======...'
	next;
    }
}
close(COORDS);


# Determine number of bases with each coverage value
for ($i=1; $i<=$genomeLength; $i++) {
    $numBasesWithCoverage{$hitCoverage[$i]}++;
}


# Calculate total number of bases covered (thanks Nisha!)
$totCovered=0;
for($i=0; $i<=$genomeLength; $i++) {
    if($hitCoverage[$i] != 0) {
        $totCovered++;
    }
}


# Print output
print "$numHits hits\n\n";
print "Cov\tNum bases\n";
print "---\t---------\n";
foreach $cov (sort(keys(%numBasesWithCoverage))) {
    print "$cov\t", $numBasesWithCoverage{$cov}, "\n";
}

print sprintf("\nWeighted average identity = %.2f", 100 * $sumWeightedPercId / $sumBases), 
    "% \n";
print sprintf("Percent of genome covered = %.2f", 100 * $totCovered / $genomeLength), 
    "% \n";

