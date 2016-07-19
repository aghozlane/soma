#!/usr/bin/perl -w

use Bio::SeqIO;
use Bio::Seq;

if (scalar(@ARGV) == 2) {
    $de_genome_file = shift @ARGV;
    $out_genome_file = shift @ARGV;
}
else {
    die "\nUsage: $0 double_encoded_genome_fasta_file output_file \n\n";
}

$in_stream  = Bio::SeqIO->new( -file => $de_genome_file, 
			       -format => 'fasta');
$out_stream = Bio::SeqIO->new( -file => ">$out_genome_file", 
			       -format => 'fasta');

$de_genome_forward = $in_stream->next_seq();
$de_genome_forward_seq = $de_genome_forward->seq();
$de_genome_reverse_seq = reverse($de_genome_forward_seq);
$de_genome_forward->seq( $de_genome_forward_seq . $de_genome_reverse_seq );
$out_stream->write_seq($de_genome_forward);
print STDERR "Wrote output to '$out_genome_file'\n";
