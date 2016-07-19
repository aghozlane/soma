#!/usr/bin/env perl

my $exeDeNovo = $ENV{"denovo2"};
my $exeMUMmer = "${exeDeNovo}/MUMmer3.22";
my $exeUTILS  = "${exeDeNovo}/utils";


print "
Assembly Analysis Pipeline for SOLiD v.2.0
Copyright (2010) by Life Technologies 
***************************************************\n";

my $usage = "

usage: \$denovo2/analyze.pl <contig_file> [-options]

* Input: *

 contig_file - fasta file with assembled contigs/scaffolds in base-space, 2-base encoded
               (color space) or double encoded space. Double encoding (de) is an equivalent to
               color space where 0=A,1=C,2=G,3=T.

* Output: *

 n50.stats.txt      - file reporting contigs length statistics.
 cumulative.len.txt - a list of contig sizes sorted in decreasing order, and their accumulation.
 coverage.stats.txt - percentage of genome coverage. for base-space contigs max coverage is 100%.
                      for (de) contigs maximum coverage is 50%, due to reference representation.
                      this file is generated only if reference sequence is provided.
 mapview.pdf        - a plot showing alignment between contigs and reference sequence. the file
                      is generated only if reference sequence is provided.

* OPTIONS: *

 -outdir dir        outputs results and intermediate files into \"dir\" directory
                    (default \"analysis\").
 -break_scaf        for a file with scaffolds this option will report contigs stats.
 -min_length xx     minimum length of contigs to be reported in statistics.
 -ref_file  rf      fasta file with reference sequence. this option must be used with -ref_type
                    and -cont_type options to provide type of representation of the sequence.
 -ref_type  rt      indicates representation of the reference. \"rt\" can be \"nt\",\"de\",\"color\",
                    or \"de2\". \"nt\" - base-space, \"de\" - double encoded, \"color\" - color
                    space, \"de2\" - double encoded forward plus reverse.
 -cont_type ct      indicates representation of contigs. \"ct\" can be \"nt\",\"de\", or \"color\".
                    \"nt\" - base-space, \"de\" - double encoded, \"color\" - color space.

Examples of running:

   \$denovo2/analyze.pl contigs.fa -outdir anls1

   \$denovo2/analyze.pl contigs.de -cont_type de -ref_file reference.csfasta -ref_type color
     -outdir anls2

   \$denovo2/analyze.pl scaffolds.fasta -cont_type nt -ref_file reference.fasta -ref_type nt
     -outdir anls3 -break_scaf                    

\n";

if ($#ARGV < 1) { print $usage; exit(1); }

# List of options
my $cont_file = shift @ARGV;
my $ref_file = "";
my $ref_type = "";
my $cont_type = "";
my $break_scaf = 0;
my $outdir="analysis";
my $min_length = 70;

#Initializing input parameters
#*************************************************************************

while ($#ARGV >= 0) {
		$opt   = shift @ARGV;
		if ($opt =~ /^\-/) {
				if ($opt eq "-outdir") {
						$outdir = shift @ARGV;
        }
    	  elsif ($opt eq "-min_length") {
						$min_length = shift @ARGV;
        }					
        elsif ($opt eq "-ref_file") {
            $ref_file = shift @ARGV;
        }					
        elsif ($opt eq "-ref_type") {
            $ref_type = shift @ARGV;
        }					
        elsif ($opt eq "-cont_type") {
            $cont_type = shift @ARGV;
        }					
        elsif ($opt eq "-break_scaf") {
            $break_scaf = 1;
        }
        else
        { print $usage; exit(1);}					        
		} else { print $usage; exit(1);}					        
}
# Running Analysis Pipeline
#********************************************************************
mkdir($outdir,0777);

if($ref_file ne "")
{
  if($cont_type eq "" or $ref_type eq "")
  {
   print "Please provide the type (format) of both contigs and reference.\n";
   exit(1);
  }
  
   
  if($cont_type eq "nt" and $ref_type ne "nt")
  {
     print "Converting contigs to double encoded.\n";
   
     my $fasta2deCmd = "java -cp ${exeUTILS}/miniAssembler.jar com.lifetech.miniAssembler.util.FormatsTranslator fasta2de $cont_file > $outdir/contigs.de"; 
     print "$fasta2deCmd\n";
 
     my $res = system($fasta2deCmd);   
     if ($res != 0) {
      print "Conversion failed: $res\n";
      exit(1);
      }
      
     $cont_file = "$outdir/contigs.de"; 
     $cont_type = "de";
   }
   elsif($cont_type eq "color")
   {
     print "Converting contigs to double encoded.\n";
   
     my $color2deCmd = "java -cp ${exeUTILS}/miniAssembler.jar com.lifetech.miniAssembler.util.FormatsTranslator color2de $cont_file > $outdir/contigs.de"; 
     print "$fasta2deCmd\n";
 
     my $res = system($fasta2deCmd);   
     if ($res != 0) {
      print "Conversion failed: $res\n";
      exit(1);
      }
      
     $cont_file = "$outdir/contigs.de"; 
     $cont_type = "de";
   }
   
 if($ref_type eq "nt" and $cont_type ne "nt")
    {
     print "Converting reference to double encoded reference.\n";
   
     my $fasta2deCmd = "java -cp ${exeUTILS}/miniAssembler.jar com.lifetech.miniAssembler.util.FormatsTranslator fasta2de $ref_file > $outdir/reference.de"; 
     print "$fasta2deCmd\n";
 
     my $res = system($fasta2deCmd);   
     if ($res != 0) {
      print "Conversion failed: $res\n";
      exit(1);
      }
      
     $ref_file = "$outdir/reference.de"; 
     $ref_type = "de";
    }
    elsif($ref_type eq "color")
    {
     print "Converting reference to double encoded reference.\n";
   
     my $color2deCmd = "java -cp ${exeUTILS}/miniAssembler.jar com.lifetech.miniAssembler.util.FormatsTranslator color2de $ref_file > $outdir/reference.de"; 
     print "$color2deCmd\n";
 
     my $res = system($color2deCmd);   
     if ($res != 0) {
      print "Conversion failed: $res\n";
      exit(1);
      }
      
     $ref_file = "$outdir/reference.de"; 
     $ref_type = "de";
    }
  
   if($ref_type eq "de")
    {
     print "Reversing and merging reference.\n";
     my $refConvertCmd = "${exeUTILS}/reverse_and_concatenate_de_genome.pl $ref_file $outdir/reference.de2";
     print "$refConvertCmd\n";
   
     my $res = system($refConvertCmd);   
     if ($res != 0) {
     print "Reversing and merging reference into a new reference failed: $res\n";
     exit(1);
     }
     $ref_file = "$outdir/reference.de2";
    }
 }

 if($break_scaf == 1)
 {
 print "Break Scafflods. \n";
 mkdir("$outdir",0777);
 my $breakScaffoldsCmd = "${exeUTILS}/break_contigs_at_N.pl -i $cont_file -o $outdir/contigs.N.split.fa -n 1";
 print "$scaffoldStatCmd\n";
 
 my $res = system($breakScaffoldsCmd);   
 if ($res != 0) {
  print "Scaffold breaking failed: $res\n";
 }
 $cont_file = "$outdir/contigs.N.split.fa";
 }
 
 print "Run Analysis. \n";
 &generate_stats($cont_file, $ref_file, "$outdir");

#********************************************************************
sub generate_stats {
    my ($contigs_file, $reference_file, $outdir) = @_;

 # Generate N50 Stats
 my $scaffoldStatCmd = "mkdir -p $outdir;${exeUTILS}/assembly_stats.pl ${contigs_file} > $outdir/n50.stats.txt";
 print "$scaffoldStatCmd\n";
 
 my $res = system($scaffoldStatCmd);   
 if ($res != 0) {
  print "Generation of N50 stats failed: $res\n";
 }

 # Generate Cumulative List
 my $scaffoldCumCmd = "${exeUTILS}/cumlength ${contigs_file} $min_length > $outdir/cumulative.len.txt";
 print "$scaffoldCumCmd\n";
 
 $res = system($scaffoldCumCmd);   
 if ($res != 0) {
  print "Generation of cumulative contigs length failed: $res\n";
 }

 if($reference_file ne "")
 {

 # Alligning to reference genome
 if ($res == 0) {
 my $scaffoldNucmerCmd = "${exeMUMmer}/nucmer $reference_file ${contigs_file} -prefix $outdir/out";
 $res = system($scaffoldNucmerCmd);
  print "$scaffoldNucmerCmd\n";
 }   
 if ($res == 0) {
 my $scaffoldNucmerCmd = "${exeMUMmer}/delta-filter -1 $outdir/out.delta > $outdir/out_1-filter.delta";
 $res = system($scaffoldNucmerCmd);
  print "$scaffoldNucmerCmd\n";  
 }
 if ($res == 0) {
 my $scaffoldNucmerCmd = "${exeMUMmer}/show-coords -rcl $outdir/out_1-filter.delta > $outdir/out_1-filter.coords";
 $res = system($scaffoldNucmerCmd); 
  print "$scaffoldNucmerCmd\n";
 }
 if ($res == 0) {
 my $scaffoldNucmerCmd = "${exeUTILS}/parse_mummer_coord_file.pl $outdir/out_1-filter.coords > $outdir/coverage.stats.txt";
 $res = system($scaffoldNucmerCmd);
  print "$scaffoldNucmerCmd\n"; 
 } 
 if ($res == 0) {
 my $scaffoldNucmerCmd = "${exeMUMmer}/mummerplot --layout $outdir/out_1-filter.delta -p $outdir/cov --terminal postscript";
 $res = system($scaffoldNucmerCmd); 
  print "$scaffoldNucmerCmd\n";
 }
  
 if ($res != 0) {
 print "Nucmer alignment failed: $res\n";
 }

 # generate plot 
 my $ps2pdfCmd = "ps2pdf $outdir/cov.ps $outdir/mapview.pdf;";
 print "$ps2pdfCmd \n";
 
 $res = system($ps2pdfCmd);   
 if ($res != 0) {
  print "Generation of pdf file failed: $res\n";
 }
 }
}

#*******************************************************************





