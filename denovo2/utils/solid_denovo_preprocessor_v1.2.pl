#!/usr/bin/perl -w

####################################################################################
#
# Program : solid_denovo_preprocessor.pl
#
# Authors : Vrunda Sheth and Craig Cummings
#
# Purpose : Convert SOLiD colorspace read data in csfasta format into a format that 
#           can be read by the velvet assembler.  This involves three steps, with
#           steps 2 and 3 applying only to SOLiD paired-end (or mate-pair) data
#
#           1) Double-encode the reads such that 0123->ACGT and first base and color 
#              are removed
#           2) Arrange reads such that F5/R3 read follows F3 read for each bead, e.g.
#              >123_456_789_F3
#              F3_sequence...
#              >123_456_789_(F5/R3)
#              R3_sequence...
#           3) Reverse the F3 read so that F3 and R3 are in expected orientation
#
# Usage   : See usage subroutine below or run with -h option
#
####################################################################################

use strict;
use Getopt::Long;
use File::Basename;

# Declare command-line options as globals
use vars qw($opt_help);
use vars qw($opt_f3);
use vars qw($opt_r3);
use vars qw($opt_m);
use vars qw($opt_o);
use vars qw($opt_run);
use vars qw($opt_debug);
use vars qw($opt_version);

use vars qw($PROG);
$PROG = basename($0);

use vars qw($version);
$version = 1.2;

my $read_id;
my $f3_seq;
my $r3_seq;
my $seq;
my %f3_sequence_of_bead;


&GetOptions( "help|h"              => \$opt_help,
             "f3_file|f=s"         => \$opt_f3,
	     "r3_file|r=s"         => \$opt_r3,
	     "mixed_tag_file|m=s"  => \$opt_m,
	     "output|o=s"          => \$opt_o,
	     "run_type|t=s"        => \$opt_run,
             "debug|d"             => \$opt_debug,
             "version|v"           => \$opt_version,
	     ) or &PrintUsage('No options provided');

my $f3_file    = $opt_f3;
my $r3_file    = $opt_r3;
my $output_dir = $opt_o;
my @readFiles  = &ParseCommandLineOptions();
my ($cs_output, $de_output) = &CreateOutputFiles($output_dir);

if ($opt_run eq 'fragment') {

    if (scalar(@readFiles) == 0) {
	die "No read files\n";
    }
    elsif (scalar(@readFiles) > 1) {
	die "Multiple read files not permitted as input\n";
    }
    else {
	my $file = shift(@readFiles);
	open (READ_FILE, $file) or die "Cannot open the input file $file: $!\n";
	while(<READ_FILE>) {
	    chomp;
	    if ( $_ !~ /^\#/ ) {
		print $cs_output "$_\n";
		if ( $_ =~ /^>/ ) {
		    print $de_output "$_\n";
		}
		else {
		    $f3_seq = $_;
		    $f3_seq = &RemoveFirstBaseAndColor($f3_seq);
		    $f3_seq = &DoubleEncode($f3_seq);
		    print $de_output "$f3_seq\n";
		}
	    }
	}
    }
}

else {  # Paired-end or Mate-pair run

    if (scalar(@readFiles) == 2) {  # F3 and F5/R3 input files separate

	open (F3, "$f3_file") or die "Cannot open F3 tag file $!\n";
	
	while (<F3>) {
	    chomp;
	    if ( $_ !~ /^\#/ ) {
		if ( /^>/ ) {
		    if ( /_R3/ || /_F5/) {
			die "F3 tag file contains an F5/R3 tag: '$_' \n";
		    }
		    else {
			$read_id =  $_;
			$read_id =~ s/_F3//;
		    }
		}
		else {
		    $f3_seq = $_;
		    $f3_sequence_of_bead{$read_id} = $f3_seq;
		    # Would be good to put a check in here for single line FASTA
		}
	    }
	}
	close (F3);
	&ScanR3FileAndWriteOutputFiles($r3_file, $cs_output, $de_output, $opt_run);
    }

    elsif (scalar(@readFiles) == 1) { # mixed F3 and F5/R3 read file

	my $temp_R3_file = "$output_dir/R3_reads.csfasta";
	open(R3, ">$temp_R3_file") or die; # create a temporary R3 file

	$seq = '';
	open(MIXED, $opt_m) or die "Can't open mixed read file '$opt_m': $!\n";
	while (<MIXED>) {
	    chomp;
	    if ( $_ !~ /^\#/ ) {
		if ( /^>/ ) {
		    if ( $seq ne '' ) {
			if ( $read_id =~ /_F3/ ) {
			    $read_id =~ s/_F3//;
			    $f3_sequence_of_bead{$read_id} = $seq;
			}
			else {
			    print R3 "$read_id\n$seq\n";
			}
		    }
		    $read_id =  $_;
		}
		else {
		    $seq = $_;
		}
	    }
	}

	if ( $read_id =~ /_F3/ ) { # Process the last record
	    $read_id =~ s/_F3//;
	    $f3_sequence_of_bead{$read_id} = $seq;
	}
	else {
	    print R3 "$read_id\n$seq\n";
	}

	close (R3);
	&ScanR3FileAndWriteOutputFiles($temp_R3_file, $cs_output, $de_output, $opt_run);
	unlink ($temp_R3_file);
    }
    elsif (scalar(@readFiles) < 1) {
	die "No read files\n";
    }
    else {
	die "Multiple read files not permitted as input\n";
    }
}

close $cs_output;
close $de_output;
 

### SUBROUTINES ###


sub ScanR3FileAndWriteOutputFiles {
    my ($r3_file, $cs_output, $de_output, $opt_run) = @_;

    open (R3, "$r3_file") or die "Cannot open F5/R3 tag file $!\n";
    my $r3_seq = '';
    while(<R3>) {
	  chomp;

	if ( $_ !~ /^\#/ ) {        
	    if ( /^>/ ) {
	    	if ( /_F3/ ) {
		      die "F5/R3 tag file contains an F3 tag: '$_' \n";
		  }
		  else {
      if ( (exists $f3_sequence_of_bead{$read_id}) && ($r3_seq ne '') ) {

			print $cs_output $read_id . "_F3\n";
			print $cs_output $f3_sequence_of_bead{$read_id} . "\n";
		  if ($opt_run eq 'paired'){
		  print $cs_output $read_id . "_F5\n";
		  } else {
			print $cs_output $read_id . "_R3\n";
			}
			print $cs_output $r3_seq . "\n";
			
			$f3_seq = $f3_sequence_of_bead{$read_id};
			$f3_seq = &RemoveFirstBaseAndColor($f3_seq);
  		if ($opt_run eq 'mates'){
			$f3_seq = reverse($f3_seq);
			}
			$f3_seq = &DoubleEncode($f3_seq);
			$r3_seq = &RemoveFirstBaseAndColor($r3_seq);
			$r3_seq = &DoubleEncode($r3_seq);
			
			print $de_output $read_id . "_F3\n";
			print $de_output $f3_seq . "\n";
		  if ($opt_run eq 'paired'){
		  print $de_output $read_id . "_F5\n";
		  } else {
			print $de_output $read_id . "_R3\n";
			}
			print $de_output $r3_seq . "\n";
		    }
		    $read_id =  $_;
        $read_id =~ s/_F5-P2//;
        $read_id =~ s/_F5-BC//;        
        $read_id =~ s/_F5//;        
		    $read_id =~ s/_R3//;
    		}
	    }
	    else {
		$r3_seq = $_;
	    }
	}
    }
    
   #processing last read
   if ( (exists $f3_sequence_of_bead{$read_id}) && ($r3_seq ne '') ) {

			print $cs_output $read_id . "_F3\n";
			print $cs_output $f3_sequence_of_bead{$read_id} . "\n";
		  if ($opt_run eq 'paired'){
		  print $cs_output $read_id . "_F5\n";
		  } else {
			print $cs_output $read_id . "_R3\n";
			}
			print $cs_output $r3_seq . "\n";
			
			$f3_seq = $f3_sequence_of_bead{$read_id};
			$f3_seq = &RemoveFirstBaseAndColor($f3_seq);
			if ($opt_run eq 'mates'){
			$f3_seq = reverse($f3_seq);
			}
			$f3_seq = &DoubleEncode($f3_seq);
			$r3_seq = &RemoveFirstBaseAndColor($r3_seq);
			$r3_seq = &DoubleEncode($r3_seq);
			
			print $de_output $read_id . "_F3\n";
			print $de_output $f3_seq . "\n";
		  if ($opt_run eq 'paired'){
		  print $de_output $read_id . "_F5\n";
		  } else {
			print $de_output $read_id . "_R3\n";
			}
			print $de_output $r3_seq . "\n";
		    }
    
    close (R3);
}


sub RemoveFirstBaseAndColor {
    my $seq = shift;
    $seq =~ s/\D.//;
    return $seq;
}		    


sub DoubleEncode {
    my $seq = shift;
    $seq =~ tr/0123/ACGT/; 
    return $seq;
}

		       
sub CreateOutputFiles {
    my $output_dir = shift;
    my $fh1;
    my $fh2;
    my @fileHandles;

    # Create output directory and output files
    #     colorspace_input.csfasta is the csfasta file that will go as input 
    #         to the postprocessor
    #     doubleEncoded_input.de is the double encoded file that will go as 
    #         input to velvet

    mkdir("$output_dir", 0777);

    open $fh1, '>', "$output_dir/colorspace_input.csfasta"
	or die "Cannot create file $output_dir/colorspace_input.csfasta\n";
    push(@fileHandles, $fh1);

    open $fh2, '>', "$output_dir/doubleEncoded_input.de"
	or die "Cannot create file $output_dir/doubleEncoded_input.de\n";
    push(@fileHandles, $fh2);

    return @fileHandles;
}


sub ParseCommandLineOptions {
    my $errors = '';
    my @readFiles;

    if($opt_version){
	print "$PROG version $version\n";
	exit(0);
    }

    if($opt_help) {
	&PrintUsage('');
    }
    
    if ($opt_run) {
	if ($opt_run =~ /\bf\b|\bfragment\b/i) {
	    if ( ! $opt_f3 ) {
		$errors .= "An F3 read file must be provided.\n";
	    }
	    else {
		unless ( -r $opt_f3 ) {
		    $errors .= "F3 read file $opt_f3 is not readable or does not exist\n";
		}
	    }
	    if ( $opt_r3 ) {
		$errors .= "An R3 read file may not be specified for a fragment run.\n"; 
	    }
	    if ( $opt_m ) {
		$errors .= "A mixed F3/(F5/R3) read file may not be specified for a fragment run.\n"; 
	    }
	    $opt_run = 'fragment';
	    push(@readFiles, $opt_f3);
	} 
	elsif ($opt_run =~ /\bm\b|\bmates\b/) {
	    if ( $opt_m ) {
		if ( $opt_f3 || $opt_r3 ) {
		    $errors .= "Cannot specify both mixed tag and F3 or R3 file.\n";
		}
		unless ( -r $opt_m ) {
		    $errors .= "Mixed read file $opt_m is not readable\n";
		}
		push(@readFiles, $opt_m);
	    }
	    elsif ( $opt_f3 && $opt_r3 ) {
		unless ( -r $opt_f3 ) {
		    $errors .= "F3 read file $opt_f3 is not readable or does not exist\n";
		}
		unless ( -r $opt_r3 ) {
		    $errors .= "R3 read file $opt_r3 is not readable or does not exist\n";
		}
		push(@readFiles, $opt_f3, $opt_r3);
	    } 
	    else {
		$errors .= "For a mate-pair run, you must specify either 1) both F3 and R3, or 2) a mixed read file\n";
	    }
	    $opt_run = 'mates';
	}
	elsif ($opt_run =~ /\bp\b|\bpaired\b/) {
	    if ( $opt_m ) {
		if ( $opt_f3 || $opt_r3 ) {
		    $errors .= "Cannot specify both mixed tag and F3 or F5 file.\n";
		}
		unless ( -r $opt_m ) {
		    $errors .= "Mixed read file $opt_m is not readable\n";
		}
		push(@readFiles, $opt_m);
	    }
	    elsif ( $opt_f3 && $opt_r3 ) {
		unless ( -r $opt_f3 ) {
		    $errors .= "F3 read file $opt_f3 is not readable or does not exist\n";
		}
		unless ( -r $opt_r3 ) {
		    $errors .= "F5 read file $opt_r3 is not readable or does not exist\n";
		}
		push(@readFiles, $opt_f3, $opt_r3);
	    } 
	    else {
		$errors .= "For a paired-end run, you must specify either 1) both F3 and F5, or 2) a mixed read file\n";
	    }
	    $opt_run = 'paired';
	}
	else {
	    $errors .= "Run type of $opt_run is not recognized\n";
	}
	
    }
    else {
	$errors .= "Run type must be specified using the --run_type option.\n";
    }

    if(! $opt_o){
	$errors .= "Output directory must be specified using the -o or --output option\n";
    }
    
    if ($errors =~ /\S+/) {
	&PrintUsage($errors);
    }
    else {
	return @readFiles;
    }
}


sub PrintUsage {
    my $error = shift;

    if ($error =~ /\S+/) {
        print "\n***INPUT ERRORS***\n$error\n";
    }

    print <<"END";

Usage: 

Fragment run:
     $PROG --run_type fragment --output out_dir  \\
     --f3_file F3.csfasta

Paired-end (or Mate-pair) run:
     $PROG --run_type paired(mates) --output out_dir     \\
     --f3_file F3.csfasta -r3_file F5/R3.csfasta
     
     or
     
     $PROG --run_type paired(mates) --output out_dir     \\
     --mixed_tag_file F3_plus_F5/R3.csfasta

Required arguments

    -t|--run_type       : Run type:  'fragment' or 'f' for a fragment run,
                          'paired' or 'p' for paired-end run, 'mates' or 'm' for 
                          a mate-paired run (required)
    -o|--output         : output directory (required)

  At least one colorspace read file (csfasta format) is required as input.  

  For a fragment run
    - Use the -f|--f3_file option to specify the csfasta file.  

  For a mate-pair run
    - Reads may be entered as separate F3 and F5/R3 files using the -f|--f3_file 
      and -r|--r3_file options.
    - Alternatively, a single file containing all reads (F3 and F5/R3) may be 
      specified using the -m|--mixed_tag_file option

    -f|--f3_file        : The F3 tag csfasta file for a mate pair or fragment run 
    -r|--r3_file        : The F5/R3 tag csfasta file for a paired-end/mate pair run
    -m|--mixed_tag_file : A csfasta file containing both F3 and F5/R3 reads for
                          a mate-pair run
Other options

    -h|--help           : print usage and exit
    -v|--version        : print version and exit

END

    exit(0);
}

