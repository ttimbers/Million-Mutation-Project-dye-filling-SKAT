#!/usr/bin/perl -w
#################################################################
#
# script to filter variants for SKAT from master file MMP.vcf
#
#################################################################
use strict;

use Getopt::Long;

my $help = 0;
my $input;
my $output;
my $affect_protein = 0;
my $strain_filename;

GetOptions (
            "help" => \$help,
            "input=s" => \$input,
            "output=s" => \$output,
            "protein" => \$affect_protein,
            "strain=s" => \$strain_filename
            );

#
# provides help
#
if ( (!defined $input) || (!defined $output) ) {
  $help = 1;
}
if ($help) {
   write;
   exit;
}
open (INPUT, "< $input") || die "Can't open input vcf file $input\n";
open (OUTPUT, "> $output") || die "Can't open output vcf file $output";

my %strains_used = ("#CHROM", 1, "POS", 1, "ID", 1, "REF", 1, "ALT", 1, "QUAL", 1, "FILTER", 1, "INFO", 1, "FORMAT", 1);
my $filter_by_strain = 0;
if (defined $strain_filename) {
  open (STRAIN, "< $strain_filename") || die "Can't open input strain file $strain_filename\n";
  while (<STRAIN>) {
    chomp;
    $strains_used{$_} = 1;
    $filter_by_strain++;
  }
  close (STRAIN);
  unless ($filter_by_strain) {
    die "Your strain file is empty, create a proper strain file for filtering or do not use the -strain option\n";
  }
}

my @columns_to_keep = ();
while(<INPUT>) {
  chomp;
  my $input_line = $_;
  my @values_to_keep = ();
  if ($input_line =~ /^#/) {
    if ($input_line =~ /^#CHROM/) {
      unless ($input_line =~ /^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT/) {
        die "Your vcf file has the wrong header format, try using the master file MMP.vcf\n";
      }
      my @entries = split(/\t/, $input_line);
      my $current_index = 0;
      for my $label (@entries) {
        if ( (defined $strains_used{$label}) || (! $filter_by_strain) ) {
          push (@columns_to_keep, $current_index);
          push (@values_to_keep, $label);
        }
        $current_index++;
      }
      print OUTPUT join ("\t", @values_to_keep);
      print OUTPUT "\n";
    } else {
      print OUTPUT "$input_line\n";
    } 
  } else {
    if ($affect_protein) {
      unless ( ($input_line =~ /GF=intron_splicing/) || ($input_line =~ /GF=coding_exon/) || ($input_line =~ /SVTYPE=DEL/) ) {next;}
      if ( ($input_line =~ /GRANT=0/) || ($input_line =~ /CODING=NA/) ) {next};
    }
    my @entries = split(/\t/, $input_line);
    my $output_variant = 0;
    for my $current_index (@columns_to_keep) {
      push (@values_to_keep, $entries[$current_index]);
      if ($entries[$current_index] eq "1/1") {
        $output_variant = 1;
      }
    }
    if ($output_variant) {
      print OUTPUT join ("\t", @values_to_keep);
      print OUTPUT "\n";
    }
  }
}

close (INPUT);
close (OUTPUT);
#
# help message
#
format STDOUT =
Usage: ./filter_MMP_variants.pl [options] -input <input_vcf_filename> -output <output_vcf_filename>

Options:
-help, -h:
	provides this help
		
-input <input_vcf_filename>, -i <input_vcf_filename>: 
	required input vcf file
		
-output <output_vcf_filename>, -o <output_vcf_filename>
	required output file
		
-protein, -p:
	only keeps the variants potentially affecting protein AA sequences
		
-strain <strain_filename>, -s <strain_filename>:
	only keeps the strains listed in <strain_filename>

The input file has to be in the vcf format like the MMP.vcf file provided, the output will
be in the same format. The optional strain file should be a list with one strain per line.
One can a use a Unix pipe to work directly with the compressed master variant file without
having to uncompress it first, for exmaple:
gzcat MMP.vcf.gz | ./filter_MMP_variants.pl -input - -output filteredMMP.vcf -strain strain_file.txt -protein
.
