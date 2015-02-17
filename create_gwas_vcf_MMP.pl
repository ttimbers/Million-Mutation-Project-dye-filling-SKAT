#!/usr/bin/perl -w
#
use strict;

my $strain_filename = "list_VCstrains_vcf.txt";
my $allele_filename = "gknames.txt";

my %strains = ();
open(STRAIN, "$strain_filename") or die "can't open list of strains file $strain_filename\n";
while (<STRAIN>) {
  chomp;
  $strains{$_} = 1;
}
close (STRAIN);

my %alleles = ();
my %gknames = ();
open(ALLELE, "$allele_filename") or die "can't open allele/gkname file $allele_filename\n";
while (<ALLELE>) {
  chomp;
  my @entries = split(/\t/);
  $entries[1] =~ /^(VC\d+)\.(\S+)/;
  my $strain = $1;
  my $mutation = $2;
  if (defined $strains{$strain}) {
    unless (defined $alleles{$mutation}) {
      $alleles{$mutation} = [];
      $gknames{$mutation} = $entries[0];
    }
    push (@{$alleles{$mutation}}, $strain);
  }
}
close (ALLELE);

my @genotype = ();
my @sorted_strains = sort keys %strains;
my $index = 0;
my %index_of_strains = ();
for my $strain_name (@sorted_strains) {
  $index_of_strains{$strain_name} = $index;
  $index++;
  push (@genotype, '0/0');
}
my $first_file = 1;
my %valid_bases = ("A", 1, "C", 1, "G", 1, "T", 1);
for my $strain (@sorted_strains) {
  my $in_vcf_file = $strain . ".combined.vcf";
  open(IN, "$in_vcf_file") or die "can't open vcf file $in_vcf_file for input\n";
  while (<IN>) {
    chomp;
    my $line = $_;
    if ($line =~ /^#/) {
      if ($first_file) {
        if ($line =~ /^#CHROM/) {
          print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
          print "$line\tFORMAT\t";
          print join ("\t", @sorted_strains);
          print "\n";
          $first_file = 0;
        } else {
          print "$line\n";
        }
      }
      next;
    }
    my @entries = split(/\t/, $line);
    if ($entries[7] =~ /IMPRECISE|SVTYPE=BND/) {next;}
    $entries[3] =~ tr/[a-z]/[A-Z]/;
    $entries[4] =~ tr/[a-z]/[A-Z]/;
    $entries[2] =~ s/^VC\d+\.//;
    my $label = $entries[2];
    unless (defined $alleles{$label}) {
      my $chromosome = $entries[0];
      my $coordinate = $entries[1];
      my $reference = $entries[3];
      my $mutation = $entries[4];
      if ($entries[7] =~ /INDEL=(\-|\+)(\d+)/) {
        my $indel_type = $1;
        my $indel_size = $2;
        if ($indel_size >= 10) {
          if ($indel_type eq "+") {
            $mutation = $reference . "+" . $indel_size . "N";
          } else {
            $reference = $mutation . "+" . $indel_size . "N";
          }
        }
        $label = $chromosome . "_" . $coordinate . "_" . $reference . "_" . $mutation;
      }
    }
    unless (defined $alleles{$label}) {next;}
    my @change_genotype = @{$alleles{$label}};
    for my $mutated_strain (@change_genotype) {
      $genotype[$index_of_strains{$mutated_strain}] = '1/1';
    }
    $entries[2] = $gknames{$label};
    print join ("\t", @entries);
    print "\tGT\t";
    print join ("\t", @genotype);
    print "\n";
    delete $alleles{$label};
    for my $mutated_strain2 (@change_genotype) {
      $genotype[$index_of_strains{$mutated_strain2}] = '0/0';
    }
  }
}
