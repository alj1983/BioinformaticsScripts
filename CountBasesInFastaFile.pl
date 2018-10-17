#!/usr/bin/perl -w

use strict;
use warnings;
# Script to count the number of bases or amino acids in a fasta file; 
if ($ARGV[0] eq '') { die "A fasta input file is needed"};

open(IN, $ARGV[0]) or die "Can not open file $ARGV[0]";

my $seq;
while(<IN>) {
    unless (/>/){
	chomp;
	$seq .= $_;
    }
}


close IN;

my $seqlength=length($seq);

print "Number of bases or amino acids in the fasta file: $seqlength\n";
