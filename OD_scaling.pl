#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

#### Script to scale the coverage and appply log2(sample/ input)
#### Order of columns is important here -- input column must be the first one!
#### Usage: perl scaling.pl meanCoverage.bed scaleFactor outputFile


my $in = $ARGV[0];
my $scale = $ARGV[1];
my $out = $ARGV[2];


open IN, $in or die "Can't open $in\n";
open OUT, ">$out" or die "Can't create $out\n";

print OUT "chr\tstart\tend\tinput\tsample1\tsample2\tInputScaled\tScaled1\tScaled2\tlogSample1\tlogSample2\n";

while (<IN>) {
	chomp $_;
	my @tmp = split (/\s+/, $_);
	if ($tmp[3] == 0 && $tmp[4] == 0 && $tmp[5] ==0) {
		print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$tmp[5]\t0\t0\t0\t0\t0\n";
	}
	elsif($tmp[3] == 0) {
		my $sample1 = log2($tmp[4]*$scale);  
		my $sample2 = log2($tmp[5]*$scale); 
		print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$tmp[5]\t0\t$sample1\t$sample2\t$sample1\t$sample2\n";
	}
	else {
		my $inpScaled = $tmp[3]*$scale;
		my $scale1 = $tmp[4]*$scale;
		my $scale2 = $tmp[5]*$scale;
		my $sample1 = log2(($tmp[4]*$scale)/($tmp[3]*$scale));
		my $sample2 = log2(($tmp[5]*$scale)/($tmp[3]*$scale));
		print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$tmp[5]\t$inpScaled\t$scale1\t$scale2\t$sample1\t$sample2\n";
	}
}


###### SUBROUTINE #####
sub log2 { 
    my $n = shift; 
    # using pre-defined log function 
    return log($n) / log(2); 
} 