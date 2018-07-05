#!/usr/bin/perl -w -I . -I ../lib -I ../lib2
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
#########################################
##  This program uses the PAMBuilder package to build a series of PAM-like
##  alignment scoring matrices from a training set of alignments.
#########################################


use strict;
use PAMBuilder;

#########################################
##  MAIN PROGRAM
#########################################

my $A = collectCounts();
printMatrix($A,"Raw Counts");
#my
$A = addPseudocounts($A);
printMatrix($A,"Raw Counts plus Pseudocounts");
my $P = freqsFromCounts($A);
printMatrix($P,"Frequencies");
my $M = ratesFromFreqs($P);
my $PAMM = rescaleRates($M, 0.01);
for (my $i=1; $i<=256; $i+=$i)
{
    print "Mutation rate of PAM$i is ", averageMutationRate($PAMM), "\n";
    printMatrix(lodsFromFreqs(freqsFromRates($PAMM)),"PAM$i (half-bits)");
    $PAMM = multiplyRates($PAMM,$PAMM);
}
