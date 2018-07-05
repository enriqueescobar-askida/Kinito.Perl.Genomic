#!/usr/bin/perl -I . -I ../perllib
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
#########################################
##  This program generates a random test case for miniphrap.
#########################################

use strict;   ## request conservative error checking
use Util;

## Generate target string of 10000 bases.
my $len = 10000;
my @target = map {('a','c','g','t')[int rand(4)]} (0..$len-1);
my %reads;

## Now generate 80 reads with errors.
foreach (1..80)
{
    my $middle = int rand($len);
    my $halflen = 200 + int rand(100);
    my $left = max($middle - $halflen, 0) || "0" ;
    my $right = min($middle + $halflen, $len-1);
    my $qualString;
    my $currQual = 1;
    
    foreach (0 .. ($right-$left))
    {
        $currQual = (0,0,1,2,3,4,5,6,7,7)[$currQual] + rand(3);
        $qualString .= int($currQual);
    }
    $reads{$left} = 
    join('', "$left..$right\n",
        mutate(@target[$left..$right]), "\n$qualString\n");
}

## For human readability, number the reads as they should appear left to
## right in the output.
my $i=1;

foreach (sort {$a <=> $b} keys %reads)
{
    print "Q>Read$i $reads{$_}";
    $i++;
}

my $currQual = 1;
print "Q>Template\n"; print join('',@target); print "\n";

foreach (1..$len)
{
    $currQual = (1,1,1,2,3,4,5,6,7,7)[$currQual] + rand(3);
    print int($currQual);
}
print "\n";

#########################################
sub mutate
{
##  mutates a DNA sequence, given as list.
##  RETURNS: mutated sequence, as string.
#########################################
    my @mutant;
    while (@_)
    {
        my $r = rand(60);
        
        if ($r == 0)
        { ## deletion 
            pop @_;
        }
        elsif ($r == 1)
        { ## insertion
            push @mutant, ('A','C','G','T')[int rand(4)];
        }
        elsif ($r == 2)
        { ## substitution 
            push @mutant, ('A','C','G','T')[int rand(4)];
            pop @_;
        }
        else
        {
            push @mutant, pop @_;
        }
    }
    return join('', reverse @mutant);
}

