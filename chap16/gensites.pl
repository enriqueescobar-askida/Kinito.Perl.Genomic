#!/usr/bin/perl -I . -I ../perllib
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################

#########################################
##  A short program to generate test inputs for exact.pl.
#########################################

use strict;
use Util;

main();

#########################################
sub main
{
##  calls generateRandomSites and prints the results.
#########################################
    my $ranSites = generateRandomSites(20,100000);
#   my $ranSites = [0,2,4,7,10,12,16,18,20];
    print "sites: ",  join(',', @$ranSites), "\n";
    printPairwiseDistances($ranSites);
}

#########################################
sub generateRandomSites
{
##  generates a set of random sites.
##  RETURNS: reference to list of the sites generated.
#########################################
    my ($count,  ## number of random sites to generate.
    $range   ## random sites lie between 1 and $range.
    ) = @_;
    my @sites;
    push @sites,0;
    foreach (2..$count-1)
    {
        push @sites, int(rand($range));
    }
    push @sites,$range;
    my @sites = sort { $a <=> $b } @sites;
    my $j=0;
    my $last = $range+1;
    foreach (my $i=1; $i<@sites; $i++)
    {  ## remove duplicates
        $sites[++$j] = $sites[$i] unless $sites[$i]==$sites[$j];
    }
    $#sites = $j;
    return \@sites;
}

#########################################
sub printPairwiseDistances
{
##  prints all n(n-1)/2 pairwise distances between the n sites in @$sites.
##  RETURNS: nothing.
#########################################
    my ($sites) = @_;   ## reference to a sorted list of integer sites.
    foreach (my $sep=1; $sep<@$sites; $sep++)
    {
        foreach (my $i=$#$sites-$sep; $i>=0; $i--)
        {
            print $$sites[$i+$sep]-$$sites[$i], " ";
        }
        print "\n";
    }
}
