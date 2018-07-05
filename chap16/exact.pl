#!/usr/bin/perl -I . -I ../perllib
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
use strict;
use Util;

#########################################
##  This program implements a backtracking approach to solving the 
##  partial-digest problem when the input distances are exact.
#########################################

#########################################
##  GLOBAL VARIABLES
#########################################
my %distAvail;  ## A hash table whose keys are the inter-site distances
                ## of the input and whose values are the number of occurrences
                ## of the distance unclaimed by the current partial placement.
my @distances;  ## A sorted list of the keys of %distAvail.
my @solutions;  ## A list of placements found to be consistent with the input.
my @sites;      ## The current partial placement.

#########################################
##  MAIN PROGRAM
#########################################
foreach (map(split(/\D+/,$_), <STDIN>))
{
    $distAvail{$_}++;
};
@distances = sort {$a <=> $b} keys %distAvail;
my $numDists;
foreach (values %distAvail)
{
    $numDists += $_;
}
@sites = (("-") x (1+int(sqrt(2*$numDists))));
$sites[0] = 0;
$sites[-1] = $distances[-1];
$distAvail{$distances[-1]}--;
placeSite(1,$#sites-1);
foreach (@solutions)
{
    print "Solution: @$_\n";
}
exit();

#########################################
sub placeSite
{
##  attempts to place the next restriction site according to the
##  largest remaining distance.  Recursive, but indirectly through
##  verifyDistances.  
#########################################
    my ($l,$r) = @_;  ## first / last available position in global @sites
    my $i = $#distances;
    while (($i>=0) && !$distAvail{$distances[$i]})
    {
        $i--;
    }
    if ($i<0)
    {   ## All distances accounted for.
        push @solutions, [@sites];  ## Makes copy of @sites.
        return;
    }
    my $dist = $distances[$i];
    my $newSite = $sites[-1] - $dist;
    verifyDistances($l, $r, "L", $newSite);
    verifyDistances($l, $r, "R", $dist);
}

#########################################
sub verifyDistances
{
##  verifies whether available distances would allow a new restriction
##  site at $newSite.  If so, places site and calls placeSite to work on
##  next site.
#########################################
    my ($l,$r,   ## first / last available position in global @sites
    $LR,     ## attempting to place on left or right end? "L" or "R"
    $newSite ## proposed location of placement
    ) = @_;

    my $n = ($LR eq "L")? $l : $r;
    my $allAvail=1;
    foreach my $j ((0..$l-1), ($r+1..$#sites))
    {
        $allAvail = $distAvail{abs($newSite-$sites[$j])}-- && $allAvail;
    }
    if ($allAvail)
    {  ### We have found all the distances we need!
        $sites[$n] = $newSite;
        if ($LR eq "L")
        {
            placeSite($l+1, $r);
        }
        else
        {
            placeSite($l, $r-1);
        }
        $sites[$n] = "-";
    }
    foreach my $j ((0..$l-1), ($r+1..$#sites))
    {
        $distAvail{abs($newSite-$sites[$j])}++;
    }
}
