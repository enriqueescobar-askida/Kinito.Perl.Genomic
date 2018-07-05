#!/usr/bin/perl -w -I . -I ../lib -I ../lib2
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################

#########################################
##  This program searches a sequence for satellites -- short sequences that
##  are repeated in tandem (in a row) with minor variation.
##  The algorithm is based on one due to Sagot and Myers.
#########################################

use strict;
use Util;

#########################################
## Default parameters; modifiable on first line of input file.
#########################################
my $maxSatLength = 15;  ## Don't search for satellites longer than this...
my $minSatLength = 5;   ## ...or shorter than this.
my $minSatRepeats = 5;  ## It must be repeated at least this many times.
my $epsilon = 0.10;     ## In edit distance, the repetitions must differ
    ## from the model at most this much times the length of the model.


#########################################
## MAIN PROGRAM -- call main subroutine
#########################################
main();

#########################################
sub main
{
##  reads input, calls findSatellites, reports results.
#########################################
    $|=1;             ## flush all output immediately.
    my $s = <STDIN>;  ## first line gives options
    eval $s;          ## sets options maxSatLength, etc.
    my $s = join('',<STDIN>); ## remainder of file is sequence
    $s =~ s/\s//sg;   ## remove all whitespace
    my $satellites = findSatellites(\$s);
    my @sats = sort { ($$a[1]+$$a[2])<=>($$b[1]+$$b[2]) } @$satellites;
    foreach (@sats)
    {
        print "($$_[1],$$_[2])~$$_[3] $$_[4]x $$_[0]\n";
    }
}

#########################################
sub findSatellites
{
##  maintains a stack of potential satellite suffixes.  Repeatedly pops
##  stack, expands suffix, and pushes valid expansions onto stack.
##  Suffixes that are complete satellites are remembered along with their
##  locations in the target string.
##  RETURNS: a list of maximal satellite chains.
#########################################
    my ($sref) = @_;   ## reference to string in which satellites are sought.

    ### Build match set for empty-string suffix & use it to prime pump.
    my @suffix0 = ("");
    my $slen = length $$sref;
    foreach my $i (0..$slen-1)
    {
        push @suffix0, $slen-$i,0;
    }
    my @suffixes = (\@suffix0);
    my @satellites;
    
    ### Main loop:  
    while (@suffixes)
    {
        my $suffix = pop @suffixes;
        next if length $$suffix[0] >= $maxSatLength;
    #   printf "Expanding /%s/\n", $$suffix[0];
        foreach my $base ('a','c','g','t')
        {
            my $eSuffix = suffixCheck($base, $suffix, $sref);
            if ($eSuffix)
            {
                push @suffixes, $eSuffix;
                push @satellites, @{satelliteCheck($eSuffix,$sref)};
            }
        }
    }
    return \@satellites;
}

#########################################
sub suffixCheck
{
##  Determines whether valid suffix $suffix is still valid if extended
##  by prepending $base.
##  RETURNS: match set for extended suffix.
#########################################
    my ($base,   ## nucleotide by which suffix is to be extended
    $suffix, ## known valid suffix (string)
    $sref    ## reference to target sequence
    ) = @_;
    my @eSuffix = ($base.$$suffix[0]);  ## string for extended suffix
    my $eLen = length $eSuffix[0];
    my $maxErr = int($epsilon * $maxSatLength);  ## permissible error

    for (my $p=1; $p<@$suffix; $p+=2)
    {
        my $k = $$suffix[$p]-1;
        last if $k<0;
        my $eq = (substr($$sref,$k,1) eq $base)? 0 : 1;
        my $eErr = $$suffix[$p+1]+$eq;   ## pair base to base
        $eErr = $$suffix[$p+3]+1         ## put gap in sequence 
            if ($$suffix[$p+2] == $k) && ($$suffix[$p+3]+1 < $eErr);
        $eErr = $eSuffix[$#eSuffix]+1    ## put gap in satellite 
            if ($eSuffix[$#eSuffix-1] == $k+1) 
            && ($eSuffix[$#eSuffix]+1<$eErr);
        push @eSuffix,$k,$eErr if ($eErr <= $maxErr);
    }

    my @jump = (max($eLen,$minSatLength)-$maxErr..$maxSatLength+$maxErr);
    my $filteredESuffix = filterForRepeats(\@eSuffix,$maxErr,\@jump);
#   print "/$eSuffix[0]/ ",($filteredESuffix? "accepted\n" : "rejected\n");
    return $filteredESuffix;
}


#########################################
sub filterForRepeats
{
##  filters a suffix's match set according to permissible deviation between
##  suffix and its repetitions and permissible distance between repetitions.
##  RETURNS:
#########################################
    my ($pattern,  ## the match set
    $maxErr,   ## permissible edit distance between suffix and repetitions.
    $jumps     ## ref to list of permissible distances between repetitions.
    ) = @_;
    my (@repeatsBefore, @repeatsAfter);
    $#repeatsBefore = $#repeatsAfter = $$pattern[1];  ## fast allocate
    for (my $p=$#$pattern-1; $p>0; $p-=2)
    {
        next if $$pattern[$p+1] > $maxErr;
        my $i = $$pattern[$p];
        $repeatsBefore[$i] = 1 + max(map($repeatsBefore[$i-$_],@$jumps));
    }
    my @filteredPattern = ($$pattern[0]);
    for (my $p=1; $p<$#$pattern; $p+=2)
    {
        next if $$pattern[$p+1] > $maxErr;
        my $i = $$pattern[$p];
        $repeatsAfter[$i] = 1 + max(map($repeatsAfter[$i+$_],@$jumps));
        push @filteredPattern,$i,$$pattern[$p+1]
        if $repeatsBefore[$i]+$repeatsAfter[$i] > $minSatRepeats;
    }
    return $#filteredPattern ? \@filteredPattern : "";
}

#########################################
sub satelliteCheck
{
##  Determines whether $sat actually qualifies as a satellite and not merely
##  a suffix.
##  RETURNS: a list of maximal chains for the satellite (if any).
#########################################
    my ($sat,  ## candidate satellite string
    $sref  ## reference to target string
    ) = @_;
    my $satString = $$sat[0];
    my $satLength = length $satString;
    return [ ] if ($satLength < $minSatLength);

    ## We have a fixed satellite now; error does not depend on $maxSatLength
    ## Clean out positions that don't pass the more stringent criteria.
    my $maxErr = int($epsilon * $satLength);
    my @jump = ($satLength-$maxErr .. $satLength+$maxErr);
    $sat = filterForRepeats($sat, $maxErr, \@jump);
    return [ ] unless $sat;

    ## Now we must look precisely at the individual substrings matching
    ## the satellite at each position and their chaining potential.
    my (@repeatsAfter,@stringEnd,@stringErr,%span);
    $#repeatsAfter = $#stringErr = $#stringEnd = $$sat[1];  ## fast allocate

    for (my $p=1; $p<@$sat; $p+=2)
    {
        my $i = $$sat[$p];
        my $goodJ = goodEndings($satString, $sref, $i, $maxErr);
        for (my $q=0; $q<@$goodJ; $q+=2)
        {
            my $j = $$goodJ[$q];
            if ($repeatsAfter[$j] >= $repeatsAfter[$i])
            {
                $repeatsAfter[$i] = 1+$repeatsAfter[$j];
                $stringEnd[$i] = $stringEnd[$j] || $j;
                $stringErr[$i] = $stringErr[$j] + $$goodJ[$q+1];
            }
        }
        next if ($repeatsAfter[$i] < $minSatRepeats);
        my $end = $stringEnd[$i];
        $span{$end} = $i
            if !defined($span{$stringEnd[$i]})
            || ($end-$i-$stringErr[$i] ## length minus error
                > $end-$span{$end}-$stringErr[$span{$end}]);
    }
#   print "/$satString/ final check\n";

    ## We return a list of the maximal chains for the satellite.
    return [map([$satString,$span{$_},$_,
        $stringErr[$span{$_}],$repeatsAfter[$span{$_}]],
        (keys %span)
        )];
}

#########################################
sub goodEndings
{
##  called by satelliteCheck to determine edit distances between
##  a suffix and all substrings beginning at a given position $i in
##  the target string.
##  RETURNS:  an even-length list of (ending position, edit distance) pairs.
#########################################
    my ($satString, ## the satellite (string)
    $sref,      ## ref to the target string
    $i,         ## starting point in target string
    $maxErr     ## permitted edit distance between satellite, repetitions
    ) = @_;
    my $satLength = length $satString;
    my @M;
    $M[0] = [0..$maxErr+1];
    foreach my $p (1..$satLength)
    {
        $M[$p][$p+$maxErr+1] = $maxErr+1; ## too big; allocates row.
        my $lo = max(1,$p-$maxErr);
        $M[$p][$lo-1] = $p-($lo-1);
        foreach my $q ($lo..$p+$maxErr)
        {
            if (substr($$sref,$i+$q-1,1) eq substr($satString,$p-1,1)) {
                $M[$p][$q] = $M[$p-1][$q-1];
            }
            else
            { 
                $M[$p][$q] = 
                1 + min($M[$p-1][$q-1],$M[$p][$q-1],$M[$p-1][$q]);
            }
        }
    }
    my @good;
    for (my $q=$satLength+$maxErr; $q>=$satLength-$maxErr; $q--)
    {
        push @good, $i+$q, $M[$satLength][$q] 
        if $M[$satLength][$q]<=$maxErr;
    }
    return \@good;
}
