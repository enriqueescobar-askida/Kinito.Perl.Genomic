﻿#!/usr/bin/perl -I . -I ../perllib
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
#########################################
## This program is a re-implementation of the Needleman-Wunsch algorithm
## for global alignment of DNA sequences.  In the original version,
## processing an entry meant assigning its final value by "pulling"
## information from its predecessors.  In this version, processing means
## "pushing" its final value into its successors.
## 
## Only subroutine similarity differs from Chapter 3.
#########################################

use strict;
use Util;

#########################################
## GLOBAL VARIABLES
#########################################
my @M;        ## alignment matrix; filled by similarity.
my $g = -2;   ## gap penalty.


#########################################
sub p
{  
## given two amino acids or nucleotides, compares
## them and returns a match reward (+1) or mismatch
## penalty (-1).  For amino acids, it should normally
## be replaced with a more complicated function.
## RETURNS: numerical reward/penalty.
#########################################
    my ($aa1, $aa2) = @_;  ## residues/bases to be compared.
    return ($aa1 eq $aa2)?1:-1;
}

#########################################
sub max
{
## Given any positive number of numerical arguments,
## RETURNS: the largest.
#########################################
    my ($m,@l) = @_; ## numerical values.
    foreach my $x (@l)
    {
        $m = $x if ($x > $m);
    }
    return $m;
}

#########################################
sub similarity
{
##  Determines score of best alignment of strings $s and $t
##  by filling in alignment matrix @M.
##  RETURNS: nothing; fills @M.
#########################################
    my($s,$t) = @_;
    foreach my $i (0..length($s))
    { 
        foreach my $j (0..length($t))
        {
            $M[$i][$j] = -999999;
        }
    }
    $M[0][0]=0;
    foreach my $layer (0..(length($s)+length($t)))
    { 
        my $i = min($layer,length($s)); 
        my $j = $layer-$i;
        while ($j<=length($t) && $i>=0)
        {
            $M[$i+1][$j] = max($M[$i+1][$j], $M[$i][$j]+$g);
            $M[$i][$j+1] = max($M[$i][$j+1], $M[$i][$j]+$g);
            my $p =  p(substr($s,$i,1),substr($t,$j,1));
            $M[$i+1][$j+1] = max($M[$i+1][$j+1], $M[$i][$j]+$p);
            $i--; $j++;
        }
    }
    return ( $M[length($s)][length($t)] );
}

#########################################
sub getAlignment
{
##  Reconstructs best alignment of strings $s and $t using information
##  stored in alignment matrix @M by similarity.  Recursive.
##  RETURNS: list of two strings representing best alignments.
##     These strings are $s and $t with gap symbols inserted.
#########################################
    my ($s,$t) = @_;  ## sequences to be aligned.
    my ($i,$j) = (length($s), length($t));
    return ( "-"x$j, $t) if ($i==0);
    return ( $s, "-"x$i) if ($j==0);
    my ($sLast,$tLast) = (substr($s,-1),substr($t,-1));
    
    if ($M[$i][$j] == $M[$i-1][$j-1] + p($sLast,$tLast))
    { ## Case 1
        ## last letters are paired in the best alignment
        my ($sa, $ta) = getAlignment(substr($s,0,-1), substr($t,0,-1));
        return ($sa . $sLast , $ta . $tLast );
    }
    elsif ($M[$i][$j] == $M[$i-1][$j] + $g)
    { ## Case 2
        ## last letter of the first string is paired with a gap
        my ($sa, $ta) = getAlignment(substr($s,0,-1), $t);
        return ($sa . $sLast , $ta . "-");
    }
    else
    { ## Case 3: last letter of the 2nd string is paired with a gap
        my ($sa, $ta) = getAlignment($s, substr($t,0,-1));
        return ($sa . "-" , $ta . $tLast );
    }
}

#########################################
##  MAIN PROGRAM
#########################################
{
    my ($s,$t) = ("SASKATCHEWAN", "SESQUICENTENNIAL");
    print "Similarity score: ", similarity($s,$t), "\n";
    print "Alignment: \n";
    foreach my $x (getAlignment($s,$t))
    {
        print $x,"\n";
    }
    print "\n";

    ## Read two strings, fill matrix @M, reconstruct alignment, print it.
    $s = <STDIN>; chomp($s);
    $t = <STDIN>; chomp($t);
    print "Similarity score: ", similarity($s,$t), "\n";
    print "Alignment: \n";
    foreach my $x (getAlignment($s,$t))
    {
        print $x,"\n";
    }
}
