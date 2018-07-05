#!/usr/bin/perl -w -I . -I ../lib -I ../lib2
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################

#########################################
##  This program implements a direct extension of the Needleman-Wunsch 
##  algorithm to align three sequences.
#########################################

use strict;
use Util;

my @M;    ## scoring matrix
my $g = -2;   ## gap penalty

#########################################
sub p
{
##  computes reward/penalty for putting residues $aa1 and $aa2 in same column.
##  RETURNS: reward/penalty.
#########################################
    my ($aa1, $aa2) = @_;
    return (($aa1 eq $aa2)?1:-1);
}

#########################################
sub similarity
{
##  uses direct dynamic programming approach to determine the score of the
##  best possible multiple alignment of three sequences.
##  RETURNS: score of best alignment.
#########################################
    my($s1,$s2,$s3) = @_;  ## the 3 sequences being aligned (strings)

    ### Fill in edges of cube.
    foreach my $i1 (0..length($s1))
    {
        $M[$i1][0][0] = $g * $i1 * 2;
    }
    foreach my $i2 (0..length($s2))
    {
        $M[0][$i2][0] = $g * $i2 * 2;
    }
    foreach my $i3 (0..length($s3))
    {
        $M[0][0][$i3] = $g * $i3 * 2;
    }

    ### Fill in sides of cube.
    ## Side 1
    foreach my $i1 (1..length($s1))
    {
        my $aa1 = substr($s1,$i1-1,1);
        foreach my $i2 (1..length($s2))
        {
            my $aa2 = substr($s2,$i2-1,1);
            $M[$i1][$i2][0] = max($M[$i1-1][$i2][0] + $g + $g,
                    $M[$i1][$i2-1][0] + $g + $g,
                    $M[$i1-1][$i2-1][0] + $g + p($aa1,$aa2));
        }
    }

    ## Side 2
    foreach my $i1 (1..length($s1))
    {
        my $aa1 = substr($s1,$i1-1,1);
        foreach my $i3 (1..length($s3))
        {
            my $aa3 = substr($s3,$i3-1,1);
            $M[$i1][0][$i3] = max($M[$i1-1][0][$i3] + $g + $g,
                $M[$i1][0][$i3-1] + $g + $g,
                $M[$i1-1][0][$i3-1] + $g + p($aa1,$aa3));
        }
    }
    
    ## Side 3
    foreach my $i2 (1..length($s2))
    {
        my $aa2 = substr($s2,$i2-1,1);
        foreach my $i3 (1..length($s3))
        {
            my $aa3 = substr($s3,$i3-1,1);
            $M[$i2][0][$i3] = max($M[0][$i2-1][$i3] + $g + $g,
                $M[0][$i2][$i3-1] + $g + $g,
                $M[0][$i2-1][$i3-1] + $g + p($aa2,$aa3));
        }
    }
    

    ### Fill in interior of cube.
    foreach my $i1 (1..length($s1))
    {
        my $aa1 = substr($s1,$i1-1,1);
        foreach my $i2 (1..length($s2))
        {
            my $aa2 = substr($s2,$i2-1,1);
            my $p12 = p($aa1,$aa2);
            foreach my $i3 (1..length($s3))
            {
                my $aa3 = substr($s3,$i3-1,1);
                my $p13 = p($aa1,$aa3);
                my $p23 = p($aa2,$aa3);
                $M[$i1][$i2][$i3] = 
                max($M[$i1-1][$i2-1][$i3-1] + $p12 + $p13 + $p23,
                $M[$i1-1][$i2-1][$i3] + $p12 + $g + $g,
                $M[$i1-1][$i2][$i3-1] + $g + $p13 + $g,
                $M[$i1][$i2-1][$i3-1] + $g + $g + $p23,
                $M[$i1][$i2][$i3-1] + 0 + $g + $g,
                $M[$i1][$i2-1][$i3] + $g + 0 + $g,
                $M[$i1-1][$i2][$i3] + $g + $g + 0);
            }
        }
    }
    return ( $M[length($s1)][length($s2)][length($s3)] );
}


#########################################
## MAIN PROGRAM
#########################################
{
    my $s1 = <DATA>; chomp($s1);
    my $s2 = <DATA>; chomp($s2);
    my $s3 = <DATA>; chomp($s3);
    print "Similarity score: ", similarity($s1,$s2,$s3), "\n";
}


__END__
GVLTDVQVALVKSSFEEFNANIPKNTHRFFTLVLEIAPGAKDLFSFLKGSSVHLSGGEKSAVTNLWGKVNINELGGEALGRLLVVYPWTQRFFEAFGDLS
VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSVLSAADKTNVKGVFSKIGGHAEEYGAETLERMFIAYPQTKTYFPHFDLS
VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSSPLTADEASLVQSSWKAVSHNEVEILAAVFAAYPDIQNKFSQFAGK



