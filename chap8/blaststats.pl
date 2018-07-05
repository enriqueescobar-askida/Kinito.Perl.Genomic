#!/usr/bin/perl -w -I . -I ../lib -I ../lib2
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################

#########################################
##  
#########################################


use strict;
use Util;
use Recurrence;

#########################################
## MAIN PROGRAM
#########################################

## Read background probabilities and scoring matrix.
my ($p,$s) = readScoreScheme();
## Compute, print average score per pair.  (Should be negative.)
my $avg = averagePairScore($p,$s);
print "\nAverage Score per Pair is $avg\n";

## Turn probabilities and scores into a recurrence relation.
my @terms;
foreach my $i (0..$#$p)
{
    foreach my $j (0..$#$p)
    {
        push @terms, $$s[$i][$j], $$p[$i]*$$p[$j];
    }
}
my $Psm = Recurrence->new(@terms);
$Psm->printRecurrence();

## Compute and print characteristic polynomial and its root.
print "\nCharacteristic Polynomial is\n";
$Psm->printPolynomial();
my $root = $Psm->polyRoot();
print "\nRelevant root is $root\n";

## Print two tables of values.
print "\nTable of Expected-s[m] and P[s,m]\n";
printTable($Psm, [1..20], [1..22]);
print "\nTable of Expected-s[m] and P[s,m]\n";
printTable($Psm, [1..20], [map(10*$_,(1..22))]);

## Print table of limit values.
print "\nTable of P[s,220] / ($root)^s\n";
foreach my $s (1..20)
{
    printf " %6d", $s;
}
print "\n";
my $power = 1;
foreach my $s (1..20)
{
    $power *= $root;
    printf " %6.4f", $Psm->value($s,220)/$power;
}
print "\n";


#########################################
sub averagePairScore
{     
## RETURNS: weighted average score per residue pair.
#########################################
    my ($p,   ## reference to an array of background probabilities
    $s)   ## reference to a scoring matix as an array.
    = @_;
    my $avg;
    foreach my $i (0..$#$p)
    {
        foreach my $j (0..$#$p)
        {
            $avg += $$p[$i]*$$p[$j]*$$s[$i][$j];
        }
    }
    return $avg;
}    

#########################################
sub readScoreScheme
{
## reads background probabilities and scoring matrix from STDIN. 
## RETURNS: list of references to two arrays: probs and score matrix.
#########################################
    ## Read background frequencies into @p.
    print "Background Frequencies (input)\n";
    $_ = <STDIN>;
    my @p = split;
    print;
    
    ## Read scoring matrix into @s; remember largest entry.
    print "\nScoring Matrix (input)\n";
    my @s;
    while (<STDIN>)
    {
        print;
        push @s, [split];
    }
    return (\@p, \@s);
}

#########################################
sub printTable
{
##  prints table of specific values of the function defined by the recurrence.
##  RETURNS: nothing.
#########################################
    my ($rec,     ## reference to a Recurrence object.
    $sbounds, ## reference to a list of scores to include in table.
    $mbounds) ## reference to a list of sequence lengths to include in table.
    = @_;
    print "\ns\\m";
    foreach my $m (@$mbounds)
    { 
        printf " %6d", $m;
    }
    print " ratio\n E:";
    foreach my $m (@$mbounds)
    {
        my $sum=0;
        for (my $s=1; $rec->value($s,$m)>0; $s++)
        {
            $sum += $rec->value($s,$m);
        }
        printf " %6.4f", $sum;
    }
    foreach my $s (@$sbounds)
    {
        printf "\n%2d:", $s;
        foreach my $m (@$mbounds)
        {
            printf " %6.4f", $rec->value($s,$m);
        }
        my $m = $$mbounds[-1];
        printf " %6.4f", 
                $rec->value($s,$m) / $rec->value($s-1,$m)
            if $s>1;
    }
    print "\n";
}
