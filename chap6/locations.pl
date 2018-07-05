#!/usr/bin/perl -w -I . -I ../lib -I ../lib2
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################

#########################################
##  This program illustrates how the subroutine extract_feature operates
##  by printing out the successive stages of the transformation of a 
##  GenBank feature location as it is transformed into a Perl regular 
##  expression.
#########################################
use strict;

#########################################
sub complement
{
##  RETURNS: the reverse complement of its string argument.
#########################################
    my ($s) = @_;
    $s =~ tr/acgtACGT/tgcaTGCA/;
    return reverse $s;
}

#########################################
sub ident
{
##  RETURNS: its argument without modification.
#########################################
    my ($s) = @_;
    return $s;
}

#########################################
sub extract_feature
{
##  RETURNS: the DNA sequence extracted from $s at location $loc.
#########################################
    my ($loc,   ## a GenBank feature location
    $s)     ## a sequence from which the feature is to be extracted.
    = @_;
    print @_, "\n";
    ## We don't want to create lots of copies of a long string.
    unless (ref $s)
    {
        my $t = $s; $t =~ s/\s//g; $s = \$t;
    }
    $loc = lc $loc;
    $loc =~ s/,/\./g;
    print "$loc\n";
    $loc =~ s/(\d+)/N$1/g;
    print "$loc\n";
    $loc =~ s/join/ident/g;
    print "$loc\n";
    $loc =~ s/N(\d+)\.\.N(\d+)/substr(\$\$s,$1-1,$2+1-$1)/g;
    print "$loc\n";
    $loc =~ s/N(\d+)/substr(\$\$s,$1-1,1)/g;
    print "$loc\n";
    print eval($loc),"\n";
    print $@, "\n";
    return eval($loc);
}



#########################################
##  MAIN PROGRAM
#########################################
my $seq = "ourfatherwhoartinheavenhallowedbethyname";
my $location = "join(18..24,complement(join(5..10,6..11)),27,40..45)";

extract_feature($location,$seq);

