#!/usr/bin/perl -I . -I ../perllib
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################

use strict;
use Prosite;
use SeqReader;

#########################################
##  MAIN PROGRAM
##  
##  Driver to produce a list of PROSITE motifs found in each of 
##  a number of input sequences found in a protein file.
#########################################
my $seqFile = new SeqReader("turkey.fasta");

my $protein;
while ($protein = $seqFile->readSeq())
{
    print "\n$protein\n";
    my $prosite = new Prosite("prosite.dat");
    while (my $motif = $prosite->readMotif())
    {
        listMotifs($protein,$motif);
    }
    $prosite->close();
}
$seqFile->close();


#########################################
sub listMotifs
{
##  prints one or more occurences of a single motif in a single sequence.
#########################################
    my ($protein,  ## query protein
    $motif)    ## single motif
    = @_;
    my $pattern = $$motif{PA};
    return unless $pattern;
    return if $pattern =~ /,/;   ## skip for fair comparison to protree.
    my $description = $$motif{DE};
    $description =~ s/\n//gs;   ## remove newlines
    $pattern =~ s/\n//gs;   ## remove newlines
    $pattern = lc Prosite->perlizePattern($pattern);
    while ($protein =~ /($pattern)/g) {   ## case-insensitive
    print "$description seen as $1 at position ";
    print 1+index($protein, $1);
    print "\n";
    }
}


