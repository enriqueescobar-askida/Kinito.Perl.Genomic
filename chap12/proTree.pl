#!/usr/bin/perl -I . -I ../perllib
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################

#########################################
##  Driver to produce a list of PROSITE motifs found in each of 
##  a number of input sequences found in a protein file.
##  Uses suffix trees.
#########################################

use strict;
use Prosite;
use PrositeSuffixTree;
use SeqReader;


##  Read all protein sequences and build a suffix tree.
my $seqFile = new SeqReader("turkey.fasta");
my $allSeqs;
while (my $protein = $seqFile->readSeq())
{
    $allSeqs .= "<$protein>"
}
$seqFile->close();

my $ptree = new PrositeSuffixTree($allSeqs."=");


##  Now read each motif in PROSITE one by one, and search for occurrences
##  of each in the big suffix tree.
my $prosite = new Prosite("prosite.dat");
while (my $motif = $prosite->readMotif())
{
    $ptree->listMotifs($motif);
}
$prosite->close();


