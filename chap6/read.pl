#!/usr/bin/perl -w -I . -I ../lib -I ../lib2
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
#########################################
##  This program is a simple driver to illustrate the use
##  of the SeqReader packages.
#########################################
use strict;
use diagnostics										;
use SeqReader;
#########################################
##  Create readers for three files, one each GenBank, FASTA, and Simple.
#########################################
my @readers;
push @readers, new SeqReader("file.gb");
push @readers, new SeqReader("turkey.fasta");
push @readers, new SeqReader("file.txt");
#########################################
##  Read two sequences from each file, then closes them.
#########################################
foreach my $r (@readers) {
	print "\n\nFILE: ", $r->fileName(), ":  ", $r->fileType(), "\n";
	my $seq = $r->readSeq();
	my $id = $r->seqId();
	print "ID=$id\nSEQ=$seq\n\n";
	my $seq = $r->readSeq();
	my $id = $r->seqId();
	print "ID=$id\nSEQ=$seq\n\n";
	$r->close();
}
#########################################
##  This time, require the files to be GenBank files.
#########################################
push @readers, new GenBankReader("file.gb");     ## succeeds
push @readers, new GenBankReader("file.fasta");  ## fails; aborts program
