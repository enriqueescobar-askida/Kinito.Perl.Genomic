#!/usr/bin/perl -w -I . -I ../lib -I ../lib2
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################

#########################################
## Simple driver for SimpleReader.
#########################################

use strict;
use SimpleReader;


my $r = SimpleReader->new("file.txt");
while (my $seq = $r->readSeq()) {
	my $id = $r->seqId();
	print "\n\n", $r->fileName(), ", $id:\n";
	print "$seq\n";
}
$r->close();
