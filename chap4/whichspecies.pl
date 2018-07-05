#!/usr/bin/perl -w -I . -I ../lib -I ../lib2
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
sub entropy {
##  Computes entropy in bits of a distribution given as list of probabilities.
##  RETURNS: the entropy.
#########################################
	my $eta=0;
	foreach (@_) { $eta -= $_ * lg($_); }
	return $eta;
}

#########################################
##  MAIN PROGRAM
#########################################
##  Read name, probabilities of A,C,G,T from DATA.
my ($species1, @model1) = split /\s+/, <DATA>;    
my ($species2, @model2) = split /\s+/, <DATA>;
##  Compute, print entropy for each species.
print "Entropy for $species1 is ", entropy(@model1), " bits.\n";
print "Entropy for $species2 is ", entropy(@model2), " bits.\n";
##  Compute lod scores for A, C, G, T;
my @lod = map(lg($model1[$_]/$model2[$_]), (0..3));

##  Make a species prediction for each line of the input file.
while ( my $dna = <STDIN> ) {
	chomp($dna);
	$dna =~ tr/acgtACGT/01230123/;
	my $score = sum(map($lod[$_], split //,$dna));
	print "lod=$score bits favors ";
	print "$species1 over $species2.\n" if $score>=0;
	print "$species2 over $species1.\n" if $score<0;
}

__END__
corn     .255 .245 .245 .255
fruitfly .225 .275 .275 .225 
