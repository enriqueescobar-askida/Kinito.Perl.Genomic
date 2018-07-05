################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
#########################################
package PAMBuilder;
##  This package provides subroutines to manipulate residue frequency matrices,
##  mutation rate matrices, and PAM matrices.
#########################################
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(printMatrix perlFormatMatrix collectCounts 
		addPseudocounts freqsFromCounts ratesFromFreqs
		rescaleRates freqsFromRates lodsFromFreqs
		averageMutationRate multiplyRates);

use strict;
use Util;
my @residues = split //,"acdefghiklmnpqrstvwy";

#########################################
sub checkType {
##  Used to check whether a matrix is of the appropriate type for the operation
##  requested.  Not exported by the package.
##  RETURNS: The matrix type if it is appropriate; ABORTS if not.
#########################################
	my ($M,    ## reference to hash representing matrix.
	$type) ## pattern describing appropriate matrix type(s).
	= @_;
	$$M{type} =~ /($type)/;
	$1 or die "Matrix is of type $$M{type}; $type is required.";
}

#########################################
sub printMatrix {
##  Prints a matrix with specified heading.
##  RETURNS: nothing.
#########################################
	my ($M,      ## reference to hash representing matrix.
	$header) ## header for printing.
	= @_; 
	my %format = (frequencies=>" %4.2f",'mutation rates'=>" %4.2f",
		counts=>"%5d", 'lod scores'=>"%5d");
	my $format = $format{$$M{type}};
	print "\#\#\# $header\n\#  bkgnd";
	foreach (@residues) { printf "%5s", $_; }
	foreach my $r1 (@residues) {
	printf "\n# $r1$format", $$M{$r1};
	foreach my $r2 (@residues) {
		printf $format, $$M{$r1.$r2};
	}
	}
	print "\n\n";
}

#########################################
sub perlFormatMatrix {
##  Prints a matrix in a form that allows it to be re-read later.
##  RETURNS: nothing.
#########################################
	my ($M,        ## reference to hash representing matrix.
	$name)     ## string giving Perl variable name of matrix.
	= @_;
	print "$name=\n";
	my $s = join(",", map("$_=>'$$M{$_}'", (keys %$M)));
	print "(";
	while ($s =~ s/^(.{0,70}\,)//) { print " $1\n";}
	print " $s);\n";
}

#########################################
sub collectCounts {
## Reads alignments from STDIN.  Counts occurrences of residue and
## residue pairs.
## RETURNS: a reference to a hash representing a "counts" matrix.
#########################################
	my %nu = (type=>'counts');
	while (chomp(my $seq1 = <STDIN>)) {
	chomp(my $seq2 = <STDIN>);    
	foreach my $i (0..((length $seq1)-1)) {
		my ($r1,$r2) = (substr($seq1,$i,1), substr($seq2,$i,1));
		next if ($r1.$r2) =~ "-";
		$nu{$r1}++; $nu{$r2}++;
		$nu{$r1.$r2}++; $nu{$r2.$r1}++;
	}
	}
	return \%nu;
}

#########################################
sub addPseudocounts {
##  Adds pseudo-counts to a counts matrix.
##  RETURNS: new counts matrix including pseudo-counts.
#########################################
	my ($M,      ## reference to hash representing matrix.
	$pseudo) ## optional count to add to each entry; default 1.
	= @_;   
	$pseudo ||= 1;
	checkType($M, 'counts');
	my %nu = %$M;
	foreach my $r1 (@residues) {
	foreach my $r2 (@residues) {
		$nu{$r1} += $pseudo ; $nu{$r2} += $pseudo;
		$nu{$r1.$r2} += $pseudo; $nu{$r2.$r1} += $pseudo;
		last if ($r1 eq $r2);
	}
	}
	return \%nu;
}

#########################################
sub freqsFromCounts {
##  Derives a frequency matrix (entries are fractions of total occurrences)
##  from a counts matrix (entries are raw whole-number counts of occurrences).
##  RETURNS: new frequency matrix.
#########################################
	my ($M) = @_;        ## reference to hash representing matrix.
	checkType($M, 'counts');
	my %nu = (type=>'frequencies');
	my $total=0;
	foreach my $r (@residues) { $total += $$M{$r}; }
	foreach my $r (@residues) { $nu{$r} = $$M{$r} / $total; }
	$total=0;
	foreach my $r1 (@residues) {
	foreach my $r2 (@residues) { 
		$total += $$M{$r1.$r2};
	}
	}
	foreach my $r1 (@residues) {
	foreach my $r2 (@residues) { 
		$nu{$r1.$r2} = $$M{$r1.$r2}/$total;
	}
	}
	return \%nu;
}

#########################################
sub ratesFromFreqs {
##  Derives a mutation-rate matrix from a frequency matrix.
##  RETURNS: new mutation-rate matrix.
#########################################
	my ($M) = @_;        ## reference to hash representing matrix.
	checkType($M, 'frequencies');
	my %nu = (type=>'mutation rates');
	foreach my $r1 (@residues) {
	my $q = $nu{$r1} = $$M{$r1};
	foreach my $r2 (@residues) { $nu{$r1.$r2} =  $$M{$r1.$r2} / $q; }
	}
	return \%nu;
}

#########################################
sub freqsFromRates {
##  Derives a frequency matrix from a mutation-rate matrix.
##  RETURNS: new frequency matrix.
#########################################
	my ($M) = @_;        ## reference to hash representing matrix.
	checkType($M, 'mutation rates');
	my %nu = (type=>'frequencies');
	foreach my $r1 (@residues) {
	my $q = $nu{$r1} = $$M{$r1};
	foreach my $r2 (@residues) { $nu{$r1.$r2} =  $$M{$r1.$r2} * $q; }
	}
	return \%nu;
}

#########################################
sub averageMutationRate {
##  Derives the average mutation rate of a frequency or mutation-rate matrix.
##  RETURNS: the average mutation rate (a real number).
#########################################
	my ($M) = @_;        ## reference to hash representing matrix.
	my $type = checkType($M, 'mutation rates|frequencies');
	my $rate = 1;
	if ($type eq "frequencies") {
	foreach my $r (@residues) { $rate -= $$M{$r.$r}; }
	} else {
	foreach my $r (@residues) { $rate -= $$M{$r} * $$M{$r.$r}};
	}
	return $rate;
}

#########################################
sub rescaleRates {
##  Rescales the off-diagonal entries of a mutation-rate matrix so that
##  their ratios remain constant but the average mutation rate of the matrix
##  is adjusted to $newrate.
##  RETURNS: a new mutation-rate matrix.
#########################################
	my ($M,        ## reference to hash representing matrix.
	$newrate)  ## mutation rate of the matrix to be created.
	= @_;
	checkType($M, 'mutation rates');
	my %nu = (type=>'mutation rates');
	my $factor = $newrate/averageMutationRate($M);
	foreach my $r1 (@residues) {
	foreach my $r2 (@residues) { 
		$nu{$r1.$r2} = $$M{$r1.$r2} * $factor;
	}
	$nu{$r1.$r1} += (1-$factor);
	$nu{$r1} = $$M{$r1};
	}
	return \%nu;
}

#########################################
sub multiplyRates {
##  Derives a new mutation-rate matrix from %$M1 and %$M2.  The new matrix
##  represents a period of mutation at rates given in %$M1 followed by a period
##  of mutation at rates given in %$M2.
##  RETURNS: a new mutation-rate matrix.
#########################################
	my ($M1,$M2) = @_;        ## reference to hashes representing the matrices.
	checkType($M1, 'mutation rates');
	checkType($M2, 'mutation rates');
	my %nu = (type=>'mutation rates');
	foreach my $r1 (@residues) {
	foreach my $r2 (@residues) { 
		my $prod;
		foreach my $r (@residues) {
		$prod += $$M1{$r1.$r} * $$M2{$r.$r2};
		}
		$nu{$r1.$r2} = $prod;
	}
	$nu{$r1} = $$M1{$r1};
	}
	return \%nu;
}

#########################################
sub lodsFromFreqs {
##  Derives a matrix of log-of-odds scores from a matrix of frequencies.
##  Lod scores are given in half-bits.
##  RETURNS: a new lod-score matrix.
#########################################
	my ($M) = @_;        ## reference to hash representing matrix.
	checkType($M, 'frequencies');
	my %nu = (type=>'lod scores');
	foreach my $r1 (@residues) {
	$nu{$r1} = round(2*lg($$M{$r1}));
	foreach my $r2 (@residues) { 
		$nu{$r1.$r2} =
		round(2*lg($$M{$r1.$r2} / ($$M{$r1}*$$M{$r2})));
	}
	}
	return \%nu;
}

