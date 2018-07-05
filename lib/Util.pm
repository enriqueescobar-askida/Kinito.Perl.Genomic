################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
package	Util;

require	Exporter;
@ISA	= qw(Exporter);
@EXPORT	= qw(max min sum lg log10 pow round);

use	strict;
use	diagnostics;

sub	max
{
	my ($m,@l) = @_;
	foreach (@l) { $m=$_ if $_>$m };
	$m;
}

sub	min
{
	my ($m,@l) = @_;
	foreach (@l) { $m=$_ if $_<$m };
	$m;
}

sub	sum
{
	my ($s,@l) = @_;
	foreach (@l) { $s += $_ };
	$s;
}

sub	lg
{
	log(shift)/log(2.0);
}

sub	log10
{
	log(shift)/log(10.0);
}

sub	pow
{ 
	my ($x,$y) = @_;
	exp($y * log($x));
}

sub	round
{
	int(0.5+shift);
}
