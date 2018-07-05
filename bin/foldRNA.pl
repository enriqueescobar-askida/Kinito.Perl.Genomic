#!/usr/bin/perl -w -I . -I ../lib -I ../lib2
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
use		strict;
use		diagnostics;

#########################################
## GLOBAL VARIABLES
#########################################
my		@s;  ## holds input string, split into single characters.
my		@c;  ## holds dynamic programming table.
			## number of hydrogen bonds different base pairs can form:
my		%bonds = (GU=>1,UG=>1,AU=>2,UA=>2,CG=>3,GC=>3);
    
#########################################
##  RETURNS: the larger of its two arguments.
#########################################
sub	max
{
	my ($x,$y) = @_;
	if ($x>$y)
	{
		return	$x;
	}
	else
	{
		return	$y;
	}
}

#########################################
##  Implements the dynamic programming scheme to determine the maximum
##  number of hydrogen bonds that can be formed by folding an RNA string.
##  RETURNS: nothing; fills @c.
#########################################
sub	foldRna
{
	my ($s) = @_;  ## the RNA to be folded, in form of a string.
	my $slen = length $s;
	@s = ('X', split(//, $s));

	for (my $len = 5; $len <= $slen; $len++)
	{
		for (my $i=1; $i<=$slen-$len+1; $i++)
		{
			my $j = $i+$len-1;
			$c[$i][$j] = max($c[$i+1][$j],
							 $bonds{$s[$i].$s[$j]}+$c[$i+1][$j-1]);

			for (my $k=$i+1; $k<$j; $k++)
			{
				$c[$i][$j] = max($c[$i][$j], 
								 $c[$i][$k]+$c[$k+1][$j]);
			}
		}
	}
}

#########################################
##  Uses the contents of the dynamic programming table @c to construct
##  a string of parentheses and dots describing the optimal folding of
##  the RNA string @s[$i..$j].  Recursive.
##  RETURNS: the string of parentheses/dots.
#########################################
sub	traceBack
{
	my ($i,$j) = @_;  ## left and right boundaries of substring being folded.
	my $cij = $c[$i][$j];

	return	("."  x ($j-$i+1)) if ($cij==0);
	return	"." . traceBack($i+1,$j) 
		if ($cij == $c[$i+1][$j]);
	return	"(" . traceBack($i+1,$j-1) . ")" 
	if ($cij == $bonds{$s[$i].$s[$j]}+ $c[$i+1][$j-1]);
	for (my $k = $i+1; $k < $j; $k++)
	{
		return	traceBack($i,$k) . traceBack($k+1,$j)
		if ($cij == ($c[$i][$k]+$c[$k+1][$j]));
	}
}

#########################################
## MAIN PROGRAM
#########################################
$|=1;
my $basestring = <STDIN>;
chomp($basestring);
foldRna($basestring);
print "$basestring\n", traceBack(1, length $basestring), "\n";
