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
my		@structure;
my		@bases;

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
##  Determines the number of hydrogen bonds in the folded RNA described
##  by @bases[$l..$r] and @structure[$l..$r].  Recursive.
##  RETURNS: the number of hydrogen bonds.
#########################################
sub	evalRna
{
	my ($l,  ## the left boundary of the subarrays of @bases and @structure.
		$r)  ## the right boundary of the subarrays of @bases and @structure.
		= @_;

	my %bonds = (GU=>1,UG=>1,AU=>2,UA=>2,CG=>3,GC=>3);
	my $energy = $bonds{$bases[$l].$bases[$r]};
	my $level=0;
	
	my $ii = $l;
	
	for (my $i=$l+1; $i<=$r; $i++)
	{
		$level-- if ($structure[$i] eq ")");
		if ($level==0)
		{
			$energy += evalRna($ii,$i) if ($structure[$i] eq ")");
			$ii = $i;
		}
		$level++ if ($structure[$i] eq "(");
	}
	return	$energy;
}

#########################################
##  Determines the number of hydrogen bonds in the folded RNA described
##  by strings $basestring and $structurestring.
##  RETURNS: the number of hydrogen bonds.
#########################################
sub	evalRnaStructure
{
	my ($basestring,$structurestring) = @_;
	@bases = split(//, 5 . $basestring . 3);
	@structure = split(//, "($structurestring)");
	return	evalRna(0, $#structure);
}

#########################################
## MAIN PROGRAM
#########################################
## Read string of bases, string of parentheses.  Find number of H bonds.
my $basestring = <STDIN>;  
chomp($basestring);
my $parenstring = <STDIN>;
chomp($parenstring);
print evalRnaStructure($basestring,$parenstring), " hydrogen bonds in this structure.\n";
