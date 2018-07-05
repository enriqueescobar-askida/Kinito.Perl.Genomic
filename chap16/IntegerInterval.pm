################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
package IntegerInterval;
#########################################
## This package implements arithmetic on intervals with integer
## endpoints.  Multiplication and division are not implemented.
#########################################

use strict;
use Util;

#########################################
sub new {
##  creates a new interval object representing the largest interval
##  with integer endpoints contained in $lo..$hi.  ($lo and $hi may
##  be floating-point.
##  RETURNS:  reference to blessed LIST reference.
#########################################
	my ($this,     ##  literal "IntegerInterval" or reference to an IntegerInterval.
	$lo, $hi   ##  limits of new interval.
	) = @_;
	my $iLo = int($lo);
	$iLo++ if $iLo < $lo;
	my $iHi = int($hi);
	$iHi-- if $iHi > $hi;
	bless [$iLo, $iHi];
}

#########################################
sub intersect {
##  intersects two intervals.
##  RETURNS: new IntegerInterval object representing intersection
#########################################
	my ($this, $other) = @_;   ## the two intervals
	my ($a,$b) = @$this;
	my ($c,$d) = @$other;
	$a>$c or $a=$c;
	$b<$d or $b=$d;
	$b>=$a or return undef;
	return $this->new($a,$b);
}

#########################################
sub plus {
##  adds two intervals.
##  RETURNS: new IntegerInterval object representing sum
#########################################
	my ($this, $other) = @_;   ## the two intervals
	my ($a,$b) = @$this;
	my ($c,$d) = @$other;
	return $this->new($a+$c,$b+$d);
}

#########################################
sub minus {
##  subtracts two intervals.
##  RETURNS: new IntegerInterval object representing difference
#########################################
	my ($this, $other) = @_;   ## the two intervals
	my ($a,$b) = @$this;
	my ($c,$d) = @$other;
	return $this->new($a-$d,$b-$c);
}
	
#########################################
sub print {
##  prints the current interval, followed by any other arguments.
##  RETURNS: nothing.
#########################################
	my ($this,     ## current interval 
	@s         ## other things to print
	) = @_;
	$$this[2] ||= "$$this[0]:$$this[1]";
	print "[$$this[2]]@s";
}

#########################################
sub toString {
##  converts an interval to a printable string form, e.g., "13:34".
##  RETURNS: the string.
#########################################
	my ($this) = @_;
	$$this[2] ||= "$$this[0]:$$this[1]";
}

1;
