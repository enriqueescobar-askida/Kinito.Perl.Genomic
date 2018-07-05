################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
package Recurrence;
#########################################
##  This package implements a class of objects representing a class of 
##  recurrence relations that come up in the analysis of the statistics
##  of the scores of gap-free alignments.
#########################################

use strict;


#########################################
sub new
##  creates a new object of class Recurrence to represent a recurrence
##  relation with specified terms.
##  For example, the call
##     Recurrence->new(2,3,4,5,6,7)
##  creates an object to represent
##     m[n] = -3*m[n-2] + -5*m[n-4] + -7*m[n-6]
##  RETURNS: bless reference to hash representing recurrence.
#########################################
{
	my ($this,  ## Normally, the literal "Recurrence".
	@terms) ## Terms of the right hand side of the recurrence,
				## as described above.
	= @_;
	my %coefficients;
	while (@terms) {
	(my $d, my $c, @terms) = @terms;
	$coefficients{0-$d} += $c;
	}
	my @coefficients 
	= map([$_,$coefficients{$_}], sort {$a<=>$b} keys %coefficients);
	my @polynomial;
	my $offset = $coefficients[0][0];
	$polynomial[-$offset] = -1.0;
	foreach (@coefficients) {
	my ($d,$c) = @$_;
	$polynomial[$d-$offset] += $c;
	}
	my $hash = {coefficients => \@coefficients,
		polynomial => \@polynomial,
		memory =>  { }};
	bless $hash, (ref($this) || $this);
}

#########################################
sub printPolynomial {
##  prints the characteristic polynomial of the recurrence.
##  RETURNS: nothing.
#########################################
	my ($this) = @_;
	my $poly = $this->{polynomial};
	print "$$poly[-1]*x^$#$poly +";
	for (my $i=$#$poly-1; $i>=0; $i--) {
	next if $$poly[$i]==0;
	print " $$poly[$i]";
	if ($i>1) { print "*x^$i +"; }
	elsif ($i==1)  { print "*x +"; }
	else { print " "; }
	}
	print "= 0\n";
}

#########################################
sub printRecurrence {
##  prints the recurrence relation.
##  RETURNS: nothing.
#########################################
	my ($this) = @_;
	print "\nRecurrence Relation is\n";
	my $leader =  "P[s,m] = ";
	foreach (@{$this->{coefficients}}) {
	my ($d,$c) = @$_;
	print "$leader $c * P[s";
	printf "%+0d", $d if $d!=0;
	print ",m-1] \n";
	$leader =  "        +";
	}
}

#########################################
sub polyRoot {  
##  finds and returns the "important" root of the characteristic polynomial
##  of the recurrence relation.  This is the unique real root in the interval
##  (0,1).
##  This is not a general routine for solving polynomials!
##  It exploits special properties of the recurrences arising in
##  analyzing gap-free alignments.
##  RETURNS: The root.
#########################################
	my ($this) = @_;
	my $poly = $this->{polynomial};
	my ($neg,$pos) = (1,0);
	while ($neg-$pos > 1E-15) {
	my $x = ($neg+$pos)/2.0;
	my $value = 0;
	my $power = 1;
	foreach my $coeff (@$poly) {
		$value += $coeff * $power;
		$power *= $x;
	}
	if ($value>0) { $pos = $x; } else { $neg = $x; }
	}
	return ($neg+$pos)/2.0;
}

#########################################
sub value {
##  Computes a specific value of the function defined by the
##  recurrence relation.  Saves value in $this->{memory}.
##  RETURNS: the value.
#########################################
	my ($this,$s,$m) = @_;
	my $mem = $this->{memory};
	return $$mem{$s,$m} if defined($$mem{$s,$m});
	return 1 if $s<=0 && $m>=0;
	return 0 if $m<=0;
	my $val;
	foreach (@{$this->{coefficients}}) {
	my ($d,$c) = @$_;
	$val += $c * $this->value($s+$d,$m-1)
	}
	return $$mem{$s,$m}=$val;
}

1;
