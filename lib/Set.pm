################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################

package	Set;

use	Util;

#########################################
## Implements set operations: union, intersect, setMinus, disjoint, member,
##    newSet, subset, superset
#########################################

use	strict;
use	diagnostics;

#########################################
#########################################
sub	newSet
{
	my %set;
	foreach (@_) { $set{$_} = 1; }
	return	bless \%set, "Set";
}

#########################################
#########################################
sub	union
{
	my ($set1,$set2) = @_;
	my %set;
	foreach (keys %$set1) { $set{$_} = 1; }
	foreach (keys %$set2) { $set{$_} = 1; }
	return	bless \%set, "Set";
}

#########################################
#########################################
sub	intersect
{
	my ($set1,$set2) = @_;
	my %set;
	foreach (keys %$set2) { $set{$_} = 1 if $$set1{$_}; }
	return	bless \%set, "Set";
}

#########################################
#########################################
sub	setMinus
{
	my ($set1,$set2) = @_;
	my %set = %$set1;
	delete @set{keys %$set2};
	return	bless \%set, "Set";
}

#########################################
#########################################
sub	member
{
	my ($set,$x) = @_;
	return	$$set{$x};
}

#########################################
#########################################
sub	disjoint
{
	my ($set1,$set2) = @_;
	foreach (keys %$set2) { return 0 if $$set1{$_}; }
	return	1;
}

#########################################
#########################################
sub	subset
{
	my ($set1,$set2) = @_;
	foreach (keys %$set2) { return 0 unless $$set1{$_}; }
	return	1;
}

#########################################
#########################################
sub	superset
{
	my ($set1,$set2) = @_;
	foreach (keys %$set1) { return 0 unless $$set2{$_}; }
	return	1;
}

#########################################
#########################################
sub	setList
{
	my ($set) = @_;
	return	(keys %$set);
}

1;
