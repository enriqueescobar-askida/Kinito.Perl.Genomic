################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################

package	UnionFind;

use	Util;

#########################################
## Implements disjoint sets. 
## Uses inverted trees, union by rank, and path compression.
#########################################
use	strict;
use	diagnostics;

#########################################
## Creates a new union/find object.
#########################################
sub	new
{
	my ($this) = @_;
	return	bless [{ },{ }], (ref($this) || $this);
}

#########################################
#########################################
sub union
{
	my ($this,
	$key1, $key2) = @_;
	my ($parent,$rank) = @$this;
	$key1 = $this->find($key1);
	$key2 = $this->find($key2);
	return	$key1	if $key1 eq $key2; 
	($key1,$key2) = ($key2,$key1) if $$rank{$key1} > $$rank{$key2};
	$$rank{$key2} = max($$rank{$key2}, 1+$$rank{$key1});
	return	($$parent{$key1} = $key2);
}


#########################################
#########################################
sub find
{
	my ($this,
	$key) = @_;
	my ($parent,$rank) = @$this;
	my $parentKey = $$parent{$key};
	if (!defined($parentKey))
	{
		return	($key);
	}
	else
	{
		return	($$parent{$key} = $this->find($parentKey));
	}
}

#########################################
#########################################
sub inSameSet
{
	my ($this,$key1, $key2) = @_;
	$this->find($key1) eq $this->find($key2);
}

1;
