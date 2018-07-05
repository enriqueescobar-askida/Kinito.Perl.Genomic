################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
#########################################
package UnionFind;
#########################################
## Implements disjoint sets. 
## Uses inverted trees, union by rank, and path compression.
#########################################

use strict;
use Util;

#########################################
sub new
##  creates a new union/find object.
##  RETURNS: blessed reference to new object.
#########################################
{
    my ($this) = @_;   ## usually, the literal "UnionFind"
    return bless [{ },{ }], (ref($this) || $this);
}

#########################################
sub union
##  merges the sets containing two given keys.
##  RETURNS: representative key of merged set.
#########################################
{
    my ($this,         ## current UnionFind object
    $key1, $key2   ## the two keys
    ) = @_;
    $key1 = $this->find($key1);
    $key2 = $this->find($key2);
    return $key1 if $key1 eq $key2; 
    ($key1,$key2) = ($key2,$key1) if $$rank{$key1} > $$rank{$key2};
    $$rank{$key2} = max($$rank{$key2}, 1+$$rank{$key1});
    return ($$parent{$key1} = $key2);
}


#########################################
sub find
##  RETURNS: representative key of set containing given key.
#########################################
{
    my ($this,
    $key) = @_;
    my ($parent,$rank) = @$this;
    my $parentKey = $$parent{$key};
    if (!defined($parentKey))
    {
        return ($key);
    }
    else
    {
        return ($$parent{$key} = $this->find($parentKey));
    }
}

#########################################
sub inSameSet
##  RETURNS: true if and only if two given keys are in the same set.
#########################################
{
    my ($this,         ## current UnionFind object
    $key1, $key2   ## the two keys
    ) = @_;
    $this->find($key1) eq $this->find($key2);
}

1;
