#!/usr/bin/perl -I. -I../perllib
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################

#########################################
## This program implements a branch-and-bound algorithm
## for finding shortest common superstrings.
#########################################

use strict;
use SuffixTree;

#########################################
##  GLOBAL VARIABLES
#########################################

my %tree;   ## hash mapping input strings to their suffix trees.
my %graph;  ## hash encoding overlap graph.

#########################################
##  MAIN PROGRAM
#########################################

## Read in strings, build suffix trees, create weighted graph.
while (my $in=<STDIN>)
{
    chomp($in);
    $graph{$in} = [];
    my $inTree = SuffixTree->new($in);

    foreach (keys %tree)
    {
        addWeightedEdge($in, $_, $inTree->maxSuffix($_));
        addWeightedEdge($_, $in, $tree{$_}->maxSuffix($in));
    }
    $tree{$in} = $inTree;
}

## Throw away suffix trees.
%tree = undef;
## Sort adjacency lists to put big weights first.
foreach (keys %graph)
{
    my @temp = sort {$$b[1] <=> $$a[1]} @{$graph{$_}};
    $graph{$_} = \@temp;
}

my $bestWeight = 0;
my $bestString;
my %visited;

foreach (keys %graph)
{
    search($_, 0, $_);
}
print "\n\nShortest Common Superstring is:\n$bestString\n";

#########################################
sub addWeightedEdge
{
##  adds a weighted edge to the graph.
##  RETURNS: nothing.
#########################################
    my ($v1, $v2, $weight) = @_;   ## two vertices, weight
    push @{$graph{$v1}}, [$v2, $weight];
}


#########################################
sub prospects
{
##  Overestimates the amount of overlap that can be achieved in
##  completing the current partial superstring.
##  RETURNS: the estimate of overlap.
#########################################
    my ($newV,  ## next vertex(fragment) to add
		$viz    ## ref to list of vertices(fragments) already visited(included)
		) = @_;
    my @unviz = grep { !$$viz{$_} } (keys %graph);
    my %maxIn;

    foreach my $v (@unviz)
    {
        foreach my $edge (@{$graph{$v}})
        {
            my ($u, $edgeWeight) = @$edge;
            last if $edgeWeight == 0;
            $maxIn{$u} = $edgeWeight if $edgeWeight>$maxIn{$u};
        }
    }
    my $result = 0;
    
    foreach my $v (@unviz)
    {
        $result += $maxIn{$v};
    }
    return ($result - $maxIn{$newV});
}
        
#########################################
sub search
{
##  recursively searches for the shortest common supersequence by adding
##  adding one fragment at a time.  Each invocation is responsible for
##  selecting the next fragment.
##  RETURNS: nothing.
#########################################
    my ($v,          ##  vertex/fragment to be added.
		$pathWeight, ##  total overlap (edge weight) achieved so far.
		$pathString  ##  current partial superstring.
		) = @_;
    return if $visited{$v};
    return if ($pathWeight + prospects($v,\%visited)) <= $bestWeight;
    
    if ($pathWeight>$bestWeight)
    { 
        ($bestWeight,$bestString) = ($pathWeight,$pathString);
        print "New Best ($pathWeight): $pathString\n";
    }
    $visited{$v} = 1;
    
    foreach my $edge (@{$graph{$v}})
    {
        my ($u, $edgeWeight) = @$edge;
        search($u,$pathWeight+$edgeWeight,$pathString.substr($u,$edgeWeight)); 
    }
    $visited{$v} = 0;
}    
