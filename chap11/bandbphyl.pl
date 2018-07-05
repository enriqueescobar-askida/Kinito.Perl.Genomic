#!/usr/bin/perl -I . -I ../perllib
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
use strict; 
use Util;
use Phylogeny;

my $count;

#########################################
## MAIN PROGRAM -- isolated in a subroutine
#########################################
main();

#########################################
sub main
{
##  Driver to build approximate and optimal phylogenies
##  using various strategies.
#########################################
    $|=1;
    my @S = <STDIN>;    ## reads all lines at once.
    foreach (@S) { chomp($_); }
    my $greedy = build_tree_greedily(\@S);
    $greedy->printTable();
    my $best = build_trees_exhaustively(\@S);
    $best->printTable();
    $best = build_trees_efficiently(\@S);
    $best->printTable();
    $best = build_trees_efficiently(\@S, $greedy->countMutations());
    $best->printTable();
}


#########################################
sub build_trees_exhaustively
{
##  Explicitly evaluates every possible tree topology to find the
##  one with the fewest mutations.
##  RETURNS: reference to phylogeny object for optimal phylogeny.
#########################################
    my ($sequences) = @_;   ## reference to list of protein sequences.
    my $phyl = Phylogeny->new($sequences);
    $phyl->joinNewNeighbors(0,1);
    my $bound = 1E99;
    $count=0;
    my $best;
    build_all_trees($phyl,2, $phyl->countMutations(),\$bound,\$best);
    warn "$count calls to build_all_trees.\n";
    return $best;
}

#########################################
sub build_all_trees
{
##  Recursive subroutine to assist build_trees_exhaustively.
##  Each invocation is responsible for trying all possible insertion
##  points for one leaf, and makes recursive calls to add other leaves.
##  RETURNS: nothing.
#########################################
    my ($phyl,       ## reference to phylogeny being built
    $newLeaf,    ## index of leaf that this invocation places
    $mutations,  ## number of mutations in phylogeny so far
    $bound,      ## number of mutations in best complete phylogeny so far
    $best) = @_; ## best complete phylogeny seen so far
    $count++;
    if (!($phyl->isLeaf($newLeaf)))
    {
        if ($mutations<=$$bound)
        {
            $$bound = $mutations;
            warn "Call no. $count: Found tree with $mutations mutations.\n";
            $$best = $phyl->copy();
        }
        return;
    }
    foreach ($phyl->allEdges())
    {
        (my $v1, my $v2) = @$_;
        my $changeList = [];
        my $deltaMut = $phyl->attachLeafBetween($newLeaf,$v1,$v2,$changeList);
        build_all_trees($phyl,$newLeaf+1,$mutations+$deltaMut,$bound,$best);
        $phyl->detachLeafBetween($newLeaf, $v1, $v2, $changeList);
    }
}
    
#########################################
sub build_tree_greedily
{
##  Uses a greedy strategy to quickly find a tree topology with
##  a small, but not necessarily minimal, number of mutations.
##  Gives useful bound quickly for branch-and-bound algorithm.
##  RETURNS: reference to phylogeny object for good phylogeny.
#########################################
    my ($sequences) = @_;   ## reference to list of protein sequences.
    my $phyl = Phylogeny->new($sequences);
    $phyl->joinNewNeighbors(0,1);
    my $mutations = $phyl->countMutations();
    for (my $newLeaf=2; $phyl->isLeaf($newLeaf); $newLeaf++)
    {
        my ($bestDelta, $bestV1, $bestV2) = (1E99, -1, -1);
        foreach ($phyl->allEdges())
        {
            (my $v1, my $v2) = @$_;

            my $changeList = [];
            my $deltaMutations
            = $phyl->attachLeafBetween($newLeaf, $v1, $v2, $changeList);
            ($bestDelta, $bestV1, $bestV2) = ($deltaMutations,$v1,$v2)
            if $deltaMutations < $bestDelta;
            $phyl->detachLeafBetween($newLeaf, $v1, $v2, $changeList);
        }
        $phyl->attachLeafBetween($newLeaf, $bestV1, $bestV2, []);
        $mutations += $bestDelta;
    }
    return $phyl;
}

#########################################
sub build_trees_efficiently
{
##  Uses branch-and-bound strategy to avoid explicitly evaluating
##  every possible tree topology to find the one with the fewest mutations.
##  RETURNS: reference to phylogeny object for optimal phylogeny.
#########################################
    my ($sequences,   ## reference to list of protein sequences.
    $upperBound)  ## number of mutations in some phylogeny known 
                    ##     a priori; optional.
    = @_;
    my $phyl = Phylogeny->new($sequences);
    $phyl->joinNewNeighbors(0,1);
    $upperBound ||= 1E99;  ## assume the absurd if no bound given.
    $count=0;
    my $best;
    build_good_trees($phyl,2, $phyl->countMutations(), \$upperBound, \$best);
    warn "$count calls to build_good_trees.\n";
    return $best;
}

#########################################
sub build_good_trees
{
##  Recursive subroutine to assist build_trees_efficiently.
##  Each invocation is responsible for trying all possible insertion
##  points for one leaf, and makes recursive calls to add other leaves.
##  Search is "pruned" whenever number of mutations exceeds number in
##  best complete phylogeny seen so far.
##  RETURNS: nothing.
#########################################
    my ($phyl,       ## reference to phylogeny being built
    $newLeaf,    ## index of leaf that this invocation places
    $mutations,  ## number of mutations in phylogeny so far
    $bound,      ## number of mutations in best complete phylogeny so far
    $best) = @_; ## best complete phylogeny seen so far
    $count++;
    if ($mutations > $$bound)
    {
    #   warn "Pruning at Leaf $newLeaf...\n";
        return;
    }
    if (!($phyl->isLeaf($newLeaf)))
    {
        warn "Call no. $count: Found tree with $mutations mutations.\n";
        $$bound = $mutations;
        $$best = $phyl->copy();
        return;
    }
    foreach ($phyl->allEdges())
    {
        (my $v1, my $v2) = @$_;
        my $changeList = [];
        my $deltaMut = $phyl->attachLeafBetween($newLeaf,$v1,$v2,$changeList);
        build_good_trees($phyl,$newLeaf+1,$mutations+$deltaMut,$bound,$best);
        $phyl->detachLeafBetween($newLeaf, $v1, $v2, $changeList);
    }
}



