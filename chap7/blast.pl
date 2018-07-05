#!/usr/bin/perl -w -I . -I ../lib -I ../lib2
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################

#########################################
##  A BLAST-like program for searching a database for sequences homologous
##  to a given query sequence.
#########################################

use strict;
use Util;
use SeqReader;

my %lod;  ## matrix of lod-scores for residue pairs is store here.

#########################################
sub fillLod
{
##  Reads a scoring matrix from the specified file; stores it in %lod.
##  RETURNS: nothing.
#########################################
    my ($matrixFile) = @_;  ## name of file from which to read matrix.
    open LOD, $matrixFile;
    my ($trash,@residues) = split /\s+/, <LOD>;
    while (<LOD>)
    {
        my ($r,@scores) = split;
        foreach (0..19)
        {
            $lod{$r.$residues[$_]} = $scores[$_];
        }
    }
}

my @residues = split //,"acdefghiklmnpqrstvwy";   ## list of 20 amino acids

#########################################
sub findSimilarWords
{
## Given three-letter $word, finds all others with lod score >= 11.
## RETURNS: reference to list of words found.
#########################################
    my ($word) = @_;
    my ($w1,$w2,$w3) = split //,$word;
#   print "Similar to $word:\n";
    return [] if $lod{$w1.$w1}+$lod{$w2.$w2}+$lod{$w3.$w3}<11;
    my @similar;
    foreach my $r1 (@residues)
    {
        foreach my $r2 (@residues)
        {
            my $t = 11-$lod{$w1.$r1}-$lod{$w2.$r2};
            foreach my $r3 (@residues)
            {
                push @similar, "$r1$r2$r3" if $lod{$w3.$r3}>=$t;
            }
        }
    }
#   foreach (@similar) { print " $_"; }
#   print "\n";
    return \@similar;
}
    
my %similarWords;

#########################################
sub preprocessQuery
{
## Given query sequence, builds index (hash) mapping three-letter words 
## into list of positions in query having similar (score>11) words.
## RETURNS: reference to hash containing index.
#########################################
    my ($query) = @_;  ## query sequence to be indexed.
    my %similarPositions;
    for (my $i=0; $i<(length $query)-2; $i++)
    {
        my $word = substr($query,$i,3);
        $similarWords{$word} ||= findSimilarWords($word);
        foreach (@{$similarWords{$word}})
        {
            $similarPositions{$_} ||= [ ];
            push @{$similarPositions{$_}}, $i+1;
        }
    }
#   print scalar keys %similarPositions, " keys:\n";
#   foreach (sort keys %similarPositions) {
#       print "$_: ", join(",", @{$similarPositions{$_}}), "\n";
#   }
    return \%similarPositions;
}


#########################################
sub extendHit
{
##  extends a three-letter "hit" between query and target to a maximal-scoring
##  gapless alignment.
##  RETURNS: reference to a three-element list:
##           [ alignment score, aligned part of query, aligned part of target ]
#########################################
    my ($query,  ## query sequence (string)
    $qPos,   ## position of hit in query sequence
    $qId,    ## identifier (database label) of query
    $target, ## target sequence (string)
    $tPos,   ## position of hit in target sequence
    $tId)    ## identifier (database label) of query
    = @_;
    my @target = ("-", split //,$target);
    my @query = ("-", split //,$query);
    my ($lo,$hi) = (0,2);

    my $maxscore = my $score =
    $lod{$target[$tPos].$query[$qPos]}
        + $lod{$target[$tPos+1].$query[$qPos+1]}
        + $lod{$target[$tPos+2].$query[$qPos+2]};

    ### Try to grow gap-free alignment to right.
    my $ilim = min(@target - $tPos, @query - $qPos);
    for (my $i = 3; $i<$ilim; $i++)
    {
        $score += $lod{$target[$tPos+$i].$query[$qPos+$i]};
        last if $score < 0;
        ($maxscore,$hi) = ($score,$i) if $score > $maxscore;
    }

    ### Try to grow gap-free alignment to left.
    my $score = $maxscore;
    my $ilim = min($tPos, $qPos);
    for (my $i = 1; $i<$ilim; $i++)
    {
        $score += $lod{$target[$tPos-$i].$query[$qPos-$i]};
        last if $score < 0;
        ($maxscore,$lo) = ($score,$i) if $score > $maxscore;
    }

    ### Return best alignment as triple.
    my $len = $hi+$lo+1;
    return 
    [$maxscore, 
    substr($query, $qPos-$lo-1, $len)." at ".($qPos-$lo-1)." in $qId",
    substr($target, $tPos-$lo-1, $len)." at ".($tPos-$lo-1)." in $tId"];
}   
    

#########################################
sub scanTarget
{
##  scans target to find "hits" with query sequence.  For each, calls
##  extendHit to get an alignment.
##  Trims list of alignments to 25 best.
##  RETURNS: reference to list of best alignments.
#########################################
    my ($target,    ## target sequence (string)
    $tId,       ## target identifier
    $query,     ## query sequence (string)
    $qId,       ## query identifier
    $queryIndex)## index of query produced by preprocessQuery
    = @_;
    my @alignments;
    for (my $i=0; $i<(length $target)-2; $i++)
    {
        my $word = substr($target,$i,3);
        my $tPos = $i+1;
        foreach my $qPos (@{$$queryIndex{$word}})
        {
    #       print "Hit $word at $tPos in target and $qPos in query.\n";
            my $alignment = extendHit($query,$qPos,$qId,$target,$tPos,$tId);
    #       print(join("\n", @$alignment), "\n\n");
            push @alignments, $alignment;
        }
    }
    @alignments = sort { $$b[0] <=> $$a[0] } @alignments;
    $#alignments = min( $#alignments, 24);

#   foreach (@alignments) { print "$$_[0] "; }
#   print "\n";
    return @alignments;
}

#########################################
### MAIN PROGRAM
#########################################
## Read in scoring matrix specified on command line -- default blosum62.
fillLod($ARGV[2] || "blosum62");
my $queryFile = new SeqReader $ARGV[0];
while (my $query=$queryFile->readSeq())
{  ## for each query in query file.
    my $qId = $queryFile->seqId();
    my $queryIndex = preprocessQuery($query);  ## build index
    my $targetFile = new SeqReader $ARGV[1];
    my @bestAlignments;
    while (my $target = $targetFile->readSeq())
    {  ## for each seq in database
        my $tId = $targetFile->seqId();
        print "\nQuery: $qId\nTarget: $tId\n";

        my @newAlignments = 
            scanTarget($target,$tId,$query,$qId,$queryIndex);  ## scan seq
        ### Add new alignments to sorted list.
        push @bestAlignments, @newAlignments;
        @bestAlignments = 
            sort { $$b[0] <=> $$a[0] } @bestAlignments;
        ### Cut off list at 25 alignments.
        $#bestAlignments = min($#bestAlignments,24);
    }
    $targetFile->close();
    print "********** 25 Best Alignments for $qId:\n\n";
    foreach (@bestAlignments)
    {
        print(join("\n", @$_), "\n\n");
    }
}
$queryFile->close();

