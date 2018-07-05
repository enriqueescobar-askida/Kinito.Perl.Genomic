#!/usr/bin/perl -w -I . -I ../lib -I ../lib2
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################

#########################################
##  This program illustrates two greedy strategies for multiple alignment
##  based on computing similarities of alignments and merging alignments.
#########################################

use strict;
use Util;

my $gapPenalty = -2;

#########################################
sub planAlignmentMerge
{
##  fills two dynamic programming tables, @M and @how, which then tell
##  the score of the best possible merge of alignments $S and $T and how
##  to add gaps to $S and $T to form the best merge.
##  RETURNS: best score, and reference @how array.
#########################################
    my($S,  ## an alignment = reference to list of sequences with gaps
    $T)  ## another alignment
    = @_;
    my ($sLen, $tLen) = (length($$S[0]), length($$T[0]));
    my @Sgaps = (("-")x@$S);
    my @Tgaps = (("-")x@$T);
    my @M;
    my @how;
    $M[0][0] = 0;
    foreach my $i (1..$sLen)
    {
        $how[$i][0]="|";
        my @Scol = map { substr($_,$i-1,1) } @$S;
        $M[$i][0] = $M[$i-1][0] + scoreColumn(@Scol, @Tgaps);
    }
    foreach my $j (1..$tLen)
    {
        $how[0][$j]="-";
        my @Tcol = map { substr($_,$j-1,1) } @$T;
        $M[0][$j] = $M[0][$j-1] + scoreColumn(@Sgaps, @Tcol);
    }
    foreach my $i (1..$sLen)
    {
        foreach my $j (1..$tLen)
        {
            my @Scol = map { substr($_,$i-1,1) } @$S;
            my @Tcol = map { substr($_,$j-1,1) } @$T;
            $M[$i][$j] = $M[$i-1][$j-1] + scoreColumn(@Scol,@Tcol);
            $how[$i][$j] = "\\";
            my $left = $M[$i][$j-1] + scoreColumn(@Sgaps, @Tcol);
            if ($left>$M[$i][$j])
            {
                $M[$i][$j] = $left;
                $how[$i][$j] = "-";
            }
            my $up = $M[$i-1][$j] + scoreColumn(@Scol, @Tgaps);
            if ($up>$M[$i][$j])
            {
                $M[$i][$j] = $up;
                $how[$i][$j] = "|";
            }
        }
    }
    return ( $M[$sLen][$tLen], \@how );
}

#########################################
sub prependColumnToAlignment
{
##  adds a new column to the left end of an existing alignment.
##  RETURNS: nothing.
#########################################
    my ($A,   ## alignment (reference to list of strings)
    @col) ## list of symbols to be added to each string of alignment
    = @_;
    foreach (@$A)
    {
        $_ = (shift @col).$_
    };
}

#########################################
sub mergeAlignments
{
##  merges two alignment to form the highest-scoring merged alignment possible.
##  RETURNS: resulting alignment (reference to list of strings).
#########################################
    my ($S,$T) = @_;
    my ($score, $how) = planAlignmentMerge($S,$T);
    my @result = (("") x (@$S+@$T));
    my ($i,$j) = (length($$S[0]), length($$T[0]));
    my @Sgaps = (("-")x@$S);
    my @Tgaps = (("-")x@$T);
    while ($i>0 || $j>0)
    {
        if ($$how[$i][$j] eq "\\")
        {
            my @Scol = map { substr($_,$i-1,1) } @$S;
            my @Tcol = map { substr($_,$j-1,1) } @$T;
            prependColumnToAlignment(\@result, @Scol, @Tcol);
            $i--; $j--;
        }
        elsif ($$how[$i][$j] eq "|")
        {
            my @Scol = map { substr($_,$i-1,1) } @$S;
            prependColumnToAlignment(\@result, @Scol, @Tgaps);
            $i--;
        }
        elsif ($$how[$i][$j] eq "-")
        {
            my @Tcol = map { substr($_,$j-1,1) } @$T;
            prependColumnToAlignment(\@result, @Sgaps, @Tcol);
            $j--;
        }
    }
    return \@result;
}

#########################################
sub scoreColumn
{
##  given a list of symbols in a column,
##  RETURNS: the sum-of-pairs score for the column.
#########################################
    my @col = @_;  ## the symbols in the column.
    my ($gaps,$aas,$score) = (0,0,0);
    while (@col)
    {
        my $aa = shift @col;
        ($gaps++, next) if $aa eq "-";
        $aas++;
        foreach my $aa1 (@col)
        {
            next if $aa1 eq "-";
    #       $score += $blosum62{$aa,$aa1};
            $score += ($aa eq $aa1) ? +1 : -1;
        }
    }
    return $score + ($gapPenalty * $gaps * $aas);
}

#########################################
sub scoreMultipleAlignment
{ 
##  sums up the sum-of-pairs scores of all columns of an alignment.
##  RETURNS: the sum-of-pairs score of the entire alignment.
#########################################
    my ($alignment) = @_;
    my $score;
    foreach my $i (0..length($$alignment[0])-1)
    {
        $score += scoreColumn(map {substr($_,$i,1)} @$alignment);
    }
    return $score;
}
    
#########################################
sub greedyStrategyA
{  
##  given a list of sequences, applies a greedy merging strategy to derive
##  a good, though not necessarily optimal, alignment.  The strategy consists
##  of repeatedly maximizing immediate gain by merging any two existing 
##  alignments.
##  RETURNS: the final multiple alignment.
#########################################
    my @aligns = map { [$_] } @_;  ## Change sequences to singleton alignments.
    
    while (@aligns>1) {
    ### Find two best ones to merge.
    my @alignScores = map { scoreMultipleAlignment($_) } @aligns;
    my ($bestI, $bestJ, $bestDelta) = (-1,-1,-999999999);
    foreach my $i (0..$#aligns)
    {
        foreach my $j ($i+1..$#aligns)
        {
            my ($score,$how) = planAlignmentMerge($aligns[$i],$aligns[$j]);
            my $delta = $score - $alignScores[$i] - $alignScores[$j];
    #       print "$i,$j: $delta\n";
            ($bestI,$bestJ,$bestDelta) = ($i,$j,$delta) 
                if $delta>$bestDelta;
        }
    }
    ### Replace them with their merge.
    $aligns[$bestI] = mergeAlignments($aligns[$bestI], $aligns[$bestJ]);
    splice(@aligns, $bestJ, 1);
    ### Print trace.
#   print "***** ALIGNMENT LIST:\n";
#   for (my $i=0; $i<@aligns; $i++) {
#       my $score = scoreMultipleAlignment($aligns[$i]);
#       print "$i:  score $score\n";
#       foreach (@{$aligns[$i]}) { print "$_\n"; }
#   }
    }

    return $aligns[0];
}

#########################################
sub greedyStrategyB
{
##  given a list of sequences, applies a greedy merging strategy to derive
##  a good, though not necessarily optimal, alignment.  The strategy consists
##  of repeatedly maximizing immediate gain by merging any single sequence
##  into the existing alignments.
##  RETURNS: the final multiple alignment.
#########################################
    my @ss = @_;     ## Argument list is a list of sequences.
    ### Find best pair to merge first.
    my ($bestI, $bestJ, $bestScore) = (-1,-1,-999999999);
    foreach my $i (0..$#ss)
    {
        foreach my $j ($i+1..$#ss)
        {
            my ($score,$how) = planAlignmentMerge([$ss[$i]],[$ss[$j]]);
            ($bestI,$bestJ,$bestScore) = ($i,$j,$score)
            if $score>$bestScore;
        }
    }
    ### Merge them.
    my $alignment = mergeAlignments([$ss[$bestI]],[$ss[$bestJ]]);
    splice(@ss, $bestJ, 1);
    splice(@ss, $bestI, 1);
    
    while (@ss)
    {
        ### Find best sequence to add to alignment.
        my ($bestI, $bestScore) = (-1,-999999999);
        foreach my $i (0..$#ss)
        {
            my ($score,$how) = planAlignmentMerge([$ss[$i]], $alignment);
            ($bestI,$bestScore) = ($i,$score) if $score>$bestScore;
        }
        ### Add best sequence to alignment; remove from @ss.
        $alignment = mergeAlignments([$ss[$bestI]],$alignment);
        splice(@ss, $bestI, 1);
    }
    return $alignment;
}

#########################################
sub printAlignment
{
##  prints a multiple alignment with a header.
##  RETURNS: nothing.
#########################################
    my ($title,      ## heading for alignment
    $alignment)  ## alignment : reference to list of strings.
    = @_;
    print "\n***** $title\n";
    foreach (@$alignment)
    {
        print "$_\n";
    }
    my $score = scoreMultipleAlignment($alignment);
    print "score $score\n\n";
}    


#########################################
##  MAIN PROGRAM
#########################################
my @ss;
while (my $s = <STDIN>)
{
    chomp($s);
    push @ss, $s;
}
printAlignment("Greedy Strategy A", greedyStrategyA(@ss));
printAlignment("Greedy Strategy B", greedyStrategyB(@ss));
