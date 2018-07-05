################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
package ContigBuilder;
#########################################
##  This package assigns DNA reads to contigs, based on alignments
##  between pairs of reads collected by JudiciousAligner and qualities
##  computed by QualityReviser.
#########################################
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(assignReadsToContigs);

use strict; 
use Util;

my $CutoffFactor = -60;  # -20

#########################################
sub assignReadsToContigs
{
##  assigns DNA reads to contigs using a greedy algorithm.
##  The only routine exported by this package.
##  RETURNS:  A list of non-empty contigs.
#########################################
	my ($readList, $alignmentsByPair) = @_;
	my @allContigs;
	
	foreach my $read (@$readList)
	{
		my $contig = Contig->new([$read], 1, 0, $read->Length()-1);
		$read->setContig($contig);
		$read->setContigStart(0);
		$read->setContigEnd($read->Length()-1);
		push @allContigs, $contig;
	}

	## make read layout, using greedy algorithm
	my $sortedAlignments = findBestAlignments($alignmentsByPair);

	## Process each pair from highest to lowest weighted score.
	foreach my $alignment (@$sortedAlignments)
	{
		mergeByAlignment($alignment);
	} ## foreach my $alignment
	return [grep { @{$_->Reads() } > 0 } @allContigs];
}

#########################################
sub findBestAlignments
{
##  eliminates overlapping alignments from its argument, keeping those with
##  higher scores.  Produces a score-sorted list of all remaining alignments
##  for use by greedy contig-building algorithm.
##  RETURNS:  ref to score-sorted list.
#########################################
	my ($alignmentsByPair) = @_;  ## ref to hash mapping pairs of reads to alignments
	my @bestAlignments;
	
	foreach my $pair (keys %$alignmentsByPair)
	{
		my ($r1,$r2) = split($;, $pair);
		next if $r1 > $r2;  ## Don't look at the same list twice.
		my @alignments = 
						sort {$a->WeightedScore() <=> $b->WeightedScore() || 
							$a->RawScore() <=> $b->RawScore() 
							}  @{$$alignmentsByPair{$pair}};
		ILOOP:
		foreach my $i (0 .. $#alignments)
		{
			my $cand = $alignments[$i];
			my ($cStart1,$cEnd1) = ($cand->Start1(), $cand->End1());
			my ($cStart2,$cEnd2) = ($cand->Start2(), $cand->End2());
			
			foreach my $j ($i+1 .. $#alignments)
			{
				my $other = $alignments[$j];
				my ($oStart1,$oEnd1)=($other->Start1(),$other->End1());
				my ($oStart2,$oEnd2)=($other->Start2(),$other->End2());
				## Dismiss candidate if it is overlapped by higher-scoring
				## alignment in either read.
				next ILOOP unless ($cStart1 > $oEnd1) || ($cEnd1 < $oStart1);
				next ILOOP unless ($cStart2 > $oEnd2) || ($cEnd2 < $oStart2);
			}
			push @bestAlignments, $cand;
		}
	}	
	## Sort out the survivors by score for greedy algorithm.
	@bestAlignments = sort {$b->WeightedScore() <=> $a->WeightedScore() || 
							$b->RawScore() <=> $a->RawScore()
							} @bestAlignments;
	return \@bestAlignments;
}



#########################################
sub mergeByAlignment
{
##  tries to merge the contigs containing the reads aligned by the given
##  alignment.
##  RETURNS:  nothing.
#########################################
	my ($alignment) = @_;  ## the alignment
	my ($read1, $read2) = ($alignment->Read1(), $alignment->Read2());
	return if $read1 == $read2;
	
	## find the contigs whereof the two reads are currently members
	my ($contig1, $contig2) = ($read1->Contig(), $read2->Contig());
	
	## find the offset between the two contigs implied by the alignment.
	my $offset = findAlignmentOffset($alignment);

	if ($contig1 == $contig2)
	{ ## if already in same contig
		## Does this alignment just bend the contig, or break it?
		if (abs($offset) < 25)
		{ ## It just bends it. Use match
			## Further support for this contig.
			$alignment->setSupportedContig($contig1);
			push @{$contig1->AlignedPairs()}, $alignment;
		}
	}
	else
	{ ## different contigs; we will try to merge the contigs.
		if (testContigMerge($alignment,$offset))
		{
			mergeContigs($contig1, $contig2, $offset);
			$alignment->setSupportedContig($contig1);
			push @{$contig1->AlignedPairs()}, $alignment;
		}
	} ## else different contigs 
}


#########################################
sub findAlignmentOffset
{
##  finds appropriate offset for two contigs implied by an 
##  alignment between two reads, one from each contig.
##  RETURNS: the offset.
#########################################
	my ($alignment) = @_;   ## the alignment
	my ($read1,$read2) = ($alignment->Read1(),$alignment->Read2());
	my ($start1,$start2) = ($read1->ContigStart(),$read2->ContigStart());
	my ($end1,$end2) = ($read1->ContigEnd(), $read2->ContigEnd());
	my $offset = $alignment->MeanOffset();
	my ($length1,$length2) = ($read1->Length(), $read2->Length());
	return ($start1 - $start2 - $offset);
}

#########################################
sub mergeContigs
{
##  Merge contig1 and contig2 with offset given.
##  All reads are assigned to contig1 and contig2 is emptied.
##  RETURNS: nothing.
#########################################
	my ($contig1, $contig2,   ## the two contigs to merge.
	$offset) = @_;        ## the offset between their origins.

	## combined linked list of reads.  
	foreach my $read (@{$contig2->Reads()})
	{
		$read->setContig($contig1);
		$read->setContigStart($offset + $read->ContigStart());
		$read->setContigEnd($offset + $read->ContigEnd());
	}

	$contig1->setReads([@{$contig1->Reads()}, @{$contig2->Reads()}]);
	$contig1->setAlignedPairs([@{$contig1->AlignedPairs()},
				@{$contig2->AlignedPairs()}]);

	## Adjust indices of ends of $contig1 w.r.t. contig origin.
	$contig1->setFirstStart(min($contig1->FirstStart(),
				$contig2->FirstStart() + $offset));
	$contig1->setLastEnd(max($contig1->LastEnd(), 
				$contig2->LastEnd() + $offset));
	
	$contig2->setReads([]);
	$contig2->setFirstStart(undef);
	$contig2->setLastEnd(undef);
	
}

#########################################
sub testContigMerge
{
##  
##  RETURNS: 1 if merge is OK; 0 if not.
#########################################
	my ($alignment, ## Alignment upon which proposed merge is to be based.
		$offset)    ## Offset between contigs implied by alignment.
		= @_;
	
	my ($aRead1, $aRead2) = ($alignment->Read1(), $alignment->Read2());
	## find the contigs whereof the two reads are currently members
	my ($contig1, $contig2) = ($aRead1->Contig(), $aRead2->Contig());

	my $alignmentScore = $alignment->WeightedScore();
	my $rejectCutoff = $alignmentScore - 0.20 * abs($alignmentScore);
	
	## for each read in contig2
	foreach my $c2Read (@{$contig2->Reads()})
	{
		## for each alignment to the read 
		foreach my $align2 (@{$c2Read->Alignments()})
		{
			## $c1Read is a read in contig 1.
			my $c1Read = $align2->Read2();
			my $swap = ($c1Read == $c2Read);
			$c1Read = $align2->Read1() if $swap;
			## We need alignments between pairs in the contigs being merged.
			next if $c1Read->Contig() != $contig1;
			
			## What offset between contigs is implied by this alignment?
			## How different is this from previous estimates?
				## If a lot, including the both $alignment and $align2 in the 
			## contig will put a "strain" on the contig.
			my $alignOffset = ($swap?-1:1)*findAlignmentOffset($align2);
			my $strain = abs($offset + $alignOffset);
			my $align2score = $align2->WeightedScore();
	
			## If the strain is large, and $align2 is not grossly less
			## reliable than $alignment, we cannot safely merge.
			return 0 if ($strain > 25) && ($align2score >= $rejectCutoff);
		} ## for my $align2
	} ## for my $read2
	return 1;  ## merge is OK.
}

1;
