################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
package QualityReviser;
#########################################
##  This package provides routines for revising the quality values of reads
##  based on the qualities of corresponding bases in good alignments of pairs 
##  of pairs of reads.
#########################################
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(reviseQualities);

use strict;
use Util;

#########################################
sub reviseQualities
{
##  computes revised qualities for all given reads based on given
##  alignments.  The only exported routine of this package.
##  RETURNS: nothing; new qualities are saved in DnaRead objects.
#########################################
	my ($readList,          ## ref to list of DnaRead objects
		$alignmentsByPair,  ## ref to hash mapping pairs of reads to alignments
		$pass               ## integer indicated which pass; fewer alignments
		) = @_;             ## are considered on later passes
	
	## Create adjusted qualities for each read.
	foreach my $read (@$readList)
	{
		adjustQuality($read, 
				collectMatchedQualities($read,$pass));
	}
	## Compute Log-Likelihood Ratios for alignments.
	foreach (values %$alignmentsByPair)
	{
		foreach (@$_)
		{
			computeWeightedScore($_); 
		}
	}
}


#########################################
sub findBaseRange
{
##  finds a maximal substring of given sequence consisting of repetition
##  of given base and containg given position of sequence.
#########################################
	my ($seq,   ## sequence string
		$base,  ## base sought
		$index  ## position in index
		) = @_;

	## Break the string into two pieces at the index position,
	## and count the number of occurrences of $base on either side
	## of the break.
	my ($left, $bases) = ($seq =~ /^(.{$index})($base*)/);
	my $hi = $index + length($bases) - 1;
	my $n = length($seq) - $index;
	($bases, my $right) = ($seq =~ /($base*)(.{$n})$/);
	my $lo = $index - length($bases);
	return ($lo, $hi);
}

#########################################
sub collectMatchedQualities
{
##  examines all alignments to given read and collects, for each position
##  of the read, 1) best score of corresponding matching base in any alignment
##  and 2) best score of corresponding nonmatching base in any alignment.
##  RETURNS: list of two references to lists for (1) and (2).
#########################################
	my ($read,       ## the read to be revised.
		$pass        ## 0=any weighted score is OK; 
					## 1= weighted score must be positive.
					## 2=alignment must be within single contig.
		) = @_;

	my @bestAlternative;
	my @bestMatch = @{$read->OrigQual()};  ## COPIES!!
	
	foreach my $alignment (@{$read->Alignments()})
	{
	#	print "\nAdjust $read in pass $pass for "; $alignment->summarize();
		## Look at each alignment involving this read.
		## In second pass, ignore alignments with negative likelihood scores. 
		## In third pass, ignore alignments not used to form contigs.
		next if $pass == 1 && $alignment->WeightedScore() < 0;
		next if $pass == 2 && !($alignment->SupportedContig());
	
		## We want "read1" to refer to the read being adjusted.
		my $swap = ($read == $alignment->Read2());
		my $read1 = $read;  ## Just to give it a name with a number.
		my $read2 = $swap ? $alignment->Read1() : $alignment->Read2();
	
		my $path = $alignment->Path();
		$path =~ tr/ID/DI/ if $swap;
		my ($seq1,$seq2) = ($read1->Bases(), $read2->Bases());
	
		my $origQual2 = $read2->OrigQual();
		
		## Work through the local alignment column by column.
		## compute indices of first bases aligned.
		my ($siteIn1,$siteIn2) 
			= ($alignment->Start1()-1, $alignment->Start2()-1);
		($siteIn1,$siteIn2) = ($siteIn2,$siteIn1) if swap;
	
		foreach my $column (split(//,$path))
		{
			
			$siteIn1++ unless $column eq "I";
			$siteIn2++ unless $column eq "D";
			
			if ($column eq "M")
			{ ## match
				$bestMatch[$siteIn1] = max($bestMatch[$siteIn1],
							$$origQual2[$siteIn2]);
			}
			elsif ($column eq "S")
			{ ## substitution
				$bestAlternative[$siteIn1] 
				= max($bestAlternative[$siteIn1],
				$$origQual2[$siteIn2]);
			}
			else
			{  ## Insertion or Deletion
				my $base = ($column eq "D")? 
					substr($seq1,$siteIn1,1):substr($seq2,$siteIn2,1);
				my ($lo2,$hi2) = findBaseRange($seq2,$base,$siteIn2);
				next if ($hi2 < $lo2);
				my $q2 = min(@$origQual2[$lo2 .. $hi2]);
		#		my $q2 = ($hi2>$lo2) ? 
		#		    min(@$origQual2[$lo2 .. $hi2]) : $$origQual2[$lo2];
				my ($lo1,$hi1) = findBaseRange($seq1,$base,$siteIn1);
				
				foreach ($lo1 .. $hi1)
				{
					$bestAlternative[$_] = max($bestAlternative[$_], $q2);
				}
			}
		} ## foreach my $column
	} ## foreach my $alignment
	## Finished looking at each alignment for current sequence. 
	return (\@bestMatch, \@bestAlternative); 
}
	
#########################################
sub adjustQuality
{
## called once the quality information has been collected from all the
## alignments against given read's sequence.  Iterates through the
## positions of sequence and to derive new quality
## values based on this information.
## RETURNS: nothing; updates DnaRead object.
#########################################
	my ($read,            ## reference to DnaRead object being adjusted
		$bestMatch,       ## ref to list of best qualities of matching bases
		$bestAlternative  ## ref to list of best qualities of mismatch bases
		) = @_;
	
	my $len = $read->Length();
	my @adjQual = @{$read->AdjQual()};
	
	for (my $i=0; $i < $len; $i++)
	{
		my $q = $$bestMatch[$i] - $$bestAlternative[$i];
		$adjQual[$i] = min(9,max(1, $q));
	}
	
	$read->setAdjQual(\@adjQual);
}

#########################################
sub computeWeightedScore
##  computes a quality-weighted score of an alignment.
##  RETURNS: the weighted score.
#########################################
{
	my ($alignment) = @_;  ## The alignment whose score is computed.

	my ($read1,$read2) = ($alignment->Read1(), $alignment->Read2());
	my ($seq1,$seq2) = ($read1->Bases(), $read2->Bases());
	my ($adjQual1,$adjQual2) = ($read1->AdjQual(),$read2->AdjQual());
	
	## Work through the local alignment column by column.

	my ($siteIn1,$siteIn2) 
	= ($alignment->Start1()-1, $alignment->Start2()-1);
	my $path = $alignment->Path();
	my $score;

	foreach my $column (split('',$path))
	{
		$siteIn1++ unless $column eq "I"; ## gap in #1 = insertion
		$siteIn2++ unless $column eq "D"; ## gap in #2 = deletion
	
		my $minQ = min($$adjQual1[$siteIn1], $$adjQual2[$siteIn2]);
	
		if ($column eq "M")
		{ ## match
			$score += $minQ;
		}
		elsif ($column eq "S")
		{ ## mismatch
			$score += $minQ * JudiciousAligner->SUBSTPENALTY;
		}
		else
		{
			$score += $minQ * JudiciousAligner->INDELPENALTY;
		}
	} ## foreach my $column
	### Finished looking at each alignment for current sequence. 
	### Save result in object and return.
	$alignment->setWeightedScore($score);
	return $score;
}
