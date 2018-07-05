################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
package ConsensusBuilder;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(constructConsensusSequences);

use strict; 
use UnionFind;
use Util;

#########################################
sub constructConsensusSequences
{
##  constructs a consensus sequence for each of a list of contigs.
#########################################
	my ($contigs) = @_;  ## reference to list of Contigs
	
	foreach (@$contigs)
	{
		if (@{$_->Reads()} > 1)
		{
			constructConsensusFromGraph($_);
		}
		else
		{  ## singleton contig
			my ($read) = @{$_->Reads()};
			$_->setLength($read->Length());
			$_->setBases($read->Bases());
			$_->setQuals(join("", @{$read->AdjQual()}));
			$_->setScore(sum(@{$read->AdjQual()}));
		}
	}
}

#########################################
sub findEquivalentReadSites
{
##  gathers sets of positions in individual reads that correspond to
##  the same position on the contig consensus.
##  RETURNS: references to three hashes encoding information found.
#########################################
	my ($contig) = @_;  ## ref to the Contig to be processed.
	my %readSitePairs;
	my $uf = UnionFind->new(); ## maps readSites to (unplaced) contigSites.
	
	foreach my $alignment (@{$contig->AlignedPairs()})
	{
		my ($site1,$site2)=($alignment->Start1(),$alignment->Start2());
		my $serial1 = $alignment->Read1()->SerialNumber();
		my $serial2 = $alignment->Read2()->SerialNumber();
		my $path = $alignment->Path();
		
		while ($path)
		{
			my $column = substr($path,0,1);
			$site1++ unless $column eq "I";
			$site2++ unless $column eq "D";
			
			if ($path =~ /^M{8}/)
			{
				while ($path =~ s/^M{4}//)
				{
					my $rsPair1 = $serial1 . $; . $site1;
					my $rsPair2 = $serial2 . $; . $site2;
					$readSitePairs{$rsPair1}++;
					$readSitePairs{$rsPair2}++;
					$uf			->union($rsPair1, $rsPair2);
					$site1 += 4; $site2 += 4;
				}
			}
			else
			{
				$path =~ s/^.//;
			}
		}
	}

	## Build one hash whose values are equivalence classes, and
	## another containing lists of sites on the same read.
	my (%equivClasses, %sitesOnRead);
	
	foreach my $rsPair (keys %readSitePairs)
	{
		my $root = $uf->find($rsPair);
		$equivClasses{$root} ||= [ ];
		push @{$equivClasses{$root}}, $rsPair;
		my ($read,$site) = split($;, $rsPair);
		$sitesOnRead{$read} ||= [ ];
		push @{$sitesOnRead{$read}}, $site;
	}
	return (\%equivClasses, \%sitesOnRead, $uf);
}

#########################################
sub constructConsensusFromGraph
{
##  constructs a consensus sequence for a single contig.
##  RETURNS: nothing; consensus stored in Contig object.
#########################################
	my ($contig) = @_;   ## ref to the contig
	my ($equivClasses,$sitesOnRead,$classNames) 
					= findEquivalentReadSites($contig);
	my ($rsSuccessors, $rsPredecessors)
					= findReadSiteOrder($equivClasses,$sitesOnRead);
	my %completion; ## gives the best possible nucleotide sequence.
	my %quality; ## gives individual base qualities for above.
	my %score; ## gives total score (sum of qualities) of best sequence.
	my %pending; ## number of successors with still-unknown completion.
	my @ready;

	foreach my $eqClass (keys %$equivClasses)
	{
		$pending{$eqClass} = @{$$equivClasses{$eqClass}};
	}
	
	foreach (keys %$rsPredecessors)
	{
		push @ready, $_ unless defined $$rsSuccessors{$_};
	}

	while (@ready)
	{
		my $eqClass = pop @ready;
		my ($bestString, $bestScore);
		
		foreach my $readSite (@{$$equivClasses{$eqClass}})
		{
			my ($readNum,$site) = split($;, $readSite);
	
			my $predReadSite = $$rsPredecessors{$readSite};
			my $predContigSite = $classNames->find($predReadSite);
			$pending{$predContigSite}--;
			push @ready, $predContigSite if $pending{$predContigSite} == 0;
			my $succReadSite = $$rsSuccessors{$readSite};
			next unless $succReadSite;  ## dummy sinks have no successor.
			my $succContigSite = $classNames->find($succReadSite);
			my $succSite = (split($;, $succReadSite))[1];
			my $readEntry = DnaRead->GetBySerial($readNum);
	
			my $score = $score{$classNames->find($succReadSite)};
			my @quals = @{$readEntry->AdjQual()}[$site..($succSite-1)];
			
			foreach (@quals)
			{
				$score += $_;
			}
	
			next unless ($score > $score{$eqClass});
	
			$score{$eqClass} = $score;
			my $snippet = substr($readEntry->Bases(), $site, $succSite-$site);
			$completion{$eqClass} = $snippet . $completion{$succContigSite};
			$quality{$eqClass} = join('', @quals, $quality{$succContigSite});
		}
	}

	my ($bestKey, $bestScore);
	
	foreach (@{$$equivClasses{SOURCE}})
	{
		my $succ = $$rsSuccessors{$_};
		($bestKey, $bestScore) = ($succ, $score{$succ}) 
			if $score{$succ} > $bestScore;
	}
	$contig->setLength(length($completion{$bestKey}));
	$contig->setBases($completion{$bestKey});
	$contig->setQuals($quality{$bestKey});
	$contig->setScore($bestScore);
}

#########################################
sub findReadSiteOrder
{
##  determines partial order of read sites of various reads in contig.
##  RETURNS: refs to two hashes giving successors and predecessors of
##  each read site.
#########################################
	my ($equivClasses,## UnionFind relating equivalent readsites
		$sitesOnRead  ## ref to hash mapping reads to lists of positions
		) = @_;
	## Construct a hash of successors for each read-site pair.
	my %rsSuccessors; ## maps each read-site pair to next pair on same read.
	my %rsPredecessors; ## maps read-site pair to previous pair on same read.
	## dummy equiv. class; contains dummy site for each read in contig.
	$$equivClasses{"SOURCE"} ||= [ ];
	
	foreach my $read (keys %$sitesOnRead)
	{
		my @sites = sort {$a <=> $b} @{$$sitesOnRead{$read}};
		$$sitesOnRead{$read} = \@sites;
		my $lastPair = $read . $; . 0;
		push @{$$equivClasses{"SOURCE"}}, $lastPair;
		
		foreach (@sites)
		{
			my $rsPair = $read . $; . $_;
			$rsSuccessors{$lastPair} = $rsPair;
			$rsPredecessors{$rsPair} = $lastPair;
			$lastPair = $rsPair;
		}
		my $rsPair = $read.$;.DnaRead->GetBySerial($read)->Length();
		$rsSuccessors{$lastPair} = $rsPair;
		$rsPredecessors{$rsPair} = $lastPair;
		$$equivClasses{$rsPair} = [ $rsPair ];  
	}
	return (\%rsSuccessors, \%rsPredecessors);
}

1;
