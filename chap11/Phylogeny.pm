################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
package Phylogeny;

#########################################
##  Objects of this class represent unoriented bifurcating phylogenies.
#########################################

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(new allEdges isLeaf neighbors printTable
		attachLeafBetween removeLeafBetween joinNewNeighbors
		countMutations);

use strict;

#########################################
sub new {
##  creates a new object to represent a phylogeny for given 
##  leaf sequences (with names).
##  RETURNS: blessed reference to hash representing new phylogeny.
#########################################
	my ($this,  ## usually, the literal "Phylogeny"
	$seqs,  ## reference to a list of sequences (strings)
	$names) ## if present, reference to a list of names corresponding to
			## the sequences.
	= @_;
	my @characters;
	for (my $i=0; $i<@$seqs; $i++) {
	my @states = split(//, $$seqs[$i]);
	for (my $j=0; $j<@states; $j++) {
		$characters[$j][$i] = $states[$j];
	}
	}
	$this =
	bless {characters=>\@characters,
		names=>($names||[]),
		numTraits=> scalar @characters,
		numLeaves => scalar @{$seqs},
		nextBranch => scalar @$seqs,
		tree=>[]},
	(ref($this) || $this);
	return $this;
}

#########################################
sub copy {  
## makes a "deep copy" of the phylogeny.
## RETURNS: reference to new hash representing copy of phylogeny.
#########################################
	my ($this) = @_;    ## this object -- the one to be copied.
	my @characters = map(\@$_, @{$this->{characters}});
	my @tree = map(\@$_, @{$this->{tree}});
	my $copy =
	bless {characters=> \@characters,
		names     => \@{$this->{names}},
		numTraits => $this->{numTraits},
		numLeaves => $this->{numLeaves},
		nextBranch=> $this->{nextBranch},
		tree      => \@tree},
	ref($this);
	return $copy;
}

#########################################
sub neighbors {   
## RETURNS: the list of neighbors of its arguement (a vertex)
#########################################
	my ($this,  ## reference to phylogeny object
	$v)     ## a vertex (species) in the phylogeny
	= @_;
	return @{${$this->{tree}}[$v]};
}

#########################################
sub joinNewNeighbors {    
##  modifies phylogeny to include an edge between v1 and v2
##  RETURNS: nothing
#########################################
	my ($this,     ## reference to phylogeny object to be modified.
	$v1, $v2)  ## two vertices (species) in the phylogeny.
	= @_;
	push @{$this->{tree}[$v1]}, $v2;
	push @{$this->{tree}[$v2]}, $v1;
}

#########################################
sub removeNeighbors {    
## removes v's neighbor list from phylogeny.
## RETURNS: list of removed neighbors. 
#########################################
	my ($this,     ## reference to phylogeny object to be modified.
	$v)        ## vertex whose neighbors are to be expunged.
	= @_;
	my $r = $this->{tree}[$v];
	$this->{tree}[$v] = undef;
	return @$r;
}


#########################################
sub allEdges {  
## RETURNS:  a list of all edges in the phylogenies
#########################################
	my ($this) = @_;  ## reference to phylogeny object.
	my $tree = $this->{tree};
	my @result;
	for (my $i=0; $i<@$tree; $i++) {
	next unless defined($$tree[$i]);
	foreach (@{$$tree[$i]}) {
		push @result, [$i,$_] if $i<$_;
	}
	}
	return @result;
}

#########################################
sub replaceNeighbor {
##  modifies a phylogeny to replace an edge between $v and $oldNeighbor
##  with an edge between $v and $newNeighbor.
##  RETURNS: nothing.
#########################################
	my ($this,         ## reference to phylogeny object to be modified
	$v,            ## vertex whose neighborhood is to be modified
	$oldNeighbor,  ## neighbor to be removed
	$newNeighbor)  ## neighbor to be added
	= @_;
	my $r = ${$this->{tree}}[$v];
	for (my $j=0; $j<@$r; $j++) {
	$$r[$j] = $newNeighbor if $$r[$j] == $oldNeighbor; 
	}
}

#########################################
sub isLeaf {  
## RETURNS: true if v is a leaf in the phylogeny; false if an internal node.
#########################################
	my ($this,  ## reference to a phylogeny object
	$v)     ## vertex to be tested
	= @_;
	return ($v < $this->{numLeaves});
}

#########################################
sub printTable {  
## prints the phylogeny in a tabular form.
## RETURNS: nothing. 
#########################################
	my ($this) = @_;  ## reference to a phylogeny object
	print "Phylogeny:";
#   foreach (keys %$this) { print "\n  $_ = $this->{$_}"; }
	my $tree = $this->{tree};
	my $chars = $this->{characters};
	for (my $i=0; $i<@$tree; $i++) {
	printf("\n%5d: (%2s,%2s,%2s) ", $i, @{$$tree[$i]||[]});
	for (my $j=0; $j<$this->{numTraits}; $j++) { print $$chars[$j][$i]; }
	}
	print "\n\n";
}

#########################################
sub attachLeafBetween {  
## modifies this phylogeny.
## creates a new branch node with neighbors leaf, old1, and old2.
## removes edge between old1 and old2.
#########################################
	my ($this,      ## reference to a phylogeny object
	$leaf,      ## new leaf node
	$old1,      ## endpoint of edge broken by new branch node
	$old2,      ## endpoint of edge broken by new branch node
	$changeList ## accumulates changes to ancestral sequences
					##    caused by attaching leaf, for later reversal.
	) = @_;
	my $newBranch = $this->{nextBranch}++;
	$this->replaceNeighbor($old1,$old2,$newBranch);
	$this->replaceNeighbor($old2,$old1,$newBranch);
	${$this->{tree}}[$newBranch] = [$leaf,$old1,$old2];
	${$this->{tree}}[$leaf] = [$newBranch];
	return ($this->updateSeqs($leaf,$changeList));
}


#########################################
sub detachLeafBetween {
## modifies this phylogeny.
## removes the branch node with neighbors leaf, old1, and old2.
## restores the edge between old1 and old2.
#########################################
	my ($this,      ## reference to a phylogeny object
	$leaf,      ## leaf node removed
	$old1,      ## endpoint of edge restored by removal of branch node
	$old2,      ## endpoint of edge restored by removal of branch node
	$changeList ## accumulated changes to ancestral sequences
					##    caused by attaching leaf, to be undone now.
	) = @_;
	foreach (@$changeList) { ## undo changes to character states of branches
	my ($ref, $state) = @$_;
	$$ref = $state;
	}
	my $oldBranch = ($this->removeNeighbors($leaf))[0];
	$this->removeNeighbors($oldBranch);
	$this->{nextBranch}--;
	$this->replaceNeighbor($old1,$oldBranch,$old2);
	$this->replaceNeighbor($old2,$oldBranch,$old1);
}

#########################################
sub countMutations {	      
## counts the total number of mutations for all characters in the phylogeny.
## normally, we update counts incrementally, but this method can be used for
## a final check.
#########################################
	my ($this) = @_;       ## [in] reference to current phylogeny
	my $chars = $this->{characters};
	my $n = $this->{numTraits};
	my $mutations = 0;
	foreach ($this->allEdges()) {
	my ($v1, $v2) = @$_;
	for (my $j=0; $j<$n; $j++) {
		$mutations++ if $$chars[$j][$v1] ne $$chars[$j][$v2];
	}
	}
	return $mutations;
}

#########################################
sub updateSeqs {
## updates the sequences of branch nodes to account for the addition of
## the leaf newLeaf.
## RETURNS: the total change in number of mutations brought about by
## additon of new leaf.
#########################################
	my ($this,      ## [in] reference to current phylogeny
	$newLeaf,   ## [in] index of vertex to be added to phylogeny
	$changeList  ## [out] list of all changes 
	) = @_;
	my $totalDelta = 0;
	my $newBranch = ($this->neighbors($newLeaf))[0];
	my ($leaf, $old1, $old2) = $this->neighbors($newBranch);
	
	foreach my $column (@{$this->{characters}}) {
	$$column[$newBranch] = $$column[$old1];
	my $newState = $$column[$leaf];
	next if	$$column[$newBranch] eq $newState;
	my @bestChanges;
	my $traitDelta = $this->findBestChange($column,$leaf,$newBranch,$newState,
						\@bestChanges);
	$this->makeBestChange($column,$leaf,$newBranch,$newState,
				\@bestChanges,$changeList);
	$totalDelta += ($traitDelta+1);
	}
return $totalDelta;
}

#########################################
sub makeBestChange {
## Carries out the most advantageous change found by findBestChange.
## Constructs a list of the changes made for later reversal during backtrack.
## RETURNS: nothing; results put into $changeList.
#########################################
	my ($this,       ## [in] ref to current phylogeny
	$column,     ## [in] ref to list of states of current character
	$prev,       ## [in] last tree node considered (parent of curr)
	$curr,       ## [in] current tree node
	$newState,   ## [in] state of the character in the new leaf
	$bestChanges,## [in] list of which nodes to change (if parent changes)
	$changeList  ## [out] list of all changes 
	) = @_;
	my ($next1,$next2) = grep($_ ne $prev, $this->neighbors($curr));
	return unless $$bestChanges[$curr];
	my $traitRef = \$$column[$curr];
	my $currState = $$traitRef;
#    print "push @$changeList, [$traitRef,$currState]\n";
	push @$changeList, [$traitRef,$currState];
	$$traitRef = $newState;
	$this->makeBestChange($column,$curr,$next1,$newState,$bestChanges,$changeList);
	$this->makeBestChange($column,$curr,$next2,$newState,$bestChanges,$changeList);
}

#########################################
sub findBestChange {
## Investigates whether it is advantageous to change branch node $curr's 
## character to $newState, given that neighbor $prev's state is to be changed.
## To determine this, recursive calls to the other two neighbors are necessary.
## Return value is the total reduction in mutation count achievable in the
## subtree rooted at $curr by changing $curr's state to $newState.
##
## In calculating the change in number of mutations, we include the tree
## edge from $curr to $prev as well as all edges in the subtree.
## 
## RETURNS: best possible improvement to mutation count.
#########################################
	my ($this,       ## [in] ref to current phylogeny
	$column,     ## [in] ref to list of states of current character
	$prev,       ## [in] last tree node considered (parent of curr)
	$curr,       ## [in] current tree node
	$newState,   ## [in] state of the character in the new leaf
	$bestChanges ## [out] list of which nodes to change (if parent changes)
	) = @_;
	
	my $currState = $$column[$curr];
	my $mutBefore = ($currState eq $$column[$prev]) ? 0 : 1;
	$$bestChanges[$curr] = 0;

	## Option 0:  don't change this vertex's state in subtree.
	my $mutAfter  = ($currState eq $newState) ? 0 : 1;
	my $delta0 = $mutAfter-$mutBefore;
	# changing is not always an option...
	return $delta0 if ($this->isLeaf($curr) || ($currState eq $newState));

	## Option 1:  Change this vertex's state in subtree.
#   my $mutAfter  = ($newState eq $newState) ? 0 : 1;
	my ($next1,$next2) = grep($_ ne $prev, $this->neighbors($curr));
	my $b1 = $this->findBestChange($column,$curr, $next1, $newState,$bestChanges);
	my $b2 = $this->findBestChange($column,$curr, $next2, $newState,$bestChanges);
	my $delta1 = $b1 + $b2 - $mutBefore;

	## Return, save best option.
	return $delta0 if $delta0 <= $delta1;

	$$bestChanges[$curr] = 1;
	return $delta1;
}

