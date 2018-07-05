################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
#########################################
package Contig;
#########################################
## This package defines a class of objects representing contigs during the
## sequence assembly process.
#########################################
use strict;

#########################################
## PACKAGE VARIABLES
#########################################

my $SerialNumber; ## Unique integer identifying next contig to be created.

my @ImmutableFields =  ## List of fields that cannot change
			("SerialNumber");  ## Records order in which Contig objects are created.

my @MutableFields = ## List of fields that can change
			("Bases",       ## The consensus sequence of the contig; a string.
			"Quals",       ## Ref to list of qualities caller assigned to each base.
			"Reads",       ## Ref to list of included reads.
			"AlignedPairs",## Ref to list of alignments pertaining to this contig.
			"Length",      ## Estimated length of contig.
			"FirstStart",  ## Where leftmost read starts, wrt contig's 0. (Often <0).
			"LastEnd",     ## Where rightmost read end, relative to contig's 0.
			"Score");      ## Sum of nucleotide qualities for this contig.

#########################################
## PACKAGE INITIALIZATION -- define access subroutines
#########################################
foreach (@MutableFields)
{
	eval "sub set$_ { \$_[0]->{$_} = \$_[1] } ;" ;
}

foreach (@MutableFields, @ImmutableFields)
{
	eval "sub $_ { \$_[0]->{$_} } ;" ;
}


#########################################
sub new
##  creates a new Contig object.
##  RETURNS: blessed reference to the new Contig
#########################################
{
	my ($this,       ## literal "Contig" or ref to an existing Contig
		$reads,      ## ref to list of DnaRead objects
		$numReads,   ## number of reads in the contig
		$firstStart, ## position in contig where leftmost read starts 
		$lastEnd,    ## position in contig where rightmost read ends
		) = @_;

	return bless
	{
		Reads => $reads,
		NumReads => $numReads,
		FirstStart => $firstStart,
		LastEnd => $lastEnd,
		AlignedPairs => [ ],
		SerialNumber => $SerialNumber++
	}, ($this || ref $this);
}


#########################################
sub dump
##  dumps the Contig in human readable form.
##  RETURNS: nothing.
#########################################
{
	my ($this) = @_;  ## ref to an existing Contig.
	print "****************************************************************\n";
	print "Contig $$this{SerialNumber} at $this\n";
	print "****************************************************************\n";
	
	foreach ((sort @ImmutableFields),(sort @MutableFields))
	{
		print "   $_ = $$this{$_}\n";
	}
	print "This Contig contains the following reads:\n";
	
	foreach my $read (@{$this->Reads()})
	{
		$read->summarize();
	}
}

1;
