################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
package DnaRead;
#########################################
## This package defines a class of objects representing DNA reads
## (sequence plus qualities).
#########################################
use strict;
use diagnostics										;
use Util;
#########################################
## PACKAGE VARIABLES
#########################################
my @AllReads;  ## List of refs to all DnaRead objects created so far.
my @ImmutableFields = 
	("Bases",       ## The nucleotide sequence of the read; a string.
	"Length",      ## Number of nucleotides in the read.
	"Name",        ## The name of the sequence; a string.
	"OrigQual",    ## Ref to list of qualities caller assigned to each base.
	"SerialNumber");  ## Records index of object in @AllReads.
my @MutableFields = 
	("Alignments",  ## Ref to list of alignments to this read.
	"AdjQual",     ## Ref to list of current quality estimates.
	"Contig",      ## Ref to contig to which read is currently assigned.
	"ContigStart", ## Approx. interval of contig covered by this read.
	"ContigEnd");   ## 
#########################################
## PACKAGE INITIALIZATION -- define access subroutines
#########################################
foreach (@MutableFields) {
	eval "sub set$_ { \$_[0]->{$_} = \$_[1] } ;" ;
}
foreach (@MutableFields, @ImmutableFields) {
	eval "sub $_ { \$_[0]->{$_} } ;" ;
}
#########################################
sub new {
##  creates a new DnaRead object.
##  RETURNS: blessed reference to the new DnaRead
#########################################
	my ($this,     ## literal "DnaRead" or ref to an existing DnaRead
	$name,     ## identifier of this sequences
	$bases,    ## the bases
	$qualities ## qualities corresponding to bases
	) = @_;
	my $this = bless {Bases => $bases, 
		Length => length($bases),
		Name => $name,
		OrigQual => $qualities,
		AdjQual => \@$qualities,  ## Makes a copy of the array.
		SerialNumber => scalar(@AllReads)}, ($this || ref $this);
	push @AllReads, $this;
	return $this;
}
#########################################
sub GetBySerial {
##  retrieves the DnaRead object with given serial number.
##  RETURNS: ref to the DnaRead object.
#########################################
	my ($this,  ## literal "DnaRead" or some DnaRead object (ignored).
	$num)   ## serial number of desired DnaRead.
	= @_;
	return $AllReads[$num];
}
#########################################
sub dump {
##  dumps the read in human readable form.
##  RETURNS: nothing.
#########################################
	my ($this) = @_;
	
	print "================ Dump of $this =======================\n";
	foreach ((sort @ImmutableFields),(sort @MutableFields)) {
	print "   $_ = $$this{$_}\n";
	}
	print "   OrigQual = ", join('',@{$this->{OrigQual}}), "\n";
	print "    AdjQual = ", join('',@{$this->{AdjQual}}), "\n";
	print "=======================================================\n";
}
#########################################
sub summarize {
##  summarizes the read in human readable form; no sequence or quality.
##  RETURNS: nothing.
#########################################
	my ($this) = @_;
	printf("**> read %3d at (%6d,%6d) in contig %3d ( %s )\n",
	@$this{"SerialNumber","ContigStart","ContigEnd"},
	$this->Contig()->SerialNumber(), 
	$this->Name());
}
1;
