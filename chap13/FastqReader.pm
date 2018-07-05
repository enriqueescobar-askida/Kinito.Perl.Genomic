################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
package FastqReader;
#########################################
## This package builds on SeqReader.  Each object of class
## FastqReader is attached to a sequence file in "FASTQ" format,
## containing sequences with one-digit quality values.
## The package provides transparent methods for reading these files
## sequence by sequence.
#########################################

require SeqReader;
@ISA = qw(SeqReader);
use strict;

#########################################
sub verifyFormat
{
##  verifies that first line of the file being opened is consistent with
##  FASTA format.
##  RETURNS: $hash blessed into class FastaReader if consistent;
##           empty string otherwise.
#########################################
	my ($this, ## normally, the literal "FastqReader".
		$hash  ## reference to a hash containing data for file attachment.
		) = @_;
	return (bless $hash) if $$hash{buff} =~ m/^[qQ]\>/;
	return "";
}

#########################################
sub readSeq
{
##  Reads the next sequence from the GenBank file.
##  Saves seqId and seqQuality in hash.
##  RETURNS: the sequence of bases/residues as a string.
#########################################
	my ($this) = @_;
	return () unless $this->{buff};
	my $id = $this->{buff};
	$id =~ s/^\>\s*//;
	$this->{seqId} = $id;
	my ($seq, $tbuff);
	my $fh = $this->{fh};
	
	while (($tbuff = <$fh>) && ($tbuff !~ /^[Qq]\>/))
	{
		chomp $tbuff;
		$tbuff = lc $tbuff;
		$tbuff =~ s/[^a-z0-9]//g;
		$seq .= $tbuff;
	}
	chomp $tbuff;
	my $quals = $seq;
	$seq =~ s/[0-9]//g;
	$quals =~ s/[a-z]//g;
	warn "Incorrect number of qualities for sequence $id in $$this{fileName}"
	if (length $quals > 0 && (length $quals != length $seq));
##    warn "No qualities for sequence $id in $$this{fileName}; using defaults"
##	if length $quals == 0;
	$quals = "1" x (length $seq) if (length $quals != length $seq);
	$this->{seqQuality} = [map {int $_} split(//, $quals)];
	$this->{buff} = $tbuff;
	return $seq;
}


#########################################
sub seqQuality
{
##  RETURNS: the string of quality information for most
##  recently read sequence.
#########################################
	my ($this) = @_;
	$this->{seqQuality};
}

1;
