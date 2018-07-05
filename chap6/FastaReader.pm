################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
package FastaReader;
#########################################
## This package builds on SeqReader.  Each object of class
## FastaReader is attached to a sequence file in FASTA format.
## The package provides transparent methods for reading these files
## sequence by sequence.
#########################################

require SeqReader;
@ISA = qw(SeqReader);
use strict;

#########################################
sub verifyFormat {
##  verifies that first line of the file being opened is consistent with
##  FASTA format.
##  RETURNS: $hash blessed into class FastaReader if consistent;
##           empty string otherwise.
#########################################
	my ($this, ## normally, the literal "FastaReader".
	$hash  ## reference to a hash containing data for file attachment.
	) = @_;
	return (bless $hash) if $$hash{buff} =~ m/^\>/;
	return "";
}

#########################################
sub readSeq {
##  Reads the next sequence from the GenBank file.
##  Saves seqId in hash.
##  RETURNS: the sequence of bases/residues as a string.
#########################################
	my ($this) = @_;
	return () unless $this->{buff};
	my $id = $this->{buff};
	$id =~ s/^\>\s*//;
	$this->{seqId} = $id;
	my ($seq, $tbuff);
	my $fh = $this->{fh};
	while (($tbuff = <$fh>) && ($tbuff !~ /^\>/)) {
	chomp $tbuff;
	$tbuff = lc $tbuff;
	$tbuff =~ s/[^a-z]//g;
	$seq .= $tbuff;
	}
	chomp $tbuff;
	$this->{buff} = $tbuff;
	return $seq;
}



