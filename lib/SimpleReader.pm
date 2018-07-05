################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################

package SimpleReader;

#########################################
## This package builds on SeqReader.  Each object of class
## SimpleReader is attached to a sequence file in "simple" format:
## one sequence per line, possibly preceded by label+colon.
## The package provides transparent methods for reading these files
## sequence by sequence.
#########################################

require	SeqReader;
@ISA	= qw(SeqReader);

use	strict;
use	diagnostics;
#########################################
##  verifies that first line of the file being opened is consistent with
##  simple format.
##  RETURNS: $hash blessed into class SimpleReader if consistent;
##           empty string otherwise.
#########################################
sub	verifyFormat
{
	my ($this, ## normally, the literal "SimpleReader".
		$hash  ## reference to a hash containing data for file attachment.
		) = @_;
	return	(bless $hash);
}

#########################################
##  Reads the next sequence from the simple file.
##  Saves seqId in hash.
##  RETURNS: the sequence of bases/residues as a string.
#########################################
sub	readSeq
{
	my ($this) = @_;   ##  reference to this SimpleReader
	return	()	unless $this->{buff};
	my ($id, $seq) = ($this->{buff} =~ m/(.*):(.*)$/);
	$this->{seqId} = $id;
	$seq ||= $this->{buff};
	$seq =~ s/[^a-z]//g;
	my $fh = $this->{fh};
	$this->{buff} = <$fh>;
	chomp $this->{buff};
	return	$seq;
}

1;
