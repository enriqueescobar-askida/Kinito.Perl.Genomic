################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
package GenBankReader;

#########################################
## This package builds on SeqReader.  Each object of class
## GenBankReader is attached to a sequence file in GenBank format.
## The package provides transparent methods for reading these files
## sequence by sequence.
#########################################

require SeqReader;
@ISA = qw(SeqReader);
use strict;

#########################################
sub verifyFormat {
##  verifies that first line of the file being opened is consistent with
##  GenBank format.
##  RETURNS: $hash blessed into class GenBankReader if consistent;
##           empty string otherwise.
#########################################
	my ($this, ## normally, the literal "GenBankReader".
	$hash  ## reference to a hash containing data for file attachment.
	) = @_;
	return (bless $hash) if $$hash{buff} =~ m/^LOCUS/;
	return "";
}

#########################################
sub readSeq {
##  Reads the next sequence from the GenBank file.
##  Saves seqId and seqNotes in hash.
##  RETURNS: the sequence of bases/residues as a string.
#########################################
	my ($this) = @_;   ##  reference to this GenBankReader
	return () unless $this->{buff};
	my $id = $this->{buff};
	$id =~ s/^LOCUS\s\s*//;
	$this->{seqId} = $id;
	my ($seq, $notes, $tbuff);
	my $fh = $this->{fh};
	while (($tbuff = <$fh>) && ($tbuff !~ /^ORIGIN/)) {
	$notes .= $tbuff;
	}
	$notes .= $tbuff;
	$this->{seqNotes} = $id;
	while (($tbuff = <$fh>) && ($tbuff !~ /^LOCUS/)) {
	chomp $tbuff;
	$tbuff = lc $tbuff;
	$tbuff =~ s/[^a-z]//g;
	$seq .= $tbuff;
	}
	$this->{buff} = $tbuff;
	return $seq;
}

#########################################
sub complement {
##  RETURNS: the reverse complement of its string argument.
#########################################
	my ($s) = @_;
	$s =~ tr/acgtACGT/tgcaTGCA/;
	return reverse $s;
}

#########################################
sub ident {
##  RETURNS: its argument without modification.
#########################################
	my ($s) = @_;
	return $s;
}


#########################################
sub extract_feature {
##  RETURNS: the DNA sequence extracted from $s at location $location.
#########################################
	my ($this,      ## reference to this GenBankReader object.
	$location,  ## a GenBank feature location (string).
	$s)         ## a sequence from which to extract the string
					##   -- a string or a reference to a string.
	= @_;
	print @_, "\n";
	## Allow $s to be either a string or a reference to a string.
	unless (ref $s) {my $t = $s; $s = \$t; }
	$location = lc $location;
	$location =~ s/,/\./g;
	$location =~ s/join/ident/g;
	$location =~ s/(\d+)/N$1/g;
	$location =~ s/N(\d+)\.\.N(\d+)/substr(\$\$s,$1-1,$2+1-$1)/g;
	$location =~ s/N(\d+)/substr(\$\$s,$1-1,1)/g;
	return eval($location);
}





