################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
package	FastaReader;

require	SeqReader;
@ISA	= qw(SeqReader);

use	strict;
use	diagnostics;

sub	verifyFormat
{
	my ($this, $hash) = @_;
	return	(bless $hash)	if $$hash{buff} =~ m/^\>/;
	return	"";
}

#########################################
## Reads, returns sequence and id from fasta file.
#########################################
sub	readSeq
{
	my ($this) = @_;
	return	()	unless $this->{buff};
	my $id = $this->{buff};
	$id =~ s/^\>\s*//;
	$this->{seqId} = $id;
	my ($seq, $tbuff);
	my $fh = $this->{fh};
	while (($tbuff = <$fh>) && ($tbuff !~ /^\>/))
	{
		chomp $tbuff;
		$tbuff = lc $tbuff;
		$tbuff =~ s/[^a-z]//g;
		$seq .= $tbuff;
	}
	chomp $tbuff;
	$this->{buff} = $tbuff;
	return	$seq;
}



