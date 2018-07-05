################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################

package	SeqReader;

#########################################
## This package allows DNA or protein sequences to be transparently read from
## files in fasta, GenBank, or one-sequence-per-line format.
#########################################
use	strict;
use	diagnostics;
use	FileHandle;
use	FastaReader;
use	FastqReader;
use	GenBankReader;
use	SimpleReader;

#########################################
## Creates a new "attachment" to the named sequence file.
## Determine the file type - fasta, GenBank, or one-per-line - for future reads.
#########################################
sub	new
{
	my ($this,     ## literal"SeqReader" or ref to an existing SeqReader
		$fileName) ## string giving operating system's name for file
					= @_;

	my	$fh = *STDIN;
	if	($fileName ne "STDIN")
	{
		$fh = new FileHandle "<$fileName"	or die	"can't open $fileName ";
	}
	my $buff = <$fh>;
	chomp $buff;
	my $hash = {fh=>$fh,                ## save filehandle
				buff=>$buff,            ## save line read for next read.
				fileName=>$fileName};   ## save filename
	my $reader = $this->verifyFormat($hash)
	or die "can't open $fileName with a " . (ref($this)||$this);
	return	$reader;
}


sub	verifyFormat
{
	my ($this, $hash) = @_;
	return(FastaReader->verifyFormat($hash)
	or GenBankReader->verifyFormat($hash)
	or FastqReader->verifyFormat($hash)
	or SimpleReader->verifyFormat($hash));
}

sub	fileName
{
	my ($this) = @_;
	$this->{fileName};
}

sub	seqId
{
	my ($this) = @_;
	$this->{seqId};
}

sub	fileType
{
	my ($this) = @_;
	my ($t) = (ref($this) =~ /^(.*)Reader/);
	$t;
}

sub	close
{
	my ($this) = @_;
	$this->{fh}->close();
	delete @{$this}{keys %$this};
}

1;
