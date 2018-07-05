################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
#########################################
package SeqReader;
#########################################
## This package allows DNA or protein sequences to be transparently read from
## files in fasta, GenBank, or one-sequence-per-line format.
#########################################
use strict;
use FileHandle;   ## This is a standard Perl package.
use FastaReader;
use GenBankReader;
use SimpleReader;

#########################################
sub new
## Creates a new "attachment" to the named sequence file.
## Determine the file type - fasta, GenBank, or one-per-line - for future reads.
## RETURNS: a reference to the hash representing the new SeqReader.
#########################################
{
	my ($this,     ## literal"SeqReader" or ref to an existing SeqReader
	$fileName) ## string giving operating system's name for file
	= @_;

	my $fh = *STDIN;
	if ($fileName ne "STDIN") {
	$fh = new FileHandle "<$fileName" 
		or die "can't open $fileName ";
	}
	my $buff = <$fh>;
	chomp $buff;
	my $hash = {fh=>$fh,                ## save filehandle
		buff=>$buff,            ## save line read for next read.
		fileName=>$fileName};   ## save filename
	my $reader = $this->verifyFormat($hash)
	or die "can't open $fileName with a " . (ref($this)||$this);
	return $reader;
	
}

#########################################
sub verifyFormat {
##  determines, from first line, the format of the file being opened.
##  RETURNS: $hash blessed into class indicated by first line;
##           empty string if none.
#########################################
	my ($this, ## normally, the literal "SeqReader".
	$hash  ## reference to a hash containing data for file attachment.
	) = @_;
	return(FastaReader->verifyFormat($hash)
	or GenBankReader->verifyFormat($hash)
	or SimpleReader->verifyFormat($hash));
}

#########################################
sub fileName {
##  RETURNS: name of file from which this object reads.
#########################################
	my ($this) = @_;
	$this->{fileName};
}

#########################################
sub seqId {
##  RETURNS: ID of most recently read sequence.
#########################################
	my ($this) = @_;
	$this->{seqId};
}

#########################################
sub fileType {
##  RETURNS: format of this file.
#########################################
	my ($this) = @_;
	my ($t) = (ref($this) =~ /^(.*)Reader/);
	$t;		   
}

#########################################
sub close {
##  closes file to which this object is attached and cleans out hash.
##  RETURNS: nothing.
#########################################
	my ($this) = @_;
	$this->{fh}->close();
	delete @{$this}{keys %$this};
}

1;


