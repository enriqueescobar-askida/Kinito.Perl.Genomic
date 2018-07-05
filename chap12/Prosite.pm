################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
#########################################
package Prosite;
use strict;
use FileHandle;

#########################################
sub new
## Creates a new "attachment" to a PROSITE database.
#########################################
{
	my ($this,     ## literal "Prosite" or ref to an existing Prosite
	$fileName) ## string giving operating system's name for file
	= @_;

	my $fh = *STDIN;
	if ($fileName ne "STDIN") {
	$fh = new FileHandle "<$fileName" 
		or die "can't open $fileName ";
	}
	bless {fh=>$fh}, (ref($this) || $this);
}


#########################################
sub readMotif {
##  Reads the next PROSITE motif from the database.
##  RETURNS: reference to a hash containing the motif ID, etc., as labeled
##  in the first two columns of each line.
#########################################
	my ($this) = @_;  ## reference to this PROSITE reader
	my $fh = $this->{fh};
	my $line;
	while (($line = <$fh>) && ($line !~ /^ID/)) { };
	return undef unless $line;
	my %entry;
	$entry{ID} = substr($line,5);
	while (($line = <$fh>) && ($line !~ m|^//|)) {
	$line =~ s/(..)...//;
	$entry{$1} .= $line;
	}
	return \%entry;
}

#########################################
sub close {
##  Releases this attachment to the PROSITE database and frees space.
##  RETURNS: nothing.
#########################################
	my ($this) = @_;  ## reference to this PROSITE reader
	$this->{fh}->close();
	delete @{$this}{keys %$this};
}

#########################################
sub perlizePattern {
##  Translates a PROSITE pattern into a Perl regular expression.
##  RETURNS: the Perl regular expression (string).
#########################################
	my ($this,          ## reference to this PROSITE reader
	$pattern) = @_; ## a PROSITE pattern (string)
	$pattern =~ s/\.$//;
	my ($left,$right);
	my $left = '^' if ($pattern =~ s/^\<//);
	my $right = '$' if ($pattern =~ s/\>$//);
	my $newpattern = join('', map(perlizeElement($_), split('-',$pattern)));
	return $left . $newpattern . $right;
}

#########################################
sub perlizeElement {
##  Translates a single element of a PROSITE pattern to Perl.
##  Auxiliary function for perlizePattern.
##  Designed for internal use only -- not a method.
#########################################
	my ($elem) = @_;  ## single element of a PROSITE pattern.
	my ($residue,$repetition) = ($elem =~ /^([^\(]*)(\(.*\))?$/);
	$repetition =~ tr|()|{}|;
	$residue =~ tr/xX/../;
	return "[^$1]$repetition" if $residue =~ /\{(.*)\}/;
	return "$residue$repetition";
}

1;
