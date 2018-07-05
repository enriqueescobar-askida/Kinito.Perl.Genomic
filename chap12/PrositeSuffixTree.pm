################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
package PrositeSuffixTree;

#########################################
##  This package extends SuffixTree to allow the target string to be
##  searched for substrings that match a PROSITE pattern in addition 
##  to simple substrings.
#########################################

require SuffixTree;
@ISA = qw(SuffixTree);
use strict;

#########################################
sub listMotifs {
##  
#########################################
	my ($this,   ## reference to PrositeSuffixTree for target to be matched
	$motif)  ## reference to a hash containing a PROSITE motif.
	= @_;
	my $pattern = lc $$motif{PA};
	return unless $pattern;
	my $description = $$motif{DE};
	$description =~ s/\n//gs;   ## remove newlines
	return "Can't handle $description\n" if $pattern =~ /\,/;
	foreach ($this->searchForPattern($pattern)) {
	print "$description seen at position $_\n";
	}
}

#########################################
sub pairForm {
##  Breaks a pattern element in Perl regular expression form into
##  its character-set and repetition parts, and returns these as a pair
##  of strings.
##  RETURNS: the pair of strings.
#########################################
	my ($perlPatElem) = @_;  ## a PROSITE pattern element in Perl r.e. form
	my ($spec, $rep) = ($perlPatElem =~ /^([^\{]*)(\{.*\})?$/);
	$rep =~ s/[\{\}]//g;
	$rep ||= 1;
	return ($spec, $rep);
}

#########################################
sub searchForPattern { 
##  Searches for all occurrences of given PROSITE pattern in target matching
##  represented by PrositeSuffixTree.
##  RETURNS: a list of all starting points of motifs in the target.
#########################################
	my ($this,            ## reference to PrositeSuffixTree for target.
	$pattern) = @_;   ## PROSITE motif in PROSITE form (string).
	$pattern =~ s/\n//gs;   ## remove newlines
	$pattern =~ s/\.$//;
	$pattern =~ s/^</<-/;
	$pattern =~ s/>$/->/;
	my $patLinks = [];
	foreach (reverse split("-",$pattern)) {
	$patLinks = [pairForm(Prosite::perlizeElement($_)), $patLinks];
	}
	return ($this->searchHelp($this->{root}, $patLinks));
}

#########################################
sub searchHelp {
##  Recursive auxiliary routine to searchForPattern.
##  RETURNS: list of positions in target where matches to 
##  pattern can be found.
#########################################
	my ($this,     ## reference to PrositeSuffixTree.
	$node,     ## reference to current node in tree.
	$patLinks) ## reference to current position in pattern.
	= @_;
	my $ssLen = $$node{len};
	my $ss = substr($this->{target},$$node{off},$ssLen);
#   print "At node $node representing $ss \n";
	my ($spec, $rep);
	while ($ss && @$patLinks) {
	($spec, $rep, $patLinks) = @$patLinks;
	my $pat = "$spec\{$rep\}";
#	print "Comparing $ss to $pat \n";
	if ($ss !~ s/^$pat//) {
		return () if ($ssLen > $rep);
		my $pat = "$spec\{$ssLen\}";
		return () if ($ss !~ s/^$pat//);
	}
	$ssLen -= $rep;
	}
	## Reached end of substring for this node.
	$patLinks = [$spec, -$ssLen, $patLinks] if $ssLen<0;
	if (@$patLinks==0) { return ($this->gatherStartPoints($node)); }
	$spec = $$patLinks[0];
	if ($spec =~ /[\^\.]/) {
#	my @goodKeys = grep(/^$spec$/, (keys %$node));
#	print "Spec $spec; recursive calls for keys @goodKeys\n";
	return map(($this->searchHelp($$node{$_}, $patLinks)),
		grep(/^$spec$/, (keys %$node)));
	} 
	else {
#	my @goodKeys = grep($$node{$_},(split //,$spec));
#	print "Spec $spec; recursive calls for keys @goodKeys\n";
	return map(($this->searchHelp($$node{$_}, $patLinks)),
		grep($$node{$_},(split //,$spec)));
	}
}
1;
