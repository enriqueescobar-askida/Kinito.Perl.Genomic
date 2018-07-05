#!/usr/bin/perl -w -I . -I ../lib -I ../lib2
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################

#########################################
##  This program implements compact (space-saving) strategies for
##  computing similarity scores and alignments of sequences.
#########################################
use		strict;
use		diagnostics;

use		Util;

#########################################
##  determines the similarity score of two sequences in space
##  proportional to the length of the longer string rather than 
##  proportional to the product of the lengths of the two strings.
#########################################
sub compactSimilarity
{
	my ($s,$t) = @_;  ## the two sequences to be scored
	my ($sLen, $tLen) = (length($s),length($t));
	my @a;
	$#a = $tLen;  ## fast allocate.
	foreach my $j (0..$tLen)
	{
		$a[$j] = -2*$j;
	}
	foreach my $i (1..$sLen)
	{
		unshift @a, -2*$i;
		foreach my $j (1..$tLen)
		{
			## Loop Invariant:
			## a[k] = M[i-1][k-1]  for k>=j; a[k] = M[i][k] for k<j;
			my $m =  (substr($s,$i-1,1) eq substr($t,$j-1,1)) ? +1 : -1;
			$a[$j] = max($a[$j]+$m, $a[$j-1]-2, $a[$j+1]-2);
		}
		pop @a;
	}
	return	(@a);
}

#########################################
##  determines the actual optimal alignment of two sequences in space
##  proportional to the length of the longer string rather than 
##  proportional to the product of the lengths of the two strings.
##  Uses a divide-and-conquer strategy.
#########################################
sub compactAlignment
{
	my ($s,$t) = @_;  ## the two sequences to be aligned.
	my ($sLen, $tLen) = (length($s),length($t));
	if ($sLen == 1)
	{
		if ($tLen==0)
		{
			return	"$s\n-\n-2";
		}
		elsif ($t =~ /^(.*)$s(.*)$/)
		{
			my $ss = "-"x length($1) . $s . "-"x length($2);
			my $score = 3 - 2*$tLen;
			return	"$ss\n$t\n$score";
		}
		else
		{
			my $score = 1 - 2*$tLen;
			return	"-"x($tLen-1) . "$s\n$t\n$score";
		}
	}
	else
	{
		my $mid = int($sLen/2);
		my $sBeg = substr($s,0,$mid);
		my $sEnd = substr($s,$mid, $sLen-$mid);
		my @aBeg = compactSimilarity($sBeg, $t);
		my @aEnd = reverse(compactSimilarity((scalar reverse $sEnd), 
						(scalar reverse $t)));
		my ($kMax, $kMaxVal) = (0, $aBeg[0]+$aEnd[0]);
		foreach my $k (1..$tLen)
		{ 
			($kMax, $kMaxVal) = ($k, $aBeg[$k]+$aEnd[$k])
			if ($aBeg[$k]+$aEnd[$k]) > $kMaxVal;
		}
		my ($ssBeg,$ttBeg,$scoreBeg) = 
			split("\n", compactAlignment($sBeg, substr($t,0,$kMax)));
		my ($ssEnd,$ttEnd,$scoreEnd) = 
			split("\n", compactAlignment($sEnd, substr($t,$kMax,$tLen-$kMax)));
			my $score = $scoreBeg + $scoreEnd;
		return	"$ssBeg$ssEnd\n$ttBeg$ttEnd\n$score";
	}
}

#########################################
##  MAIN PROGRAM
##  read examples from bottom of file and process.
#########################################
{
	while (my $s = <DATA>)
	{
		chomp($s);
		my $t = <DATA>; chomp($t);
		print "\n", compactAlignment($s,$t), "\n";
	}
}


__END__
SASKATCHEWAN
SESQUICENTENNIAL
RECONSTRUCTION
UNCONDITIONALLY
THANKSGIVING
CHANGELING
ACURSEDFIENDBRINGSGRIEFANDPAIN
ABLESSEDFRIENDBRINGSRELIEFAGAIN
B
BC
AB
ABC
ABC
AB
AC
ABC
ABC
AC
