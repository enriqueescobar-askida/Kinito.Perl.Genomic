#!/usr/bin/perl -w -I . -I ../lib -I ../lib2
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################
use		strict;		## request conservative error checking
use		diagnostics;

#########################################
## GLOBAL VARIABLES
#########################################
my		%codonMap;	## declare a hash table to hold the residue abbreviations
					## corresponding to each possible codon.

#########################################
##  Software transcriptase: turns a DNA string into an RNA string.
##  RETURNS: the RNA string.
#########################################
sub	transcribe 
{
	my ($dna) = @_;  ## the DNA string to be transcribed
	my $rna = scalar reverse $dna;
	$rna =~ tr/ACGT/UGCA/;
	return	$rna;
}


#########################################
##  Software ribosome: translates mRNA strings to proteins
##  RETURNS: the protein string. (3-letter abbreviations.)
#########################################
sub	translate 
{
	my ($mrna) = @_;  ## the mRNA string to be translated.
	my $pro = "";
	while ( $mrna =~ s/(...)// )
	{
		$pro = $pro . $codonMap{$1};
	}
	return	$pro;
}

#########################################
## MAIN PROGRAM
#########################################
## Construct hash that maps codons to amino acids by reading table
## from DATA lines at the end of the program, which have the form:
## Residue Codon1 Codon2 ...

while (my $in = <DATA> )
{  ## assigns next line of DATA to $in; fails if none
	chomp($in);             ## remove line feed at end of $in
	my @codons = split " ",$in;
	my $residue = shift @codons;  ##remove first item from @codons and assign
	foreach my $nnn (@codons)
	{
		$codonMap{$nnn} = $residue;
	}
}

## Now read DNA strands from input <STDIN> and print translations in all six
## possible reading frames
while ( my $dna = <STDIN> )
{
	chomp($dna);
	print "DNA: ", $dna, "\n";
	my $rna = transcribe($dna);
	print "RNA: ", $rna, "\n";
	my $protein = translate($rna);
	print "RF1: ", $protein, "\n";
	$rna =~ s/.//;
	$protein = translate($rna);
	print "RF2:  ", $protein, "\n";
	$rna =~ s/.//;
	$protein = translate($rna);
	print "RF3:   ", $protein, "\n\n";
}
    
## The lines below are not perl statements and are not executed as part of the 
## program.  Instead, they are available to be read as input by the program
## using the file handle DATA.
__END__
Ala GCU GCC GCA GCG
Arg CGU CGC CGA CGG AGA AGG
Asn AAU AAC
Asp GAU GAC 
Cys UGU UGC
Gln CAA CAG
Glu GAA GAG
Gly GGU GGC GGA GGG
His CAU CAC
Ile AUU AUC AUA
Leu UUA UUG CUU CUC CUA CUG
Lys AAA AAG
Met AUG
Phe UUU UUC
Pro CCU CCC CCA CCG
Ser UCU UCC UCA UCG AGU AGC
Thr ACU ACC ACA ACG
Trp UGG
Tyr UAU UAC
Val GUU GUC GUA GUG
... UAA UAG UGA
