#!/usr/bin/perl -w -I . -I ../lib -I ../lib2
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################

use strict; 

#########################################
##  GLOBAL VARIABLES
#########################################
my @structure;  ## list of dots and parens describing input structure.
my @bases;      ## list of bases from input string.
### Largest and smallest x- and y-coordinates used in the drawing;
### ultimately used to scale drawing onto one page.
my ($minx, $maxx, $miny, $maxy);


#########################################
sub max
{
##  RETURNS: the larger of its two arguments.
#########################################
    my ($x,$y) = @_;
    if ($x>$y)
    {
        return $x;
    }
    else
    {
        return $y
    };
}

#########################################
sub begin_PostScriptPicture
{
##  issues PostScript commands needed at beginning of .ps file.
##  RETURNS: nothing.
#########################################
    print "/rnaline { pt setlinewidth moveto lineto stroke} def\n";
    print "/rnabase { moveto -0.12 -0.12 rmoveto show} def\n";
    print "/rnapicture {\n\n";
}

#########################################
sub end_PostScriptPicture
{
##  issues PostScript commands needed at end of .ps file.
##  RETURNS: nothing.
#########################################
    print "\n} def\n\n";
    my $scale = 0.95 * (8.5*72) / max($maxx-$minx, $maxy-$miny);
    my $xorigin = ((8.5*72) * -$minx / ($maxx-$minx))  || "0";
    my $yorigin = ((8.5*72) * -$miny / ($maxy-$miny)) || "0";
    print "/pt {$scale div} def\n\n";
    print "/Helvetica findfont ", 8, " pt scalefont setfont\n";
    print "$xorigin $yorigin translate\n";
    print "$scale dup scale\n";
    print "rnapicture showpage\n";
}


#########################################
sub drawLine
{
##  issues PostScript commands needed to draw a line.
##  RETURNS: nothing.
#########################################
    my ($x1, $y1,  ## coordinates of one endpoint
    $x2, $y2,  ## coordinates of other endpoint
    $thick     ## thickness of line
    ) = @_;
    ($x1, $y1, $x2, $y2) = (0.8 * $x1 + 0.2 * $x2,
                0.8 * $y1 + 0.2 * $y2,
                0.2 * $x1 + 0.8 * $x2,
                0.2 * $y1 + 0.8 * $y2);
    print "$x1 $y1 $x2 $y2 $thick rnaline\n";
    ## record this info to know overall size of picture.
    ($x1,$x2) = ($x2,$x1) if ($x1>$x2);
    $maxx = $x2 if $x2>$maxx;
    $minx = $x1 if $x1<$minx;
    ($y1,$y2) = ($y2,$y1) if ($y1>$y2);
    $maxy = $y2 if $y2>$maxy;
    $miny = $y1 if $y1<$miny;
}

#########################################
sub drawBase
{
##  issues PostScript commands needed to print letter for a base.
##  RETURNS: nothing.
#########################################
    my ($x, $y, ## coordinates of base
    $b      ## letter (character) to print; can be 5 or 3
    ) = @_;
    $b .= "'" if ("53" =~ $b);
    print "($b) $x $y rnabase\n"
}

#########################################
sub drawRna
{
##  assists drawRnaStructure to create PostScript drawing of RNA.
##  Invokes itself recusively for each ring of the structure.
##  RETURNS: nothing.
#########################################
    my ($l, $r,   ## range of @structure and @bases to be drawn.
    $lx, $ly, ## coordinates of first position of ring.
    $rx, $ry  ## coordinates of last position of ring.
    ) = @_;
    my $level=0;
    my $count=2;
    for (my $i=$l+1; $i<$r; $i++)
    {
        $level-- if ($structure[$i] eq ")");
        $count++ if $level==0;
        $level++ if ($structure[$i] eq "(");
    }

    my $theta = 2 * 3.14159 / $count;
    my $rad = 1 / (2*sin($theta/2));            ## radius 
    my $h = $rad * cos($theta/2);
    my ($cx, $cy) = ((($lx+$rx)/2.0)+$h*($ly-$ry),    ## center of circle
            (($ly+$ry)/2.0)+$h*($rx-$lx)); 
    my $alpha = atan2($ly-$cy,$lx-$cx);

    my ($ii,$xx,$yy) = ($l,$lx,$ly);
    
    for (my $i=$l+1; $i<=$r; $i++)
    {
        $level-- if ($structure[$i] eq ")");
        if ($level==0)
        {
            $alpha -= $theta;
            my ($x,$y) = ($cx+$rad*cos($alpha), $cy+$rad*sin($alpha));
            drawLine($xx,$yy,$x,$y,0);
            drawRna($ii,$i,$xx,$yy, $x, $y) if ($structure[$i] eq ")");
            drawBase($xx,$yy, $bases[$ii]);
            ($xx,$yy)=($x,$y);
            $ii = $i;
        }
        $level++ if ($structure[$i] eq "(");
    }
    drawLine($xx,$yy,$rx,$ry,0);
    drawBase($xx,$yy, $bases[$r-1]);
    drawBase($rx,$ry, $bases[$r]);
    my %bonds = (GU=>1,UG=>1,AU=>2,UA=>2,CG=>3,GC=>3);
    drawLine($lx,$ly,$rx,$ry,$bonds{$bases[$l].$bases[$r]}+0)
            unless ($lx==0 || $ly==0);   ## 3-5 pair               
}

#########################################
sub drawRnaStructure
{
##  creates PostScript drawing of RNA structure.
##  Most of the work is delegated to the recursive subroutine drawRna.
##  RETURNS: nothing.
#########################################
    my ($basestring,$structurestring) = @_;
    @bases = split(//, 5 . $basestring . 3);
    @structure = split(//, "($structurestring)");
    begin_PostScriptPicture();
    drawRna(0, $#structure, 0,0, 1,0);
    end_PostScriptPicture();
}


#########################################
## MAIN PROGRAM
#########################################
my $basestring = <STDIN>;
chomp($basestring);
my $parenstring = <STDIN>;
chomp($parenstring);
drawRnaStructure($basestring,$parenstring);
