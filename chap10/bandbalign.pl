.          0 �4�5252  ��~2R    ..          �4�5252  4�52}�    �OMPACT PL   �4�52F�  4�52��                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   @_;  ## reference to list of strings
    my $score;
    foreach my $i (0..length($$alignment[0])-1)
    {
        $score += scoreColumn(map {substr($_,$i,1)} @$alignment);
    }
    return $score;
}

#########################################
sub printAlignment
{
##  prints a multiple alignment with a header.
##  RETURNS: nothing.
#########################################
    my ($title,      ## heading for alignment
    $alignment)  ## reference to list of strings.
    = @_;
    print "\n***** $title\n";
    foreach (@$alignment)
    {
        print "$_\n";
    }
    my $score = scoreMultipleAlignment($alignment);
    print "score $score\n\n";
}    

#########################################
##  MAIN PROGRAM
#########################################
my @ss;
while (my $s = <DATA>)
{
    chomp($s);
    push @ss, $s if $s;  ## tolerate blank lines
}

my ($tunnelScore, $tunnelAlignment) = tunnelAlign(2,@ss);
printAlignment("Tunnel Strategy", $tunnelAlignment);

my ($lcScore, $lcAlignment) = branchBoundAlign($tunnelScore, @ss);
printAlignment("Branch & Bound Strategy", $lcAlignment);



__END__
GVLTDVQVALVKSSFEEFNANIPKNTHRFFTLVLEIAPGAKDLFSFLKGSS
SPLTADEASLVQSSWKAVSHNEVEILAAVFAAYPDIQNKFSQFAGK
VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLS
VHLSGGEKSAVTNLWGKVNINELGGEALGRLLVVYPWTQRFFEAFGDLS
VLSAADKTNVKGVFSKIGGHAEEYGAETLERMFIAYPQTKTYFPHFDLS
