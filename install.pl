
### Create the directory where you would like to install the Genomic Perl
### software, and copy this script into that directory.
### Then, edit the next line to give the mount point of the CD ROM drive.
### Finally, execute the script with the command   "perl install.pl"

my $cdrom = "/mnt/cdrom";    ## typical

foreach (1..17)
{
    install_chapter("chap$_");
}
foreach ("appa", "appb", "appc", "perllib")
{
    install_chapter($_);
}

sub install_chapter
{
    my $chap = shift;
    system_or_die("mkdir $chap");
    system_or_die("cp $cdrom/unix/$chap/* $chap");
}

system_or_die("mv appc/unionfind.pm appc/UnionFind.pm");
system_or_die("mv chap11/phylogeny.pm chap11/Phylogeny.pm");
system_or_die("mv chap12/prosite.pm chap12/Prosite.pm");
system_or_die("mv chap12/prositesuffixtree.pm chap12/PrositeSuffixTree.pm");
system_or_die("mv chap12/suffixtree.pm chap12/SuffixTree.pm");
system_or_die("mv chap13/alignment.pm chap13/Alignment.pm");
system_or_die("mv chap13/consensusbuilder.pm chap13/ConsensusBuilder.pm");
system_or_die("mv chap13/contigbuilder.pm chap13/ContigBuilder.pm");
system_or_die("mv chap13/contig.pm chap13/Contig.pm");
system_or_die("mv chap13/dnaread.pm chap13/DnaRead.pm");
system_or_die("mv chap13/fastqreader.pm chap13/FastqReader.pm");
system_or_die("mv chap13/judiciousaligner.pm chap13/JudiciousAligner.pm");
system_or_die("mv chap13/qualityreviser.pm chap13/QualityReviser.pm");
system_or_die("mv chap13/seqreader.pm chap13/SeqReader.pm");
system_or_die("mv chap14/cspredictor.pm chap14/CSPredictor.pm");
system_or_die("mv chap16/integerinterval.pm chap16/IntegerInterval.pm");
system_or_die("mv chap4/util.pm chap4/Util.pm");
system_or_die("mv chap5/pambuilder.pm chap5/PAMBuilder.pm");
system_or_die("mv chap6/fastareader.pm chap6/FastaReader.pm");
system_or_die("mv chap6/genbankreader.pm chap6/GenBankReader.pm");
system_or_die("mv chap6/seqreader.pm chap6/SeqReader.pm");
system_or_die("mv chap6/simplereader.pm chap6/SimpleReader.pm");
system_or_die("mv chap8/recurrence.pm chap8/Recurrence.pm");
system_or_die("mv perllib/fastareader.pm lib/FastaReader.pm");
system_or_die("mv perllib/fastqreader.pm lib/FastqReader.pm");
system_or_die("mv perllib/genbankreader.pm lib/GenBankReader.pm");
system_or_die("mv perllib/seqreader.pm lib/SeqReader.pm");
system_or_die("mv perllib/set.pm lib/Set.pm");
system_or_die("mv perllib/simplereader.pm lib/SimpleReader.pm");
system_or_die("mv perllib/suffixtree.pm lib/SuffixTree.pm");
system_or_die("mv perllib/unionfind.pm lib/UnionFind.pm");
system_or_die("mv perllib/util.pm lib/Util.pm");


system_or_die("chmod go-rwx */*");
system_or_die("chmod u+rw */*");
system_or_die("chmod u-x */*");
system_or_die("chmod u+x */*.pl");
system_or_die("chmod u+x */*.pm");


sub system_or_die
{
    my ($cmd) = @_;
    !system($cmd) or die "couldn't $cmd ";
}
