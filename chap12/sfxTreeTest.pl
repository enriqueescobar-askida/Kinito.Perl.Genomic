#!/usr/bin/perl -I . -I ../perllib
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################

################################################################
## A small driver program to create and dump a few suffix trees.
################################################################

use strict;
use SuffixTree;

$|=1;

my $t = SuffixTree->new("banana");
$t->dump();

my $t = 
bless
{
    target=>'banana$',
    root=>
    {
            off=>0, len=>0,
            'a'=>{off=>1, len=>1,    
                'n'=>{off=>2, len=>2,
                        'n'=>{off=>4, len=>3},
                        '$'=>{off=>6, len=>1}
                    },
                '$'=>{off=>6, len=>1}
                },
            'b'=>{off=>0, len=>7},
            'n'=>{off=>2, len=>2,
                'n'=>{off=>4, len=>3},
                '$'=>{off=>6, len=>1}
                },
            '$'=>{off=>6, len=>1}
    }
},
"SuffixTree";
$t->dump();

my $t = SuffixTree->new("abracadabra");
$t->dump();
my $t = SuffixTree->new("bookkeeper");
$t->dump();
my $t = SuffixTree->new("<bookkeeper><abracadabra>");
$t->dump();
my $t = SuffixTree->new("baaaaaaaaaaaaaaa");
$t->dump();
my $t = SuffixTree->new("misssisssisssisssippippi");
$t->dump();







