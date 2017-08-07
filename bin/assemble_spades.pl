#!/usr/bin/env perl

use warnings;
use strict;

my $i;
for($i=0; $i<@ARGV; $i++){
    last if $ARGV[$i] eq "-o";
}
my $pre = $ARGV[$i+1];
open(LOG, ">", $pre.".log") or die "$!: ",$pre,".log\n";

open(SPADES , "spades.py --only-assembler @ARGV |") or die "$!: Running SPAdes\n";
print LOG $_ while <SPADES>;
close SPADES;
close LOG;


my $file_target = $pre."/scaffolds.fasta";
my $link_name = $pre.".fa";
my $symlink_exists = eval { symlink($file_target, $link_name); 1 };
unless ($symlink_exists)
{
    die "Something went wrong during creation of the file link: $file_target --> $link_name\n";
}
