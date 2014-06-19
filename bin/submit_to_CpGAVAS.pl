#!/usr/bin/perl -w
use strict;
use WWW::Mechanize::Firefox;
use WWW::Mechanize::FormFiller;
use URI::URL;

use File::Spec;

# USAGE perl submit_to_CpGAVAS.pl <chloro.fa> <species>
my $usage = "USAGE perl submit_to_CpGAVAS.pl <chloro.fa> <species>";
# Returns: <Species>\t<File>\t<Weblink to results>

my $agent = WWW::Mechanize::Firefox->new( autocheck => 1 );
my $formfiller = WWW::Mechanize::FormFiller->new();

die "No species name given\n$usage\n" unless($ARGV[1]);

my ($species, $project, $file) = ($ARGV[1], "chloroExtractor", $ARGV[0]);

die "File $file not found\n$usage\n" unless(-e $file);

$file=File::Spec->rel2abs($file); # absolute path required

$agent->get('http://www.herbalgenomics.org/0506/cpgavas/analyzer/annotate');
$agent->form_number(1) if $agent->forms and scalar @{$agent->forms};
$agent->field('projectName' => $project);
$agent->field('kingdom' => 'plants', [], ['change']);
$agent->field('speciesName' => $species);
$agent->field('file1' => $file);
$agent->field('imageFormat' => 'png' );
$agent->field('emailAddress' => '' );
$agent->field('blastn_evalue' => '1e-10' );
$agent->field('blastx_evalue' => '1e-10' );
$agent->field('numOfHits' => '10' );
$agent->field('runRepeatMasker' => 'no' );
$agent->field('refdb_par1' => '0' );
$agent->field('tRNAscan_par1' => 'C' );
$agent->field('tRNAscan_par2' => 'O' );
$agent->field('tRNAscan_par3' => '15' );
$agent->field('tRNAscan_par4' => '116' );
$agent->field('codontable_par1' => '');
$agent->field('species' => 'all' );

#sleep 10000;
$agent->click_button(value => 'Submit');

print "$species\t$file\t",$agent->base(),"\n";

#print Dumper($agent); use Data::Dumper;

