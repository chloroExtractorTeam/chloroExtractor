#!/usr/bin/perl -w
use strict;
use WWW::Mechanize;
use WWW::Mechanize::FormFiller;
use URI::URL;

# USAGE perl submit_to_CpGAVAS.pl <chloro.fa> <species>
my $usage = "USAGE perl submit_to_CpGAVAS.pl <chloro.fa> <species>";
# Returns: <Species>\t<File>\t<Weblink to results>

my $agent = WWW::Mechanize->new( autocheck => 1 );
my $formfiller = WWW::Mechanize::FormFiller->new();
$agent->env_proxy();

die "No species name given\n$usage\n" unless($ARGV[1]);

my ($species, $project, $file) = ($ARGV[1], "chloroExtractor", $ARGV[0]);

die "File $file not found\n$usage\n" unless(-e $file);

  $agent->get('http://www.herbalgenomics.org/0506/cpgavas/analyzer/annotate');
  $agent->form_number(1) if $agent->forms and scalar @{$agent->forms};
  $formfiller->add_filler( 'projectName' => Fixed => $project );
  $formfiller->add_filler( 'kingdom' => Fixed => 'plants' );
  $formfiller->add_filler( 'speciesName' => Fixed => $species );
  $formfiller->add_filler( 'file1' => Fixed => $file );
  $formfiller->add_filler( 'imageFormat' => Fixed => 'png' );
  $formfiller->add_filler( 'emailAddress' => Fixed => '' );
  $formfiller->add_filler( 'blastn_evalue' => Fixed => '1e-10' );
  $formfiller->add_filler( 'blastx_evalue' => Fixed => '1e-10' );
  $formfiller->add_filler( 'numOfHits' => Fixed => '10' );
  $formfiller->add_filler( 'runRepeatMasker' => Fixed => 'no' );
  $formfiller->add_filler( 'refdb_par1' => Fixed => '0' );
  $formfiller->add_filler( 'species' => Fixed => 'all' );
  $formfiller->add_filler( 'tRNAscan_par1' => Fixed => 'C' );
  $formfiller->add_filler( 'tRNAscan_par2' => Fixed => 'O' );
  $formfiller->add_filler( 'tRNAscan_par3' => Fixed => '15' );
  $formfiller->add_filler( 'tRNAscan_par4' => Fixed => '116' );
  $formfiller->add_filler( 'codontable_par1' => Fixed => '' );
  $formfiller->add_filler( 'Clear' => Fixed => '' );
  $formfiller->fill_form($agent->current_form);
  $agent->submit();

print "$species\t$file\t",$agent->base(),"\n";

#print Dumper($agent); use Data::Dumper;

