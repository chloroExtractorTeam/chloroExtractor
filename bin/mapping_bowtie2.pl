#!/usr/bin/env perl

=head1 NAME


=head1 DESCRIPTION


=head1 SYNOPSIS


=head1 OPTIONS

=over

=item -s|--sam

SAM file. REQUIRED.

=item --create-config

Create a config file with default settings for user customization.

=item -c|--config

Use user customized config file. Superseeds default config.

=item -V|--version

Display version.

=item -h|--help

Display this help.

=back

=head1 CHANGELOG

see git log.

=head1 TODO

=head1 CODE

=cut

#-----------------------------------------------------------------------------#
# Modules

# core
use strict;
use warnings;
no warnings 'qw';

use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);
use Pod::Usage;
use Log::Log4perl qw(:no_extra_logdie_message);
use Log::Log4perl::Level;

use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

use FindBin qw($RealBin);
use lib "$RealBin/../lib/";

use File::Basename;
use File::Copy;

# additional modules
use Cfg;

use Sam::Parser;
use Sam::Seq;
use Sam::Alignment;

#enable in script mapping
use Bowtie2;


#-----------------------------------------------------------------------------#
# Globals

our $VERSION = 0.01;

# get a logger
my $L = Log::Log4perl::get_logger();
Log::Log4perl->init( \q(
	log4perl.rootLogger                     = INFO, Screen
	log4perl.appender.Screen                = Log::Log4perl::Appender::Screen
	log4perl.appender.Screen.stderr         = 1
	log4perl.appender.Screen.layout         = PatternLayout
	log4perl.appender.Screen.layout.ConversionPattern = [%d{yy-MM-dd HH:mm:ss}] [%C] %m%n
));


#-----------------------------------------------------------------------------#
# Config

# core
my $core_cfg = "$RealBin/../chloroExtractor.cfg";
my %opt = Cfg->Read_Cfg($core_cfg); 

# user defaults and overwrite core
my $user_cfg;
for(my $i=0; $i<@ARGV; $i++){
        if($ARGV[$i] =~ /-c$|--config$/){
                $user_cfg = $ARGV[$i+1];
                last;
        }
}

%opt = (%opt, Cfg->Read_Cfg($user_cfg)) if $user_cfg; # simple overwrite


#TODO: custom config


#-----------------------------------------------------------------------------#
# GetOptions

GetOptions( # use %opt (Cfg) as defaults
	\%opt, qw(
		version|V!
		debug|D!
		help|h!
		config|c=s
		create_config|create-config
                input_reads|1=s
                input_mates|2=s
                target_coverage|coverage=i
	)
) or $L->logcroak('Failed to "GetOptions"');

# help
$opt{help} && pod2usage(1);

# version
if($opt{version}){
	print "$VERSION\n"; 
	exit 0;
}

# create template for user cfg
if($opt{create_config}){
	pod2usage(-msg => 'To many arguments', -exitval=>1) if @ARGV > 1;
	my $user_cfg = @ARGV ? $ARGV[0] : basename($core_cfg);
	copy($core_cfg, $user_cfg) or $L->logdie("Creatring config failed: $!");
	$L->info("Created config file: $user_cfg");
	exit 0;
}


# debug level
$L->level($DEBUG) if $opt{debug};
$L->debug('Verbose level set to DEBUG');

$L->debug(Dumper(\%opt));



#-----------------------------------------------------------------------------#
# MAIN


my $bowtie2 = Bowtie2->new(
  path => $opt{bowtie2_path},
  log => $opt{bowtie2_log},

);

$bowtie2->bowtie2(
  %{$opt{bowtie2_params}},
  "-U" => $opt{input_reads},
  "-x" => $opt{bowtie2_DB}
);

# WARNING: Bowtie required "-U" twice for paired; not possible in option hash
# Therefore: if input mates option exists, set singleorpair to 2 to map only reads
# and multiply current coverage by two

my $singleorpair;
if ($opt{input_mates}) {
  $singleorpair = 2;
}
else {
  $singleorpair = 1;
}


# read output on the fly
my $sp = Sam::Parser->new(
  fh => $bowtie2->stdout
);

my %h; 
my $ss;

my %seqs;

while(%h = $sp->next_header_line('@SQ')){
  $seqs{$h{'SN'}} = Sam::Seq->new(
    id => $h{'SN'},
    len => $h{'LN'},
  );  
  print '@SQ'."\tSN:$h{'SN'}\tLN:$h{'LN'}\n" if ($opt{debug});
}

my $current_cov;
my $last_id;

my $c; 
while(my $aln = $sp->next_aln){
  my $id = $aln->rname();
  $seqs{$id}->add_aln($aln);
  print "$aln" if ($opt{debug});
  $c++;
  unless ($c % $opt{coverage_check_interval}){
    $current_cov = estimate_coverage(\%seqs)*$singleorpair;
    $L->debug("Coverage: ",$current_cov);
    $last_id = $aln->qname;
    if ($current_cov >= $opt{target_coverage}) { 
      last;
    }
  }
}

$bowtie2->cancel("Coverage = $current_cov\n");

$last_id =~ s?/[12]$??;

my $sed_cmd1 = 'sed "1~4 {/^\@'.$last_id.'\(\/[12]\)*\s/{N;N;N;q}}" '.$opt{input_reads}." >".$opt{reads_out};

$L->debug($sed_cmd1);
qx($sed_cmd1);

if ($opt{input_mates}) {

  my $sed_cmd2 = 'sed "1~4 {/^\@'.$last_id.'\(\/[12]\)*\s/{N;N;N;q}}" '.$opt{input_mates}." >".$opt{mates_out};

  $L->debug($sed_cmd2);
  qx($sed_cmd2);

}










# ----- SUBS ------ #



sub pairwise_sum {
  my @array1 = @{$_[0]};
  my @array2 = @{$_[1]};

  my @len;
  push @len, scalar @array1;
  push @len, scalar @array2;

  my @sort = sort {$a <=> $b} @len;

  my @added;

  for (my $i=0;$i<=$sort[-1];$i++) {

    if (exists $array1[$i] && exists $array2[$i]) {
      $added[$i] = $array1[$i] + $array2[$i];
    }
    elsif (exists $array1[$i]) {
      $added[$i] = $array1[$i];
    }
    elsif (exists $array2[$i]) {
      $added[$i] = $array2[$i];
    }
  }
  return @added;
}

sub median {
  my @array = @{$_[0]};
  my $median = (sort{$a<=>$b}@array)[@array/2];
  return $median;
}


sub estimate_coverage {

  # get coverages for complete sam file
  # and build hash with protein name as key
  # and protein wise sum of coverage in
  # array

  my %protein_wise_coverage;
  for my $ss(values %{$_[0]}){
    my @covs=$ss->coverage();
    my $protein = (split /_/, $ss->id())[-1];
    @{$protein_wise_coverage{$protein}} = pairwise_sum(\@{$protein_wise_coverage{$protein}}, \@covs);
  }



  # print proteinwise coverage // omitting 
  # zeros and empty proteins // also trim
  # both ends

  my %medians;
  my @median_array;

  for (keys%protein_wise_coverage) {
    @{$protein_wise_coverage{$_}} = grep { $_ != 0 } @{$protein_wise_coverage{$_}};
    if (@{$protein_wise_coverage{$_}}-100>200) {
      @{$protein_wise_coverage{$_}} = @{$protein_wise_coverage{$_}}[100..@{$protein_wise_coverage{$_}}-100];
    }
    else {
      @{$protein_wise_coverage{$_}} = ();
    }
    if (@{$protein_wise_coverage{$_}}) {
      my $prot_median = median(\@{$protein_wise_coverage{$_}});
      $medians{$_} = $prot_median;
      push @median_array, $prot_median if ($prot_median > $opt{min_single_coverage});
    }
  }

  $L->debug("Coverage Array ","@median_array"," -- ",scalar @median_array);

  if (@median_array > $opt{min_covered_CDS}) {
    return median(\@median_array);
  }
  else {
    return -1;
  }

}



#-----------------------------------------------------------------------------#

=head1 AUTHOR

Clemens Weiss S<clemens.weiss@stud-mail.uni-wuerzburg.de>

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



