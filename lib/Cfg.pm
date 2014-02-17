package Cfg;

=head1 NAME

Cfg.pm

=head1 DESCRIPTION

Read a simple config in perl hash style from file. Inspired from
 L<http://www.perlmonks.org/?node_id=464358>.

=head1 SYNOPSIS

  use Cfg;
  
  Cfg->Read_Cfg("path/to/my.cfg") || die "Couldn't read *.cfg";
  
  my $log_file = $Cfg::Cfg{log_file};

=head1 CHANGELOG

see git log.

=head1 TODO

=cut

#-----------------------------------------------------------------------------#

=head1 Include

=cut


use warnings;
no warnings 'qw';
use strict;

use Carp;
use Log::Log4perl qw(:easy :no_extra_logdie_message);

#use Data::Dumper;
#$Data::Dumper::Sortkeys = 1;



#-----------------------------------------------------------------------------#

=head2 Globals

=cut


my $L = Log::Log4perl::get_logger();


#-----------------------------------------------------------------------------#

=head2 Class Attributes

=cut

our $VERSION = 0.01;

our %Cfg;

#-----------------------------------------------------------------------------#

=head2 Class Methods

Read config from simple perl syntax config file (LIST context). Returns config
 , dies on error.

Simple *.cfg format:

  #--------------------------------#
  log => "path/to/somewhere",
  
  servers => [
    "foo", "bar"
  ],
  
  more => {
    complex => [qw( s t u f f )]
  },
  #--------------------------------#

=cut

sub Read_Cfg{
	my ($class,$cfg_file) = @_;

	unless(-f $cfg_file){
		$L->logdie("Cannot find config file '$cfg_file'");
	}

	%Cfg = do($cfg_file);
	
        if ($@) {
            $L->logdie("Failed to read config '$cfg_file' - $@");
        }

	return %Cfg;
}

=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut

1;


















