package __TEaMPLATE__;

=head1 NAME


=head1 DESCRIPTION


=head1 SYNOPSIS


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

our $Some_Thing = '';

#-----------------------------------------------------------------------------#

=head2 Class Methods

=cut

sub Do_sth{
	my ($class) = @_;
}



#-----------------------------------------------------------------------------#

=head2 Constructor

=cut

sub new{
        my $proto = shift;
        my $self;
        my $class;
        
        # object method -> clone + overwrite
        if($class = ref $proto){ 
                return bless ({%$proto, @_}, $class);
		# NOTE: Deep copy required for complex objects, e.g.:
		#  use Storable qw(dclone);
		#    +
		#  return bless( {%{dclone($proto)}}, $class )
        }
	
	
	my $self = {
		# defaults
		
		
		# custom
		@_,
		
		# protected
		#id => $class->Auto_Id(),
	};

	bless $self;

	return $self;
}


#-----------------------------------------------------------------------------#

=head2 Object Methods

=cut

sub do_sth{
	my ($self) = @_;
}



#-----------------------------------------------------------------------------#

=head2 Object Accessors

Generic get/set methods to access object attributes.

=cut

=head2 id

Get/Set id. 

=cut

sub id{
	my ($self, $value, $force) = @_;
	if(defined $value || $force){
		$self->{id} = $value;
	}
	return $self->{id}
}


#-----------------------------------------------------------------------------#

=head2 Private Methods

=cut

sub _do_sth_private{

}







#-----------------------------------------------------------------------------#

=head1 AUTHOR

NAME S<MAIL>

=cut

1;


















