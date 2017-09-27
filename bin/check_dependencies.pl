#!/usr/bin/env perl
use strict;
use warnings;

use File::Which;

my %modules2check = (

    );

my %programs2check = (

    );

foreach my $module (keys %modules2check)
{
    eval { require $module };

    # assuming no error
    my ($error, $msg) = (0, "");

    if ($@)
    {
	($error, $msg) = (1, $@);
    }

    $modules2check{$module}{error} = $error;
    $modules2check{$module}{errormsg} = $msg;
}

foreach my $program (keys %programs2check)
{
    my $path = which($program);

    # assuming no error
    my $error = 0;

    unless (defined $path)
    {
	$error = 1;
    }

    $programs2check{$program}{error} = $error;
    $programs2check{$program}{path} = $path;
}

=head1 NAME

check_dependencies.pl

=head1 DESCRIPTION

Checks all dependencies for chloroExtactor. Returns without error, if
all dependencies are fullfilled. Will list missing dependencies
otherwise.

=head1 SYNOPSIS

  $ ./check_dependencies.pl

=head1 CHANGELOG

see git log.

=head1 CODE

=cut
