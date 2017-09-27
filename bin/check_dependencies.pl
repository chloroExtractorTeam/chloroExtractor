#!/usr/bin/env perl
use strict;
use warnings;

eval { require File::Which; };

my $problems_with_file_which = $@;

my $missing_which = ! -e "which";

if ($problems_with_file_which && $missing_which)
{
    die "Module 'File::Which' or the little program 'which' are required for determination of executable files. Please install it\n";
}

unless ($missing_which)
{
    warn "Helper program 'which' was found and will be used instead of perl module 'File::Which'\n";
}

sub dep_which
{
    my ($prog) = @_;
    my $val = undef;

    unless ($problems_with_file_which)
    {
	$val = File::Which::which($prog);
    } else {
	my $cmd = "which $prog";
	$val = qx($cmd);
	chomp($val);
    }

    return $val;
}

my %modules2check = (

    'Moose' => {},
    'Log::Log4Perl' => {},
    'Term::ProgressBar' => {},
    'Graph' => {},
    'IPC::Run' => {},
    'File::Which' => {}

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
    my $path = dep_which($program);

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
