#!/usr/bin/perl

use strict;
use warnings;

use lib '/star/u/dkap7827/Tools2/Tools/PerlScripts';
use StarFileList;
use Cwd qw(abs_path);
use Getopt::Long qw(GetOptions);
use Getopt::Long qw(HelpMessage);
#use Getopt::Long qw(VesrionMessage);
use Pod::Usage;

=pod

=head1 NAME

MakeRunList - Handy script for generating the run list to use for Fcs Run 22 Analysis. For this reason certain serach conditions are fixed

=head1 SYNOPSIS

MakeRunList.pl

  --keys, -k     Keys to return from search (default path,filename,events)
  --file, -f     File to write output to (default is 'test.list')
  --limit, -l    Number of files to return (default is 10)
  --print, -p    Just print to terminal the commands or file list, will not create a file  (1=command, 2=command and list)
  --verbose, -v  Set the printout level (default is 0)
  --skip-status-check  Proceed even if RunLog.runStatus is unreachable, skipping
                        the shift-crew bad-run filter instead of dying. Default
                        is to die on an unreachable DB, since that filter is a
                        correctness gate, not a best-effort lookup.
  --help, -h     Print this help

=head1 VERSION

0.3

=head1 DESCRIPTION

This program will generate a run list of MuDst files to use for Run 22 analysis. It will use "StarFileList" to call the neccessary methods to do this. Certain search conditions will be fixed since this makes it easier to 

=cut

=begin comment
@[December 15, 2023](David Kapukchyan)
> Added some more command line options and cleaned up a little bit (v0.2)

@[December 7, 2023](David Kapukchyan)
> First instance (v0.1)
=cut


my $KEYS = "path,filename,events";
#my $CONDITIONS = "";
my $FILE = "test.list";
my $LIMIT = 10;
my $PRINTLEVEL = 0;
my $VERBOSE = 0;
my $SKIP_STATUS_CHECK = 0;

GetOptions(
    'keys|k=s'    => \$KEYS,
    'file|f=s'    => \$FILE,
    'limit|l=i'   => \$LIMIT,
    'print|p=i'   => \$PRINTLEVEL,
    'verbose|v=i' => \$VERBOSE,
    'skip-status-check' => \$SKIP_STATUS_CHECK,
    'help|h'      => sub { HelpMessage() }
    ) or HelpMessage(1);


# Runs excluded from the wave-1 window below for known data-quality
# reasons (per the 2026-10 investigation writeup):
#   23003063-23004033: FEE low-voltage problem (below the wave-1 start,
#                       kept here only so this list stays correct if the
#                       start is ever widened back below 23005044)
#   23007067-23007076: beam and luminosity tests
#   23027046-23027048: 56-bunch and calibration runs
# ECAL is only calibrated from run 23005043 (the LED reference run)
# onward, hence the wave-1 start at 23005044.
my @EXCLUDED_RANGES = (
    [23003063, 23004033],
    [23007067, 23007076],
    [23027046, 23027048],
);

sub IsExcludedRun
{
    my $run = shift;
    foreach my $range (@EXCLUDED_RANGES){
        return 1 if ($run >= $range->[0] && $run <= $range->[1]);
    }
    return 0;
}

# RunLog.runStatus's shift-crew assessment, queried live so this covers
# whatever run range is actually requested (not just wave 1) without
# needing a hardcoded list maintained by hand. Good-run rule matches
# STAR's own production tooling (packages/scripts/RunDAQ.pm):
# rtsStatus=0 AND shiftLeaderStatus<=1. Checked directly against wave 1
# (23005044-23027047): 517/519 runs pass, only 23023038/23023039 fail
# (both already low-file-count/incomplete runs, not otherwise flagged) --
# this is a genuinely separate signal from the EXCLUDED_RANGES above, not
# a duplicate of it.
#
# Unlike RunNumberToDateTime() elsewhere in this project (which falls
# back silently to day-of-year arithmetic), an unreachable DB here is
# fatal by default: this filter is a correctness gate (dropping runs the
# shift crew marked bad), not a best-effort date lookup, so silently
# proceeding without it would risk quietly shipping a run list that
# includes known-bad runs. Pass --skip-status-check to proceed anyway
# (e.g. from a worker node with no DB access) -- the run list will then
# NOT have this filter applied, which is loudly reported either way.
sub GetBadRunsFromRunLog
{
    my ($runMin, $runMax, $skipCheck) = @_;
    my %badRuns = ();
    my $cmd = "mysql -h db04.star.bnl.gov --port=3421 --connect-timeout=10 -N -s -e " .
        "\"SELECT runNumber, rtsStatus, shiftLeaderStatus FROM runStatus " .
        "WHERE runNumber BETWEEN $runMin AND $runMax AND deactive=0 AND flavor='ofl'\" RunLog 2>/dev/null";
    my @lines = `$cmd`;
    if( $? != 0 || !@lines ){
        if( $skipCheck ){
            print "WARNING: MakeRunList.pl::GetBadRunsFromRunLog - could not query RunLog.runStatus " .
                  "(unreachable?) -- --skip-status-check was given, proceeding WITHOUT the " .
                  "shift-crew run-status filter. No runs dropped by this check (it did not run).\n";
            return (\%badRuns, 1);
        }
        else{
            die "ERROR: MakeRunList.pl::GetBadRunsFromRunLog - could not query RunLog.runStatus " .
                "(unreachable?) -- refusing to generate a run list without the shift-crew " .
                "run-status filter. Pass --skip-status-check to proceed anyway.\n";
        }
    }
    foreach my $line (@lines){
        chomp $line;
        my ($run, $rts, $sl) = split(/\t/, $line);
        next unless defined $run;
        if( !($rts == 0 && $sl <= 1) ){ $badRuns{$run} = 1; }
    }
    return (\%badRuns, 0);
}

my $FileList = new StarFileList();
$FileList->Keys($KEYS);
$FileList->Storage("=local");
$FileList->Production("=P24ia");
$FileList->Trgsetupname("=production_pp500_2022");
$FileList->Filetype("=daq_reco_mudst");
$FileList->Sname2("=st_fwd");
#$FileList->Runnumber(">23005043");
#$FileList->Runnumber("=23007011");
$FileList->AddExtraOpts("-limit $LIMIT");

# Wave 1 only (23005044-23027047); wave 2 (23027048-23048036) needs its
# own deliberate change here once that production is actually run.
$FileList->Runnumber("[]23005044-23027047");

#Get Time
#my $epochtime = time();            #UNIX time (seconds from Jan 1, 1970)
#my $localtime = localtime();       #Human readaable time
if( $PRINTLEVEL>2 ){ die "ERROR:Invalid print level:'$PRINTLEVEL'\n"; }
if($PRINTLEVEL==1 ){  $FileList->Command(1); exit; }

if( $VERBOSE>=1 ){
    print "Keys:",$FileList->Keys(),"\n";
    print "Conditions:",$FileList->Conditions(),"\n";
}

my $FileArrRef = $FileList->RetrieveModFiles();

my ($BadRunLogRunsRef, $StatusCheckSkipped) = GetBadRunsFromRunLog(23005044, 23027047, $SKIP_STATUS_CHECK);
my %BadRunLogRuns = %{$BadRunLogRunsRef};
if( $StatusCheckSkipped ){
    print "RunLog shift-crew bad runs dropped: check skipped (DB unreachable, --skip-status-check given)\n";
}
elsif( %BadRunLogRuns ){
    print "RunLog shift-crew bad runs dropped (", scalar(keys %BadRunLogRuns), "): ",
          join(",", sort keys %BadRunLogRuns), "\n";
}
else{
    print "RunLog shift-crew bad runs dropped: none\n";
}

# Post-filter the excluded sub-ranges above, plus RunLog's shift-crew bad
# runs -- StarFileList's Runnumber() only holds one range/condition at a
# time, so neither is expressible in the catalog query itself.
my @KeptFiles = ();
foreach my $file ( @{$FileArrRef} ){
    if( $file =~ m/(\d{8})_\w*_\d{7}\.MuDst\.root/ ){
        my $run = $1;
        next if IsExcludedRun($run);
        next if exists $BadRunLogRuns{$run};
    }
    push @KeptFiles, $file;
}
$FileList->Files(\@KeptFiles);
$FileArrRef = $FileList->Files();

$FileList->SortByRunNumber();

if($PRINTLEVEL==2 ){ $FileList->PrintFiles(); exit; }

my $filesize = @{$FileArrRef};
if( $VERBOSE>=1 ){ $FileList->PrintFiles(); }

$FileList->WriteFiles($FILE);

