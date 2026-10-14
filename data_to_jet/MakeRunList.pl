#!/usr/bin/perl

use strict;
use warnings;

# StarFileList.pm lives in another user's area; update if it moves
use lib '/star/u/dkap7827/Tools2/Tools/PerlScripts';
use StarFileList;
use Cwd qw(abs_path);
use Getopt::Long qw(GetOptions);
use Getopt::Long qw(HelpMessage);
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


# Runs excluded for data quality:
#   23003063-23004033: FEE low-voltage problem (below the wave-1 start)
#   23007067-23007076: beam and luminosity tests
#   23027046-23027048: 56-bunch and calibration runs
# ECAL is calibrated from run 23005043 (LED reference) onward, hence the
# wave-1 start at 23005044.
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

# Shift-crew good-run flag from RunLog.runStatus, queried live, with STAR
# production's rule (RunDAQ.pm): rtsStatus = 0 and shiftLeaderStatus <= 1.
# An unreachable DB is fatal unless --skip-status-check is given, in which
# case the list is made without this filter (and says so).
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
$FileList->AddExtraOpts("-limit $LIMIT");

# Wave 1 only (23005044-23027047); wave 2 is 23027048-23048036
$FileList->Runnumber("[]23005044-23027047");

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

