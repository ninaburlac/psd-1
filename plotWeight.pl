#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use YAML::Tiny 'LoadFile';
use YAML::Tiny;
use Data::Dumper;
use POSIX qw(ceil floor);

use constant false => 0;
use constant true  => 1;

my $localdir = $ENV{PWD};

my $len = @ARGV;

#Open the configuration file
#my $configfile = YAML::Tiny::LoadFile( $localdir . "/config.yml");

my $firstchn = $ARGV[0];#$configfile->{FirstChannel};
my $nchn = $ARGV[1];#$configfile->{NChannels};
my $lastchn = $firstchn + $nchn;

my $myUser = $ENV{'USER'};

# save current umask
my $old_umask = umask;
umask 0000;

my $resdir = $localdir . "/plot_weight";
mkdir "$resdir", 0770 unless -d "$resdir";

for(my $chn = $firstchn; $chn < $lastchn; $chn++){
    if ( $chn == 8 || $chn == 9 || $chn == 10 || $chn == 27 || $chn == 28 || $chn == 29 ) {next;}
    my $thisdir = $resdir . "/chn" . $chn;
    mkdir "$thisdir", 0770 unless -d "$thisdir";
    
    my $chnfile = $thisdir . "/PSA_results_chn" . $chn . ".txt";
    my $catCmd = "cat " . $localdir . "/analysis_*/chn" . $chn . "/PSA_results.txt > " . $chnfile;
    print $catCmd."\n";
    system($catCmd);
    
    my $plotcmd = $localdir . "/plotSF_AoE_weight " . $chnfile . " " . $thisdir . " " . $chn . " > " . $thisdir . "/plot.out"; 
    print $plotcmd."\n";
    system($plotcmd);
}
