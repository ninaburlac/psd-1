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
my $tierfile = $ARGV[2];
my $dplms = 0;
if ( $len > 3 ) { $dplms = $ARGV[3];}
my $weight = 0;
if ( $len > 4 ) { $weight = $ARGV[4];}

my $myUser = $ENV{'USER'};

# save current umask
my $old_umask = umask;
umask 0000;

my $resdir;
if ( $dplms ) {$resdir = $localdir . "/analysis_dplms_" . $weight;}
else { $resdir = $localdir . "/analysis_standard";}
mkdir "$resdir", 0770 unless -d "$resdir";

my @cmd = ( "!psdcommand!");


for(my $chn = $firstchn; $chn < $lastchn; $chn++){
    #my $thisdir = $localdir . "/chn" . $chn;
    my $thisdir = $resdir . "/chn" . $chn;
    mkdir "$thisdir", 0770 unless -d "$thisdir";
    
    my $psdprogram = $localdir . "/processPSD";
    my $psdout = $thisdir . "/psd_chn" . $chn . ".out";
    my $psdcommand = $psdprogram . " " . $thisdir . " " . $tierfile . " " . $chn . " " . $dplms . " " . $weight . " > " . $psdout; 
    
    ### submit noise script
    
    my $scripttobecopied = $localdir . "/script.template.sh";
    open IN, $scripttobecopied or die "Can't read source file $scripttobecopied: $!\n";
    
    my $currscript = $thisdir . "/script_chn" . $chn . ".sh";
    open OUT, ">$currscript" or die "Can't write on file $currscript: $!\n";
    
    while (<IN>) {
	s/$cmd[0]/$psdcommand/g;
	print OUT $_;
    }
    close IN;
    close OUT;
    
    ### submit PSD job to queue
    my $scriptlog = $thisdir . "/psd_" . $chn . ".out";
    my $scripterr = $thisdir . "/psd_" . $chn . ".err";
    my $jobName;
    if ( $dplms ) {$jobName = "pd_" . $chn;}
    else {$jobName = "ps_" . $chn;}
    
    my $QUEUEcmd = "qsub -N " . $jobName . " -q gerda -V -d " . $thisdir . " -m abe -e localhost:". $scripterr . " -o localhost:" . $scriptlog . " -l mem=4000mb " . $currscript;
    system($QUEUEcmd);
    print $QUEUEcmd . "\n";
    
}

my $countJob;
if ( $dplms ) { $countJob = "qstat -u " . $myUser . " |  grep pd_ | wc -l ";}
else { $countJob = "qstat -u " . $myUser . " |  grep ps_ | wc -l ";}

print $countJob . "\n";
my $inQueue = 1;

while ($inQueue) {
    sleep 30;
    my $actualJob = `$countJob`;
    print $actualJob;
    if($actualJob==0) {
	$inQueue=0;
	
	my $resultsfile;
	if ( $dplms ) { $resultsfile = $resdir . "/PSA_results_dplms.txt";}
	else { $resultsfile = $resdir . "/PSA_results_stand.txt";}
	my $catCmd = "cat ";
	for(my $chn = $firstchn; $chn < $lastchn; $chn++){
	    if ( $chn == 8 || $chn == 9 || $chn == 10 || $chn == 27 || $chn == 28 || $chn == 29 ) {next;}
	    my $chnfile = $resdir . "/chn" . $chn . "/PSA_results.txt";
	    if (-e $chnfile){$catCmd = $catCmd . $chnfile . " ";}
	}
	$catCmd = $catCmd . " > " . $resultsfile;
	print $catCmd;
	system($catCmd);
	
	my $text;
	if ( $dplms ) { $text = "DPLMS-" . $weight;}
	else { $text = "standard";}
	my $plotcmd = $localdir . "/plotSF_AoE " . $resultsfile . " " . $resdir . " " . $text . " > " . $resdir . "/plot.out"; 
	print $plotcmd;
	system($plotcmd);
	
        # Reset umask to old value
	umask $old_umask;
	
	my $from = 'psdtest';
	my $to = 'valerio.dandrea\@lngs.infn.it';
	my $subject = "psd";
	#my $fwhmFile = $scandir . "/FWHM_chn" . $chn . ".txt";
	my $body = "Hi, \n this is an automatic message from the test of the PSD";
	my $sendMail = "echo \"". $body . "\" | mailx -s \"". $subject . "\" -r " . $from . " " . $to;
	system($sendMail);
    }
}


