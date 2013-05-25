#!/usr/bin/perl
#  Usage: ./runDataExport.pl simBatch nSimulations
#
#Author: Jonathan Karr, jkarr@stanford.edu
#Affiliation: Covert Lab, Department of Bioengineering, Stanford University
#Last updated: 8/22/2012

use HTML::Template;
use strict;
require 'library.pl';

#options
my $conditionSetTimeStamp = $ARGV[0];
my $nSimulations = $ARGV[1];

my %config = getConfiguration();
my $linuxUser = $config{'fileUserName'};
my $pathToRunTime = $config{'mcrPath'};
my $emailAddress = $config{'email'};
my $baseDir = $config{'simulationPath'};
my $outDir = "$baseDir/output/runSimulation";

#compile runDataExport project
`./build.sh runDataExport`;

#data export job
my $submitJobs = '';
my $jobFileName = '';
my $template = HTML::Template->new(filename => 'job.dataexport.sh.tmpl');
$template->param(baseDir => $baseDir);
$template->param(conditionSetTimeStamp => $conditionSetTimeStamp);
$template->param(linuxUser => $linuxUser);
$template->param(emailAddress => $emailAddress);
$template->param(outDir => $outDir);
$template->param(pathToRunTime => $pathToRunTime);

for (my $i = 1; $i <= $nSimulations; $i++) {
	my $simOutDir = sprintf("%s/%s/%d/json", $outDir, $conditionSetTimeStamp, $i);
	`mkdir -p $simOutDir`;
	
	$jobFileName = sprintf("%s/job.dataexport.sh", $simOutDir);
	$template->param(iSim => $i);
	open(FH, '>', $jobFileName) or die $!;
	print FH $template->output;
	close (FH);
	
	`chmod 775 $simOutDir`;
	`chmod 775 $jobFileName`;
	`sudo chown $linuxUser:$linuxUser $simOutDir`;
	`sudo chown $linuxUser:$linuxUser $jobFileName`;
	`sudo qsub $jobFileName`;
}

#print status message with total number of jobs submitted
print "Data export jobs queued.\n";
