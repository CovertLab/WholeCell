#!/usr/bin/perl
#Simulates high-throughput experiments, averages results, and stores averages in a new
#subdirectory of /path/to/WholeCell/output/.
#
#Inptus:
# 1) Simulation name
# 2) Path to .mat file containing parameter values
# 3) Number of cells to simulate
#
#Usage:
# 1) cd /path/to/WholeCell/
# 2) ./build.sh simulateHighthroughputExperiments
# 3) ./build.sh averageHighthroughputExperiments
# 4) ./runSimHtExpts.pl simName parameterValsPath.{mat|xml} nSim
#
#Author: Jonathan Karr, jkarr@stanford.edu
#Affiliation: Covert Lab, Department of Bioengineering, Stanford University
#Last updated: 6/2/2013

use Date::Format;
use Date::Language;
use HTML::Template;
use MIME::Lite;
use strict;
require 'library.pl';

#options
my %config = getConfiguration();

my $linuxUser = $config{'fileUserName'};
my $linuxRunUser = $config{'execUserName'};
my $storageServer = $config{'simulationHostName'};
my $nodeTmpDir = $config{'nodeTmpDir'};
my $pathToRunTime = $config{'mcrPath'};
my $emailAddress = $config{'email'};
my $baseDir = $config{'simulationPath'};

my $outDir = "$baseDir/output/simHtExpts";
my $lang = Date::Language->new('English');

my $simName = $ARGV[0];
my $parametersFile = $ARGV[1];
my $nSim = $ARGV[2];

my $jobFileName = '';

#make output directory
unless(-d $outDir){
    mkdir "$outDir" or die "Unable to make directory $outDir/$simName: $!";
}
mkdir "$outDir/$simName" or die "Unable to make directory $outDir/$simName: $!";
`cp "$parametersFile" "$outDir/$simName/parameters.mat"`;

#compile MATLAB code
#`./build.sh simulateHighthroughputExperiments`;
#`./build.sh averageHighthroughputExperiments`;

#copy executables
mkdir "$outDir/$simName/bin" or die "Unable to make directory $outDir/$simName/bin: $!";;
`cp -R bin/simulateHighthroughputExperiments $outDir/$simName/bin`;
`cp -R bin/averageHighthroughputExperiments $outDir/$simName/bin`;

#submit jobs for each condition
my $template = HTML::Template->new(filename => 'job.simHtExpts.sh.tmpl');
$template->param(simName => $simName);
$template->param(linuxRunUser => $linuxRunUser);
$template->param(emailAddress => $emailAddress);
$template->param(outDir => $outDir);
$template->param(storageServer => $storageServer);
$template->param(pathToRunTime => $pathToRunTime);
$template->param(nodeTmpDir => $nodeTmpDir);

my $submitJobs = '';
for (my $n = 1; $n <= $nSim; $n++){
	$jobFileName = sprintf("%s/%s/job.sim-%d.sh", $outDir, $simName, $n);
	open(FH, '>', $jobFileName) or die $!;
	$template->param(n => $n);
    print FH $template->output;
	close (FH);
	$submitJobs .= "sudo qsub $jobFileName; ";
}

#set permissions and run jobs
`sudo chmod -R 775 $outDir/$simName`;
`sudo chown -R $linuxUser:$linuxUser $outDir/$simName`;
`$submitJobs`;

#get job id
my $simulationId = `qstat | tail -n 1`;
$simulationId =~ /^(\d+)/;
my $simulationIdx = $1 - $nSim + 1;

#averaging job
my $afterany = '';
for (my $n = $simulationIdx; $n < $simulationIdx + $nSim; $n++){
	$afterany.=":$n";
}

$template = HTML::Template->new(filename => 'job.avgHtExpts.sh.tmpl');
$template->param(simName => $simName);
$template->param(linuxRunUser => $linuxRunUser);
$template->param(emailAddress => $emailAddress);
$template->param(outDir => $outDir);
$template->param(storageServer => $storageServer);
$template->param(pathToRunTime => $pathToRunTime);
$template->param(nodeTmpDir => $nodeTmpDir);
$template->param(afterany => $afterany);

$jobFileName = sprintf("%s/%s/job.avg.sh", $outDir, $simName);
open(FH, '>', $jobFileName) or die $!;
print FH $template->output;
close (FH);

#set permissions and submit job
`sudo chmod -R 775 $outDir/$simName`;
`sudo chown -R $linuxUser:$linuxUser $outDir/$simName`;
`sudo qsub $jobFileName`;

#print status message with total number of jobs submitted
print "$nSim simulations queued. Results will be stored at output/simHtExpts/$simName.\n";