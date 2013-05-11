#!/usr/bin/perl

use HTML::Template;
use strict;
require 'library.pl';

#options
my %config = getConfiguration();
my $linuxUser = $config{'fileUserName'};
my $pathToRunTime = $config{'mcrPath'};
my $emailAddress = $config{'email'};
my $baseDir = $config{'simulationPath'};

my $outDir = "$baseDir/documentation/paper/figures/figure5";
my $conditionSetTimeStamp = '2011_10_25_01_45_44';
my $nJobs = 128;

#output directory
`mkdir -p "$outDir"`;
`sudo chmod -R 775 "$outDir"`;
`sudo chown -R $linuxUser:$linuxUser "$outDir"`;

#compile project
#`./build.sh runDisplacementCalculation`;

#analysis job
my $template = HTML::Template->new(filename => 'job.displacementCalculation.sh.tmpl');
$template->param(baseDir => $baseDir);
$template->param(conditionSetTimeStamp => $conditionSetTimeStamp);
$template->param(outDir => $outDir);
$template->param(linuxUser => $linuxUser);
$template->param(emailAddress => $emailAddress);
$template->param(pathToRunTime => $pathToRunTime);

#frame rendering
for (my $i = 1; $i <= $nJobs; $i++) {
	my $jobFileName = sprintf("%s/job.displacementCalculation_%s_%d.sh", $outDir, $conditionSetTimeStamp, $i);
	$template->param(iJob => $i);
	open(FH, '>', $jobFileName) or die $!;
	print FH $template->output;
	close (FH);
	
	`sudo chmod 775 $jobFileName`;
	`sudo chown -R $linuxUser:$linuxUser $jobFileName`;
	`sudo qsub $jobFileName`;
}

#print status message with total number of jobs submitted
print "Jobs queued.\n";