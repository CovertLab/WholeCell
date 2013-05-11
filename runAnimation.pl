#!/usr/bin/perl
#
#mv /home/projects/WholeCell/simulation/output/runSimulation/2011_10_03_17_17_32/movie /home/projects/WholeCell/simulation/output/runSimulation/2011_10_03_17_17_32/movie1
#./runAnimation.pl 2011_10_03_17_17_32
#cd /home/projects/WholeCell/simulation/output/runSimulation/2011_10_03_17_17_32/movie
#mencoder "mf://frame_%d.png" -mf type=png:fps=30.000000:w=2222:h=1666 -vf scale=1200:900 -o video4.avi -ovc lavc -lavcopts vcodec=mpeg4:vbitrate=1200:vpass=1:aspect=4/3:mbd=2:turbo -info name="Whole-Cell Animation":subject="Animation of Mycoplasma genitalium whole-cell simulations.":artist="Jonathan Karr, Covert Lab, Department of Bioengineering, Stanford University":copyright="Copyright Jonathan Karr, Covert Lab, Department of Bioengineering, Stanford University 2012" 
#mencoder "mf://frame_%d.png" -mf type=png:fps=30.000000:w=800:h=600 -vf scale=800:600 -o video_h264.avi -ovc x264 -x264encopts pass=2:turbo:bitrate=400:frameref=5:bframes=1: -info name="Whole-Cell Animation":subject="Animation of Mycoplasma genitalium whole-cell simulations.":artist="Jonathan Karr, Covert Lab, Department of Bioengineering, Stanford University":copyright="Copyright Jonathan Karr, Covert Lab, Department of Bioengineering, Stanford University 2012" 
#ffmpeg -i video4.avi -vcodec wmv2 -r 30 -b 5000k -s 1200x900 -an -y video4.wmv
#ffmpeg -f image2 -i "frame_%d.png" -vcodec wmv2 -r 30 -qscale 4 -pass 1 -s 1200x900 -an -y video1.wmv
#ffmpeg -f image2 -i "frame_%d.png" -vcodec libx264 -r 30 -b 1500k -vpre medium -s 800x600 -an -y movie.mp4
#mencoder "mf://frame_%d.png" -mf type=png:fps=30.000000:w=800:h=600 -vf scale=800:600 -o video2.mp4 -ovc x264 -x264encopts pass=2:turbo:bitrate=400:frameref=5:bframes=1: -info name="Whole-Cell Animation":subject="Animation of Mycoplasma genitalium whole-cell simulations.":artist="Jonathan Karr, Covert Lab, Department of Bioengineering, Stanford University":copyright="Copyright Jonathan Karr, Covert Lab, Department of Bioengineering, Stanford University 2012" 
#
#Author: Jonathan Karr, jkarr@stanford.edu
#Affiliation: Covert Lab, Department of Bioengineering, Stanford University
#Last updated: 10/6/2011

use Cwd;
use HTML::Template;
use strict;
require 'library.pl';

#options
my %config = getConfiguration();
my $linuxUser = $config{'fileUserName'};
my $pathToRunTime = $config{'mcrPath'};
my $emailAddress = $config{'email'};
my $baseDir = $config{'simulationPath'};

my $outDir = "$baseDir/output/runSimulation";
my $conditionSetTimeStamp = $ARGV[0];
my $nFrameRenderJobs = 32;
my $movieSubFolder = 'movie';

#output directory
`mkdir -p "$outDir/$conditionSetTimeStamp/$movieSubFolder"`;
`sudo chown -R $linuxUser:$linuxUser "$outDir/$conditionSetTimeStamp/$movieSubFolder"`;
`chmod -R 775 "$outDir/$conditionSetTimeStamp/$movieSubFolder"`;

#compile runAnimation project
`./build.sh runAnimation`;

#analysis job
my $afterany = '';
my $submitJobs = '';
my $jobFileName = '';
my $template = HTML::Template->new(filename => 'job.animation.sh.tmpl');
$template->param(baseDir => $baseDir);
$template->param(conditionSetTimeStamp => $conditionSetTimeStamp);
$template->param(movieSubFolder => $movieSubFolder);
$template->param(linuxUser => $linuxUser);
$template->param(emailAddress => $emailAddress);
$template->param(outDir => $outDir);
$template->param(pathToRunTime => $pathToRunTime);
$template->param(nFrameRenderJobs => $nFrameRenderJobs);
$template->param(afterany => "");
$template->param(renderFrames => 1);

#frame rendering
for (my $i = 1; $i <= $nFrameRenderJobs; $i++) {
	$jobFileName = sprintf("%s/%s/%s/job.animation_%d.sh", $outDir, $conditionSetTimeStamp, $movieSubFolder, $i);
	$template->param(iJob => $i);
	open(FH, '>', $jobFileName) or die $!;
	print FH $template->output;
	close (FH);
	
	`chmod 775 $jobFileName`;
	`sudo chown -R $linuxUser:$linuxUser $jobFileName`;
	`sudo qsub $jobFileName`;
	
	my $jobId = `qstat | tail -n 1`;
	$jobId =~ /^(\d+)/;
	$afterany.=":$1";
}

#movie rendering
$template->param(iJob => $nFrameRenderJobs+1);
$template->param(afterany => $afterany);
$template->param(renderFrames => 0);
$jobFileName = sprintf("%s/%s/%s/job.animation_%d.sh", $outDir, $conditionSetTimeStamp, $movieSubFolder, $nFrameRenderJobs+1);
open(FH, '>', $jobFileName) or die $!;
print FH $template->output;
close (FH);

`chmod 775 $jobFileName`;
`sudo chown -R $linuxUser:$linuxUser $jobFileName`;
`sudo qsub $jobFileName`;

#print status message with total number of jobs submitted
print "Animation jobs queued.\n";
