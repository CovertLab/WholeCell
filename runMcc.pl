#!/usr/bin/perl
#runMcc.pl
# 
# Runs MATLAB compiler (mcc) with settings defined in .prj file. This replace some of the
# functionality of the deploytool which replaces the deprecated mcc -F functionality. In 
# particular this script can be used in environments without DISPLAY whereas deploytool cannot.
# This is particularly useful for running mcc under hudson.
# 
# Requirements: XML-DOM, Date-Format, Date-Language
# Usage: ./runMcc.pl runSmallTests.prj
# Author: Jonathan Karr, jkarr@stanford.edu
# Affiliation: Covert Lab, Department of Bioengineering, Stanford University
# Last updated: 2/24/2011

use XML::DOM;
use Switch;
use strict;
use Date::Format;
use Date::Language;
use Sys::Hostname;

#parse xml document
my $parser = new XML::DOM::Parser;
my $doc = $parser->parsefile ($ARGV[0]);
my $deploymentProject = $doc->getElementsByTagName('deployment-project')->item(0);
my $configuration = $deploymentProject->getElementsByTagName('configuration',0)->item(0);

#app name
my $appname = $configuration->getElementsByTagName('param.appname',0)->item(0)->getFirstChild()->getNodeValue();

#output
my $output = $configuration->getElementsByTagName('param.output',0)->item(0)->getFirstChild()->getNodeValue();
$output =~ s/\${PROJECT_ROOT}\///;

#toolboxes
my $toolboxes = $configuration->getElementsByTagName('param.tbx_on_path',0)->item(0)->getElementsByTagName('item',0);
my $toolboxOptions="";
for (my $i=0; $i < $toolboxes->getLength(); $i++){
	$toolboxOptions .= " -p ".$toolboxes->item($i)->getFirstChild()->getNodeValue();
}

#log file
my $outputLogfile = $configuration->getElementsByTagName('param.create.log',0)->item(0)->getFirstChild()->getNodeValue();
my $logfile = $configuration->getElementsByTagName('param.log.file',0)->item(0)->getFirstChild()->getNodeValue();

#main file
my $main = $configuration->getElementsByTagName('fileset.main',0)->item(0)->getElementsByTagName('file',0)->item(0)->getFirstChild()->getNodeValue();
$main =~ s/\${PROJECT_ROOT}\///;

#fileset
my $fileSet = $configuration->getElementsByTagName('fileset.resources',0)->item(0);
my $files = $fileSet->getElementsByTagName('file',0);
my $filesetOptions = '';
for (my $i=0; $i < $files->getLength(); $i++){
	my $file = $files->item($i)->getFirstChild()->getNodeValue();
	$file =~ s/\${PROJECT_ROOT}\///;
	$filesetOptions .= " -a $file";
}

#override .prj file settings
if ($#ARGV + 1 >= 2) {
	$appname = $ARGV[1];
	$output = 'bin/'.$ARGV[1];
	$logfile = 'output/'.$ARGV[1].'/output.log';	
	$main = $ARGV[1];
}

#run mcc
my $cmd="mcc";
$cmd .= " -o $appname";
$cmd .= " -W main:$appname";
$cmd .= " -T link:exe";
$cmd .= " -d $output";
$cmd .= " $toolboxOptions";
$cmd .= " -N";
$cmd .= " -w enable:specified_file_mismatch";
$cmd .= " -w enable:repeated_file";
$cmd .= " -w enable:switch_ignored";
$cmd .= " -w enable:missing_lib_sentinel";
$cmd .= " -w enable:demo_license";
$cmd .= " -R -nodisplay";
if ($outputLogfile eq 'true'){
	$cmd .= " -R '-logfile,$logfile'";
}
$cmd .= " -v $main";
$cmd .= " $filesetOptions";

#sbin is added to path to avoid no ifconfig warning
#- http://hpc.ucalgary.ca/quickstart/terminus/software/matlab 
#
#Export display is set to empty string to prevent mcc from requiring X display
#- http://www.mathworks.com/support/solutions/en/data/1-37K72B/index.html
#
#Note this bundles matlab.prf into the compile application. Occassionally this leads to errors. 
#There are a few suggestions to work around this by seting the prefdir environment variable to empty directory:
#- http://www.mathworks.com/support/solutions/en/data/1-1UOGW7/?solution=1-1UOGW7
#- http://www.mathworks.com/support/solutions/en/data/1-37KDWN/?solution=1-37KDWN
#- http://wiki.rcs.manchester.ac.uk/community/MatlabWithCondor
#- http://svn.broadinstitute.org/CellProfiler/trunk/CellProfiler_old_matlab/BuildCellProfiler.m
#I didn't have success following these suggestions.
#
#MCR_CACHE_ROOT needs to be set to a directory to prevent hagging from sharing network resources
#-http://www.mathworks.com/support/solutions/en/data/1-7RH0IV/index.html?product=CO&solution=1-7RH0IV
`export PATH="\$PATH:/sbin"; export DISPLAY=""; export MATLAB_PREFDIR="/tmp/emptydir"; export MCR_CACHE_ROOT="/tmp/emptydir"; export XAPPLRESDIR="/usr/local/bin/MATLAB/MATLAB_Compiler_Runtime/v716/X11/app-defaults"; export LD_LIBRARY_PATH="\$LD_LIBRARY_PATH:/usr/local/bin/MATLAB/MATLAB_Compiler_Runtime/v716/runtime/glnxa64:/usr/local/bin/MATLAB/MATLAB_Compiler_Runtime/v716/sys/os/glnxa64:/usr/local/bin/MATLAB/MATLAB_Compiler_Runtime/v716/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/usr/local/bin/MATLAB/MATLAB_Compiler_Runtime/v716/sys/java/jre/glnxa64/jre/lib/amd64/server:/usr/local/bin/MATLAB/MATLAB_Compiler_Runtime/v716/sys/java/jre/glnxa64/jre/lib/amd64"; $cmd`;

if ($? == 0){
    exit 0;
}
exit 1;

