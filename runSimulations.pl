#!/usr/bin/perl
#
#Usage:
# 1) ./build.sh runSimulation
# 2) ./build.sh runReindexing
# 3) ./build.sh runAnalysis
# 4) ./runSimulations.pl runSimulations.xml
#
#Author: Jonathan Karr, jkarr@stanford.edu
#Affiliation: Covert Lab, Department of Bioengineering, Stanford University
#Last updated: 1/9/2011

use Cwd;
use Date::Format;
use Date::Language;
use HTML::Entities;
use HTML::Template;
use MIME::Lite;
use strict;
use Switch;
use Sys::Hostname;
use XML::DOM;
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
my $baseURL = 'http://'.$config{'URL'}.'/simulation';

my $webSVNURL = $config{'webSVNURL'};

my $xmlns = 'http://'.$config{'hostName'};
my $schemaLocation = 'http://'.$config{'hostName'}.' conditions.xsd';

my $outDir = "$baseDir/output/runSimulation";
my $lang = Date::Language->new('English');
my $conditionSetTimeStamp = $lang->time2str("%Y_%m_%d_%H_%M_%S", time);
my @tmp = split('_', $conditionSetTimeStamp, 6);
my $timeStamp = $tmp[1].'/'.$tmp[2].'/'.$tmp[0].' '.
				$tmp[3].':'.$tmp[4].':'.$tmp[5];
my $useSeedOffset = 1;

#check output directory exists
if (-d $outDir){}
else {die "Output directory doesn't exist";}

#make output directory
mkdir "$outDir/$conditionSetTimeStamp" or die "Unable to make directory $outDir/$conditionSetTimeStamp: $!";;

#compile runSimulation, runReindexing, runAnalysis projects
`./build.sh runSimulation`;
`./build.sh runReindexing`;
`./build.sh runAnalysis`;

#copy executables
mkdir "$outDir/$conditionSetTimeStamp/bin" or die "Unable to make directory $outDir/$conditionSetTimeStamp/bin: $!";;
`cp -R bin/runSimulation $outDir/$conditionSetTimeStamp/bin`;
`cp -R bin/runReindexing $outDir/$conditionSetTimeStamp/bin`;
`cp -R bin/runAnalysis $outDir/$conditionSetTimeStamp/bin`;

#parse xml document
my $conditionsXMLFile = $ARGV[0];
my $parser = new XML::DOM::Parser;
my $doc = $parser->parsefile ($conditionsXMLFile);
`cp $conditionsXMLFile $outDir/$conditionSetTimeStamp/conditions.xml`;

#error check document
my $conditions;
my $hasConditions = 0;
for (my $i=0; $i < $doc->getChildNodes()->getLength(); $i++){
	my $child = $doc->getChildNodes()->item($i);
	switch ($child->getNodeName()) {
		case 'conditions' { $hasConditions++;}
		case '#comment' {}
		else{ die "Unsupported tag ".$child->getNodeName()."\n";}
	}
}
if ($hasConditions!=1) {die "Invalid XML\n";}

#error check conditions -- has non empty name, description
$conditions = $doc->getElementsByTagName('conditions', 0)->item(0);
my %metadata = validate_metadata($conditions);
my $firstName = $metadata{'firstName'};
my $lastName = $metadata{'lastName'};
my $email = $metadata{'email'};
my $affiliation = $metadata{'affiliation'};
my $userName = $metadata{'userName'};
my $hostName = $metadata{'hostName'};
my $ipAddress = $metadata{'ipAddress'};
my $revision = $metadata{'revision'};
my $differencesFromRevision = $metadata{'differencesFromRevision'};

#error check conditions
my $conditionsArr = $conditions->getElementsByTagName('condition', 0);
my $conditionsHTML = validate_conditions($conditionsArr);

#submit jobs for each condition
my $nJobs = 0;
my $N;
my $simulationId;
my $simulationIdx;
my $condition = '';
my $submitJobs = '';
my $afterany = '';
my $afteranyAnalysis = '';
my $template;
my $allSingleGeneDeletions = 1;
my $seedOffset = ($useSeedOffset ? time : 0);
for (my $i=0; $i < $conditionsArr->getLength(); $i++){
	$condition = $conditionsArr->item($i);
	my $replicates = $condition->getElementsByTagName("replicates", 0);
	if ($replicates->getLength() > 0) {
		$N = $replicates->item(0)->getFirstChild()->getNodeValue();
		$condition->removeChild($replicates->item(0));
	}else{
		$N = 1;
	}
	my $options = $condition->getElementsByTagName("options", 0);
	if ($options->getLength() == 0) {
		$options = $condition->appendChild($doc->createElement("options"));
	}else{
          $options = $options->item(0);
        }
	my $tmpOptions = $options->getElementsByTagName("option", 0);
	my $seed = "";
	for (my $i=0; $i< $tmpOptions->getLength(); $i++){
		my $option = $tmpOptions->item($i);
		my $tf1 = 0;
		my $tf2 = 1;
		my $tf3 = 1;
		for (my $j = 0; $j < $option->getAttributes()->getLength(); $j++){
			my $child = $option->getAttributes()->item($j);
			if ($child->getNodeName() == 'seed'){
				$tf1 = 1;
			}
			if ($child->getNodeName() == 'state'){
				$tf2 = 0;
			}
			if ($child->getNodeName() == 'process'){
				$tf3 = 0;
			}
		}
		if ($tf1==1 && !$tf2==0 && !$tf3==0){
			$seed = $option;
		}
	}
	if ($seed == ""){
		$seed = $options->appendChild($doc->createElement("option"));
		$seed->setAttribute("name", "seed");
	}
	
	my $nGeneDeletions = 0;
	my $perturbations = $condition->getElementsByTagName("perturbations", 0);
	if ($perturbations->getLength > 0){
		$perturbations = $perturbations->item(0)->getElementsByTagName("perturbation", 0);
	}
	for (my $n = 0; $n < $perturbations->getLength(); $n++){
		my $perturbation = $perturbations->item($n);
		for (my $j = 0; $j < $perturbation->getAttributes()->getLength(); $j++){
			my $child = $perturbation->getAttributes()->item($j);
			if ($child->getNodeName() eq 'type' && $child->getFirstChild()->getNodeValue() eq 'geneticKnockout'){
				$nGeneDeletions++;
			}
		}
	}
	$allSingleGeneDeletions = $allSingleGeneDeletions && ($nGeneDeletions == 1);
	
	#submit job
	for (my $n = 0; $n < $N; $n++){
		$nJobs++;

		#create output file, setup output file names
		my $dirName = sprintf("%s/%s/%d", $outDir, $conditionSetTimeStamp, $nJobs);
		mkdir $dirName or die "Unable to make directory $dirName: $!";
		
		my $conditionFileName = sprintf("%s/conditions.xml", $dirName);
		$seed->setAttribute("value", $seedOffset + $nJobs);
		
		#output condition xml file
		open(FH, '>', $conditionFileName) or die $!;
		print FH "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
		print FH "<!--\n";
		print FH "Autogenerated condition set.\n";
		print FH "\n";
		print FH "-->\n";
		print FH "<conditions\n";
		print FH "\txmlns=\"$xmlns\"\n";
		print FH "\txmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n";
		print FH "\txsi:schemaLocation=\"$schemaLocation\">\n";
		print FH sprintf("\t<firstName><![CDATA[%s]]></firstName>\n", $firstName);
		print FH sprintf("\t<lastName><![CDATA[%s]]></lastName>\n", $lastName);
		print FH sprintf("\t<email><![CDATA[%s]]></email>\n", $email);
		print FH sprintf("\t<affiliation><![CDATA[%s]]></affiliation>\n", $affiliation);
		print FH sprintf("\t<userName><![CDATA[%s]]></userName>\n", $userName);
		print FH sprintf("\t<hostName><![CDATA[%s]]></hostName>\n", $hostName);
		print FH sprintf("\t<ipAddress><![CDATA[%s]]></ipAddress>\n", $ipAddress);
		print FH sprintf("\t<revision><![CDATA[%s]]></revision>\n", $revision);
		print FH sprintf("\t<differencesFromRevision><![CDATA[%s]]></differencesFromRevision>\n", encode_entities($differencesFromRevision));
		print FH "\t".$condition->toString()."\n";
		print FH "</conditions>\n";
		close (FH);
	}
}

$submitJobs = '';
$template = HTML::Template->new(filename => 'job.simulation.sh.tmpl');
$template->param(conditionSetTimeStamp => $conditionSetTimeStamp);
$template->param(linuxRunUser => $linuxRunUser);
$template->param(emailAddress => $emailAddress);
$template->param(outDir => $outDir);
$template->param(storageServer => $storageServer);
$template->param(pathToRunTime => $pathToRunTime);
$template->param(nodeTmpDir => $nodeTmpDir);
for (my $n = 1; $n <= $nJobs; $n++){
	my $jobFileName = sprintf("%s/%s/%d/job.simulation.sh", $outDir, $conditionSetTimeStamp, $n);
	open(FH, '>', $jobFileName) or die $!;
	$template->param(n => $n);
    print FH $template->output;
	close (FH);
	
	$submitJobs .= "sudo qsub $jobFileName; ";
}

#set permissions and run jobs
`sudo chmod -R 775 $outDir/$conditionSetTimeStamp`;
`sudo chown -R $linuxUser:$linuxUser $outDir/$conditionSetTimeStamp`;
`$submitJobs`;

#get job id
$simulationId = `qstat | tail -n 1`;
$simulationId =~ /^(\d+)/;
$simulationIdx = $1-$nJobs+1;

#save PBS id
open(FH, '>', "$outDir/$conditionSetTimeStamp/pbsid") or die $!;
print FH "$simulationIdx\n";
close (FH);

#reindexing job
$submitJobs = '';
$template = HTML::Template->new(filename => 'job.reindexing.sh.tmpl');
$template->param(conditionSetTimeStamp => $conditionSetTimeStamp);
$template->param(linuxRunUser => $linuxRunUser);
$template->param(emailAddress => $emailAddress);
$template->param(outDir => $outDir);
$template->param(storageServer => $storageServer);
$template->param(pathToRunTime => $pathToRunTime);
$template->param(nodeTmpDir => $nodeTmpDir);
for (my $n = 1; $n <= $nJobs; $n++){
	$afterany = $simulationIdx+$n-1;
	my $jobFileName2 = sprintf("%s/%s/%d/job.reindexing.sh", $outDir, $conditionSetTimeStamp, $n);
	open(FH, '>', $jobFileName2) or die $!;
	$template->param(n => $n);
	$template->param(afterany => $afterany);
	print FH $template->output;
	close (FH);
	
	$submitJobs .= "sudo qsub $jobFileName2; ";
}

#set permissions and run jobs
`sudo chmod -R 775 $outDir/$conditionSetTimeStamp`;
`sudo chown -R $linuxUser:$linuxUser $outDir/$conditionSetTimeStamp`;
`$submitJobs`;


#analysis job
$afterany = '';
for (my $n = $simulationIdx + $nJobs; $n < $simulationIdx + 2 * $nJobs; $n++){
	$afterany.=":$n";
}

my $nFirstStageAnalysisJobs = ($allSingleGeneDeletions ? 1 : 8);
$submitJobs = '';
$afteranyAnalysis = '';
$template = HTML::Template->new(filename => 'job.analysis.sh.tmpl');
$template->param(baseDir => $baseDir);
$template->param(conditionSetTimeStamp => $conditionSetTimeStamp);
$template->param(linuxUser => $linuxUser);
$template->param(emailAddress => $emailAddress);
$template->param(outDir => $outDir);
$template->param(pathToRunTime => $pathToRunTime);
$template->param(simulationIdx => $simulationIdx);
$template->param(nFirstStageJobs => $nFirstStageAnalysisJobs);
for (my $i = 1; $i <= $nFirstStageAnalysisJobs + 1; $i++) {
	my $jobFileName3 = sprintf("%s/%s/job.analysis_%d.sh", $outDir, $conditionSetTimeStamp, $i);
	$template->param(iJob => $i);
	if ($i <= $nFirstStageAnalysisJobs) {
		$template->param(afterany => $afterany);
		$template->param(firstStageAnalysis => 1);
		$afteranyAnalysis.=":".($simulationIdx + 2 * $nJobs + $i - 1);
	} else {
		$template->param(afterany => $afteranyAnalysis);
		$template->param(firstStageAnalysis => 0);
	}
	open(FH, '>', $jobFileName3) or die $!;
	print FH $template->output;
	close (FH);
	
	$submitJobs .= "sudo qsub $jobFileName3; ";
}

#set permissions and submit job
`sudo chmod -R 775 $outDir/$conditionSetTimeStamp`;
`sudo chown -R $linuxUser:$linuxUser $outDir/$conditionSetTimeStamp`;
`$submitJobs`;

#print status message with total number of jobs submitted
print "\Simulation set #$simulationIdx queued with $nJobs simulations.\n";

#mail message that job submitted
my $webSVNURLRevision = sprintf($webSVNURL, $revision);
my $emailBody = <<HTML;
<html>
<style type="text/css">
body, table {
	font-family: Arial;
	font-size:10pt;
}
th {
	font-weight:bold;
	text-align:left;
}
th, td {
	padding-left:10px;
}
th:first-child, td:first-child{
	padding-left:0px;
}
</style>
<body>
Simulation set <a href="$baseURL/viewSimulationSet.php?id=$conditionSetTimeStamp">#$simulationIdx</a> has been queued with $nJobs simulations. You will receive a summary email upon completion of this simulation set.<br/>
<br/>
<table cellspacing="0" cellpadding="0">
	<tr><th>Researcher</th><td>$firstName $lastName, <a href="mailto:$email">$email</a></td></tr>
	<tr><th>Affiliation</th><td>$affiliation</td></tr>
	<tr><th>Hostname</th><td>$hostName, $ipAddress</td></tr>
	<tr><th>Revision</th><td><a href="$webSVNURLRevision">$revision</a></td></tr>
	<tr><th>Differences</th><td style="font-style:italic;">See attached</td></tr>
	$conditionsHTML
</table>
</html>
HTML

my $msg = MIME::Lite->new(
	From    => 'PBS Daemon <no-reply@stanford.edu>',
	To      => $emailAddress,
	Subject => "Simulation set #$simulationIdx $timeStamp queued with $nJobs simulations",
	Type    => 'multipart/mixed'
);

$msg->attach(
	Type     => 'text/html',
	Data     => $emailBody
);

$msg->attach(
	Type     => 'text/plain',
	Data     => $differencesFromRevision,
	Filename => "diff.txt"
);

$msg->send;
