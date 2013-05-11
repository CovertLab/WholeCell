#!/usr/bin/perl

#Usage: ./mailSimulationsSummary.pl 259 2011_06_23_01_01_01
#Author: Jonathan Karr, jkarr@stanford.edu
#Affiliation: Covert Lab, Department of Bioengineering, Stanford University
#Last updated: 6/23/2011

use Switch;
use strict;
use Cwd;
use MIME::Lite;
require "/home/projects/WholeCell/simulation/library.pl";

my $simulationIdx = $ARGV[0];
my $simulationDateTime = $ARGV[1];

my %config = getConfiguration();
my $emailAddress = $config{'email'};
my $baseDir = $config{'simulationPath'};
my $outDir = "$baseDir/output/runSimulation";
my $baseURL = 'http://'.$config{'URL'}.'/simulation';
my $webSVNURL = $config{'webSVNURL'};
my @timeStamp = split('_', $simulationDateTime, 6);

chdir($baseDir) or die "$!";

my $parser = new XML::DOM::Parser;
my $doc = $parser->parsefile ("$outDir/$simulationDateTime/conditions.xml");
my $conditions = $doc->getElementsByTagName('conditions', 0)->item(0);
my %metadata = validate_metadata($conditions);
my $firstName = $metadata{'firstName'};
my $lastName = $metadata{'lastName'};
my $email = $metadata{'email'};
my $affiliation = $metadata{'affiliation'};
my $userName = $metadata{'userName'};
my $hostName = $metadata{'hostName'};
my $ipAddress = $metadata{'ipAddress'};
my $revision = $metadata{'revision'};
my $differencesFromRevision = $metadata{'$differencesFromRevision'};
my $conditionsArr = $conditions->getElementsByTagName('condition', 0);
my $conditionsHTML = validate_conditions($conditionsArr);
my $nJobs = count_simulations($conditionsArr);
my $webSVNURLRevision = sprintf($webSVNURL, $revision);
my $errors = '';
unless ((-e "$outDir/$simulationDateTime/summary.html") && (-e "$outDir/$simulationDateTime/population-Growth.pdf")) {
	$errors = "Simulation summary not available. See error logs.";
}

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
td {
	padding-left:10px;
}
</style>
<body>
Simulation set <a href="$baseURL/viewSimulationSet.php?id=$simulationDateTime">#$simulationIdx</a> has finished with $nJobs simulations.<br/>
<br/>
<table cellspacing="0" cellpadding="0">
	<tr><th>Researcher</th><td>$firstName $lastName, <a href="mailto:$email">$email</a></td></tr>
	<tr><th>Affiliation</th><td>$affiliation</td></tr>
	<tr><th>Hostname</th><td>$hostName, $ipAddress</td></tr>
	<tr><th>Revision</th><td><a href="$webSVNURLRevision">$revision</a></td></tr>	
	<tr><th>Differences</th><td style="font-style:italic;">See initial status email</td></tr>	
	$conditionsHTML
	<tr><th>Errors</th><td style="color: red">$errors</td></th>
</table>
</html>
HTML

#mail summary using MIME::LITE
#- http://www.perlfect.com/articles/sendmail.shtml
#- http://www.gettingclever.com/2009/01/sending-mail-with-multiple-attachments.html
my $msg = MIME::Lite->new(		
	From    => 'PBS Daemon <no-reply@stanford.edu>',
	To      => $emailAddress,
	Subject => "Simulation set #$simulationIdx ".
		$timeStamp[1].'/'.$timeStamp[2].'/'.$timeStamp[0].' '.
		$timeStamp[3].':'.$timeStamp[4].':'.$timeStamp[5].' '.
		"finished",
	Type    => 'multipart/mixed'
);

$msg->attach(
	Type     => 'text/html',
	Data     => $emailBody
);

if (-e "$outDir/$simulationDateTime/summary.html") {
	$msg->attach(
		Type        => 'text/html',
		Disposition => 'attachment',
		Filename    => 'summary.html',
		Path        => "$outDir/$simulationDateTime/summary.html"
	);
}

if (-e "$outDir/$simulationDateTime/population-Growth.pdf") {
	$msg->attach(
		Type        => 'application/pdf',
		Disposition => 'attachment',
		FileName    => 'populationSummary.pdf',
		Path        => "$outDir/$simulationDateTime/population-Growth.pdf"
	);
}

$msg->send;