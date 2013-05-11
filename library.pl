#!/usr/bin/perl

#Author: Jonathan Karr, jkarr@stanford.edu
#Affiliation: Covert Lab, Department of Bioengineering, Stanford University
#Last updated: 7/21/2011

use HTML::Template;
use POSIX qw(ceil floor);
use Switch;
use XML::DOM;

sub getConfiguration(){
	return (
		'hostName' => 'covertlab.stanford.edu',
		'URL' => 'covertlab.stanford.edu/projects/WholeCell',
		'schema' => 'wholecell',
		'dbUserName' => 'wholecell',
		'dbPassword' => 'wholecell',
		
		'email' => 'jkarr@stanford.edu',
		
		'execUserName' => 'jkarr',
		'fileUserName' => 'WholeCell',
		'simulationHostName' => 'covertlab',
		'simulationPath' => '/home/projects/WholeCell/simulation',
		'mcrPath' => '/usr/local/bin/MATLAB/MATLAB_Compiler_Runtime/v716',
		'nodeTmpDir' => '/state/partition1',
		
		'webSVNURL' => 'http://covertlab.stanford.edu/websvn/wsvn/WholeCell?op=revision&rev=%d',
		'gangliaJobURL' => 'http://covertlab-cluster.stanford.edu/ganglia/addons/rocks/job.php?c=covertlab-cluster&id=%d-%d',
		'gangliaQueueURL' => 'http://covertlab-cluster.stanford.edu/ganglia/addons/rocks/queue.php?c=covertlab-cluster',
		'checkJobURL' => 'http://covertlab-cluster.stanford.edu/cgi-bin/checkjob.cgi?jobid=%d-%d',
	);
}

sub validate_metadata{
	my $conditions = $_[0];
	
	my $firstName;
	my $lastName;
	my $email;
	my $affiliation;
	my $userName = "";
	my $hostName = "";
	my $ipAddress = "";
	my $revision = "";
	my $differencesFromRevision = "";
	my $hasFirstName = 0;
	my $hasLastName = 0;
	my $hasEmail = 0;
	my $hasAffiliation = 0;
	for (my $i=0; $i < $conditions->getChildNodes()->getLength(); $i++){
		my $child = $conditions->getChildNodes()->item($i);
		switch ($child->getNodeName()) {
			case 'firstName' {$hasFirstName++; if (!$child->getFirstChild() || !length($child->getFirstChild()->getNodeValue())) { die "Invalid name"; } $firstName = $child->getFirstChild()->getNodeValue();}
			case 'lastName' {$hasLastName++; if (!$child->getFirstChild() || !length($child->getFirstChild()->getNodeValue())) { die "Invalid name"; } $lastName = $child->getFirstChild()->getNodeValue();}
			case 'email' {$hasEmail++; if (!$child->getFirstChild() || !length($child->getFirstChild()->getNodeValue())) { die "Invalid email"; } $email = $child->getFirstChild()->getNodeValue();}
			case 'affiliation' {$hasAffiliation++; if (!$child->getFirstChild() || !length($child->getFirstChild()->getNodeValue())) { die "Invalid affiliation"; } $affiliation = $child->getFirstChild()->getNodeValue();}
			case 'userName' {if (!$child->getFirstChild() || !length($child->getFirstChild()->getNodeValue())) { die "Invalid user name"; } $userName = $child->getFirstChild()->getNodeValue();}
			case 'hostName' {if (!$child->getFirstChild() || !length($child->getFirstChild()->getNodeValue())) { die "Invalid host name"; } $hostName = $child->getFirstChild()->getNodeValue();}
			case 'ipAddress' {if (!$child->getFirstChild() || !length($child->getFirstChild()->getNodeValue())) { die "Invalid IP Address"; } $ipAddress = $child->getFirstChild()->getNodeValue();}
			case 'revision' {if (!$child->getFirstChild() || !length($child->getFirstChild()->getNodeValue())) { die "Invalid revision"; } $revision = $child->getFirstChild()->getNodeValue();}
			case 'differencesFromRevision' {if ($child->getFirstChild()) {$differencesFromRevision = $child->getFirstChild()->getNodeValue();}}
			case 'condition' {}
			case '#comment' {}
			case '#text' {}
			else{ die "Unsupported tag ".$child->getNodeName()."\n";}
		}
	}
	if ($hasFirstName!=1) {die "Condition set must have 1 first name\n";}
	if ($hasLastName!=1) {die "Condition set must have 1 last name\n";}
	if ($hasEmail!=1) {die "Condition set must have 1 email\n";}
	if ($hasAffiliation!=1) {die "Condition set must have 1 affiliation\n";}

	if ($userName eq ""){
	  $userName = getlogin();
	}
	if ($hostName eq ""){
	  $hostName = hostname();
	}
	if ($ipAddress eq ""){
	  $ipAddress = `/sbin/ifconfig | grep 'inet addr:'| grep -v '127.0.0.1' | grep -v '10.1.1.1' | cut -d: -f2 | cut "-d " -f1`;
	  $ipAddress =~ s/\s+$//;
	}
	if ($revision eq ""){
	  $revision = `svn info | grep 'Revision:' | cut "-d:" -f2`;
	  $revision =~ s/^\s+//g;
	  $revision =~ s/\s+$//g;
	}
	if ($differencesFromRevision eq ""){
	  $differencesFromRevision = `svn diff`;
	  $differencesFromRevision =~ s/\]\]>//g;
	  $differencesFromRevision =~ s/<!\[CDATA\[//g;
	}
	
	return (
		firstName => $firstName, 
		lastName => $lastName,
		email => $email,
		affiliation => $affiliation,
		userName => $userName,
		hostName => $hostName,
		ipAddress => $ipAddress,
		revision => $revision,
		differencesFromRevision => $differencesFromRevision
	);
}

sub validate_conditions{
	my $conditionsArr = $_[0];
	
	my $conditionsHTML = '';	
	for (my $i=0; $i < $conditionsArr->getLength(); $i++){
		my $condition = $conditionsArr->item($i);

		my $name;
		my $description;
		my $replicates;
		my @options;
		my @parameters;
		my @perturbations;
		my $hasName = 0;
		my $hasDescription = 0;	
		for (my $j = 0; $j < $condition->getChildNodes()->getLength(); $j++){
			my $child = $condition->getChildNodes()->item($j);	
			switch ($child->getNodeName()) {
				case 'shortDescription' {
					$hasName++; 
					if (!$child->getFirstChild() || !length($child->getFirstChild()->getNodeValue())) {
						die "Invalid name"; 
					} 
					$name = $child->getFirstChild()->getNodeValue();
				}
				case 'longDescription' {
					$hasDescription++; 
					if (!$child->getFirstChild() || !length($child->getFirstChild()->getNodeValue())) {
						die "Invalid description"; 
					}
					$description = $child->getFirstChild()->getNodeValue();
				}
				case 'replicates' {				
					if (!$child->getFirstChild() || !($child->getFirstChild()->getNodeValue()=~m/^\d+$/)) { 
						die "Replicates must be a non-negative integer\n";
					}
					$replicates = $child->getFirstChild()->getNodeValue();
				}
				case 'options' { @options = validate_options($child);}
				case 'parameters'{ @parameters = validate_parameters($child);}
				case 'perturbations'{ @perturbations = validate_perturbations($child);}
				case '#comments' {}
				case '#text' {}
				else { die "Unsupported tag ".$child->getNodeName()."\n";}
			}
		}
		if ($hasName != 1) {die "Condition must have 1 name\n";}
		if ($hasDescription != 1) {die "Condition must have 1 description\n";}
		
		my $idx=$i+1;
		my $descriptionHTML = $description;
		$descriptionHTML =~ s/\n/<br\/>/g;
		$conditionsHTML .= <<HTML;
		<tr><th colspan="2" style="padding-top:8px;">Condition Set #$idx</th>
		<tr><th style="padding-left:10px;">Name</th><td>$name</td></tr>
		<tr><th style="padding-left:10px;vertical-align:top;">Description</th><td>$descriptionHTML</td></tr>
		<tr><th style="padding-left:10px;">Replicates</th><td>$replicates</td></tr>
HTML
		if (@options > 0){
			my $optionsHTML = '';
			for (my $j=0; $j< scalar(@options); $j++){
				$optionsHTML.="<tr>\n";
				$optionsHTML.="<td>".($options[$j]{'state'} || $options[$j]{'process'} ? $options[$j]{'state'}.$options[$j]{'process'} : "&nbsp;")."</td>\n";
				$optionsHTML.="<td>".($options[$j]{'name'} ? $options[$j]{'name'} : "&nbsp;")."</td>\n";
				$optionsHTML.="<td>".($options[$j]{'value'} ? $options[$j]{'value'} : "&nbsp;")."</td>\n";
				$optionsHTML.="</tr>\n";
			}
			$conditionsHTML .= <<HTML;
		<tr><th style="padding-left:10px;">Options</th><td>&nbsp;</td></tr>
		<tr><td colspan="2" style="padding-left:20px;">
			<table cellspacing="0" cellpadding="0">
			<thead>
				<tr>
					<th>State/Process</th>
					<th>Name</th>
					<th>Value</th>
				</tr>
			</thead>
			<tbody>$optionsHTML</tbody>
			</table></td></tr>
HTML
		}
		if (@parameters > 0){
			my $parametersHTML = '';
			for (my $j=0; $j< scalar(@parameters); $j++){
				$parametersHTML.="<tr>\n";
				$parametersHTML.="<td>".($parameters[$j]{'state'} || $parameters[$j]{'process'} ? $parameters[$j]{'state'}.$parameters[$j]{'process'} : "&nbsp;")."</td>\n";
				$parametersHTML.="<td>".($parameters[$j]{'name'} ? $parameters[$j]{'name'} : "&nbsp;")."</td>\n";
				$parametersHTML.="<td>".($parameters[$j]{'index'} ? $parameters[$j]{'index'} : "&nbsp;")."</td>\n";
				$parametersHTML.="<td>".($parameters[$j]{'value'} ? $parameters[$j]{'value'} : "&nbsp;")."</td>\n";
				$parametersHTML.="</tr>\n";
			}
			$conditionsHTML .= <<HTML;
		<tr><th style="padding-left:10px;">Parameters</th><td>&nbsp;</td></tr>
		<tr><td colspan="2" style="padding-left:20px;">
			<table cellspacing="0" cellpadding="0">
			<thead>
				<tr>
					<th>State/Process</th>
					<th>Name</th>				
					<th>Index</th>
					<th>Value</th>
				</tr>
			</thead>
			<tbody>$parametersHTML</tbody>
			</table></td></tr>
HTML
		}	
		if (@perturbations > 0){
			my $perturbationsHTML = '';
			for (my $j=0; $j< scalar(@perturbations); $j++){
				$perturbationsHTML.="<tr>\n";
				$perturbationsHTML.="<td>".($perturbations[$j]{'type'} ? $perturbations[$j]{'type'} : "&nbsp;")."</td>\n";
				$perturbationsHTML.="<td>".($perturbations[$j]{'component'} ? $perturbations[$j]{'component'} : "&nbsp;")."</td>\n";
				$perturbationsHTML.="<td>".($perturbations[$j]{'compartment'} ? $perturbations[$j]{'compartment'} : "&nbsp;")."</td>\n";
				$perturbationsHTML.="<td>".($perturbations[$j]{'initialTime'} ? $perturbations[$j]{'initialTime'} : "&nbsp;")."</td>\n";
				$perturbationsHTML.="<td>".($perturbations[$j]{'finalTime'} ? $perturbations[$j]{'finalTime'} : "&nbsp;")."</td>\n";
				$perturbationsHTML.="<td>".($perturbations[$j]{'value'} ? $perturbations[$j]{'value'} : "&nbsp;")."</td>\n";
				$perturbationsHTML.="</tr>\n";
			}
			$conditionsHTML .= <<HTML;
		<tr><th style="padding-left:10px;">Perturbations</th><td>&nbsp;</td></tr>
		<tr><td colspan="2" style="padding-left:20px;">
			<table cellspacing="0" cellpadding="0">
			<thead>
				<tr>
					<th>Type</th>
					<th>Component</th>				
					<th>Compartment</th>
					<th>Initial Time</th>
					<th>Final Time</th>
					<th>Value</th>
				</tr>
			</thead>
			<tbody>$perturbationsHTML</tbody>
			</table></td></tr>
HTML
		}
	}
	
	return $conditionsHTML;
}

sub count_simulations{
	my $conditionsArr = $_[0];
	
	my $nSimulations = 0;
	for (my $i=0; $i < $conditionsArr->getLength(); $i++){
		my $condition = $conditionsArr->item($i);		
		for (my $j = 0; $j < $condition->getChildNodes()->getLength(); $j++){
			my $child = $condition->getChildNodes()->item($j);	
			switch ($child->getNodeName()) {
				case 'shortDescription' {}
				case 'longDescription' {}
				case 'replicates' {				
					$nSimulations += $child->getFirstChild()->getNodeValue();
				}
				case 'options' {}
				case 'parameters'{}
				case 'perturbations'{}
				case '#comments' {}
				case '#text' {}
				else { die "Unsupported tag ".$child->getNodeName()."\n";}
			}
		}
	}
	
	return $nSimulations;
}

sub validate_options{
	my $optionsXML = $_[0];	

	for (my $i=0; $i < $optionsXML->getChildNodes()->getLength(); $i++){
		my $child = $optionsXML->getChildNodes()->item($i);	
		switch ($child->getNodeName()) {
			case 'option' {}
			case '#comments' {}
			case '#text' {}
			else { die "Unsupported tag ".$child->getNodeName()."\n";}
		}
	}
	
	my @options = ();	
	my $optionsXMLArr = $optionsXML->getElementsByTagName('option', 0);
	for (my $i = 0; $i < $optionsXMLArr->getLength(); $i++){
		my $option = $optionsXMLArr->item($i);
		my $hasName = 0;
		my $name;
		my $state;
		my $process;
		my $value;
		for (my $j = 0; $j < $option->getAttributes()->getLength(); $j++){
			my $child = $option->getAttributes()->item($j);	
			switch ($child->getNodeName()) {
				case 'name' {$hasName++; if (!$child->getFirstChild() || !length($child->getFirstChild()->getNodeValue())) { die "Invalid name"; } $name = $child->getFirstChild()->getNodeValue();}
				case 'state' {if ($child->getFirstChild() && length($child->getFirstChild()->getNodeValue())) {$state = $child->getFirstChild()->getNodeValue();}}
				case 'process' {if ($child->getFirstChild() && length($child->getFirstChild()->getNodeValue())) {$process = $child->getFirstChild()->getNodeValue();}}
				case 'value' {if ($child->getFirstChild() && length($child->getFirstChild()->getNodeValue())) {$value = $child->getFirstChild()->getNodeValue();}}				
				case '#comments' {}
				case '#text' {}
				else { die "Unsupported tag ".$child->getNodeName()."\n";}
			}
		}
		if ($hasName != 1) {die "Option must have 1 name\n";}
		push(@options, {name=>$name, state=>$state, process=>$process, value=>$value});
	}
	
	return @options;
}

sub validate_parameters{
	my $parametersXML = $_[0];

	for (my $i=0; $i < $parametersXML->getChildNodes()->getLength(); $i++){
		my $child = $parametersXML->getChildNodes()->item($i);	
		switch ($child->getNodeName()) {
			case 'parameter' {}
			case '#comments' {}
			case '#text' {}
			else { die "Unsupported tag ".$child->getNodeName()."\n";}
		}
	}
	
	my @parameters = ();
	my $parametersArrXML = $parametersXML->getElementsByTagName('parameter', 0);
	for (my $i=0; $i < $parametersArrXML->getLength(); $i++){
		my $parameter = $parametersArrXML->item($i);
		my $hasName = 0;
		my $hasState = 0;
		my $hasProcess = 0;
		my $name;
		my $index;
		my $state;
		my $process;
		my $value;
		for (my $j = 0; $j < $parameter->getAttributes()->getLength(); $j++){
			my $child = $parameter->getAttributes()->item($j);	
			switch ($child->getNodeName()) {
				case 'name' {$hasName++; if (!$child->getFirstChild() || !length($child->getFirstChild()->getNodeValue())) { die "Invalid name"; } $name = $child->getFirstChild()->getNodeValue();}
                case 'index' {if ($child->getFirstChild() && length($child->getFirstChild()->getNodeValue())) {$index = $child->getFirstChild()->getNodeValue();}}
				case 'state' {$hasState++; if (!$child->getFirstChild() || !length($child->getFirstChild()->getNodeValue())) { die "Invalid state"; } $state = $child->getFirstChild()->getNodeValue();}
				case 'process' {$hasProcess++; if (!$child->getFirstChild() || !length($child->getFirstChild()->getNodeValue())) { die "Invalid process"; } $process = $child->getFirstChild()->getNodeValue();}
				case 'value' {if ($child->getFirstChild() && length($child->getFirstChild()->getNodeValue())) {$value = $child->getFirstChild()->getNodeValue();}}				
				case '#comments' {}
				case '#text' {}
				else { die "Unsupported tag ".$child->getNodeName()."\n";}
			}
		}
		if ($hasName!=1) {die "Parameter must have 1 name\n";}
		if ($hasState+$hasProcess!=1) {die "Parameter must have 1 state/process\n";}
		
		push(@parameters, {name=>$name, index=>$index, state=>$state, process=>$process, value=>$value});
	}
	
	return @parameters;
}

sub validate_perturbations{
	my $perturbationsXML = $_[0];

	for (my $i=0; $i < $perturbationsXML->getChildNodes()->getLength(); $i++){
		my $child = $perturbationsXML->getChildNodes()->item($i);	
		switch ($child->getNodeName()) {
			case 'perturbation' {}
			case '#comments' {}
			case '#text' {}
			else { die "Unsupported tag ".$child->getNodeName()."\n";}
		}
	}
	
	my @perturbations = ();
	my $perturbationsArrXMl = $perturbationsXML->getElementsByTagName('perturbation', 0);
	for (my $i=0; $i < $perturbationsArrXMl->getLength(); $i++){
		my $perturbation = $perturbationsArrXMl->item($i);
		my $type;
		my $hasType = 0;
		my $hasComponent = 0;
		my $hasCompartment = 0;		
		my $hasTime = 0;
		my $hasValue = 0;
		my $component;
		my $compartment;
		my $initialTime;
		my $finalTime;
		my $value;
		for (my $j = 0; $j < $perturbation->getAttributes()->getLength(); $j++){
			my $child = $perturbation->getAttributes()->item($j);	
			switch ($child->getNodeName()) {
				case 'type' {
					$hasType++;					
					if (!$child->getFirstChild() || !length($child->getFirstChild()->getNodeValue()) || 
						!($child->getFirstChild()->getNodeValue() eq 'geneticKnockout' ||
						$child->getFirstChild()->getNodeValue() eq 'stimulus' ||
						$child->getFirstChild()->getNodeValue() eq 'media')) {
						die "Type must be geneticKnockout, stimulus, or media";
					}
					$type = $child->getFirstChild()->getNodeValue();
				}
				case 'component' {$hasComponent++; if (!$child->getFirstChild() || !length($child->getFirstChild()->getNodeValue())) { die "Invalid component"; } $component = $child->getFirstChild()->getNodeValue();}
				case 'compartment' {$hasCompartment++; if (!$child->getFirstChild() || !length($child->getFirstChild()->getNodeValue())) { die "Invalid compartment"; } $compartment = $child->getFirstChild()->getNodeValue();}
				case 'initialTime' {$hasTime++; if (!$child->getFirstChild() || !length($child->getFirstChild()->getNodeValue())) { die "Invalid initial time"; } $initialTime = $child->getFirstChild()->getNodeValue();}
				case 'finalTime' {$hasTime++; if (!$child->getFirstChild() || !length($child->getFirstChild()->getNodeValue())) { die "Invalid final time"; } $finalTime = $child->getFirstChild()->getNodeValue();}
				case 'value' {$hasValue++; if (!$child->getFirstChild() || !length($child->getFirstChild()->getNodeValue())) { die "Invalid value"; } $value = $child->getFirstChild()->getNodeValue();}
				case '#comments' {}
				case '#text' {}
				else { die "Unsupported tag ".$child->getNodeName()."\n";}
			}
		}
		if ($hasType!=1) {die "perturbation must have 1 type\n";}
		if ($hasComponent!=1) {die "perturbation must have 1 component\n";}
		if ($type eq 'geneticKnockout' && ($hasCompartment!=0 || $hasTime!=0 || $hasValue!=0)){ die "Invalid perturbation\n";}
		if ($type ne 'geneticKnockout' && ($hasCompartment!=1 || $hasValue!=1)){ die "Invalid perturbation\n";}
		
		push(@perturbations, {type=>$type, component=>$component, compartment=>$compartment, initialTime=>$initialTime, finalTime=>$finalTime, value=>$value});
	}
	
	return @perturbations;
}

sub compileFigureSummary{
	my $title = $_[0];
	my @figures = @{$_[1]};
	my $outDir = $_[2];
	my $texFileName = $_[3];
	my $verbosity = $_[4];

	my $pages = '';
	my $pageIdx = 0;
	my $nSections = 0;
	foreach (@figures) {
		my %section = %{$_};
		
		if ($section{'verbosity'} <= $verbosity){
			my $sectionPageIdx = 0;	
	
			my @sectionfigures = @{$section{'figures'}};
	
			foreach (@sectionfigures) {
				my %figure = %{$_};
				if ($figure{'verbosity'} <= $verbosity){
					my $figFileName = $figure{'fileName'};
					if (-e "$outDir/$figFileName") {
						$pageIdx++;
						$sectionPageIdx++;
						$pages .= "\\clearpage\n";
						if ($pageIdx == 1){
							$pages .= "\\setcounter{page}{1}\n";
						}
						if ($sectionPageIdx == 1){							
							$pages .= "\\section[{".$section{'name'}."}]{}\n";	
							$nSections++;
						}
						$pages .= "\\subsection[{".$figure{'name'}."}]{}\n";
						$nSections++;
						my $viewport = '';						
						$pages .= "\\begin{figure}[h!]%\n";
						$pages .= "\\centering%\n";
						if ($figure{'landscape'}){
							if ($figure{'crop'}){
								$viewport = sprintf("angle=90,viewport=%d %d %d %d,", 
									$figure{'crop'}[2], 
									$figure{'crop'}[1], 
									11 * 72 - $figure{'crop'}[3],
									8.5 * 72 - $figure{'crop'}[0]);
							}							
						}else{
							if ($figure{'crop'}){
								$viewport = sprintf("viewport=%d %d %d %d,", 
									$figure{'crop'}[2], 
									$figure{'crop'}[1], 
									8.5 * 72 - $figure{'crop'}[3],
									11 * 72 - $figure{'crop'}[0]);
							}	
						}
						$pages .= "\\includegraphics[".$viewport."width=0.90\\textwidth,height=0.9\\textheight,keepaspectratio]{$outDir/$figFileName}%\n";
						$pages .= "\\end{figure}\n";							
						$pages .= "\\fancyhf[FC]{\\LARGE ".$figure{'name'}."}\n\n";
					}
				}
			}
		}
	}

	my $template = HTML::Template->new(filename => 'data/summary.tex.tmpl');
	$template->param(title => $title);
	$template->param(pages => $pages);
    
	open(FH, '>', $texFileName) or die $!;
	print FH $template->output;
	close FH;
	
	#cleanup 
	my $texFileNameBase = substr($texFileName, 0, -4);
	`rm -rf $texFileNameBase.log`;
	`rm -rf $texFileNameBase.aux`;
	`rm -rf $texFileNameBase.out`;
	`rm -rf $texFileNameBase.toc`;
	`rm -rf $texFileNameBase.pdf`;

	#run pdflatex
	`pdflatex -output-directory $outDir $texFileName`;
		
	#rerun pdflatex
	`pdflatex -output-directory $outDir $texFileName`;
}

#Source: http://www.perlmonks.org/?node_id=406883
sub max ($$) { $_[$_[0] < $_[1]] }
sub min ($$) { $_[$_[0] > $_[1]] }

1;
