#!/usr/bin/perl

#Usage: ./compileSummary.pl 259 2011_06_23_01_01_01
#Author: Jonathan Karr, jkarr@stanford.edu
#Author: Derek Macklin, macklin@stanford.edu
#Affiliation: Covert Lab, Department of Bioengineering, Stanford University
#Last updated: 7/19/2011

use Cwd;
use Mysql;
use strict;
use Switch;
require "library.pl";

my %config = getConfiguration();
my $URL = 'http://'.$config{'URL'};

my $link = Mysql->connect($config{'hostName'}, $config{'schema'}, $config{'dbUserName'}, $config{'dbPassword'});
my $result;
my @data;
$result = $link->query("
	SELECT DataSet.WID
	FROM DataSet
	JOIN Entry ON Entry.OtherWID = DataSet.WID
	WHERE DataSet.HomeURL = '$URL'
	ORDER BY DataSet.Version DESC, Entry.InsertDate DESC, DataSet.WID DESC
	LIMIT 1");
@data = $result->fetchrow();
my $kbWID = $data[0];

my $simulationIdx = $ARGV[0];
my $simulationDateTime = $ARGV[1];
my @tmp = split('_', $simulationDateTime, 6);
my $timeStamp = $tmp[1].'/'.$tmp[2].'/'.$tmp[0].' '.
		$tmp[3].':'.$tmp[4].':'.$tmp[5];

my $baseDir = $config{'simulationPath'};
my $outDir;
my $texFileName;

#simulations
$outDir = "$baseDir/output/runSimulation/$simulationDateTime";
$texFileName = "$outDir/summary.tex";

my @figures = (
	{verbosity => 0, 'name' => 'Single Cell Dynamics', figures => [
		{verbosity => 0, name => 'Growth', fileName => 'singleCell-Growth.pdf', crop => [201, 377, 78, 20]},
		{verbosity => 0, name => 'Cellular Composition', fileName => 'singleCell-Mass.pdf', crop => [202, 200, 62, 57]},
		{verbosity => 0, name => 'Cellular Composition', fileName => 'singleCell-MassNormalized.pdf', crop => [202, 200, 62, 57]},
		{verbosity => 0, name => 'Nucleotides', fileName => 'singleCell-Metabolites.pdf', crop => [202, 200, 62, 57]},
		{verbosity => 0, name => 'Nucleotides', fileName => 'singleCell-MetabolitesNormalized.pdf', crop => [202, 200, 62, 57]},
		{verbosity => 0, name => 'Translation', fileName => 'singleCell-Ribosomes.pdf', crop => [76, 54, 202, 212], landscape => 1},
		{verbosity => 0, name => 'Translation', fileName => 'singleCell-RibosomesGenomicCoordinates.pdf', crop => [208, 202, 52, 55]},
		{verbosity => 0, name => 'Translation', fileName => 'TranslationAnalysis-MonomerSynthesis.pdf', crop => [204, 192, 49, 64]},
		{verbosity => 0, name => 'Stalled Ribosomes', fileName => 'TranslationAnalysis-Pauses.pdf', crop => [190, 195, 36, 72]},
		{verbosity => 0, name => 'Stalled Ribosomes', fileName => 'TranslationAnalysis-PauseDistribution.pdf', crop => [192, 198, 61, 64]},
		{verbosity => 0, name => 'Chromosome Organization', fileName => 'ChromosomeSpaceTimePlot-SpaceTimeOverlay.pdf', crop => [184, 197, 52, 53]},
		{verbosity => 0, name => 'Chromosome Organization', fileName => 'ChromosomeSpaceTimePlot-CircularDensity.pdf', crop => [172, 172, 34, 32]},
		{verbosity => 0, name => 'Chromosome Organization', fileName => 'ChromosomeSpaceTimePlot-Ring.pdf', crop => [208, 221, 89, 51]},
		{verbosity => 0, name => 'Chromosome Organization', fileName => 'CellState-OnChromosome.pdf', crop => [231, 225, 10, 12]},
		{verbosity => 0, name => 'Chromosome Organization', fileName => 'ChromosomePositionHistogram-1.pdf', crop => [198, 204, 41, 49]},
		{verbosity => 0, name => 'Chromosome Organization', fileName => 'ChromosomePositionHistogram-2.pdf', crop => [198, 204, 41, 49]},
		]},
	{verbosity => 0, 'name' => 'Population Dynamics', figures => [
		{verbosity => 0, name => 'Growth', fileName => 'summary-CellOverview-Lines1.pdf', crop => [199, 200, 39, 53]},
		{verbosity => 0, name => 'Growth', fileName => 'population-Growth.pdf', crop => [197, 279, 28, 210]},
		{verbosity => 0, name => 'Cellular Composition', fileName => 'population-MassDistribution.pdf', crop => [51, 56, 201, 201], landscape=>1},
		{verbosity => 0, name => 'Metabolism', fileName => 'population-Metabolism.pdf', crop => [209, 208, 50, 68]},
		{verbosity => 0, name => 'Energy Production', fileName => 'population-EnergyProduction.pdf', crop => [287, 292, 51, 31]},
		{verbosity => 0, name => 'Energy Production', fileName => 'summary-CellOverview-Lines4.pdf', crop => [202, 201, 31, 54]},
		{verbosity => 0, name => 'Nucleotide Concentrations', fileName => 'population-NucleotideDistribution.pdf', crop => [200, 206, 23, 31]},
		{verbosity => 0, name => 'Amino Acid Counts', fileName => 'population-AminoAcidCounts.pdf', crop => [209, 206, 36, 32]},
		{verbosity => 0, name => 'Process ATP Requirements', fileName => 'processMetaboliteUsage-ATP-Requirements.pdf', crop => [176, 174, 21, 26]},
		{verbosity => 0, name => 'Process ATP Allocations', fileName => 'processMetaboliteUsage-ATP-Allocations.pdf', crop => [176, 174, 21, 26]},
		{verbosity => 0, name => 'Process ATP Usages', fileName => 'processMetaboliteUsage-ATP-Usages.pdf', crop => [176, 174, 21, 26]},
		{verbosity => 0, name => 'Process GTP Requirements', fileName => 'processMetaboliteUsage-GTP-Requirements.pdf', crop => [176, 174, 21, 26]},
		{verbosity => 0, name => 'Process GTP Allocations', fileName => 'processMetaboliteUsage-GTP-Allocations.pdf', crop => [176, 174, 21, 26]},
		{verbosity => 0, name => 'Process GTP Usages', fileName => 'processMetaboliteUsage-GTP-Usages.pdf', crop => [176, 174, 21, 26]},
		{verbosity => 0, name => 'Replication Initiation', fileName => 'summary-CellOverview-Lines2.pdf', crop => [200, 201, 40, 53]},
		{verbosity => 0, name => 'Replication', fileName => 'summary-CellOverview-Lines3.pdf', crop => [196, 199, 31, 52]},
		{verbosity => 0, name => 'Replication', fileName => 'population-SSBs.pdf', crop => [47, 53, 203, 210], landscape=>1},
		{verbosity => 0, name => 'Replication', fileName => 'population-GeneCopyNumber.pdf', crop => [202, 202, 37, 28]},
		{verbosity => 0, name => 'Transcription', fileName => 'population-RnaPolymerases.pdf', crop => [199, 203, 43, 31]},
		{verbosity => 0, name => 'RNA Maturation', fileName => 'population-RnaSynthesis.pdf', crop => [200, 195, 29, 21]},
		{verbosity => 0, name => 'RNA Maturation', fileName => 'summary-CellOverview-Lines5.pdf', crop => [199, 202, 43, 52]},
		{verbosity => 0, name => 'Translation', fileName => 'population-Translation.pdf', crop => [210, 197, 12, 32]},
		{verbosity => 0, name => 'Translation', fileName => 'population-Ribosomes.pdf', crop => [195, 205, 41, 33]},
		{verbosity => 0, name => 'Protein Maturation', fileName => 'population-ProteinSynthesis.pdf', crop => [206, 196, 26, 20]},
		{verbosity => 0, name => 'Protein Maturation', fileName => 'summary-CellOverview-Lines6.pdf', crop => [200, 202, 39, 53]},
		{verbosity => 0, name => 'Supercoiling', fileName => 'summary-CellOverview-Lines8.pdf', crop => [202, 200, 43, 54]},
		{verbosity => 0, name => 'DNA Repair', fileName => 'population-DnaRepair.pdf', crop => [25, 26, 205, 201], landscape=>1},
		{verbosity => 0, name => 'Cell Division', fileName => 'summary-CellOverview-Lines7.pdf', crop => [198, 202, 41, 53]},
		{verbosity => 0, name => 'Cell Division', fileName => 'population-CellShape.pdf', crop => [38, 52, 200, 200], landscape=>1},
		{verbosity => 0, name => 'Immune Activation', fileName => 'population-ImmuneActivation.pdf', crop => [40, 23, 204, 196], landscape=>1},
		{verbosity => 0, name => 'Unsynchronized Population', fileName => 'population-UnsynchronizedPopulation.pdf', crop => [39, 53, 201, 207], landscape=>1},
		{verbosity => 0, name => 'Simulation Runtime', fileName => 'summary-CellOverview-Lines9.pdf', crop => [198, 199, 54, 53]},
		]},
	{verbosity => 0, 'name' => 'Population Statistics', figures => [
		{verbosity => 0, name => 'Biomass Composition', fileName => 'BiomassCompositionProduction.pdf', crop => [212, 173, 43, 89]},
		{verbosity => 0, name => 'Biomass Composition', fileName => 'BiomassCompositionProduction-WeightFractions.pdf', crop => [200, 232, 117, 145]},
		{verbosity => 0, name => 'Metabolite Concentrations', fileName => 'population-MetaboliteConcentrations.pdf', crop => [186, 65, 186, 51], landscape=>0},
		{verbosity => 0, name => 'dNMP Composition', fileName => 'BiomassCompositionProduction-dNMPComposition.pdf', crop => [201, 216, 112, 85]},
		{verbosity => 0, name => 'Nucleotide Composition', fileName => 'BiomassCompositionProduction-NMPComposition.pdf', crop => [201, 216, 112, 85]},
		{verbosity => 0, name => 'Nucleotide Composition', fileName => 'Constants-Experimental_Vs_Calculated_NMP_Composition.pdf', crop => [209, 195, 94, 115]},
		{verbosity => 1, name => 'Nucleotide Composition', fileName => 'Constants-Experimental_Calculated_NMP_Composition_Ratios.pdf', crop => [208, 196, 47, 68]},
		{verbosity => 0, name => 'Amino Acid Composition', fileName => 'BiomassCompositionProduction-AAComposition.pdf', crop => [201, 216, 112, 85]},
		{verbosity => 0, name => 'Amino Acid Composition', fileName => 'Constants-Experimental_Vs_Calculated_AA_Composition.pdf', crop => [209, 195, 94, 115]},
		{verbosity => 1, name => 'Amino Acid Composition', fileName => 'Constants-Experimental_Calculated_AA_Composition_Ratios.pdf', crop => [208, 196, 47, 68]},
		{verbosity => 0, name => 'RNA Composition', fileName => 'Constants-Experimental_Vs_Calculated_RNA_Weight_Fractions.pdf', crop => [209, 195, 94, 115]},
		{verbosity => 1, name => 'RNA Composition', fileName => 'Constants-Experimental_Calculated_RNA_Weight_Fraction_Ratios.pdf', crop => [208, 196, 47, 68]},
		{verbosity => 0, name => 'Gene Expression', fileName => 'Constants-Experimental_Vs_Calculated_Gene_Expression.pdf', crop => [209, 195, 94, 115]},
		{verbosity => 1, name => 'Gene Expression', fileName => 'Constants-Experimental_Calculated_Gene_Expression_Ratios.pdf', crop => [208, 196, 47, 68]},
		{verbosity => 0, name => 'tRNA Expression', fileName => 'Constants-Experimental_Vs_Calculated_TRNA_Expression.pdf', crop => [209, 195, 94, 115]},
		{verbosity => 1, name => 'tRNA Expression', fileName => 'Constants-Experimental_Calculated_TRNA_Expression_Ratios.pdf', crop => [208, 196, 47, 68]},
		{verbosity => 0, name => 'Protein Expression', fileName => 'Constants-Experimental_Vs_Calculated_Monomer_Expression.pdf', crop => [209, 195, 94, 115]},
		{verbosity => 1, name => 'Protein Expression', fileName => 'Constants-Experimental_Calculated_Monomer_Expression_Ratios.pdf', crop => [208, 196, 47, 68]},
		{verbosity => 0, name => 'Macromolecule Expression', fileName => 'population-MacromoleculeExpression.pdf', crop => [203, 223, 36, 29]},
		{verbosity => 0, name => 'Macromolecular Complexes', fileName => 'population-MacromolecularComplexes.pdf', crop => [212, 196, 72, 56]},
		{verbosity => 0, name => 'RNA Half Lives', fileName => 'Constants-Experimental_Vs_Calculated_Gene_Decay_Rates.pdf', crop => [209, 195, 94, 115]},
		{verbosity => 1, name => 'RNA Half Lives', fileName => 'Constants-Experimental_Calculated_Gene_Decay_Rate_Ratios.pdf', crop => [208, 196, 47, 68]},
		{verbosity => 0, name => 'Replication', fileName => 'population-Replication.pdf', crop => [202, 215, 69, 71]},
		{verbosity => 0, name => 'Secondary Replication Initiation', fileName => 'population-SecondaryReplicationInitiation.pdf', crop => [210, 205, 27, 30]},
		{verbosity => 0, name => 'RNA Synthesis', fileName => 'population-RnaSynthesisDuration.pdf', crop => [207, 204, 54, 70]},
		{verbosity => 0, name => 'Protein Synthesis', fileName => 'population-ProteinSynthesisDuration.pdf', crop => [207, 204, 54, 70]},
		{verbosity => 0, name => 'DNA Bound Protein Displacement', fileName => 'population-DnaBoundProteinDisplacement.pdf'},
		{verbosity => 0, name => 'Cell Cycle Phase Durations', fileName => 'cellCyclePhaseDurations.pdf', crop => [207, 200, 68, 66]},
		{verbosity => 0, name => 'Cell Cycle Phase Distribution', fileName => 'population-CellCyclePhases.pdf', crop => [198, 379, 47, 402]},
		{verbosity => 0, name => 'Cell Cycle Phase Lengths', fileName => 'cellCyclePhaseLengths.pdf', crop => [207, 200, 68, 66]},
		{verbosity => 0, name => 'Cumulative Growth Vs. Cell Cycle Phase Durations', fileName => 'cumulativeGrowthVsCellCyclePhaseDurations.pdf', crop => [207, 200, 68, 66]},
		{verbosity => 0, name => 'Cumulative Growth Vs. Cell Cycle Phase Times', fileName => 'cumulativeGrowthVsCellCyclePhaseTimes.pdf', crop => [207, 200, 68, 66]},
		{verbosity => 0, name => 'Initial Growth Vs. Cell Cycle Phase Durations', fileName => 'initialGrowthVsCellCyclePhaseDurations.pdf', crop => [207, 200, 68, 66]},
		{verbosity => 0, name => 'Initial Growth Vs. Cell Cycle Phase Times', fileName => 'initialGrowthVsCellCyclePhaseTimes.pdf', crop => [207, 200, 68, 66]},
		{verbosity => 0, name => 'Initial Vs. Cumulative Growth', fileName => 'initialVsCumulativeGrowth.pdf', crop => [207, 200, 68, 66]},
		{verbosity => 0, name => 'Initial Growth Rate Vs. End Time', fileName => 'summary-CellOverview-InitialGrowthRateVsEndTimes.pdf', crop => [190, 192, 40, 65]},
		{verbosity => 0, name => 'Final Growth Rate Vs. End Time', fileName => 'summary-CellOverview-FinalGrowthRateVsEndTimes.pdf', crop => [190, 192, 40, 65]},
		{verbosity => 0, name => 'Cell Mass Distribution', fileName => 'population-CellMassDistribution.pdf', crop => [206, 200, 60, 63] },
		{verbosity => 0, name => 'Blocked Decay Events', fileName => 'population-BlockedDecayEvents.pdf', crop => [202, 204, 45, 68]},
		{verbosity => 0, name => 'Process ATP Usage', fileName => 'processMetaboliteUsage-Expected-ATP-Usages.pdf', crop => [180, 223, 47, 28]},
		{verbosity => 0, name => 'Process GTP Usage', fileName => 'processMetaboliteUsage-Expected-GTP-Usages.pdf', crop => [180, 223, 47, 28]},
		{verbosity => 0, name => 'Warnings', fileName => 'population-Warnings.pdf', crop => [207, 194, 12, 62]},
		]},
	{verbosity => 0, 'name' => 'Simulation Structure', figures => [
		{verbosity => 0, name => 'Metabolic Model', fileName => 'FBA-NetworkReduction.pdf', crop => [183, 227, 48, 60]},
		{verbosity => 0, name => 'Shared Gene Products', fileName => 'SimulationStructure-ProcessGeneProducts.pdf', crop => [209, 193, 43, 29]},
		{verbosity => 1, name => 'Shared Metabolites', fileName => 'SimulationStructure-ProcessMetaboliteSharing.pdf', crop => [196, 178, 42, 61]},
		{verbosity => 0, name => 'Shared Metabolites', fileName => 'SimulationStructure-ProcessMetabolites.pdf', crop => [196, 180, 33, 23]},
		{verbosity => 0, name => 'Shared Metabolites', fileName => 'SimulationStructure-ProcessSharedMetabolites.pdf', crop => [196, 180, 33, 23]},
		]},
);

my $title = "Whole Cell Simulation \\#$simulationIdx \\\\$timeStamp Summary";
compileFigureSummary($title, \@figures, $outDir, $texFileName, 1);

#single-gene deletions
$outDir = "$baseDir/output/runSimulation/singleGeneDeletions";
$texFileName = "$outDir/summary.tex";
my @figures = (
	{verbosity => 0, 'name' => 'Wild-Type Strain', figures => [
		{verbosity => 0, name => 'Wild-Type Strain', fileName => 'WT.pdf', crop => [198, 204, 50, 32]},
		]},
	{verbosity => 0, 'name' => 'Single Gene Deletion Strains', figures => [
		{verbosity => 0, name => 'Summary', fileName => 'overview.pdf'},
		{verbosity => 0, name => 'Summary', fileName => 'summaryGrid.pdf'},
		{verbosity => 0, name => 'Summary', fileName => 'all.pdf', crop => [207, 192, 43, 30]},
		{verbosity => 0, name => 'Mass', fileName => 'all-01.pdf', crop => [207, 192, 43, 30]},
		{verbosity => 0, name => 'Growth', fileName => 'all-02.pdf', crop => [77, 32, 194, 207], landscape=>1},
		{verbosity => 0, name => 'ATP Usage', fileName => 'all-03.pdf', crop => [207, 192, 43, 30]},
		{verbosity => 0, name => 'GTP Usage', fileName => 'all-04.pdf', crop => [207, 192, 43, 30]},
		{verbosity => 0, name => 'Chromosome Copy Number', fileName => 'all-05.pdf', crop => [207, 192, 43, 30]},
		{verbosity => 0, name => 'Superhelicity', fileName => 'all-06.pdf', crop => [207, 192, 43, 30]},
		{verbosity => 0, name => 'RNA', fileName => 'all-07.pdf', crop => [207, 192, 43, 30]},
		{verbosity => 0, name => 'Protein', fileName => 'all-08.pdf', crop => [207, 192, 43, 30]},
		{verbosity => 0, name => 'dNTP', fileName => 'all-09.pdf', crop => [207, 192, 43, 30]},
		{verbosity => 0, name => 'NTP', fileName => 'all-10.pdf', crop => [207, 192, 43, 30]},
		{verbosity => 0, name => 'Amino Acid', fileName => 'all-11.pdf', crop => [207, 192, 43, 30]},
		{verbosity => 0, name => 'Damaged RNA', fileName => 'all-12.pdf', crop => [207, 192, 43, 30]},
		{verbosity => 0, name => 'Damaged Protein', fileName => 'all-13.pdf', crop => [207, 192, 43, 30]},
		{verbosity => 0, name => 'Antibiotics', fileName => 'all-14.pdf', crop => [207, 192, 43, 30]},
		{verbosity => 0, name => 'Replication Initiation Duration', fileName => 'all-15.pdf', crop => [46, 32, 194, 209], landscape=>1},
		{verbosity => 0, name => 'Replication Duration', fileName => 'all-16.pdf', crop => [46, 32, 194, 209], landscape=>1},
		{verbosity => 0, name => 'Cytokinesis Duration', fileName => 'all-17.pdf', crop => [46, 32, 194, 209], landscape=>1},
		{verbosity => 0, name => 'Mass Doubling Duration', fileName => 'all-18.pdf', crop => [207, 193, 30, 30]},
		{verbosity => 0, name => 'Cell Cycle Length', fileName => 'all-19.pdf', crop => [46, 32, 194, 209], landscape=>1},
		]},
	{verbosity => 0, 'name' => 'Single Gene Deletion Strain Classes', figures => [
		{verbosity => 0, name => 'Non-Essential', fileName => 'class-1.pdf', crop => [198, 204, 50, 32]},
		{verbosity => 0, name => 'Degrading', fileName => 'class-2.pdf', crop => [198, 204, 50, 32]},
		{verbosity => 0, name => 'Non-Growing', fileName => 'class-3.pdf', crop => [198, 204, 50, 32]},
		{verbosity => 0, name => 'Decying Growth - Non-RNA Synthesizing', fileName => 'class-4.pdf', crop => [198, 204, 50, 32]},
		{verbosity => 0, name => 'Decying Growth - Non-RNA, Non-Protein Synthesizing', fileName => 'class-5.pdf', crop => [198, 204, 50, 32]},
		{verbosity => 0, name => 'Non-Replicative', fileName => 'class-6.pdf', crop => [198, 204, 50, 32]},
		{verbosity => 0, name => 'Non-Fissive', fileName => 'class-7.pdf', crop => [198, 204, 50, 32]},
		{verbosity => 0, name => 'Slow Growing', fileName => 'class-8.pdf', crop => [198, 204, 50, 32]},
		{verbosity => 0, name => 'Toxin Accumulation', fileName => 'class-9.pdf', crop => [198, 204, 50, 32]},
		{verbosity => 0, name => 'Non-perpetuating', fileName => 'class-10.pdf', crop => [198, 204, 50, 32]},
		{verbosity => 0, name => 'No Terminal Organelle', fileName => 'class-11.pdf', crop => [198, 204, 50, 32]},
		]},
	{verbosity => 0, 'name' => 'Single Gene Deletion Strain Classification', figures => [	
		{verbosity => 0, name => 'Deletion Strain Growth Rate Distribution', fileName => 'deletionGrowthRateDistribution.pdf', crop => [195, 190, 52, 68]},
		{verbosity => 0, name => 'Degrading, Non-Growing, Slow-Growing, Non-Essential', fileName => 'classification-1.pdf', crop => [195, 190, 52, 68]},
		{verbosity => 0, name => 'Decaying Growth', fileName => 'classification-2.pdf', crop => [195, 190, 52, 68]},
		{verbosity => 0, name => 'Non-Replicative', fileName => 'classification-3.pdf', crop => [195, 190, 52, 68]},
		{verbosity => 0, name => 'Non-Fissive', fileName => 'classification-4.pdf', crop => [195, 190, 52, 68]},
		{verbosity => 0, name => 'Toxin Accumulation', fileName => 'classification-5.pdf', crop => [195, 190, 52, 68]},
		]},
	{verbosity => 0, 'name' => 'Model / Experiment Comparison', figures => [
		{verbosity => 0, name => 'Doubling Time vs. Doubling Time', fileName => 'ModelDoublingTimeVsExperimentDoublingTime.pdf', crop => [206, 194, 57, 64]},
		{verbosity => 0, name => 'Cell Cycle Length vs. Doubling Time', fileName => 'ModelCellCycleLengthVsExperimentDoublingTime.pdf', crop => [206, 194, 57, 64]},		
		]},
	{verbosity => 0, 'name' => 'Specific Individual Single Gene Deletion Strains', figures => [
		]},
	{verbosity => 0, 'name' => 'Individual Single Gene Deletion Strains', figures => [
		]},
);

opendir(DIR, $outDir) || die("Cannot open directory");
my @files = readdir(DIR);
closedir(DIR);

foreach my $file (sort(@files)) {
	unless ($file eq "." || $file eq ".." || substr($file, 0, 2) ne 'MG' || substr($file, -4) ne '.pdf'){
		if (substr($file, -13) eq '-analysis.pdf'){
			my $geneID = substr($file, 0, -13);

			$result = $link->query("SELECT Gene.Name FROM Gene JOIN DBID ON DBID.OtherWID = Gene.WID WHERE DBID.XID='$geneID' && Gene.DataSetWID=$kbWID");
			@data = $result->fetchrow();

			$geneID =~ s/_/\\_/g;
			my $figure = {verbosity => 0, name => ($data[0] ? (length($data[0]) > 40 ? substr($data[0], 0, 40) ."{\\ldots}" : $data[0])." ($geneID)" : $geneID), fileName => $file, crop => [206, 203, 85, 27]};
			push(@{$figures[5]{'figures'}}, $figure);
		}else{
			my $geneID = substr($file, 0, -4);

			$result = $link->query("SELECT Gene.Name FROM Gene JOIN DBID ON DBID.OtherWID = Gene.WID WHERE DBID.XID='$geneID' && Gene.DataSetWID=$kbWID");
			@data = $result->fetchrow();

			$geneID =~ s/_/\\_/g;
			my $figure = {verbosity => 0, name => ($data[0] ? (length($data[0]) > 40 ? substr($data[0], 0, 40) ."{\\ldots}" : $data[0])." ($geneID)" : $geneID), fileName => $file, crop => [198, 204, 50, 32]};
			push(@{$figures[6]{'figures'}}, $figure);
		}
	}
}

my $title = "Whole Cell Simulation \\Single Gene Deletions Summary";
compileFigureSummary($title, \@figures, $outDir, $texFileName, 1);
