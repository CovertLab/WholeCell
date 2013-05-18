classdef SinglearaCOptimized1 < edu.stanford.covert.cell.sim.runners.SimulationRunner
	methods
		function this = SinglearaCOptimized1(varargin)
			this = this@edu.stanford.covert.cell.sim.runners.SimulationRunner(varargin{:});
		end
	end

	methods (Access = protected)
		function modifyNetworkStructure(this, kb)
			import edu.stanford.covert.cell.kb.Gene;
			import edu.stanford.covert.cell.kb.ProteinMonomer;
			import edu.stanford.covert.cell.kb.ProteinComplex;
			import edu.stanford.covert.cell.kb.Stimuli;
			import edu.stanford.covert.cell.kb.TranscriptionUnit;
			meanHL = mean([kb.mRNAGenes.halfLife]);
			expression = zeros(1, 3);
			geneA1 = Gene(kb, NaN, {'GeneA1'}, {'Gene A1'}, {'genA1'}, {''}, {'mRNA'}, 0, {''}, 580177, 879, true, {''}, meanHL, expression);
			tuA1 = TranscriptionUnit(kb, NaN, {'TuA1'}, {'Transcription unit A1'}, -35, 6, -35, 6, -1);
			monA1 = ProteinMonomer(kb, NaN, {'MonomerA1'}, {'Protein monomer A1'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			cpxAA1 = ProteinComplex(kb, NaN, {'ComplexAA1'}, {'Macromolecular complex AA1'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			genome = kb.genome;
			promoter = 'TATCTACAGAGGCCCGCTTTAAGCATCCACGACCCTACTCACTTCAAAGTGGAACCACCCGTCGACGTGTGTCTTAACCCCTGCCGCTGCAAGGTGTGAG';
			gene = 'ATGGCTGAAGCTCAAAATGATCCTTTATTACCTGGTTATTCATTTAATGCTCATTTAGTTGCTGGTTTAACTCCTATTGAAGCTAATGGTTATTTAGATTTTTTTATTGATAGACCTTTAGGTATGAAAGGTTATATTTTAAATTTAACTATTAGAGGTCAAGGTGTTGTTAAAAATCAAGGTAGAGAATTTGTTTGTAGACCTGGTGATATTTTATTATTTCCTCCTGGTGAAATTCATCATTATGGTAGACATCCTGAAGCTAGAGAATGGTATCATCAATGGGTTTATTTTAGACCTAGAGCTTATTGGCATGAATGGTTAAATTGGCCTTCAATTTTTGCTAATACTGGTTTTTTTAGACCTGATGAAGCTCATCAACCTCATTTTTCAGATTTATTTGGTCAAATTATTAATGCTGGTCAAGGTGAAGGTAGATATTCAGAATTATTAGCTATTAATTTATTAGAACAATTATTATTAAGAAGAATGGAAGCTATTAATGAATCATTACATCCTCCTATGGATAATAGAGTTAGAGAAGCTTGTCAATATATTTCAGATCATTTAGCTGATTCAAATTTTGATATTGCTTCAGTTGCTCAACATGTTTGTTTATCACCTTCAAGATTATCACATTTATTTAGACAACAATTAGGTATTTCAGTTTTATCATGGAGAGAAGATCAAAGAATTTCACAAGCTAAATTATTATTATCAACTACTAGAATGCCTATTGCTACTGTTGGTAGAAATGTTGGTTTTGATGATCAATTATATTTTTCAAGAGTTTTTAAAAAATGTACTGGTGCTTCACCTTCAGAATTTAGAGCTGGTTGTGAAGAAAAAGTTAATGATGTTGCTGTTAAATTATCATAA';
			genome.sequence = [genome.sequence repmat([promoter gene], [1, 1])];
			genome.sequenceLength = length(genome.sequence);
			cytosol = kb.compartments(kb.cytosolCompartmentIndexs);
			extracellularSpace = kb.compartments(kb.extracellularCompartmentIndexs);
			geneA1.genome = genome;
			geneA1.transcriptionUnits = tuA1;
			geneA1.proteinMonomers = monA1;
			geneA1.compartment = cytosol;

			tuA1.genome = genome;
			tuA1.genes = geneA1;
			tuA1.geneCompartments = cytosol;
			tuA1.compartment = cytosol;

			monA1.gene = geneA1;
			monA1.geneCompartments = cytosol;
			monA1.compartment = cytosol;

			cpxAA1.proteinMonomers = monA1;
			cpxAA1.proteinMonomerCompartments = cytosol;
			cpxAA1.proteinMonomerCoefficients = 4;
			cpxAA1.compartment = cytosol;
			cpxAA1.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			genome.genes = [genome.genes; geneA1];
			genome.transcriptionUnits = [genome.transcriptionUnits; tuA1];
			kb.genes = genome.genes;
			kb.transcriptionUnits = genome.transcriptionUnits;
			kb.proteinMonomers = [kb.proteinMonomers; monA1];
			kb.proteinComplexs = [kb.proteinComplexs; cpxAA1];
			this.modifyNetworkStructure@edu.stanford.covert.cell.sim.runners.SimulationRunner(kb);

		end		function modifyNetworkParameters(~, sim)
			g = sim.gene;
			time = sim.state('Time');
			rna = sim.state('Rna');
			trn = sim.process('Transcription');
			nascentMRNAIndexs = find(any(rna.nascentRNAGeneComposition(g.mRNAIndexs, :), 1));
			[~, modTuIndexs] = ismember({'TuA1'}, rna.wholeCellModelIDs(rna.nascentIndexs));
			tuBindProb = trn.transcriptionUnitBindingProbabilities;
			meanBindProb = mean(tuBindProb(setdiff(nascentMRNAIndexs, modTuIndexs)));
			rnaDecayRates = rna.decayRates(rna.matureIndexs);
			modTuBindProb = tuBindProb;
			modTuBindProb(modTuIndexs) = meanBindProb;
			modTuBindProb(nascentMRNAIndexs) = modTuBindProb(nascentMRNAIndexs) * sum(tuBindProb(nascentMRNAIndexs)) / sum(modTuBindProb(nascentMRNAIndexs));
			modRnaExp = (rna.nascentRNAMatureRNAComposition * modTuBindProb) ./ (log(2) / time.cellCycleLength + rnaDecayRates);
			modRnaExp = modRnaExp / sum(modRnaExp);
			trn.transcriptionUnitBindingProbabilities = modTuBindProb;
			rna.expression(rna.matureIndexs) = modRnaExp;
		end
	end
end