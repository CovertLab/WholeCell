% Instantiate simulation from modified knowledge base.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Bonny Jain, bonny.jain@mit.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 12/20/2012
classdef Repressilator < edu.stanford.covert.cell.sim.runners.SimulationRunner
    methods
        function this = Repressilator(varargin)
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
            
            %% create new genes, transcription units, etc.
            meanHL = mean([kb.mRNAGenes.halfLife]);
            expression = zeros(1, 3); %normal, cold shock, heat shock
            geneA = Gene(kb, NaN, {'GeneA'}, {'Gene A'}, {'genA'}, {''}, {'mRNA'}, 0, {''}, ...
                580500, 300, true, {''}, meanHL, expression);
            geneB = Gene(kb, NaN, {'GeneB'}, {'Gene B'}, {'genB'}, {''}, {'mRNA'}, 0, {''}, ...
                581000, 300, true, {''}, meanHL, expression);
            geneC = Gene(kb, NaN, {'GeneC'}, {'Gene C'}, {'genC'}, {''}, {'mRNA'}, 0, {''}, ...
                581500, 300, true, {''}, meanHL, expression);
            
            tuA = TranscriptionUnit(kb, NaN, {'TuA'}, {'Transcription unit A'},...
                -35, 6, ...
                -10, 6, ...
                -1);
            tuB = TranscriptionUnit(kb, NaN, {'TuB'}, {'Transcription unit B'},...
                -35, 6, ...
                -10, 6, ...
                -1);
            tuC = TranscriptionUnit(kb, NaN, {'TuC'}, {'Transcription unit C'},...
                -35, 6, ...
                -10, 6, ...
                -1);
            
            monA = ProteinMonomer(kb, NaN, {'MonomerA'}, {'Protein monomer A'},...
                {''}, {''}, {''}, ...
                10, {'dsDNA'}, {'dsDNA'}, ...
                {''}, {''}, {''}, {''}, {''}, {''}, ...
                false, {''}, {''}, 0, ...
                {''});
            monB = ProteinMonomer(kb, NaN, {'MonomerB'}, {'Protein monomer B'},...
                {''}, {''}, {''}, ...
                10, {'dsDNA'}, {'dsDNA'}, ...
                {''}, {''}, {''}, {''}, {''}, {''}, ...
                false, {''}, {''}, 0, ...
                {''});
            monC = ProteinMonomer(kb, NaN, {'MonomerC'}, {'Protein monomer C'},...
                {''}, {''}, {''}, ...
                10, {'dsDNA'}, {'dsDNA'}, ...
                {''}, {''}, {''}, {''}, {''}, {''}, ...
                false, {''}, {''}, 0, ...
                {''});
            
            cpxAA = ProteinComplex(kb, NaN, {'ComplexAA'}, {'Macromolecular complex AA'},...
                10, {'dsDNA'}, {'dsDNA'}, ...
                {''}, {''}, {''}, {''}, {''}, {''}, ...
                {'!StimulusA'}, {''});
            
            cpxBB = ProteinComplex(kb, NaN, {'ComplexBB'}, {'Macromolecular complex BB'},...
                10, {'dsDNA'}, {'dsDNA'}, ...
                {''}, {''}, {''}, {''}, {''}, {''}, ...
                {'!StimulusB'}, {''});
            
            cpxCC = ProteinComplex(kb, NaN, {'ComplexCC'}, {'Macromolecular complex CC'},...
                10, {'dsDNA'}, {'dsDNA'}, ...
                {''}, {''}, {''}, {''}, {''}, {''}, ...
                {''}, {''});
            
            stimA = Stimuli(kb, NaN, {'StimulusA'}, {'Stimulus A'});
            stimB = Stimuli(kb, NaN, {'StimulusB'}, {'Stimulus B'});
            
            this.createDefaultDisplacementReactions(kb, [cpxAA; cpxBB; cpxCC]);
            
            acgt = 'ACGT';
            genome = kb.genome;
            randSeq = ceil(rand(1, 582500-genome.sequenceLength) * 4);
            genome.sequence = [genome.sequence  acgt(randSeq)];
            genome.sequence(geneA.startCoordinate : geneA.endCoordinate) = acgt(ceil(rand(1, geneA.endCoordinate - geneA.startCoordinate + 1) * 4));
            genome.sequence(geneB.startCoordinate : geneB.endCoordinate) = acgt(ceil(rand(1, geneB.endCoordinate - geneB.startCoordinate + 1) * 4));
            genome.sequence(geneC.startCoordinate : geneC.endCoordinate) = acgt(ceil(rand(1, geneC.endCoordinate - geneC.startCoordinate + 1) * 4));
            genome.sequence(geneA.startCoordinate+(0:2)) = 'ATG';
            genome.sequence(geneB.startCoordinate+(0:2)) = 'ATG';
            genome.sequence(geneC.startCoordinate+(0:2)) = 'ATG';
            genome.sequence(geneA.endCoordinate+(-2:0)) = 'TAA';
            genome.sequence(geneB.endCoordinate+(-2:0)) = 'TAA';
            genome.sequence(geneC.endCoordinate+(-2:0)) = 'TAA';
            
            cytosol = kb.compartments(kb.cytosolCompartmentIndexs);
            extracellularSpace = kb.compartments(kb.extracellularCompartmentIndexs);
            
            %% create new relationships among new genes, transcription units, etc.
            geneA.genome = genome;
            geneA.transcriptionUnits = tuA;
            geneA.proteinMonomers = monA;
            geneA.compartment = cytosol;
            
            geneB.genome = genome;
            geneB.transcriptionUnits = tuB;
            geneB.proteinMonomers = monB;
            geneB.compartment = cytosol;
            
            geneC.genome = genome;
            geneC.transcriptionUnits = tuC;
            geneC.proteinMonomers = monC;
            geneC.compartment = cytosol;
            
            tuA.genome = genome;
            tuA.genes = geneA;
            tuA.geneCompartments = cytosol;
            tuA.compartment = cytosol;
            tuA.transcriptionFactorProteinComplexs = cpxBB;
            tuA.transcriptionFactorProteinComplexAffinitys = NaN; %nM
            tuA.transcriptionFactorProteinComplexActivitys = 1/10;
            tuA.transcriptionFactorProteinComplexConditions = {''};
            tuA.transcriptionFactorProteinComplexBindingSiteStartCoordinates = 580460;
            tuA.transcriptionFactorProteinComplexBindingSiteLengths = 10;
            tuA.transcriptionFactorProteinComplexBindingSiteDirections = true;
            tuA.transcriptionFactorProteinComplexCompartments = cytosol;
            
            tuB.genome = genome;
            tuB.genes = geneB;
            tuB.geneCompartments = cytosol;
            tuB.compartment = cytosol;
            tuB.transcriptionFactorProteinComplexs = cpxCC;
            tuB.transcriptionFactorProteinComplexAffinitys = NaN; %nM
            tuB.transcriptionFactorProteinComplexActivitys = 1/10;
            tuB.transcriptionFactorProteinComplexConditions = {''};
            tuB.transcriptionFactorProteinComplexBindingSiteStartCoordinates = 580960;
            tuB.transcriptionFactorProteinComplexBindingSiteLengths = 10;
            tuB.transcriptionFactorProteinComplexBindingSiteDirections = true;
            tuB.transcriptionFactorProteinComplexCompartments = cytosol;
            
            tuC.genome = genome;
            tuC.genes = geneC;
            tuC.geneCompartments = cytosol;
            tuC.compartment = cytosol;
            tuC.transcriptionFactorProteinComplexs = cpxAA;
            tuC.transcriptionFactorProteinComplexAffinitys = NaN; %nM
            tuC.transcriptionFactorProteinComplexActivitys = 1/10;
            tuC.transcriptionFactorProteinComplexConditions = {''};
            tuC.transcriptionFactorProteinComplexBindingSiteStartCoordinates = 581460;
            tuC.transcriptionFactorProteinComplexBindingSiteLengths = 10;
            tuC.transcriptionFactorProteinComplexBindingSiteDirections = true;
            tuC.transcriptionFactorProteinComplexCompartments = cytosol;
            
            monA.gene = geneA;
            monA.geneCompartments = cytosol;
            monA.compartment = cytosol;
            
            monB.gene = geneB;
            monB.geneCompartments = cytosol;
            monB.compartment = cytosol;
            
            monC.gene = geneC;
            monC.geneCompartments = cytosol;
            monC.compartment = cytosol;
            
            cpxAA.proteinMonomers = monA;
            cpxAA.proteinMonomerCompartments = cytosol;
            cpxAA.proteinMonomerCoefficients = 2;
            cpxAA.regulatedTranscriptionUnits = tuC;
            cpxAA.compartment = cytosol;
            cpxAA.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');
            cpxAA.stimuliRegulators = stimA;
            
            cpxBB.proteinMonomers = monB;
            cpxBB.proteinMonomerCompartments = cytosol;
            cpxBB.proteinMonomerCoefficients = 2;
            cpxBB.regulatedTranscriptionUnits = tuA;
            cpxBB.compartment = cytosol;
            cpxBB.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');
            cpxBB.stimuliRegulators = stimB;
            
            cpxCC.proteinMonomers = monC;
            cpxCC.proteinMonomerCompartments = cytosol;
            cpxCC.proteinMonomerCoefficients = 2;
            cpxCC.regulatedTranscriptionUnits = tuB;
            cpxCC.compartment = cytosol;
            cpxCC.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');
            
            %specificy dynamics of stimuli piece-wise
            stimA.compartments = [extracellularSpace; extracellularSpace];
            stimA.values = [true; false];
            stimA.initialTimes = [0; 10];
            stimA.finalTimes = [10; inf];
            stimA.regulatedProteinComplexs = cpxAA;
            
            stimB.compartments = extracellularSpace;
            stimB.values = true;
            stimB.initialTimes = 0;
            stimB.finalTimes = inf;
            stimB.regulatedProteinComplexs = cpxBB;
            
            genome.genes = [genome.genes; geneA; geneB; geneC];
            genome.transcriptionUnits = [genome.transcriptionUnits; tuA; tuB; tuC];
            
            kb.genes = genome.genes;
            kb.transcriptionUnits = genome.transcriptionUnits;
            kb.proteinMonomers = [kb.proteinMonomers; monA; monB; monC];
            kb.proteinComplexs = [kb.proteinComplexs; cpxAA; cpxBB; cpxCC];
            kb.stimulis = [kb.stimulis; stimA; stimB];
            
            %% class super class method
            this.modifyNetworkStructure@edu.stanford.covert.cell.sim.runners.SimulationRunner(kb);
        end
        
        function modifyNetworkParameters(~, sim)
            %% get handles
            g = sim.gene;
            time = sim.state('Time');
            rna = sim.state('Rna');
            trn = sim.process('Transcription');
            
            %% get constants
            nascentMRNAIndexs = find(any(rna.nascentRNAGeneComposition(g.mRNAIndexs, :), 1));
            [~, modTuIndexs] = ismember({'TuA'; 'TuB'; 'TuC'}, rna.wholeCellModelIDs(rna.nascentIndexs));
            
            %% get parameter values
            tuBindProb = trn.transcriptionUnitBindingProbabilities;
            meanBindProb = mean(tuBindProb(setdiff(nascentMRNAIndexs, modTuIndexs)));
            rnaDecayRates = rna.decayRates(rna.matureIndexs);
            
            %% calculate modified parameter values
            %set binding probability of new TUs
            modTuBindProb = tuBindProb;
            modTuBindProb(modTuIndexs) = meanBindProb;
            
            %renormalize
            modTuBindProb(nascentMRNAIndexs) = modTuBindProb(nascentMRNAIndexs) * sum(tuBindProb(nascentMRNAIndexs)) / sum(modTuBindProb(nascentMRNAIndexs));
            
            %update expression
            modRnaExp = (rna.nascentRNAMatureRNAComposition * modTuBindProb) ./ (log(2) / time.cellCycleLength + rnaDecayRates);
            modRnaExp = modRnaExp / sum(modRnaExp);
            
            %% update parameter values
            trn.transcriptionUnitBindingProbabilities = modTuBindProb;
            rna.expression(rna.matureIndexs) = modRnaExp;
        end
    end
end