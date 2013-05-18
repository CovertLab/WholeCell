% Instantiate simulation from modified knowledge base.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Bonny Jain, bonny.jain@mit.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 12/20/2012
classdef SingleGeneOptimized < edu.stanford.covert.cell.sim.runners.SimulationRunner
    methods
        function this = SingleGeneOptimized(varargin)
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
                580177, 1083, true, {''}, meanHL, expression);
            
            tuA = TranscriptionUnit(kb, NaN, {'TuA'}, {'Transcription unit A'},...
                -35, 6, ...
                -10, 6, ...
                -1);
            
            monA = ProteinMonomer(kb, NaN, {'MonomerA'}, {'Protein monomer A'},...
                {''}, {''}, {''}, ...
                10, {'dsDNA'}, {'dsDNA'}, ...
                {''}, {''}, {''}, {''}, {''}, {''}, ...
                false, {''}, {''}, 0, ...
                {''});
            
            cpxAA = ProteinComplex(kb, NaN, {'ComplexAA'}, {'Macromolecular complex AA'},...
                10, {'dsDNA'}, {'dsDNA'}, ...
                {''}, {''}, {''}, {''}, {''}, {''}, ...
                {''}, {''});

            genome = kb.genome;
            promoter = 'TATCTACAGAGGCCCGCTTTAAGCATCCACGACCCTACTCACTTCAAAGTGGAACCACCCGTCGACGTGTGTCTTAACCCCTGCCGCTGCAAGGTGTGAG';
            lacI = 'ATGAAACCTGTTACTTTATATGATGTTGCTGAATATGCTGGTGTTTCATATCAAACTGTTTCAAGAGTTGTTAATCAAGCTTCACATGTTTCAGCTAAAACTAGAGAAAAAGTTGAAGCTGCTATGGCTGAATTAAATTATATTCCTAATAGAGTTGCTCAACAATTAGCTGGTAAACAATCATTATTAATTGGTGTTGCTACTTCATCATTAGCTTTACATGCTCCTTCACAAATTGTTGCTGCTATTAAATCAAGAGCTGATCAATTAGGTGCTTCAGTTGTTGTTTCAATGGTTGAAAGATCAGGTGTTGAAGCTTGTAAAGCTGCTGTTCATAATTTATTAGCTCAAAGAGTTTCAGGTTTAATTATTAATTATCCTTTAGATGATCAAGATGCTATTGCTGTTGAAGCTGCTTGTACTAATGTTCCTGCTTTATTTTTAGATGTTTCAGATCAAACTCCTATTAATTCAATTATTTTTTCACATGAAGATGGTACTAGATTAGGTGTTGAACATTTAGTTGCTTTAGGTCATCAACAAATTGCTTTATTAGCTGGTCCTTTATCATCAGTTTCAGCTAGATTAAGATTAGCTGGTTGGCATAAATATTTAACTAGAAATCAAATTCAACCTATTGCTGAAAGAGAAGGTGATTGGTCAGCTATGTCAGGTTTTCAACAAACTATGCAAATGTTAAATGAAGGTATTGTTCCTACTGCTATGTTAGTTGCTAATGATCAAATGGCTTTAGGTGCTATGAGAGCTATTACTGAATCAGGTTTAAGAGTTGGTGCTGATATTTCAGTTGTTGGTTATGATGATACTGAAGATTCATCATGTTATATTCCTCCTTCAACTACTATTAAACAAGATTTTAGATTATTAGGTCAAACTTCAGTTGATAGATTATTACAATTATCACAAGGTCAAGCTGTTAAAGGTAATCAATTATTACCTGTTTCATTAGTTAAAAGAAAAACTACTTTAGCTCCTAATACTCAAACTGCTTCACCTAGAGCTTTAGCTGATTCATTAATGCAATTAGCTAGACAAGTTTCAAGATTAGAATCAGGTCAATAA';
            genome.sequence = [genome.sequence promoter lacI];
            genome.sequenceLength = length(genome.sequence);
            
            cytosol = kb.compartments(kb.cytosolCompartmentIndexs);
            extracellularSpace = kb.compartments(kb.extracellularCompartmentIndexs);
            
            %% create new relationships among new genes, transcription units, etc.
            geneA.genome = genome;
            geneA.transcriptionUnits = tuA;
            geneA.proteinMonomers = monA;
            geneA.compartment = cytosol;
            
            tuA.genome = genome;
            tuA.genes = geneA;
            tuA.geneCompartments = cytosol;
            tuA.compartment = cytosol;
            
            monA.gene = geneA;
            monA.geneCompartments = cytosol;
            monA.compartment = cytosol;

            cpxAA.proteinMonomers = monA;
            cpxAA.proteinMonomerCompartments = cytosol;
            cpxAA.proteinMonomerCoefficients = 4;
            cpxAA.compartment = cytosol;
            cpxAA.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');
            
            genome.genes = [genome.genes; geneA];
            genome.transcriptionUnits = [genome.transcriptionUnits; tuA];
            
            kb.genes = genome.genes;
            kb.transcriptionUnits = genome.transcriptionUnits;
            kb.proteinMonomers = [kb.proteinMonomers; monA];
            kb.proteinComplexs = [kb.proteinComplexs; cpxAA];
                        
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
            [~, modTuIndexs] = ismember({'TuA'}, rna.wholeCellModelIDs(rna.nascentIndexs));
            
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