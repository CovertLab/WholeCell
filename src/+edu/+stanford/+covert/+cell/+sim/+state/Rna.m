%Rna
%1. nascent
%2. processed
%3. intergenic segments
%4. modified
%5. bound
%6. misfolded
%7. damaged
%8. aminoacylated
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 1/5/2011
classdef Rna < edu.stanford.covert.cell.sim.MoleculeCountState
    properties (Constant)
        optionNames             = {   %names of properties that are options
            'verbosity';
            'seed';
            };
        fixedConstantNames      = {   %names of process properties that are considered fixed constants
            'molecularWeights';
            'baseCounts';
            'lengths';            
            'compartments'            
            'weightFractionMRNA'
            'weightFractionRRNA5S'
            'weightFractionRRNA16S'
            'weightFractionRRNA23S'
            'weightFractionTRNA'
            'weightFractionSRNA'
            'nascentRNAGeneComposition'
            'nascentRNAMatureRNAComposition'
            'matureRNAGeneComposition'
            'intergenicRNAMatrix'
            'geneExpressionRobustness'
			'minTRnaCnt'
            }
        fittedConstantNames     = {   %names of process properties that are considered fitted constants, and should be stored with the simulation as such
            'expectedGeneExpression'
            'expectedGeneHalfLives'
            'expression'
            'halfLives';
            };
        stateNames              = {   %names of properties which are part of the simulation's state
            'counts'
            };
        dependentStateNames     = {}; %names of properties which can be calculated from the simulation's state
    end
    
    properties (Constant)
        mRNAWeightFractionIndexs    = 1;
        rRNAWeightFractionIndexs    = (2:4)';
        sRNAWeightFractionIndexs    = 5;
        tRNAWeightFractionIndexs    = 6;
        rRNA5SWeightFractionIndexs  = 2;
        rRNA16SWeightFractionIndexs = 3;
        rRNA23SWeightFractionIndexs = 4;
    end
    
    %indices
    properties
        nascentIndexs               %index within RNAs
        processedIndexs             %index within RNAs
        intergenicIndexs            %index within RNAs
        matureIndexs                %index within RNAs
        boundIndexs                 %index within RNAs
        misfoldedIndexs             %index within RNAs
        damagedIndexs               %index within RNAs
        aminoacylatedIndexs         %index within RNAs
        
        matureMRNAIndexs            %mRNAs within mature RNAs
        matureRRNAIndexs            %rRNAs within mature RNAs
        matureSRNAIndexs            %sRNAs within mature RNAs
        matureTRNAIndexs            %tRNAs within mature RNAs
        matureRibosomalRRNAIndexs   %ribosomal RNAs within mature RNAs
        matureTMRNAIndexs           %tmRNAs within mature RNAs
    end
    
    properties
        expression
        expectedGeneExpression      %experimental gene expression
        expectedGeneHalfLives       %experimental gene half lives
        weightFractionMRNA          %mRNA fraction of RNAs, by weight
        weightFractionRRNA5S        %5S rRNA fraction of RNAs, by weight
        weightFractionRRNA16S       %16S rRNA fraction of RNAs, by weight
        weightFractionRRNA23S       %23S rRNA fraction of RNAs, by weight
        weightFractionTRNA          %tRNA fraction of RNAs, by weight
        weightFractionSRNA          %sRNA fraction of RNAs, by weight
        
        nascentRNAGeneComposition      %genes and nascent RNAs (genes X transcription units)
        nascentRNAMatureRNAComposition %processed/mature RNAs and nascent RNAs (processed RNAs X transcription units)
        matureRNAGeneComposition       %genes and processed/mature RNAs (genes X processed/mature RNAs)
        intergenicRNAMatrix            %unprocessed RNAs (transcription units) that are precursors to discarded intergenic RNA segments
        
        geneExpressionRobustness     %how robust to make growth rate against stochastic gene expression        
		minTRnaCnt                   %minimum mean count for each tRNA species
    end
    
    properties (Dependent)
        geneHalfLives               %half lifes of each gene (s)
        geneDecayRates              %decay rates of RNAs of each gene (molecules/s)
        weightFractions             %fractions of RNA weight (mRNA, rRNA 5S, rRNA 16S, rRNA 23S, s/tRNA)
        geneExpression              %mol fractions of genes
        expectedGeneDecayRates      %experimentally determined RNA decay rates (s)
        expectedWeightFractions     %experimentally determined RNA weight fractions
    end
    
    %object references
    properties
        ribosome
    end
    
    %constructor
    methods
        function this = Rna(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.MoleculeCountState(wholeCellModelID, name);
        end
    end
    
    methods
        function storeObjectReferences(this, simulation)
            this.storeObjectReferences@edu.stanford.covert.cell.sim.MoleculeCountState(simulation);
            
            this.ribosome = simulation.state('Ribosome');
        end
    end
    
    methods
        function initializeConstants(this, knowledgeBase, simulation)
            this.initializeConstants@edu.stanford.covert.cell.sim.MoleculeCountState(knowledgeBase, simulation);
            
            %genes
            this.expectedGeneHalfLives  = [knowledgeBase.genes.halfLife]' * 60;
            this.expectedGeneExpression = reshape([knowledgeBase.genes.expression], [], knowledgeBase.numGenes)';
            this.expectedGeneExpression = this.expectedGeneExpression ./ repmat(sum(this.expectedGeneExpression), knowledgeBase.numGenes, 1);
            
            %RNA
            this.nascentRNAGeneComposition      = sum(knowledgeBase.transcriptionUnitComposition_Genes, 3);
            this.nascentRNAMatureRNAComposition = sum(knowledgeBase.transcriptionUnitComposition_ProcessedRNAs, 3);
            this.matureRNAGeneComposition       = sum(knowledgeBase.processedRNAComposition_Genes, 3);
            this.intergenicRNAMatrix            = sum(knowledgeBase.transcriptionUnitComposition_IntergenicRNAs, 3);
                        
            nascentRNAs       = knowledgeBase.nascentRNAs;
            processedRNAs     = knowledgeBase.processedRNAs;
            intergenicRNAs    = knowledgeBase.intergenicRNAs;
            matureRNAs        = knowledgeBase.matureRNAs;
            aminoacylatedRNAs = knowledgeBase.aminoacylatedRNAs;
            
            this.nascentIndexs       = (1:numel(nascentRNAs.wholeCellModelIDs))';
            this.processedIndexs     = (1:numel(processedRNAs.wholeCellModelIDs))'     + this.nascentIndexs(end);
            this.intergenicIndexs    = (1:numel(intergenicRNAs.wholeCellModelIDs))'    + this.processedIndexs(end);
            this.matureIndexs        = (1:numel(matureRNAs.wholeCellModelIDs))'        + this.intergenicIndexs(end);
            this.boundIndexs         = (1:numel(matureRNAs.wholeCellModelIDs))'        + this.matureIndexs(end);
            this.misfoldedIndexs     = (1:numel(matureRNAs.wholeCellModelIDs))'        + this.boundIndexs(end);
            this.damagedIndexs       = (1:numel(matureRNAs.wholeCellModelIDs))'        + this.misfoldedIndexs(end);
            this.aminoacylatedIndexs = (1:numel(aminoacylatedRNAs.wholeCellModelIDs))' + this.damagedIndexs(end);
            
            this.wholeCellModelIDs = [
                nascentRNAs.wholeCellModelIDs;
                processedRNAs.wholeCellModelIDs;
                intergenicRNAs.wholeCellModelIDs;
                repmat(matureRNAs.wholeCellModelIDs, 4, 1);
                aminoacylatedRNAs.wholeCellModelIDs];
            this.expression = [
                zeros(length(this.nascentIndexs), 1);
                zeros(length(this.processedIndexs), 1);
                edu.stanford.covert.math.unit(...
                this.matureRNAGeneComposition' * this.expectedGeneExpression(:, 1), 1);
                zeros(length(this.boundIndexs), 1);
                zeros(length(this.misfoldedIndexs), 1);
                zeros(length(this.damagedIndexs), 1);
                zeros(length(this.aminoacylatedIndexs), 1)];
            this.names = [
                nascentRNAs.names;
                processedRNAs.names;
                intergenicRNAs.names;
                repmat(matureRNAs.names, 4, 1);
                aminoacylatedRNAs.names];
            this.lengths = [
                nascentRNAs.lengths;
                processedRNAs.lengths;
                intergenicRNAs.lengths;
                repmat(matureRNAs.lengths, 4, 1);
                aminoacylatedRNAs.lengths];
            this.baseCounts =[
                nascentRNAs.baseCounts;
                processedRNAs.baseCounts;
                intergenicRNAs.baseCounts;
                repmat(matureRNAs.baseCounts, 4, 1);
                aminoacylatedRNAs.baseCounts];
            this.molecularWeights = [
                nascentRNAs.molecularWeights;
                processedRNAs.molecularWeights;
                intergenicRNAs.molecularWeights;
                repmat(matureRNAs.molecularWeights, 4, 1);
                aminoacylatedRNAs.molecularWeights];
            
            this.matureMRNAIndexs = find(sum(this.matureRNAGeneComposition([knowledgeBase.mRNAGenes.idx], :), 1))';
            this.matureRRNAIndexs = find(sum(this.matureRNAGeneComposition([knowledgeBase.rRNAGenes.idx], :), 1))';
            this.matureSRNAIndexs = find(sum(this.matureRNAGeneComposition([knowledgeBase.sRNAGenes.idx], :), 1))';
            this.matureTRNAIndexs = find(sum(this.matureRNAGeneComposition([knowledgeBase.tRNAGenes.idx], :), 1))';
            this.matureRibosomalRRNAIndexs = zeros(size(knowledgeBase.ribosomalRRNAIndexs));
            for i = 1:length(knowledgeBase.ribosomalRRNAIndexs)
                this.matureRibosomalRRNAIndexs(i) = find(this.matureRNAGeneComposition(knowledgeBase.ribosomalRRNAIndexs(i), :));
            end
            this.matureTMRNAIndexs = this.getIndexs('MG_0004');
        end
    end
    
    methods
        function notUpdatedRnas = updateExternalState(this, deltaRnas, rnaIsDegraded) %#ok<INUSD>
            r = this.ribosome;
            
            notUpdatedRnas = zeros(size(deltaRnas));
            deltaMatureRnas = deltaRnas(this.matureIndexs, this.compartment.cytosolIndexs);
            deltaBoundRnas = deltaRnas(this.boundIndexs, this.compartment.cytosolIndexs);
            
            %mRNA or tmRNA
            if any(deltaMatureRnas(this.matureMRNAIndexs)) || deltaBoundRnas(this.matureTMRNAIndexs) < 0
                r.releaseRibosome(0, -deltaBoundRnas(this.matureTMRNAIndexs), -deltaMatureRnas(this.matureMRNAIndexs));
            end
        end
    end
    
    %helper methods
    methods
        function value = getIndexs(this, wholeCellModelIDs)
            [~, value] = ismember(wholeCellModelIDs, this.wholeCellModelIDs(this.matureIndexs));
        end
    end
    
    %getters
    methods
        function value = get.geneHalfLives(this)
            value = this.matureRNAGeneComposition * ...
                this.halfLives(this.matureIndexs);
        end
        
        %decay rates of RNAs (moleculeces/s)
        function value = get.geneDecayRates(this)
            value = log(2) ./ this.geneHalfLives;
        end
        
        %experimental decay rates
        function value = get.expectedGeneDecayRates(this)
            value = log(2) ./ this.expectedGeneHalfLives;
        end
        
        %RNA weight fractions
        function value = get.weightFractions(this)
            value = [
                this.molecularWeights(this.matureIndexs(this.matureMRNAIndexs))'          * this.expression(this.matureIndexs(this.matureMRNAIndexs));
                this.molecularWeights(this.matureIndexs(this.matureRibosomalRRNAIndexs)) .* this.expression(this.matureIndexs(this.matureRibosomalRRNAIndexs));
                this.molecularWeights(this.matureIndexs(this.matureSRNAIndexs))'          * this.expression(this.matureIndexs(this.matureSRNAIndexs));
                this.molecularWeights(this.matureIndexs(this.matureTRNAIndexs))'          * this.expression(this.matureIndexs(this.matureTRNAIndexs));
                ];
            value = value / sum(value);
        end
        
        %experimental RNA weight fractions
        function value = get.expectedWeightFractions(this)
            value = [
                this.weightFractionMRNA
                this.weightFractionRRNA5S
                this.weightFractionRRNA16S
                this.weightFractionRRNA23S
                this.weightFractionSRNA
                this.weightFractionTRNA];
            value = value / sum(value);
        end
        
        %calculated expression of RNAs, monomers,
        function value = get.geneExpression(this)
            value = this.matureRNAGeneComposition * this.expression(this.matureIndexs);
            value = value / sum(value);
        end
    end
end
