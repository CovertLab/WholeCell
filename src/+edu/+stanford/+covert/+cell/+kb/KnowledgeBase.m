% Defines a knowledge base
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 9/17/2009
classdef KnowledgeBase < edu.stanford.covert.cell.kb.KnowledgeBaseObject
    properties
        version
        investigator
        taxonomy
        translationTable

        parameters         = edu.stanford.covert.cell.kb.Parameter.empty(0,0);
        processes          = edu.stanford.covert.cell.kb.Process.empty(0,0);
        states             = edu.stanford.covert.cell.kb.State.empty(0,0);
        compartments       = edu.stanford.covert.cell.kb.Compartment.empty(0,0);
        metabolites        = edu.stanford.covert.cell.kb.Metabolite.empty(0,0);
        genome             = edu.stanford.covert.cell.kb.Genome.empty(0,0);
        genes              = edu.stanford.covert.cell.kb.Gene.empty(0,0);
        transcriptionUnits = edu.stanford.covert.cell.kb.TranscriptionUnit.empty(0,0);
        genomeFeatures     = edu.stanford.covert.cell.kb.GenomeFeature.empty(0,0);
        proteinMonomers    = edu.stanford.covert.cell.kb.ProteinMonomer.empty(0,0);
        proteinComplexs    = edu.stanford.covert.cell.kb.ProteinComplex.empty(0,0);
        reactions          = edu.stanford.covert.cell.kb.Reaction.empty(0,0);
        pathways           = edu.stanford.covert.cell.kb.Pathway.empty(0,0);
        stimulis           = edu.stanford.covert.cell.kb.Stimuli.empty(0,0);
        notes              = edu.stanford.covert.cell.kb.Note.empty(0,0);

        maxAvgMetConc
        maxAvgRxnFlux
    end

    properties (Constant = true)
        cytosolCompartmentWholeCellModelIDs                   = {'c'};
        chromosomeCompartmentWholeCellModelIDs                = {'d'};
        membraneCompartmentWholeCellModelIDs                  = {'m'};
        terminalOrganelleCytosolCompartmentWholeCellModelIDs  = {'tc'};
        terminalOrganelleMembraneCompartmentWholeCellModelIDs = {'tm'};
        extracellularCompartmentWholeCellModelIDs             = {'e'};
    end

    %computed and cached properties
    properties %(SetAccess = protected)
        numProcesses
        numStates
        numParameters
        numCompartments
        numMetabolites
        numGenes
        numTranscriptionUnits
        numGenomeFeatures
        numProteinMonomers
        numProteinComplexs
        numReactions
        numPathways
        numStimulis
        numNotes
        numReferences

        cytosolCompartmentIndexs
        chromosomeCompartmentIndexs
        membraneCompartmentIndexs
        terminalOrganelleCytosolCompartmentIndexs
        terminalOrganelleMembraneCompartmentIndexs
        extracellularCompartmentIndexs
        cellularCompartmentIndexs

        waterIndexs
        hydrogenIndexs
        dnmpIndexs
        nmpIndexs
        modifiedNMPIndexs
        aminoAcidIndexs
        cysteineIndexs
        methionineIndexs
        fmethionineIndexs
        modifiedAminoAcidIndexs
        diacylglycerolCysteineIndexs

        mRNAGenes
        rRNAGenes
        sRNAGenes
        tRNAGenes
        
        ribosomalRRNAIndexs
        
        metaboliteMolecularWeights
        
        biomassComposition
        
        geneComposition
        probRNAPolBinding
        
        geneticCodeCodons
        geneticCodeTRNAs
        geneticCodeAminoAcids
        geneticCodeStartCodons
        geneticCodeStartTRNAs
        geneticCodeStartAminoAcids
        geneticCodeTRNAAminoacylation
        
        nascentRNAs
        intergenicRNAs
        processedRNAs
        matureRNAs
        aminoacylatedRNAs
        transcriptionUnitComposition_Genes
        transcriptionUnitComposition_IntergenicRNAs
        transcriptionUnitComposition_ProcessedRNAs
        processedRNAComposition_Genes
        
        proteinMonomerTRNASequences
        proteinMonomerTRNACounts
        proteinMonomerCompartments
        proteinMonomerComposition
        proteinMonomerNTerminalMethionineCleavages
        proteinMonomerSignalSequenceBaseCounts
        proteinMonomerSignalSequenceMolecularWeights
        proteinProstheticGroupMatrix
        proteinChaperoneMatrix
        proteinComplexRNAComposition
        proteinComplexMonomerComposition
        proteinComplexComplexComposition
        proteinComplexAllRNAComposition
        proteinComplexAllMonomerComposition
        proteinComplexCompartments

        reactionStimuliStoichiometryMatrix
        reactionMetaboliteStoichiometryMatrix
        reactionRNAStoichiometryMatrix
        reactionProteinMonomerStoichiometryMatrix
        reactionProteinComplexStoichiometryMatrix
        reactionModifiedProteinStoichiometryMatrix
        reactionModificationMatrix
        reactionBounds
        enzymeBounds
        reactionCatalysisMatrix
        reactionCoenzymeMatrix

        maxGeneLength
        maxTranscriptionUnitLength
        maxGeneExpression
        maxTranscriptionUnitExpression
    end
    
    properties %(SetAccess = protected)
        areLinksSerialized = true
    end

    methods
        function this = KnowledgeBase(database, wid)
            if ~exist('wid', 'var')
                wid = edu.stanford.covert.cell.kb.KnowledgeBaseUtil.selectLatestKnowledgeBase(database);
            end
            this.idx = 1;            
            this.wid = wid;
            this.knowledgeBase = this;
            this.invalidate();
            this.loadFromKnowledgeBase(database);
            this.computeDependentProperties();
        end

        %% biomass composition
        %metabolite composition of biomass
        function value = get.biomassComposition(this)
            %retrieve
            if ~isempty(this.biomassComposition)
                value = this.biomassComposition;
                return;
            end
            
            %calculate
            value = zeros(this.numMetabolites,this.numCompartments);
            for i = 1:this.numMetabolites
                if ~isempty(this.metabolites(i).biomassCompartments)
                    value(i,[this.metabolites(i).biomassCompartments.idx]) = ...
                        this.metabolites(i).biomassCoefficients;
                end
            end
            value = -value;
            
            %store
            this.biomassComposition = value;
        end

        %% genes
        %metabolite composition of genes
        function value = get.geneComposition(this)
            %retrieve
            if ~isempty(this.geneComposition)
                value = this.geneComposition;
                return;
            end
            
            %calculate
            value = reshape([this.genes.baseCount],[],length(this.genes))';
            
            %store
            this.geneComposition = value;
        end       

        %% transcription
        %probability of RNA polymerase binding to transcription units
        function value = get.probRNAPolBinding(this)
            %retrieve
            if ~isempty(this.probRNAPolBinding)
                value = this.probRNAPolBinding;
                return;
            end
            
            %calculate
            synthesisRates = reshape([this.transcriptionUnits.synthesisRate], [], this.numTranscriptionUnits)';
            value = synthesisRates ./ repmat(sum(synthesisRates), this.numTranscriptionUnits, 1);
            
            %store
            this.probRNAPolBinding = value;
        end

        %% genetic code
        function value = get.geneticCodeCodons(this)
            %retrieve
            if ~isempty(this.geneticCodeCodons)
                value = this.geneticCodeCodons;
                return;
            end
            
            %calculate
            value = {};
            for i = 1:length(this.tRNAGenes)
                if this.tRNAGenes(i).startCodon; continue; end;
                value = [value this.tRNAGenes(i).codons]; %#ok<AGROW>
            end
            
            %store
            this.geneticCodeCodons = value;
        end

        function value = get.geneticCodeTRNAs(this)
            %retrieve
            if ~isempty(this.geneticCodeTRNAs)
                value = this.geneticCodeTRNAs;
                return;
            end
            
            %calculate
            value=[];
            for i=1:length(this.tRNAGenes)
                if this.tRNAGenes(i).startCodon; continue; end;
                value = [value; repmat(i,length(this.tRNAGenes(i).codons), 1)]; %#ok<AGROW>
            end
            
            %store
            this.geneticCodeTRNAs = value;
        end

        function value = get.geneticCodeAminoAcids(this)
            %retrieve
            if ~isempty(this.geneticCodeAminoAcids)
                value = this.geneticCodeAminoAcids;
                return;
            end
            
            %calculate
            value={};
            for i=1:length(this.tRNAGenes)
                if this.tRNAGenes(i).startCodon; continue; end;
                for j=1:length(this.tRNAGenes(i).codons)
                    value{end + 1} = this.tRNAGenes(i).aminoAcid.wholeCellModelID; %#ok<AGROW>
                end
            end
            
            %store
            this.geneticCodeAminoAcids = value;
        end

        function value = get.geneticCodeStartCodons(this)
            %retrieve
            if ~isempty(this.geneticCodeStartCodons)
                value = this.geneticCodeStartCodons;
                return;
            end
            
            %calculate
            value={};
            for i=1:length(this.tRNAGenes)
                if ~this.tRNAGenes(i).startCodon; continue; end;
                value = [value this.tRNAGenes(i).codons]; %#ok<AGROW>
            end
            
            %store
            this.geneticCodeStartCodons = value;
        end

        function value = get.geneticCodeStartTRNAs(this)
            %retrieve
            if ~isempty(this.geneticCodeStartTRNAs)
                value = this.geneticCodeStartTRNAs;
                return;
            end
            
            %calculate
            value=[];
            for i=1:length(this.tRNAGenes)
                if ~this.tRNAGenes(i).startCodon; continue; end;
                value = [value; repmat(i, length(this.tRNAGenes(i).codons), 1)]; %#ok<AGROW>
            end
            
            %store
            this.geneticCodeStartTRNAs = value;
        end

        function value = get.geneticCodeStartAminoAcids(this)
            %retrieve
            if ~isempty(this.geneticCodeStartAminoAcids)
                value = this.geneticCodeStartAminoAcids;
                return;
            end
            
            %calculate
            value={};
            for i=1:length(this.tRNAGenes)
                if ~this.tRNAGenes(i).startCodon; continue; end;
                for j=1:length(this.tRNAGenes(i).codons)
                    value{end + 1} = this.tRNAGenes(i).aminoAcid.wholeCellModelID; %#ok<AGROW>
                end
            end
            
            %store
            this.geneticCodeStartAminoAcids = value;
        end

        function value = get.geneticCodeTRNAAminoacylation(this)
            %retrieve
            if ~isempty(this.geneticCodeTRNAAminoacylation)
                value = this.geneticCodeTRNAAminoacylation;
                return;
            end
            
            %calculate
            value = zeros(length(this.aminoAcidIndexs), length(this.tRNAGenes), this.numCompartments);
            for i = 1:length(this.tRNAGenes)
                value(this.tRNAGenes(i).aminoAcid.idx == this.aminoAcidIndexs, i, this.cytosolCompartmentIndexs) = 1;
            end
            
            %store
            this.geneticCodeTRNAAminoacylation = value;
        end
        
        %% RNAs - nascent, intergenic, processed, mature, aminoacylated
        function value = get.nascentRNAs(this)
            %retrieve
            if ~isempty(this.nascentRNAs)
                value = this.nascentRNAs;
                return;
            end
            
            %calculate
            value = struct();
            value.wholeCellModelIDs = {this.transcriptionUnits.wholeCellModelID}';
            value.names             = {this.transcriptionUnits.name}';
            value.sequences         = {this.transcriptionUnits.sequence}';
            value.lengths           = [this.transcriptionUnits.sequenceLength]';
            value.baseCounts        = reshape([this.transcriptionUnits.baseCount], [], this.numTranscriptionUnits)';
            value.decayReactions    = reshape([this.transcriptionUnits.decayReaction], [], this.numTranscriptionUnits)';
            value.molecularWeights  = [this.transcriptionUnits.molecularWeight]';
            value.halfLives         = [this.transcriptionUnits.halfLife]';
            
            %store
            this.nascentRNAs = value;
        end
        
        function value = get.intergenicRNAs(this)
            %retrieve
            if ~isempty(this.intergenicRNAs)
                value = this.intergenicRNAs;
                return;
            end
            
            %calculate
            value = struct();
            value.wholeCellModelIDs = cell(0, 1);
            value.names             = cell(0, 1);
            value.sequences         = cell(0, 1);
            value.lengths           = zeros(0, 1);
            value.baseCounts        = zeros(0, this.numMetabolites);
            value.decayReactions    = zeros(0, this.numMetabolites);
            value.molecularWeights  = zeros(0, 1);
            value.halfLives         = zeros(0, 1);
            
            for i=1:this.numTranscriptionUnits
                numIntergenicSequences = numel(this.transcriptionUnits(i).intergenicSequences);
                
                value.wholeCellModelIDs = [value.wholeCellModelIDs; repmat({this.transcriptionUnits(i).wholeCellModelID}, numIntergenicSequences, 1)];
                value.names             = [value.names;             repmat({this.transcriptionUnits(i).name}, numIntergenicSequences, 1)];
                value.sequences         = [value.sequences;         this.transcriptionUnits(i).intergenicSequences];
                value.lengths           = [value.lengths;           this.transcriptionUnits(i).intergenicSequenceLengths];
                value.baseCounts        = [value.baseCounts;        this.transcriptionUnits(i).intergenicSequenceBaseCounts];
                value.decayReactions    = [value.decayReactions;    this.transcriptionUnits(i).intergenicSequenceDecayReactions];
                value.molecularWeights  = [value.molecularWeights;  this.transcriptionUnits(i).intergenicSequenceMolecularWeights];
                value.halfLives         = [value.halfLives;         this.transcriptionUnits(i).intergenicSequenceHalfLives];
            end
            
            %store
            this.intergenicRNAs = value;
        end
        
        function value = get.processedRNAs(this)
            %retrieve
            if ~isempty(this.processedRNAs)
                value = this.processedRNAs;
                return;
            end
            
            %calculate
            value = struct();
            value.wholeCellModelIDs = cell(0, 1);
            value.names             = cell(0, 1);
            value.sequences         = cell(0, 1);
            value.lengths           = zeros(0, 1);
            value.baseCounts        = zeros(0, this.numMetabolites);
            value.decayReactions    = zeros(0, this.numMetabolites);
            value.molecularWeights  = zeros(0, 1);
            value.halfLives         = zeros(0, 1);
            
            for i=1:this.numTranscriptionUnits
                if strcmp(this.transcriptionUnits(i).type,'mRNA')
                    value.wholeCellModelIDs = [value.wholeCellModelIDs; this.transcriptionUnits(i).wholeCellModelID];
                    value.names             = [value.names;             this.transcriptionUnits(i).name];
                    value.sequences         = [value.sequences;         this.transcriptionUnits(i).sequence];
                    value.lengths           = [value.lengths;           this.transcriptionUnits(i).sequenceLength];
                    value.baseCounts        = [value.baseCounts;        this.transcriptionUnits(i).baseCount];
                    value.decayReactions    = [value.decayReactions;    this.transcriptionUnits(i).decayReaction];
                    value.molecularWeights  = [value.molecularWeights;  this.transcriptionUnits(i).molecularWeight];
                    value.halfLives         = [value.halfLives;         this.transcriptionUnits(i).halfLife];
                else
                    numProcessedSequences = numel(this.transcriptionUnits(i).genes);
                    
                    value.wholeCellModelIDs = [value.wholeCellModelIDs; {this.transcriptionUnits(i).genes.wholeCellModelID}'];
                    value.names             = [value.names;             {this.transcriptionUnits(i).genes.name}'];
                    value.sequences         = [value.sequences;         {this.transcriptionUnits(i).genes.sequence}'];
                    value.lengths           = [value.lengths;           [this.transcriptionUnits(i).genes.sequenceLength]'];
                    value.baseCounts        = [value.baseCounts;        reshape([this.transcriptionUnits(i).genes.baseCount], [], numProcessedSequences)'];
                    value.decayReactions    = [value.decayReactions;    reshape([this.transcriptionUnits(i).genes.decayReaction], [], numProcessedSequences)'];
                    value.molecularWeights  = [value.molecularWeights;  [this.transcriptionUnits(i).genes.molecularWeight]'];
                    value.halfLives         = [value.halfLives;         [this.transcriptionUnits(i).genes.halfLife]'];
                end
            end
            
            %store
            this.processedRNAs = value;
        end
        
        function value = get.matureRNAs(this)
            %retrieve
            if ~isempty(this.matureRNAs)
                value = this.matureRNAs;
                return;
            end
            
            %calculate
            value = struct();
            value.wholeCellModelIDs = cell(0, 1);
            value.names             = cell(0, 1);
            value.sequences         = cell(0, 1);
            value.lengths           = zeros(0, 1);
            value.baseCounts        = zeros(0, this.numMetabolites);
            value.decayReactions    = zeros(0, this.numMetabolites);
            value.molecularWeights  = zeros(0, 1);
            value.halfLives         = zeros(0, 1);
            
            for i=1:this.numTranscriptionUnits
                if strcmp(this.transcriptionUnits(i).type,'mRNA')
                    value.wholeCellModelIDs = [value.wholeCellModelIDs; this.transcriptionUnits(i).wholeCellModelID];
                    value.names             = [value.names;             this.transcriptionUnits(i).name];
                    value.sequences         = [value.sequences;         this.transcriptionUnits(i).sequence];
                    value.lengths           = [value.lengths;           this.transcriptionUnits(i).sequenceLength];
                    value.baseCounts        = [value.baseCounts;        this.transcriptionUnits(i).baseCount];
                    value.decayReactions    = [value.decayReactions;    this.transcriptionUnits(i).decayReaction];
                    value.molecularWeights  = [value.molecularWeights;  this.transcriptionUnits(i).molecularWeight];
                    value.halfLives         = [value.halfLives;         this.transcriptionUnits(i).halfLife];
                else
                    numProcessedSequences = numel(this.transcriptionUnits(i).genes);
                    
                    value.wholeCellModelIDs = [value.wholeCellModelIDs; {this.transcriptionUnits(i).genes.wholeCellModelID}'];
                    value.names             = [value.names;             {this.transcriptionUnits(i).genes.name}'];
                    value.sequences         = [value.sequences;         {this.transcriptionUnits(i).genes.sequence}'];
                    value.lengths           = [value.lengths;           [this.transcriptionUnits(i).genes.sequenceLength]'];
                    value.baseCounts        = [value.baseCounts;        reshape([this.transcriptionUnits(i).genes.matureBaseCount], [], numProcessedSequences)'];
                    value.decayReactions    = [value.decayReactions;    reshape([this.transcriptionUnits(i).genes.matureDecayReaction], [], numProcessedSequences)'];
                    value.molecularWeights  = [value.molecularWeights;  [this.transcriptionUnits(i).genes.matureMolecularWeight]'];
                    value.halfLives         = [value.halfLives;         [this.transcriptionUnits(i).genes.halfLife]'];
                end
            end
            
            %store
            this.matureRNAs = value;
        end
        
        function value = get.aminoacylatedRNAs(this)
            %retrieve
            if ~isempty(this.aminoacylatedRNAs)
                value = this.aminoacylatedRNAs;
                return;
            end
            
            %calculate
            value = struct();
            value.wholeCellModelIDs = cell(0, 1);
            value.names             = cell(0, 1);
            value.sequences         = cell(0, 1);
            value.lengths           = zeros(0, 1);
            value.baseCounts        = zeros(0, this.numMetabolites);
            value.decayReactions    = zeros(0, this.numMetabolites);
            value.molecularWeights  = zeros(0, 1);
            value.halfLives         = zeros(0, 1);
            
            for i=1:this.numTranscriptionUnits
                if strcmp(this.transcriptionUnits(i).type,'mRNA')
                    value.wholeCellModelIDs = [value.wholeCellModelIDs; this.transcriptionUnits(i).wholeCellModelID];
                    value.names             = [value.names;             this.transcriptionUnits(i).name];
                    value.sequences         = [value.sequences;         this.transcriptionUnits(i).sequence];
                    value.lengths           = [value.lengths;           this.transcriptionUnits(i).sequenceLength];
                    value.baseCounts        = [value.baseCounts;        this.transcriptionUnits(i).baseCount];
                    value.decayReactions    = [value.decayReactions;    this.transcriptionUnits(i).decayReaction];
                    value.molecularWeights  = [value.molecularWeights;  this.transcriptionUnits(i).molecularWeight];
                    value.halfLives         = [value.halfLives;         this.transcriptionUnits(i).halfLife];
                else
                    numProcessedSequences = numel(this.transcriptionUnits(i).genes);
                    
                    value.wholeCellModelIDs = [value.wholeCellModelIDs; {this.transcriptionUnits(i).genes.wholeCellModelID}'];
                    value.names             = [value.names;             {this.transcriptionUnits(i).genes.name}'];
                    value.sequences         = [value.sequences;         {this.transcriptionUnits(i).genes.sequence}'];
                    value.lengths           = [value.lengths;           [this.transcriptionUnits(i).genes.sequenceLength]'];
                    value.baseCounts        = [value.baseCounts;        reshape([this.transcriptionUnits(i).genes.aminoacylatedBaseCount], [], numProcessedSequences)'];
                    value.decayReactions    = [value.decayReactions;    reshape([this.transcriptionUnits(i).genes.aminoacylatedDecayReaction], [], numProcessedSequences)'];
                    value.molecularWeights  = [value.molecularWeights;  [this.transcriptionUnits(i).genes.aminoacylatedMolecularWeight]'];
                    value.halfLives         = [value.halfLives;         [this.transcriptionUnits(i).genes.halfLife]'];
                end
            end
            
            %store
            this.aminoacylatedRNAs = value;
        end
        
        %transcription unit composition
        function value = get.transcriptionUnitComposition_Genes(this)
            %retrieve
            if ~isempty(this.transcriptionUnitComposition_Genes)
                value = this.transcriptionUnitComposition_Genes;
                return;
            end
            
            %calculate
            value = zeros(this.numGenes,this.numTranscriptionUnits,this.numCompartments);
            for i = 1:length(this.transcriptionUnits)
                genes = [this.transcriptionUnits(i).genes];
                compartments = [this.transcriptionUnits(i).geneCompartments];
                for j=1:length(genes)
                    value(genes(j).idx, i, compartments(j).idx)=1;
                end
            end
            
            %store
            this.transcriptionUnitComposition_Genes = value;
        end
        
        function value = get.transcriptionUnitComposition_IntergenicRNAs(this)
            %retrieve
            if ~isempty(this.transcriptionUnitComposition_IntergenicRNAs)
                value = this.transcriptionUnitComposition_IntergenicRNAs;
                return;
            end
            
            %calculate
            value = zeros(0, this.numTranscriptionUnits, this.numCompartments);
            
            for i=1:this.numTranscriptionUnits
                numIntergenicSequences = numel(this.transcriptionUnits(i).intergenicSequences);
                
                value = [value; zeros(numIntergenicSequences, this.numTranscriptionUnits, this.numCompartments)]; %#ok<AGROW>
                value(end-numIntergenicSequences+1:end, ...
                    this.transcriptionUnits(i).idx, ...
                    this.transcriptionUnits(i).compartment.idx) = 1;
            end
            
            %store
            this.transcriptionUnitComposition_IntergenicRNAs = value;
        end
        
        function value = get.transcriptionUnitComposition_ProcessedRNAs(this)
            %retrieve
            if ~isempty(this.transcriptionUnitComposition_ProcessedRNAs)
                value = this.transcriptionUnitComposition_ProcessedRNAs;
                return;
            end
            
            %calculate
            value = zeros(0, this.numTranscriptionUnits, this.numCompartments);
            
            for i=1:this.numTranscriptionUnits
                if strcmp(this.transcriptionUnits(i).type,'mRNA')
                    numProcessedSequences = 1;
                else
                    numProcessedSequences = numel(this.transcriptionUnits(i).genes);
                end
                
                value = [value; zeros(numProcessedSequences, this.numTranscriptionUnits, this.numCompartments)]; %#ok<AGROW>
                value(end-numProcessedSequences+1:end, ...
                    this.transcriptionUnits(i).idx, ...
                    this.transcriptionUnits(i).compartment.idx) = 1;
            end
            
            %store
            this.transcriptionUnitComposition_ProcessedRNAs = value;
        end
        
        function value = get.processedRNAComposition_Genes(this)
            %retrieve
            if ~isempty(this.processedRNAComposition_Genes)
                value = this.processedRNAComposition_Genes;
                return;
            end
            
            %calculate
            value = zeros(this.numGenes, 0, this.numCompartments);
            
            for i=1:this.numTranscriptionUnits
                if strcmp(this.transcriptionUnits(i).type,'mRNA')
                    value = [value zeros(this.numGenes, 1, this.numCompartments)]; %#ok<AGROW>
                    value([this.transcriptionUnits(i).genes.idx], end, this.transcriptionUnits(i).compartment.idx) = 1;
                else
                    numProcessedSequences = numel(this.transcriptionUnits(i).genes);
                    
                    value = [value zeros(this.numGenes, numProcessedSequences, this.numCompartments)]; %#ok<AGROW>
                    value([this.transcriptionUnits(i).genes.idx], end-numProcessedSequences+1:end, this.transcriptionUnits(i).compartment.idx) = eye(numProcessedSequences);
                end
            end
            
            %store
            this.processedRNAComposition_Genes = value;
        end

        %% proteins
        % Converts DNA sequences to tRNA sequences using genetic code.
        % For each DNA sequence
        % 1. Identify first start codon, record correspoding tRNA
        % 2. Store tRNA of successive codons until stop codon or end of sequence is
        %    reached
        function value = get.proteinMonomerTRNASequences(this)
            %retrieve
            if ~isempty(this.proteinMonomerTRNASequences)
                value = this.proteinMonomerTRNASequences;
                return;
            end
            
            %calculate
            genes = [this.proteinMonomers.gene]; %#ok<*PROP>
            dnaSequences = {genes.dnaSequence}'; %#ok<*PROP>
            value = this.computeTRNASequences(dnaSequences);
            
            %store
            this.proteinMonomerTRNASequences = value;
        end

        function [tRNASequences, aaSequences] = computeTRNASequences(this, dnaSequences, findStartCodon)
            %options
            if ~exist('findStartCodon', 'var')
                findStartCodon = true;
            end

            %tRNA sequence, protein lengths
            tRNASequences = cell(size(dnaSequences));        %tRNA sequences of each protein-coding gene
            aaSequences   = cell(size(dnaSequences));        %aa sequences of each protein-coding gene

            geneticCodeCodons          = this.geneticCodeCodons;          %DNA sequence of codons
            geneticCodeTRNAs           = this.geneticCodeTRNAs;           %index of tRNA within tRNAGenes corresponding to each codon
            geneticCodeAminoAcids      = this.geneticCodeAminoAcids;      %index of AA within tRNAGenes corresponding to each start codon
            geneticCodeStartCodons     = this.geneticCodeStartCodons;     %DNA sequence of start codons
            geneticCodeStartTRNAs      = this.geneticCodeStartTRNAs;      %index of tRNA within tRNAGenes corresponding to each start codon
            geneticCodeStartAminoAcids = this.geneticCodeStartAminoAcids; %index of AA within tRNAGenes corresponding to each start codon

            %convert nucleic acids sequences to tRNA sequences
            for i = 1:length(dnaSequences)
                sequence = dnaSequences{i};
                tRNASequence = zeros(floor(length(sequence)/3),1);
                aaSequence = cell(floor(length(sequence)/3),1);

                startCodonPosition = -2; %position of first base of start codon within sequence
                if findStartCodon
                    %Start codon
                    %1. Step along sequence until first start codon
                    %2. Store tRNA corresponding to first start codon
                    for j=1:length(sequence)-2
                        tRNA=geneticCodeStartTRNAs(strcmp(geneticCodeStartCodons, sequence(j:j+2)));  %index of tRNA within tRNAGenes
                        if ~isempty(tRNA)
                            tRNASequence(1) = tRNA;
                            aaSequence{1} = geneticCodeStartAminoAcids{strcmp(geneticCodeStartCodons, sequence(j:j+2))};
                            startCodonPosition=j;
                            break;
                        end
                    end

                    if startCodonPosition == -2
                        throw(MException('KnowledgeBase:computeTRNASequences','No start codon found'));
                    end
                end

                %Non-start codons
                %1. Step along sequence until stop codon or the end of the sequence
                %2. Store tRNA corresponding to each codon
                endCodonPosition = startCodonPosition; %position of first base of last protein-coding codon within sequence
                for j = startCodonPosition+3:3:length(sequence)-2
                    tRNA = geneticCodeTRNAs(strcmp(geneticCodeCodons, sequence(j:j+2)));
                    if ~isempty(tRNA)
                        endCodonPosition=j;
                    else
                        break;
                    end
                    tRNASequence((j-startCodonPosition)/3+findStartCodon) = tRNA;
                    aaSequence{(j-startCodonPosition)/3+findStartCodon} = geneticCodeAminoAcids{strcmp(geneticCodeCodons, sequence(j:j+2))};
                end
                proteinLength = (endCodonPosition-startCodonPosition)/3+findStartCodon;
                tRNASequences{i} = tRNASequence(1:proteinLength);
                aaSequences{i} = aaSequence(1:proteinLength);
            end
        end

        function value = get.proteinMonomerTRNACounts(this)
            %retrieve
            if ~isempty(this.proteinMonomerTRNACounts)
                value = this.proteinMonomerTRNACounts;
                return;
            end
            
            %calculate
            monomerTRNASequences = this.proteinMonomerTRNASequences;
            value = zeros(length(this.proteinMonomers), length(this.tRNAGenes));

            for i = 1:length(this.proteinMonomers)
                value(i, :) = histc(monomerTRNASequences{i}, 1:length(this.tRNAGenes));
            end
            
            %store
            this.proteinMonomerTRNACounts = value;
        end

        function value = get.proteinMonomerCompartments(this)
            %retrieve
            if ~isempty(this.proteinMonomerCompartments)
                value = this.proteinMonomerCompartments;
                return;
            end
            
            %calculate
            monomerCompartments = [this.proteinMonomers.compartment];
            value = [monomerCompartments.idx]';
            
            %store
            this.proteinMonomerCompartments = value;
        end

        %base composition of protein monomers (amino acids)
        function value = get.proteinMonomerComposition(this)
            %retrieve
            if ~isempty(this.proteinMonomerComposition)
                value = this.proteinMonomerComposition;
                return;
            end
            
            %calculate
            value = reshape([this.proteinMonomers.baseCount], [], this.numProteinMonomers)';
            
            %store
            this.proteinMonomerComposition = value;
        end

        function value = get.proteinMonomerNTerminalMethionineCleavages(this)
            %retrieve
            if ~isempty(this.proteinMonomerNTerminalMethionineCleavages)
                value = this.proteinMonomerNTerminalMethionineCleavages;
                return;
            end
            
            %calculate
            value = [this.proteinMonomers.nTerminalMethionineCleavage]';
            
            %store
            this.proteinMonomerNTerminalMethionineCleavages = value;
        end

        function value = get.proteinMonomerSignalSequenceBaseCounts(this)
            %retrieve
            if ~isempty(this.proteinMonomerSignalSequenceBaseCounts)
                value = this.proteinMonomerSignalSequenceBaseCounts;
                return;
            end
            
            %calculate
            value = reshape([this.proteinMonomers.signalSequenceBaseCount], [], this.numProteinMonomers)';
            
            %store
            this.proteinMonomerSignalSequenceBaseCounts = value;
        end

        function value = get.proteinMonomerSignalSequenceMolecularWeights(this)
            %retrieve
            if ~isempty(this.proteinMonomerSignalSequenceMolecularWeights)
                value = this.proteinMonomerSignalSequenceMolecularWeights;
                return;
            end
            
            %calculate
            value = [this.proteinMonomers.signalSequenceMolecularWeight]';
            
            %store
            this.proteinMonomerSignalSequenceMolecularWeights = value;
        end

        function value = get.proteinProstheticGroupMatrix(this)
            %retrieve
            if ~isempty(this.proteinProstheticGroupMatrix)
                value = this.proteinProstheticGroupMatrix;
                return;
            end
            
            %calculate
            value=zeros(this.numProteinMonomers + this.numProteinComplexs, this.numMetabolites, this.numCompartments);

            for i=1:this.numProteinMonomers
                prostheticGroups = this.proteinMonomers(i).prostheticGroups;
                compartments     = this.proteinMonomers(i).prostheticGroupCompartments;
                coefficients     = [this.proteinMonomers(i).prostheticGroupCoefficients];
                for j=1:length(prostheticGroups)
                    value(i,prostheticGroups(j).idx,compartments(j).idx) = coefficients(j);
                end
            end

            for i=1:this.numProteinComplexs
                prostheticGroups = this.proteinComplexs(i).prostheticGroups;
                compartments     = this.proteinComplexs(i).prostheticGroupCompartments;
                coefficients     = [this.proteinComplexs(i).prostheticGroupCoefficients];
                for j=1:length(prostheticGroups)
                    value(i + this.numProteinMonomers,prostheticGroups(j).idx,compartments(j).idx) = coefficients(j);
                end
            end
            
            %store
            this.proteinProstheticGroupMatrix = value;
        end

        function value = get.proteinChaperoneMatrix(this)
            %retrieve
            if ~isempty(this.proteinChaperoneMatrix)
                value = this.proteinChaperoneMatrix;
                return;
            end
            
            %calculate
            value = zeros(this.numProteinMonomers + this.numProteinComplexs, ...
                this.numProteinMonomers + this.numProteinComplexs, ...
                this.numCompartments);

            for i=1:this.numProteinMonomers
                substrates   = this.proteinMonomers(i).chaperoneSubstrates;
                compartments = this.proteinMonomers(i).chaperoneCompartments;
                coefficients = [this.proteinMonomers(i).chaperoneCoefficients];
                for j=1:length(substrates)
                    if isa(substrates(j),'edu.stanford.covert.cell.kb.ProteinMonomer')
                        value(i,substrates(j).idx,compartments(j).idx) = coefficients(j);
                    else
                        value(i,substrates(j).idx+this.numProteinMonomers,compartments(j).idx) = coefficients(j);
                    end
                end
            end

            for i=1:this.numProteinComplexs
                substrates   = this.proteinComplexs(i).chaperoneSubstrates;
                compartments = this.proteinComplexs(i).chaperoneCompartments;
                coefficients = [this.proteinComplexs(i).chaperoneCoefficients];
                for j=1:length(substrates)
                    if isa(substrates(j),'edu.stanford.covert.cell.kb.ProteinMonomer')
                        value(i+this.numProteinMonomers,substrates(j).idx,compartments(j).idx) = coefficients(j);
                    else
                        value(i+this.numProteinMonomers,substrates(j).idx+this.numProteinMonomers,compartments(j).idx) = coefficients(j);
                    end
                end
            end
            
            %store
            this.proteinChaperoneMatrix = value;
        end

        %protein complex RNA composition (immediate children on
        %composition tree)
        function value = get.proteinComplexRNAComposition(this)
            %retrieve
            if ~isempty(this.proteinComplexRNAComposition)
                value = this.proteinComplexRNAComposition;
                return;
            end
            
            %calculate
            value=zeros(this.numGenes,this.numProteinComplexs,this.numCompartments);
            for i=1:this.numProteinComplexs
                if isempty(this.proteinComplexs(i).rnas); continue; end;
                genes        = [this.proteinComplexs(i).rnas];
                compartments = [this.proteinComplexs(i).rnaCompartments.idx];
                coefficients = this.proteinComplexs(i).rnaCoefficients;

                for j=1:length(genes)
                    value(genes(j).idx, i, compartments(j)) = coefficients(j);
                end
            end
            
            %store
            this.proteinComplexRNAComposition = value;
        end

        %protein complex protein monomer composition (immediate children on
        %composition tree)
        function value = get.proteinComplexMonomerComposition(this)
            %retrieve
            if ~isempty(this.proteinComplexMonomerComposition)
                value = this.proteinComplexMonomerComposition;
                return;
            end
            
            %calculate
            value = zeros(this.numGenes,this.numProteinComplexs,this.numCompartments);
            for i = 1:this.numProteinComplexs
                if isempty(this.proteinComplexs(i).proteinMonomers); continue; end;
                genes        = [this.proteinComplexs(i).proteinMonomers.gene];
                compartments = [this.proteinComplexs(i).proteinMonomerCompartments.idx];
                coefficients = this.proteinComplexs(i).proteinMonomerCoefficients;

                for j = 1:length(genes)
                    value(genes(j).idx, i, compartments(j)) = coefficients(j);
                end
            end
            
            %store
            this.proteinComplexMonomerComposition = value;
        end

        %protein complex protein complex composition (immediate children on
        %composition tree)
        function value = get.proteinComplexComplexComposition(this)
            %retrieve
            if ~isempty(this.proteinComplexComplexComposition)
                value = this.proteinComplexComplexComposition;
                return;
            end
            
            %calculate
            value = zeros(this.numProteinComplexs, this.numProteinComplexs, this.numCompartments);
            for i = 1:this.numProteinComplexs
                if isempty(this.proteinComplexs(i).proteinComplexs); continue; end;
                complexs     = [this.proteinComplexs(i).proteinComplexs];
                compartments = [this.proteinComplexs(i).proteinComplexCompartments.idx];
                coefficients = this.proteinComplexs(i).proteinComplexCoefficients;

                for j=1:length(complexs)
                    value(complexs(j).idx, i, compartments(j)) = coefficients(j);
                end
            end
            
            %store
            this.proteinComplexComplexComposition = value;
        end

        %protein complex RNA composition (leaves of composition tree)
        function value = get.proteinComplexAllRNAComposition(this)
            %retrieve
            if ~isempty(this.proteinComplexAllRNAComposition)
                value = this.proteinComplexAllRNAComposition;
                return;
            end
            
            %calculate
            value = reshape([this.proteinComplexs.rnaComposition], this.numGenes, this.numProteinComplexs, this.numCompartments);
            
            %store
            this.proteinComplexAllRNAComposition = value;
        end

        %protein complex protein monomer composition (leaves of composition
        %tree)
        function value = get.proteinComplexAllMonomerComposition(this)
            %retrieve
            if ~isempty(this.proteinComplexAllMonomerComposition)
                value = this.proteinComplexAllMonomerComposition;
                return;
            end
            
            %calculate
            value = reshape([this.proteinComplexs.proteinMonomerComposition], ...
                this.numProteinMonomers, this.numProteinComplexs, this.numCompartments);
                        
            %store
            this.proteinComplexAllMonomerComposition = value;
        end

        function value = get.proteinComplexCompartments(this)
            %retrieve
            if ~isempty(this.proteinComplexCompartments)
                value = this.proteinComplexCompartments;
                return;
            end
            
            %calculate
            complexCompartments = [this.proteinComplexs.compartment];
            value = [complexCompartments.idx]';
            
            %store
            this.proteinComplexCompartments = value;
        end

        %% reactions
        function value = get.reactionStimuliStoichiometryMatrix(this)
            %retrieve
            if ~isempty(this.reactionStimuliStoichiometryMatrix)
                value = this.reactionStimuliStoichiometryMatrix;
                return;
            end
            
            %calculate
            value = zeros(this.numStimulis, this.numReactions, this.numCompartments);

            for i = 1:this.numReactions
                if isempty(this.reactions(i).stimulis); continue; end;

                stimulis  = [this.reactions(i).stimulis.idx];
                compartments = [this.reactions(i).stimuliCompartments.idx];
                coefficients = this.reactions(i).stimuliCoefficients;

                for j = 1:length(stimulis);
                    value(stimulis(j), i, compartments(j)) = coefficients(j);
                end
            end
            
            %store
            this.reactionStimuliStoichiometryMatrix = value;
        end

        function value = get.reactionMetaboliteStoichiometryMatrix(this)
            %retrieve
            if ~isempty(this.reactionMetaboliteStoichiometryMatrix)
                value = this.reactionMetaboliteStoichiometryMatrix;
                return;
            end
            
            %calculate
            value = zeros(this.numMetabolites, this.numReactions, this.numCompartments);

            for i = 1:this.numReactions
                if isempty(this.reactions(i).metabolites); continue; end;

                metabolites  = [this.reactions(i).metabolites.idx];
                compartments = [this.reactions(i).metaboliteCompartments.idx];
                coefficients = this.reactions(i).metaboliteCoefficients;

                for j = 1:length(metabolites);
                    value(metabolites(j), i, compartments(j)) = coefficients(j);
                end
            end
            
            %store
            this.reactionMetaboliteStoichiometryMatrix = value;
        end

        function value = get.reactionRNAStoichiometryMatrix(this)
            %retrieve
            if ~isempty(this.reactionRNAStoichiometryMatrix)
                value = this.reactionRNAStoichiometryMatrix;
                return;
            end
            
            %calculate
            value = zeros(this.numTranscriptionUnits, this.numReactions, this.numCompartments);

            for i=1:this.numReactions
                if isempty(this.reactions(i).rnas); continue; end;

                rnas  = [this.reactions(i).rnas.idx];
                compartments = [this.reactions(i).rnaCompartments.idx];
                coefficients = this.reactions(i).rnaCoefficients;

                for j=1:length(rnas);
                    value(rnas(j), i, compartments(j))=coefficients(j);
                end
            end
            
            %store
            this.reactionRNAStoichiometryMatrix = value;
        end

        function value = get.reactionProteinMonomerStoichiometryMatrix(this)
            %retrieve
            if ~isempty(this.reactionProteinMonomerStoichiometryMatrix)
                value = this.reactionProteinMonomerStoichiometryMatrix;
                return;
            end
            
            %calculate
            value = zeros(this.numProteinMonomers, this.numReactions, this.numCompartments);

            for i = 1:this.numReactions
                if isempty(this.reactions(i).proteinMonomers); continue; end;

                proteinMonomers  = [this.reactions(i).proteinMonomers.idx];
                compartments = [this.reactions(i).proteinMonomerCompartments.idx];
                coefficients = this.reactions(i).proteinMonomerCoefficients;

                for j = 1:length(proteinMonomers);
                    value(proteinMonomers(j), i, compartments(j)) = coefficients(j);
                end
            end
            
            %store
            this.reactionProteinMonomerStoichiometryMatrix = value;
        end

        function value = get.reactionProteinComplexStoichiometryMatrix(this)
            %retrieve
            if ~isempty(this.reactionProteinComplexStoichiometryMatrix)
                value = this.reactionProteinComplexStoichiometryMatrix;
                return;
            end
            
            %calculate
            value = zeros(this.numProteinComplexs, this.numReactions, this.numCompartments);

            for i = 1:this.numReactions
                if isempty(this.reactions(i).proteinComplexs); continue; end;

                proteinComplexs  = [this.reactions(i).proteinComplexs.idx];
                compartments = [this.reactions(i).proteinComplexCompartments.idx];
                coefficients = this.reactions(i).proteinComplexCoefficients;

                for j=1:length(proteinComplexs);
                    value(proteinComplexs(j), i, compartments(j)) = coefficients(j);
                end
            end
            
            %store
            this.reactionProteinComplexStoichiometryMatrix = value;
        end

        function value = get.reactionModificationMatrix(this)
            %retrieve
            if ~isempty(this.reactionModificationMatrix)
                value = this.reactionModificationMatrix;
                return;
            end
            
            %calculate
            value = zeros(this.numReactions, this.numGenes, this.numCompartments);
            for i = 1:this.numReactions
                modifications = this.reactions(i).stableModifications;
                compartments = this.reactions(i).stableModificationCompartments;
                positions = this.reactions(i).stableModificationPositions;
                for j=1:length(modifications)
                    if isa(modifications(j),'edu.stanford.covert.cell.kb.ProteinMonomer')
                        value(i, modifications(j).gene.idx, compartments(j).idx) = positions(j);
                    else
                        value(i, modifications(j).idx, compartments(j).idx) = positions(j);
                    end
                end
            end
            
            %store
            this.reactionModificationMatrix = value;
        end

        %converts units to one of:
        % - reactions/(gram dry biomass)/s
        % - reactions/s
        function value = get.reactionBounds(this)
            %retrieve
            if ~isempty(this.reactionBounds)
                value = this.reactionBounds;
                return;
            end
            
            %calculate
            import edu.stanford.covert.util.ConstantUtil;
            value = [[this.reactions.lowerBound]' [this.reactions.upperBound]'];

            boundUnits = {this.reactions.boundUnits};
            value(strcmp(boundUnits,'mmol/gDCW/h'),:)   = value(strcmp(boundUnits,'mmol/gDCW/h'),:) / ConstantUtil.secondsPerHour * (ConstantUtil.nAvogadro/1000);
            value(strcmp(boundUnits,'dimensionless'),:) = 0;
            
            %store
            this.reactionBounds = value;
        end

        %reactions/enzyme/s
        function value = get.enzymeBounds(this)
            import edu.stanford.covert.util.ConstantUtil;
            
            %retrieve
            if ~isempty(this.enzymeBounds)
                value = this.enzymeBounds;
                return;
            end
            
            %calculate

            value = repmat([-Inf Inf],length(this.reactions),1);

            for i = 1:length(this.reactions)
                reaction = this.reactions(i);
                enzymes = reaction.enzymes;
                if isempty(enzymes); continue; end;
                switch this.reactions(i).direction
                    case 'F'
                        value(i,2) = this.enzymeBound(reaction.vmaxForward, reaction.vmaxUnitsForward, enzymes.molecularWeight);
                    case 'B'
                        value(i,1) = -this.enzymeBound(reaction.vmaxBackward, reaction.vmaxUnitsBackward, enzymes.molecularWeight);
                    case 'R'
                        value(i,1) = -this.enzymeBound(reaction.vmaxBackward, reaction.vmaxUnitsBackward, enzymes.molecularWeight);
                        value(i,2) = this.enzymeBound(reaction.vmaxForward, reaction.vmaxUnitsForward, enzymes.molecularWeight);
                end
            end

            value = value / ConstantUtil.secondsPerMinute;
            
            %store
            this.enzymeBounds = value;
        end

        %units: 1/min
        function bound = enzymeBound(~, vmax, vmaxUnits, enzymeMolecularWeight)
            if isnan(vmax) || vmax == 0
                bound=Inf;
                return;
            end

            switch vmaxUnits
                case '1/min'
                    bound = vmax;
                case 'U/mg'
                    bound = vmax * enzymeMolecularWeight * 1e-3;
                otherwise
                    fprintf('Error: incompatible units %s\n', vmaxUnits);
                    bound = Inf;
            end
        end

        function value = get.reactionCatalysisMatrix(this)
            %retrieve
            if ~isempty(this.reactionCatalysisMatrix)
                value = this.reactionCatalysisMatrix;
                return;
            end
            
            %calculate
            value = zeros(length(this.reactions), this.numProteinMonomers + this.numProteinComplexs, this.numCompartments);
            for i=1:length(this.reactions)
                enzyme = this.reactions(i).enzymes;
                if isempty(enzyme); continue; end;
                if isa(enzyme,'edu.stanford.covert.cell.kb.ProteinMonomer');
                    value(i,enzyme.idx, this.reactions(i).enzymeCompartments.idx) = 1;
                else
                    value(i,enzyme.idx + this.numProteinMonomers, this.reactions(i).enzymeCompartments.idx) = 1;
                end
            end
            
            %store
            this.reactionCatalysisMatrix = value;
        end

        function value = get.reactionCoenzymeMatrix(this)
            %retrieve
            if ~isempty(this.reactionCoenzymeMatrix)
                value = this.reactionCoenzymeMatrix;
                return;
            end
            
            %calculate
            value = zeros(this.numReactions,this.numMetabolites,this.numCompartments);
            for i = 1:this.numReactions
                coenzymes = this.reactions(i).coenzymes;
                coenzymeCompartments = this.reactions(i).coenzymeCompartments;
                coenzymeCoefficients = this.reactions(i).coenzymeCoefficients;
                for j=1:length(coenzymes)
                    value(i,coenzymes(j).idx,coenzymeCompartments(j).idx) = coenzymeCoefficients(j);
                end
            end
            
            %store
            this.reactionCoenzymeMatrix = value;
        end

        %% statistics
        function value = get.maxGeneLength(this)
            %retrieve
            if ~isempty(this.maxGeneLength)
                value = this.maxGeneLength;
                return;
            end
            
            %calculate
            value = max([this.genes.sequenceLength]);
            
            %store
            this.maxGeneLength = value;
        end

        function value = get.maxTranscriptionUnitLength(this)
            %retrieve
            if ~isempty(this.maxTranscriptionUnitLength)
                value = this.maxTranscriptionUnitLength;
                return;
            end
            
            %calculate
            value = max([this.transcriptionUnits.sequenceLength]);
            
            %store
            this.maxTranscriptionUnitLength = value;
        end

        function value = get.maxGeneExpression(this)
            %retrieve
            if ~isempty(this.maxGeneExpression)
                value = this.maxGeneExpression;
                return;
            end
            
            %calculate
            value = reshape([this.genes.expression],[],this.numGenes)';
            
            %store
            this.maxGeneExpression = value;
        end

        function value = get.maxTranscriptionUnitExpression(this)
            %retrieve
            if ~isempty(this.maxTranscriptionUnitExpression)
                value = this.maxTranscriptionUnitExpression;
                return;
            end
            
            %calculate
            value = reshape([this.transcriptionUnits.expression], [], this.numTranscriptionUnits)';
            
            %store
            this.maxTranscriptionUnitExpression = value;
        end
    end

    methods (Access = protected)
        function loadFromKnowledgeBase(this, database)                      
            %summary
            database.setNullValue(0);
            database.prepareStatement('CALL get_knowledgebase("{Si}")', this.wid);
            data = database.query();
            this.wholeCellModelID = data.WholeCellModelID{1};
            this.version = data.Version{1};
            this.investigator = data.Investigator{1};
            this.name = data.Name{1};
            this.taxonomy = data.Taxonomy{1};
            this.translationTable = data.TranslationTable;

            %statistics
            this.maxAvgMetConc = data.MaxAvgMetConc;
            this.maxAvgRxnFlux = data.MaxAvgRxnFlux;

            %% create objects

            %processes
            database.setNullValue(0);
            database.prepareStatement('CALL get_processs("{Si}",null)', this.wid);
            data = database.query();
            this.processes = edu.stanford.covert.cell.kb.Process(...
                this, data.WID, data.WholeCellModelID, data.Name, ...
                data.InitializationOrder, data.EvaluationOrder, data.Class);
            
            %states
            database.setNullValue(0);
            database.prepareStatement('CALL get_states("{Si}",null)', this.wid);
            data = database.query();
            this.states = edu.stanford.covert.cell.kb.State(...
                this, data.WID, data.WholeCellModelID, data.Name, ...
                data.Class);

            %parameters
            database.setNullValue(0);
            database.prepareStatement('CALL get_parameters("{Si}",null)', this.wid);
            data = database.query();
            this.parameters = edu.stanford.covert.cell.kb.Parameter(...
                this, data.WID, data.WholeCellModelID, data.Name, ...
                data.Index, data.DefaultValue, data.Units, strcmp(data.ExperimentallyConstrained,'Y'));

            %compartments
            database.setNullValue(0);
            database.prepareStatement('CALL get_compartments("{Si}",null)', this.wid);
            data = database.query();
            this.compartments = edu.stanford.covert.cell.kb.Compartment(...
                this, data.WID, data.WholeCellModelID, data.Name);

            %metabolites
            database.setNullValue(NaN);
            database.prepareStatement('CALL get_metabolites("{Si}",null)', this.wid);
            data = database.query();
            this.metabolites = edu.stanford.covert.cell.kb.Metabolite(...
                this, data.WID, data.WholeCellModelID, data.Name, data.TraditionalName, data.IUPACName, ...
                data.Category, data.Subcategory, ...
                data.EmpiricalFormula, data.Smiles, data.Charge, strcmp(data.Hydrophobic, 'Y'), ...
                data.pKa, data.pI, data.logP, data.logD, data.Volume, data.MolecularWeightCalc, ...
                data.ExchangeLowerBound, data.ExchangeUpperBound);

            %genome
            database.setNullValue(0);
            database.prepareStatement('CALL get_genome("{Si}",null)', this.wid);
            data = database.query();
            this.genome = edu.stanford.covert.cell.kb.Genome(...
                this, data.WID,[],[],...
                data.Topology,data.Sequence);

            %genes
            database.setNullValue(0);
            database.prepareStatement('CALL get_genes("{Si}",null)', this.wid);
            data = database.query();
            this.genes = edu.stanford.covert.cell.kb.Gene(...
                this, data.WID, data.WholeCellModelID, data.Name, ...
                data.Symbol, data.Synonyms, data.Type, strcmp(data.StartCodon,'Y'), data.Codons, data.Coordinate, ...
                data.Length, strcmp(data.Direction,'forward'), ...
                data.Essential, data.HalfLifeExp, ...
                [data.Expression, data.ExpressionColdShock, data.ExpressionHeatShock]);

            %transcriptionUnits
            database.setNullValue(0);
            database.prepareStatement('CALL get_transcriptionunits("{Si}",null)', this.wid);
            data = database.query();
            this.transcriptionUnits = edu.stanford.covert.cell.kb.TranscriptionUnit(...
                this, data.WID, data.WholeCellModelID, data.Name, ...
                data.Promoter35Coordinate, data.Promoter35Length,...
                data.Promoter10Coordinate, data.Promoter10Length,...
                data.TSSCoordinate);

            %genomeFeatures
            database.setNullValue(0);
            database.prepareStatement('CALL get_genomefeatures("{Si}",null)', this.wid);
            data = database.query();
            this.genomeFeatures = edu.stanford.covert.cell.kb.GenomeFeature(...
                this, data.WID, data.WholeCellModelID, data.Name, ...
                data.Type, data.Subtype, data.Coordinate, data.Length, strcmp(data.Direction,'forward'));

            %proteinMonomers
            database.setNullValue(0);
            database.prepareStatement('CALL get_proteinmonomers("{Si}",null)', this.wid);
            data = database.query();
            this.proteinMonomers = edu.stanford.covert.cell.kb.ProteinMonomer(...
                this, data.WID, data.WholeCellModelID, data.Name, ...
                data.Topology, data.ActiveSite, data.MetalBindingSite, data.DNAFootprint, data.DNAFootprintBindingStrandedness, data.DNAFootprintRegionStrandedness, data.MolecularInteraction, ...
                data.ChemicalRegulation, data.Subsystem, data.GeneralClassification, data.ProteaseClassification, ...
                data.TransporterClassification, ...
                strcmp(data.NTerminalMethionineCleavage,'Y'), data.SignalSequenceType, data.SignalSequenceLocation, data.SignalSequenceLength, ...
                data.ActivationRule);

            %proteinComplexs
            database.setNullValue(0);
            database.prepareStatement('CALL get_proteincomplexs("{Si}", null)', this.wid);
            data = database.query();
            this.proteinComplexs = edu.stanford.covert.cell.kb.ProteinComplex(...
                this, data.WID, data.WholeCellModelID, data.Name, ...
                data.DNAFootprint, data.DNAFootprintBindingStrandedness, data.DNAFootprintRegionStrandedness, data.MolecularInteraction, data.ChemicalRegulation, ...
                data.Subsystem, data.GeneralClassification, data.ProteaseClassification, ...
                data.TransporterClassification, ...
                data.ActivationRule, data.DisulfideBonds);

            %reactions
            database.setNullValue(0);
            database.prepareStatement('CALL get_reactions("{Si}", null)', this.wid);
            data = database.query();
            this.reactions = edu.stanford.covert.cell.kb.Reaction(...
                this, data.WID, data.WholeCellModelID, data.Name, data.Type, ...
                data.ECNumber, data.Spontaneous, data.Direction, data.DeltaG, data.Keq, ...
                data.RateLawForward, data.KmForward, data.VmaxExpForward, data.VmaxExpUnitForward, ...
                data.RateLawBackward, data.KmBackward, data.VmaxExpBackward, data.VmaxExpUnitBackward, ...
                data.OptimalpH, data.OptimalTemp, ...
                data.Activators, data.Inhibitors, ...
                data.LowerBound, data.UpperBound, data.BoundUnits);

            %pathways
            database.setNullValue(0);
            database.prepareStatement('CALL get_pathways("{Si}",null)', this.wid);
            data = database.query();
            this.pathways = edu.stanford.covert.cell.kb.Pathway(...
                this, data.WID, data.WholeCellModelID, data.Name);

            %stimuli
            database.setNullValue(0);
            database.prepareStatement('CALL get_stimulis("{Si}",null)', this.wid);
            data = database.query();
            this.stimulis = edu.stanford.covert.cell.kb.Stimuli(...
                this, data.WID, data.WholeCellModelID, data.Name);

            %notes
            database.setNullValue(0);
            database.prepareStatement('CALL get_notes("{Si}",null)', this.wid);
            data = database.query();
            this.notes = edu.stanford.covert.cell.kb.Note(...
                this, data.WID, data.WholeCellModelID, []);

            %references
            database.setNullValue(0);
            database.prepareStatement('CALL get_references("{Si}",null)', this.wid);
            data = database.query();
            this.references = edu.stanford.covert.cell.kb.Reference(...
                this, data.WID, data.WholeCellModelID, [], ...
                data.Type, data.PMID, data.ISBN, ...
                data.Authors, data.Editors, data.Year, data.Title, data.Publication, data.Volume, data.Issue, data.Pages, ...
                data.Publisher, data.URL, data.Citations);

            %% link objects together
            %parameters-processes, states, reactions, molecules
            this.linkObjects(database, 'process', 'parameter', struct(...
                'sqlProcedure', 'get_processs_parameters'));
            this.linkObjects(database, 'state', 'parameter', struct(...
                'field2', 'state'));
            this.linkObjects(database, 'reaction', 'parameter', struct(...
                'field2', 'reactions'));
            this.linkObjects(database, 'proteinMonomer', 'parameter', struct(...
                'field2', 'proteinMonomers'));
            this.linkObjects(database, 'proteinComplex', 'parameter', struct(...
                'field2', 'proteinComplexs'));

            %genome
            this.genome.genes = this.genes;
            this.genome.transcriptionUnits = this.transcriptionUnits;
            this.genome.features = this.genomeFeatures;
            for i = 1:size(this.genes,1)
                this.genes(i).genome = this.genome;
            end
            for i = 1:size(this.transcriptionUnits,1)
                this.transcriptionUnits(i).genome = this.genome;
            end
            for i = 1:size(this.genomeFeatures,1)
                this.genomeFeatures(i).genome = this.genome;
            end

            %compartments
            this.linkObjects(database, 'gene', 'compartment', struct(...
                'field1', 'compartment'));
            this.linkObjects(database, 'transcriptionUnit', 'compartment', struct(...
                'field1', 'compartment'));
            this.linkObjects(database, 'proteinMonomer', 'compartment', struct(...
                'field1', 'compartment'));
            this.linkObjects(database, 'proteinComplex', 'compartment', struct(...
                'field1', 'compartment'));

            %prosthetic groups
            this.linkObjects(database,'proteinMonomer','metabolite',struct(...
                'sqlProcedure','get_proteinmonomers_prostheticgroups',...
                'field1', 'prostheticGroups', ...
                'field2', '',...
                'compartmentName', 'MetaboliteCompartmentIdx',...
                'compartmentField1', 'prostheticGroupCompartments',...
                'compartmentField2', '',...
                'weightName', 'Coefficient',...
                'weightField1', 'prostheticGroupCoefficients',...
                'weightField2', '',...
                'nanValue',NaN));
            this.linkObjects(database,'proteinComplex','metabolite',struct(...
                'sqlProcedure','get_proteincomplexs_prostheticgroups',...
                'field1', 'prostheticGroups', ...
                'field2', '',...
                'compartmentName', 'MetaboliteCompartmentIdx',...
                'compartmentField1', 'prostheticGroupCompartments',...
                'compartmentField2', '',...
                'weightName', 'Coefficient',...
                'weightField1', 'prostheticGroupCoefficients',...
                'weightField2', '',...
                'nanValue',NaN));

            %chaperones
            this.linkObjects(database,'proteinMonomer','proteinMonomer',struct(...
                'sqlProcedure','get_proteinmonomers_proteinmonomerchaperones',...
                'name2','ProteinMonomerChaperoneIdx',...
                'field1','',...
                'field2','chaperoneSubstrates',...
                'compartmentName','ProteinMonomerChaperoneCompartmentIdx',...
                'compartmentField1','',...
                'compartmentField2','chaperoneCompartments',...
                'weightName','Coefficient',...
                'weightField1','',...
                'weightField2','chaperoneCoefficients'));
            %this.linkObjects(database,'proteinComplex','proteinMonomer',struct(...
            %    'sqlProcedure','get_proteincomplexs_proteinmonomerchaperones',...
            %    'name2','ProteinMonomerChaperoneIdx',...
            %    'field1',''...
            %    'field2','chaperoneSubstrates',...
            %    'compartmentName','ProteinMonomerChaperoneCompartmentIdx',...
            %    'comparatmentField1','',...
            %    'comparatmentField2','chaperoneCompartments',...
            %    'weightName','Coefficient',...
            %    'weightField1','',...
            %    'weightField2','chaperoneCoefficients'));
            this.linkObjects(database,'proteinMonomer','proteinComplex',struct(...
                'sqlProcedure','get_proteinmonomers_proteincomplexchaperones',...
                'name2','ProteinComplexChaperoneIdx',...
                'field1','',...
                'field2','chaperoneSubstrates',...
                'compartmentName','ProteinComplexChaperoneCompartmentIdx',...
                'compartmentField1','',...
                'compartmentField2','chaperoneCompartments',...
                'weightName','Coefficient',...
                'weightField1','',...
                'weightField2','chaperoneCoefficients'));
            %this.linkObjects(database,'proteinComplex','proteinComplex',struct(...
            %    'sqlProcedure','get_proteincomplexs_proteincomplexchaperones',...
            %    'name2','ProteinComplexChaperoneIdx',...
            %    'field1',''...
            %    'field2','chaperoneSubstrates',...
            %    'compartmentName','ProteinComplexChaperoneCompartmentIdx',...
            %    'comparatmentField1','',...
            %    'comparatmentField2','chaperoneCompartments',...
            %    'weightName','Coefficient',...
            %    'weightField1','',...
            %    'weightField2','chaperoneCoefficients'));

            %coenzymes
            this.linkObjects(database,'reaction','metabolite',struct(...
                'sqlProcedure','get_reactions_coenzymes',...
                'field1','coenzymes',...
                'field2','coenzymeReactions',...
                'compartmentName','MetaboliteCompartmentIdx',...
                'compartmentField1','coenzymeCompartments',...
                'compartmentField2','',...
                'weightName','Coefficient',...
                'weightField1','coenzymeCoefficients',...
                'weightField2',''));

            %stable modifications
            this.linkObjects(database,'reaction','gene',struct(...
                'sqlProcedure','get_reactions_stablernamodifications',...
                'field1','stableModifications',...
                'field2','stableModificationReactions',...
                'compartmentName','GeneCompartmentIdx',...
                'compartmentField1','stableModificationCompartments',...
                'compartmentField2','',...
                'weightName','Position',...
                'weightField1','stableModificationPositions',...
                'weightField2','',...
                'nanValue',NaN));
            this.linkObjects(database,'reaction','proteinMonomer',struct(...
                'sqlProcedure','get_reactions_stableproteinmonomermodifications',...
                'field1','stableModifications',...
                'field2','stableModificationReactions',...
                'compartmentName','ProteinMonomerCompartmentIdx',...
                'compartmentField1','stableModificationCompartments',...
                'compartmentField2','',...
                'weightName','Position',...
                'weightField1','stableModificationPositions',...
                'weightField2','',...
                'nanValue',NaN));

            %tRNA aminoacylation
            this.linkObjects(database,'gene','metabolite',struct(...
                'sqlProcedure','get_genes_aminoacids',...
                'field1','aminoAcid',...
                'field2',''));

            %transcription unit composition
            this.linkObjects(database,'gene','transcriptionUnit',struct(...
                'compartmentName','GeneCompartmentIdx',...
                'compartmentField2','geneCompartments'));

            %gene->protein monomer
            this.linkObjects(database,'gene','proteinMonomer',struct(...
                'field2','gene',...
                'compartmentName','GeneCompartmentIdx',...
                'compartmentField2','geneCompartments'));

            %protein complex composition
            this.linkObjects(database,'metabolite','proteinComplex',struct(...
                'sqlProcedure','get_metabolites_proteincomplexs',...
                'name1','MetaboliteIdx',...
                'name2','ProteinComplexIdx',...
                'field1','',...
                'field2','metabolites',...
                'compartmentName','MetaboliteCompartmentIdx',...
                'compartmentField1','',...
                'compartmentField2','metaboliteCompartments',...
                'weightName','Coefficient',...
                'weightField1','',...
                'weightField2','metaboliteCoefficients'));
            this.linkObjects(database,'gene','proteinComplex',struct(...
                'sqlProcedure','get_rnasubunits_proteincomplexs',...
                'name1','RNASubunitIdx',...
                'name2','ProteinComplexIdx',...
                'field1','',...
                'field2','rnas',...
                'weightName','Coefficient',...
                'weightField1','',...
                'weightField2','rnaCoefficients',...
                'compartmentName','RNASubunitCompartmentIdx',...
                'compartmentField1','',...
                'compartmentField2','rnaCompartments'));
            this.linkObjects(database, 'proteinMonomer', 'proteinComplex',struct(...
                'sqlProcedure', 'get_proteinmonomersubunits_proteincomplexs',...
                'weightName','Coefficient',...
                'weightField1', '',...
                'weightField2', 'proteinMonomerCoefficients',...
                'compartmentName', 'ProteinMonomerCompartmentIdx',...
                'compartmentField1', '',...
                'compartmentField2', 'proteinMonomerCompartments'));
            this.linkObjects(database,'proteinComplex','proteinComplex',struct(...
                'sqlProcedure','get_proteincomplexsubunits_proteincomplexs',...
                'name1','ProteinComplexSubunitIdx',...
                'name2','ProteinComplexIdx',...
                'field1','',...
                'field2','proteinComplexs',...
                'weightName','Coefficient',...
                'weightField1','',...
                'weightField2','proteinComplexCoefficients',...
                'compartmentName','ProteinComplexSubunitCompartmentIdx',...
                'compartmentField1','',...
                'compartmentField2','proteinComplexCompartments'));

            %reaction catalysis
            this.linkObjects(database,'proteinMonomer','reaction',struct(...
                'field2','enzymes',...
                'compartmentName','ProteinMonomerCompartmentIdx',...
                'compartmentField2','enzymeCompartments'));
            this.linkObjects(database,'proteinComplex','reaction',struct(...
                'field2','enzymes',...
                'compartmentName','ProteinComplexCompartmentIdx',...
                'compartmentField2','enzymeCompartments'));

            %reaction stoichiometry
            this.linkObjects(database,'stimuli','reaction',struct(...
                'sqlProcedure','get_stimulireactantproducts_reactions',...
                'name1','StimulusIdx',...
                'name2','ReactionIdx',...
                'field1','',...
                'field2','stimulis',...
                'weightName','Coefficient',...
                'weightField1','',...
                'weightField2','stimuliCoefficients',...
                'compartmentName','StimulusCompartmentIdx',...
                'compartmentField1','',...
                'compartmentField2','stimuliCompartments'));
            this.linkObjects(database,'metabolite','reaction',struct(...
                'sqlProcedure','get_metabolitereactantproducts_reactions',...
                'weightName','Coefficient',...
                'weightField1','',...
                'weightField2','metaboliteCoefficients',...
                'compartmentName','MetaboliteCompartmentIdx',...
                'compartmentField1','',...
                'compartmentField2','metaboliteCompartments'));
            this.linkObjects(database,'transcriptionUnit','reaction',struct(...
                'sqlProcedure','get_transcriptionunitreactantproducts_reactions',...
                'name1','TranscriptionUnitIdx',...
                'name2','ReactionIdx',...
                'field1','',...
                'field2','rnas',...
                'weightName','Coefficient',...
                'weightField1','',...
                'weightField2','rnaCoefficients',...
                'compartmentName','CompartmentIdx',...
                'compartmentField1','',...
                'compartmentField2','rnaCompartments'));
            this.linkObjects(database,'proteinMonomer','reaction',struct(...
                'sqlProcedure','get_proteinmonomerreactantproducts_reactions',...
                'name1','MonomerIdx',...
                'name2','ReactionIdx',...
                'field1','',...
                'field2','proteinMonomers',...
                'weightName','Coefficient',...
                'weightField1','',...
                'weightField2','proteinMonomerCoefficients',...
                'compartmentName','CompartmentIdx',...
                'compartmentField1','',...
                'compartmentField2','proteinMonomerCompartments'));
            this.linkObjects(database,'proteinComplex','reaction',struct(...
                'sqlProcedure','get_proteincomplexreactantproducts_reactions',...
                'name1','ComplexIdx',...
                'name2','ReactionIdx',...
                'field1','',...
                'field2','proteinComplexs',...
                'weightName','Coefficient',...
                'weightField1','',...
                'weightField2','proteinComplexCoefficients',...
                'compartmentName','CompartmentIdx',...
                'compartmentField1','',...
                'compartmentField2','proteinComplexCompartments'));

            %transcription unit - transcription factor
            this.linkObjects(database,'transcriptionUnit','proteinMonomer',struct(...
                'sqlProcedure','get_transcriptionunits_transcriptionfactorproteinmonomers',...
                'field1','transcriptionFactorProteinMonomers',...
                'field2','regulatedTranscriptionUnits',...
                'weightName','Affinity',...
                'weightField1','transcriptionFactorProteinMonomerAffinitys',...
                'weightField2','',...
                'weight2Name','Activity',...
                'weight2Field1','transcriptionFactorProteinMonomerActivitys',...
                'weight2Field2','',...
                'weight3Name','Condition',...
                'weight3Field1','transcriptionFactorProteinMonomerConditions',...
                'weight3Field2','',...
                'weight4Name','BindingSiteStartCoordinate',...
                'weight4Field1','transcriptionFactorProteinMonomerBindingSiteStartCoordinates',...
                'weight4Field2','',...
                'weight5Name','BindingSiteLength',...
                'weight5Field1','transcriptionFactorProteinMonomerBindingSiteLengths',...
                'weight5Field2','',...
                'weight6Name','BindingSiteDirection',...
                'weight6Field1','transcriptionFactorProteinMonomerBindingSiteDirections',...
                'weight6Field2','',...
                'compartmentName','ProteinMonomerCompartmentIdx',...
                'compartmentField1','transcriptionFactorProteinMonomerCompartments',...
                'compartmentField2','',...
                'nanValue',NaN));
            this.linkObjects(database,'transcriptionUnit','proteinComplex',struct(...
                'sqlProcedure','get_transcriptionunits_transcriptionfactorproteincomplexs',...
                'field1','transcriptionFactorProteinComplexs',...
                'field2','regulatedTranscriptionUnits',...
                'weightName','Affinity',...
                'weightField1','transcriptionFactorProteinComplexAffinitys',...
                'weightField2','',...
                'weight2Name','Activity',...
                'weight2Field1','transcriptionFactorProteinComplexActivitys',...
                'weight2Field2','',...
                'weight3Name','Condition',...
                'weight3Field1','transcriptionFactorProteinComplexConditions',...
                'weight3Field2','',...
                'weight4Name','BindingSiteStartCoordinate',...
                'weight4Field1','transcriptionFactorProteinComplexBindingSiteStartCoordinates',...
                'weight4Field2','',...
                'weight5Name','BindingSiteLength',...
                'weight5Field1','transcriptionFactorProteinComplexBindingSiteLengths',...
                'weight5Field2','',...
                'weight6Name','BindingSiteDirection',...
                'weight6Field1','transcriptionFactorProteinComplexBindingSiteDirections',...
                'weight6Field2','',...
                'compartmentName','ProteinComplexCompartmentIdx',...
                'compartmentField1','transcriptionFactorProteinComplexCompartments',...
                'compartmentField2','',...
                'nanValue',NaN));

            %protein - stimuli
            this.linkObjects(database,'proteinMonomer','stimuli',struct(...
                'sqlProcedure','get_proteinactivations_proteinmonomers_stimulis',...
                'name2','StimulusRegulatorIdx',...
                'field1','stimuliRegulators',...
                'field2','regulatedProteinMonomers'));
            this.linkObjects(database,'proteinComplex','stimuli',struct(...
                'sqlProcedure','get_proteinactivations_proteincomplexs_stimulis',...
                'name2','StimulusRegulatorIdx',...
                'field1','stimuliRegulators',...
                'field2','regulatedProteinComplexs'));

            this.linkObjects(database,'proteinMonomer','metabolite',struct(...
                'sqlProcedure','get_proteinactivations_proteinmonomers_metabolites',...
                'name2','MetaboliteRegulatorIdx',...
                'field1','metaboliteRegulators',...
                'field2','regulatedProteinMonomers'));
            this.linkObjects(database,'proteinComplex','metabolite',struct(...
                'sqlProcedure','get_proteinactivations_proteincomplexs_metabolites',...
                'name2','MetaboliteRegulatorIdx',...
                'field1','metaboliteRegulators',...
                'field2','regulatedProteinComplexs'));

            this.linkObjects(database,'proteinMonomer','proteinMonomer',struct(...
                'sqlProcedure','get_proteinactivations_proteinmonomers_proteinmonomers',...
                'name2','ProteinMonomerRegulatorIdx',...
                'field1','proteinMonomerRegulators',...
                'field2','regulatedProteinMonomers'));
            this.linkObjects(database,'proteinComplex','proteinMonomer',struct(...
                'sqlProcedure','get_proteinactivations_proteincomplexs_proteinmonomers',...
                'name2','ProteinMonomerRegulatorIdx',...
                'field1','proteinMonomerRegulators',...
                'field2','regulatedProteinComplexs'));

            this.linkObjects(database,'proteinMonomer','proteinComplex',struct(...
                'sqlProcedure','get_proteinactivations_proteinmonomers_proteincomplexs',...
                'name2','ProteinComplexRegulatorIdx',...
                'field1','proteinComplexRegulators',...
                'field2','regulatedProteinMonomers'));
            this.linkObjects(database,'proteinComplex','proteinComplex',struct(...
                'sqlProcedure','get_proteinactivations_proteincomplexs_proteincomplexs',...
                'name2','ProteinComplexRegulatorIdx',...
                'field1','proteinComplexRegulators',...
                'field2','regulatedProteinComplexs'));

            %biomass composition
            this.linkObjects(database,'metabolite','compartment',struct(...
                'sqlProcedure','get_metabolites_biomasscompositions',...
                'field1','biomassCompartments',...
                'field2','',...
                'weightName','Coefficient',...
                'weightField1','biomassCoefficients',...
                'weightField2',''));

            %media composition
            this.linkObjects(database,'metabolite','compartment',struct(...
                'sqlProcedure','get_metabolites_mediacomponents',...
                'field1','mediaCompartments',...
                'field2','',...
                'weightName','Concentration',...
                'weightField1','mediaConcentrations',...
                'weightField2','',...
                'initialTimeName','InitialTime',...
                'initialTimeField1','mediaInitialTimes',...
                'initialTimeField2','',...
                'finalTimeName','FinalTime',...
                'finalTimeField1','mediaFinalTimes',...
                'finalTimeField2','',...
                'nanValue',Inf));

            %stimuli values
            this.linkObjects(database,'stimuli','compartment',struct(...
                'sqlProcedure','get_stimulis_values',...
                'field1','compartments',...
                'field2','',...
                'weightName','Value',...
                'weightField1','values',...
                'weightField2','',...
                'initialTimeName','InitialTime',...
                'initialTimeField1','initialTimes',...
                'initialTimeField2','',...
                'finalTimeName','FinalTime',...
                'finalTimeField1','finalTimes',...
                'finalTimeField2','',...
                'nanValue',Inf));

            %reactions-process, state
            this.linkObjects(database,'reaction','process',struct(...
                'sqlProcedure','get_reactions_processs'));
            this.linkObjects(database,'reaction','state',struct(...
                'field1','state'));

            %complexes-process
            this.linkObjects(database,'proteinComplex','process',struct(...
                'sqlProcedure','get_proteincomplexs_processs',...
                'field1','complexFormationProcess',...
                'field2',''));

            %reactions-pathways
            this.linkObjects(database,'reaction','pathway');

            %references
            this.linkObjects(database,'parameter','reference');
            this.linkObjects(database,'process','reference', struct(...
                'sqlProcedure','get_processs_references'));
            this.linkObjects(database,'state','reference');
            this.linkObjects(database,'compartment','reference');
            this.linkObjects(database,'metabolite','reference');
            this.linkObjects(database,'gene','reference');
            this.linkObjects(database,'transcriptionUnit','reference');
            this.linkObjects(database,'genomeFeature','reference');
            this.linkObjects(database,'proteinMonomer','reference');
            this.linkObjects(database,'proteinComplex','reference');
            this.linkObjects(database,'reaction','reference');
            this.linkObjects(database,'pathway','reference');
            this.linkObjects(database,'stimuli','reference');
            this.linkObjects(database,'note','reference');
            
            this.areLinksSerialized = false;
        end

        function linkObjects(this, database, objectName1, objectName2, options)
            %process options
            if ~exist('options', 'var') || ~isstruct(options)
                options = struct;
            end

            if isfield(options, 'nanValue')
                nanValue = options.nanValue;
            else
                nanValue = 0;
            end

            if isfield(options, 'name1')
                name1 = options.name1;
            else
                name1 = [upper(objectName1(1)) objectName1(2:end) 'Idx'];
            end
            if isfield(options, 'name2')
                name2 = options.name2;
            else
                name2 = [upper(objectName2(1)) objectName2(2:end) 'Idx'];
            end
            
            if objectName1(end) == 's'
                objectName1 = [objectName1 'es'];
            else
                objectName1 = [objectName1 's'];
            end
            if objectName2(end) == 's'
                objectName2 = [objectName2 'es'];
            else
                objectName2 = [objectName2 's'];
            end

            if isfield(options, 'field1')
                field1 = options.field1;
            else
                field1 = objectName2;
            end
            if isfield(options, 'field2')
                field2 = options.field2;
            else
                field2 = objectName1;
            end

            if isfield(options, 'sqlProcedure')
                sqlProcedure = options.sqlProcedure;
            else
                sqlProcedure = sprintf('get_%s_%s', objectName1, objectName2);
            end

            database.setNullValue(nanValue);
            database.prepareStatement(sprintf('call %s(%d,null)', sqlProcedure, this.wid));
            data = database.query();
            objects1 = this.(objectName1);
            objects2 = this.(objectName2);

            data1 = data.(name1);
            data2 = data.(name2);
            if isfield(options, 'compartmentName')
                compartmentData = data.(options.compartmentName);
            end
            if isfield(options, 'weightName')
                weights = data.(options.weightName);
            end
            if isfield(options, 'weight2Name')
                weight2s = data.(options.weight2Name);
            end
            if isfield(options, 'weight3Name')
                weight3s = data.(options.weight3Name);
            end
            if isfield(options, 'weight4Name')
                weight4s = data.(options.weight4Name);
            end
            if isfield(options, 'weight5Name')
                weight5s = data.(options.weight5Name);
            end
            if isfield(options, 'weight6Name')
                weight6s = data.(options.weight6Name);
            end
            if isfield(options, 'initialTimeName')
                initialTimes = data.(options.initialTimeName);
            end
            if isfield(options, 'finalTimeName')
                finalTimes = data.(options.finalTimeName);
            end

            for j=1:length(objects1)
                if isempty(findprop(objects1(j),field1)); continue; end;
                if ~isempty(objects1(j).(field1)); continue; end;
                objects1(j).(field1) = objects2(data2(data1 == j));
                if isfield(options,'compartmentName') && isfield(options,'compartmentField1') && ~isempty(options.compartmentField1)
                    objects1(j).(options.compartmentField1) = this.compartments(compartmentData(data1 == j));
                end
                if isfield(options,'weightName') && isfield(options,'weightField1') && ~isempty(options.weightField1)
                    objects1(j).(options.weightField1) = weights(data1 == j);
                end
                if isfield(options,'weight2Name') && isfield(options,'weight2Field1') && ~isempty(options.weight2Field1)                  
                    objects1(j).(options.weight2Field1) = weight2s(data1 == j);
                end
                if isfield(options,'weight3Name') && isfield(options,'weight3Field1') && ~isempty(options.weight3Field1)       
                    objects1(j).(options.weight3Field1) = weight3s(data1 == j);
                end
                if isfield(options,'weight4Name') && isfield(options,'weight4Field1') && ~isempty(options.weight4Field1)       
                    objects1(j).(options.weight4Field1) = weight4s(data1 == j);
                end
                if isfield(options,'weight5Name') && isfield(options,'weight5Field1') && ~isempty(options.weight5Field1)       
                    objects1(j).(options.weight5Field1) = weight5s(data1 == j);
                end
                if isfield(options,'weight6Name') && isfield(options,'weight6Field1') && ~isempty(options.weight6Field1)       
                    objects1(j).(options.weight6Field1) = weight6s(data1 == j);
                end
                if isfield(options,'initialTimeName') && isfield(options,'initialTimeField1') && ~isempty(options.initialTimeField1)
                    objects1(j).(options.initialTimeField1) = initialTimes(data1 == j);
                end
                if isfield(options,'finalTimeName') && isfield(options,'finalTimeField1') && ~isempty(options.finalTimeField1)
                    objects1(j).(options.finalTimeField1) = finalTimes(data1 == j);
                end
            end

            for j = 1:length(objects2)
                if isempty(findprop(objects2(j),field2)); continue; end;
                if ~isempty(objects2(j).(field2)); continue; end;
                objects2(j).(field2) = objects1(data1(data2 == j));
                if isfield(options,'compartmentName') && isfield(options,'compartmentField2') && ~isempty(options.compartmentField2)
                    objects2(j).(options.compartmentField2) = this.compartments(compartmentData(data2 == j));
                end
                if isfield(options,'weightName') && isfield(options,'weightField2') && ~isempty(options.weightField2)
                    objects2(j).(options.weightField2) = weights(data2 == j);
                end
                if isfield(options,'weight2Name') && isfield(options,'weight2Field2') && ~isempty(options.weight2Field2)
                    objects2(j).(options.weight2Field2) = weight2s(data2 == j);
                end
                if isfield(options,'weight3Name') && isfield(options,'weight3Field2') && ~isempty(options.weight3Field2)
                    objects2(j).(options.weight3Field2) = weight3s(data2 == j);
                end
                if isfield(options,'weight4Name') && isfield(options,'weight4Field2') && ~isempty(options.weight4Field2)
                    objects2(j).(options.weight4Field2) = weight4s(data2 == j);
                end
                if isfield(options,'weight5Name') && isfield(options,'weight5Field2') && ~isempty(options.weight5Field2)
                    objects2(j).(options.weight5Field2) = weight5s(data2 == j);
                end
                if isfield(options,'weight6Name') && isfield(options,'weight6Field2') && ~isempty(options.weight6Field2)
                    objects2(j).(options.weight6Field2) = weight6s(data2 == j);
                end
                if isfield(options,'initialTimeName') && isfield(options,'initialTimeField2') && ~isempty(options.initialTimeField2)
                    objects2(j).(options.initialTimeField2) = initialTimes(data2 == j);
                end
                if isfield(options,'finalTimeName') && isfield(options,'finalTimeField2') && ~isempty(options.finalTimeField2)
                    objects2(j).(options.finalTimeField2) = finalTimes(data2 == j);
                end
            end
        end
    end

    methods
        %invalidate cached properties
        function invalidate(this)
            this.invalidate@edu.stanford.covert.cell.kb.KnowledgeBaseObject();
            
            this.parameters.invalidate();
            this.processes.invalidate();
            this.states.invalidate();
            this.compartments.invalidate();
            this.metabolites.invalidate();
            this.genome.invalidate();
            this.genes.invalidate();
            this.transcriptionUnits.invalidate();
            this.genomeFeatures.invalidate();
            this.proteinMonomers.invalidate();
            this.proteinComplexs.invalidate();
            this.reactions.invalidate();
            this.pathways.invalidate();
            this.stimulis.invalidate();
            this.notes.invalidate();
        end
        
        function calcIndices(this)
            this.calcIndices@edu.stanford.covert.cell.kb.KnowledgeBaseObject();
            
            this.parameters.calcIndices();
            this.processes.calcIndices();
            this.states.calcIndices();
            this.compartments.calcIndices();
            this.metabolites.calcIndices();
            this.genome.calcIndices();
            this.genes.calcIndices();
            this.transcriptionUnits.calcIndices();
            this.genomeFeatures.calcIndices();
            this.proteinMonomers.calcIndices();
            this.proteinComplexs.calcIndices();
            this.reactions.calcIndices();
            this.pathways.calcIndices();
            this.stimulis.calcIndices();
            this.notes.calcIndices();
        end
        
        %compute dependent properties of metabolites, RNAs, proteins
        %- indices
        %- base counts
        %- molecular weight
        function computeDependentProperties(this)            
            %% cardinalities
            this.numProcesses          = numel(this.processes);
            this.numStates             = numel(this.states);
            this.numParameters         = numel(this.parameters);
            this.numCompartments       = numel(this.compartments);
            this.numMetabolites        = numel(this.metabolites);
            this.numGenes              = numel(this.genes);
            this.numTranscriptionUnits = numel(this.transcriptionUnits);
            this.numGenomeFeatures     = numel(this.genomeFeatures);
            this.numProteinMonomers    = numel(this.proteinMonomers);
            this.numProteinComplexs    = numel(this.proteinComplexs);
            this.numReactions          = numel(this.reactions);
            this.numPathways           = numel(this.pathways);
            this.numStimulis           = numel(this.stimulis);
            this.numNotes              = numel(this.notes);
            this.numReferences         = numel(this.references);
            
            %% indices
            %compartments
            wholeCellModelIDs = {this.compartments.wholeCellModelID};
            [~,this.cytosolCompartmentIndexs]                   = ismember(this.cytosolCompartmentWholeCellModelIDs, wholeCellModelIDs);
            [~,this.chromosomeCompartmentIndexs]                = ismember(this.chromosomeCompartmentWholeCellModelIDs, wholeCellModelIDs);
            [~,this.membraneCompartmentIndexs]                  = ismember(this.membraneCompartmentWholeCellModelIDs, wholeCellModelIDs);
            [~,this.terminalOrganelleCytosolCompartmentIndexs]  = ismember(this.terminalOrganelleCytosolCompartmentWholeCellModelIDs, wholeCellModelIDs);
            [~,this.terminalOrganelleMembraneCompartmentIndexs] = ismember(this.terminalOrganelleMembraneCompartmentWholeCellModelIDs, wholeCellModelIDs);
            [~,this.extracellularCompartmentIndexs]             = ismember(this.extracellularCompartmentWholeCellModelIDs, wholeCellModelIDs);
            this.cellularCompartmentIndexs                      = setdiff((1:this.numCompartments)', this.extracellularCompartmentIndexs);

            %metabolites
            wholeCellModelIDs = {this.metabolites.wholeCellModelID}';
            [~, this.waterIndexs]                       = ismember({'H2O'}, wholeCellModelIDs);
            [~, this.hydrogenIndexs]                    = ismember({'H'}, wholeCellModelIDs);
            [~, this.dnmpIndexs]                        = ismember({'DAMP';'DCMP';'DGMP';'DTMP'}, wholeCellModelIDs);
            [~, this.nmpIndexs]                         = ismember({'AMP';'CMP';'GMP';'UMP'}, wholeCellModelIDs);
            [~, this.modifiedNMPIndexs]                 = ismember({'PSIURIMP';'m5CMP';'m62AMP';'m2GMP';'m7GMP';'GmMP';'UmMP';'s2UMP';'m1GMP';'k2CMP';'s4UMP';'cmnm5s2UMP';'m6A MP'}, wholeCellModelIDs);
            [~, this.aminoAcidIndexs]                   = ismember({'ALA'; 'ARG'; 'ASN'; 'ASP'; 'CYS'; 'GLN'; 'GLU'; 'GLY'; 'HIS'; 'ILE'; 'LEU'; 'LYS'; 'MET'; 'PHE'; 'PRO'; 'SER'; 'THR'; 'TRP'; 'TYR'; 'VAL'; 'FMET'}, wholeCellModelIDs);
            [~, this.cysteineIndexs]                    = ismember({'CYS'}, wholeCellModelIDs);
            [~, this.methionineIndexs]                  = ismember({'MET'}, wholeCellModelIDs);
            [~, this.fmethionineIndexs]                 = ismember({'FMET'}, wholeCellModelIDs);
            [~, this.modifiedAminoAcidIndexs]           = ismember({'diacylglycerolCys';'LIPOYLLYS'; 'pSER'; 'pTHR'; 'pTYR'}, wholeCellModelIDs);
            [~, this.diacylglycerolCysteineIndexs]      = ismember({'diacylglycerolCys'}, wholeCellModelIDs);

            %genes
            this.mRNAGenes = findobj(this.genes, 'type', 'mRNA');
            this.rRNAGenes = findobj(this.genes, 'type', 'rRNA');
            this.sRNAGenes = findobj(this.genes, 'type', 'sRNA');
            this.tRNAGenes = findobj(this.genes, 'type', 'tRNA');

            %RNAs
            names = {this.genes.name}';
            [~,this.ribosomalRRNAIndexs] = ismember({'5S ribosomal rRNA';'16S ribosomal rRNA';'23S ribosomal rRNA'}, names);
            
            %% other dependent properties
            %metabolites
            this.metaboliteMolecularWeights  = [this.metabolites.molecularWeight]';
        end

        function delete(this)
            %deletes pointers between objects
            %(early versions of MATLAB (eg R2010b) seem to have difficulty
            %clearing objects when these links are present. By manually
            %removing the links, we improve the run time of deleting
            %instances of this object to 10-15s)
            this.deleteLinks();

            %super class delete method
            this.delete@edu.stanford.covert.cell.kb.KnowledgeBaseObject();
        end

        %deletes pointers between objects
        function serializeLinks(this)
            if this.areLinksSerialized
               return;
            end
            
            this.areLinksSerialized = true;
            
            this.parameters.serializeLinks();
            this.processes.serializeLinks();
            this.states.serializeLinks();
            this.compartments.serializeLinks();
            this.metabolites.serializeLinks();
            this.genome.serializeLinks();
            this.genes.serializeLinks();
            this.transcriptionUnits.serializeLinks();
            this.genomeFeatures.serializeLinks();
            this.proteinMonomers.serializeLinks();
            this.proteinComplexs.serializeLinks();
            this.reactions.serializeLinks();
            this.pathways.serializeLinks();
            this.stimulis.serializeLinks();
            this.notes.serializeLinks();
            
            this.knowledgeBase = {class(this); this.idx};

            %serializeLinks@edu.stanford.covert.cell.kb.KnowledgeBaseObject(this);
        end
        
        function deserializeLinks(this)
            if ~this.areLinksSerialized
                return
            end
            
            this.areLinksSerialized = false;
            
            this.parameters.deserializeLinks(this);
            this.processes.deserializeLinks(this);
            this.states.deserializeLinks(this);
            this.compartments.deserializeLinks(this);
            this.metabolites.deserializeLinks(this);
            this.genome.deserializeLinks(this);
            this.genes.deserializeLinks(this);
            this.transcriptionUnits.deserializeLinks(this);
            this.genomeFeatures.deserializeLinks(this);
            this.proteinMonomers.deserializeLinks(this);
            this.proteinComplexs.deserializeLinks(this);
            this.reactions.deserializeLinks(this);
            this.pathways.deserializeLinks(this);
            this.stimulis.deserializeLinks(this);
            this.notes.deserializeLinks(this);
            
            this.knowledgeBase = this;
            %deserializeLinks@edu.stanford.covert.cell.kb.KnowledgeBaseObject(this);
        end
        
        function deleteLinks(this)
            this.parameters.deleteLinks();
            this.processes.deleteLinks();
            this.states.deleteLinks();
            this.compartments.deleteLinks();
            this.metabolites.deleteLinks();
            this.genome.deleteLinks();
            this.genes.deleteLinks();
            this.transcriptionUnits.deleteLinks();
            this.genomeFeatures.deleteLinks();
            this.proteinMonomers.deleteLinks();
            this.proteinComplexs.deleteLinks();
            this.reactions.deleteLinks();
            this.pathways.deleteLinks();
            this.stimulis.deleteLinks();
            this.notes.deleteLinks();
            
            deleteLinks@edu.stanford.covert.cell.kb.KnowledgeBaseObject(this);
        end
    end
end