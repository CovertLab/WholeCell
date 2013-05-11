% Defines a transcription unit
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/7/2009
classdef TranscriptionUnit < edu.stanford.covert.cell.kb.ssRNA
    properties
        genome           = edu.stanford.covert.cell.kb.Genome.empty(0,0);
        
        genes            = edu.stanford.covert.cell.kb.Gene.empty(0,0);
        geneCompartments = edu.stanford.covert.cell.kb.Compartment.empty(0,0);
        
        compartment      = edu.stanford.covert.cell.kb.Compartment.empty(0,0);
        
        transcriptionFactorProteinMonomers            = edu.stanford.covert.cell.kb.ProteinMonomer.empty(0,0);        
        transcriptionFactorProteinMonomerAffinitys    = []; %nM
        transcriptionFactorProteinMonomerActivitys    = [];
        transcriptionFactorProteinMonomerConditions   = [];
        transcriptionFactorProteinMonomerBindingSiteStartCoordinates = [];
        transcriptionFactorProteinMonomerBindingSiteLengths          = [];
        transcriptionFactorProteinMonomerBindingSiteDirections       = [];
        transcriptionFactorProteinMonomerCompartments = edu.stanford.covert.cell.kb.Compartment.empty(0,0);
        
        transcriptionFactorProteinComplexs            = edu.stanford.covert.cell.kb.ProteinComplex.empty(0,0);        
        transcriptionFactorProteinComplexAffinitys    = []; %nM
        transcriptionFactorProteinComplexActivitys    = [];
        transcriptionFactorProteinComplexConditions   = [];
        transcriptionFactorProteinComplexBindingSiteStartCoordinates = [];
        transcriptionFactorProteinComplexBindingSiteLengths          = [];
        transcriptionFactorProteinComplexBindingSiteDirections       = [];
        transcriptionFactorProteinComplexCompartments = edu.stanford.covert.cell.kb.Compartment.empty(0,0);
    end
    
    properties %(SetAccess = protected)
        halfLife
        sequence
        promoter35Coordinate
        promoter35Length
        promoter10Coordinate
        promoter10Length
        tssCoordinate
    end
    
    %computed properties
    properties %(SetAccess = protected)
        type
        dnaSequence
        startCoordinate
        endCoordinate
        fivePrimeCoordinate
        threePrimeCoordinate
        direction
        synthesisRate
        expression
        
        promoter1Sequence
        promoter10Sequence
        promoter35Sequence
        
        intergenicSequences
        intergenicSequenceLengths
        intergenicSequenceBaseCounts
        intergenicSequenceDecayReactions
        intergenicSequenceMolecularWeights
        intergenicSequenceHalfLives
        
        numGenes
        features
    end
    
    methods
        function this = TranscriptionUnit(knowledgeBase, wid,wholeCellModelID,name,...
                promoter35Coordinate, promoter35Length, ...
                promoter10Coordinate, promoter10Length, ...
                tssCoordinate, ...
                comments,crossReferences)
            
            if nargin == 0; return; end;
            
            this = edu.stanford.covert.cell.kb.TranscriptionUnit.empty(size(wid,1),0);
            this(size(wid,1),1) = edu.stanford.covert.cell.kb.TranscriptionUnit;
            for i = 1:size(wid,1)
                this(i,1).idx = i;
                this(i,1).knowledgeBase = knowledgeBase;
                this(i,1).wid = wid(i);
                this(i,1).wholeCellModelID = wholeCellModelID{i};
                this(i,1).name = name{i};
                if exist('comments','var') && ~isempty(comments); this(i,1).comments = comments{i}; end;
                if exist('crossReferences','var')
                    if size(crossReferences,1)>1
                        this(i,1).crossReferences = crossReferences(i);
                    else
                        this(i,1).crossReferences = struct;
                        fields = fieldnames(crossReferences);
                        for j = 1:size(fields,1)
                            values = crossReferences.(fields{j});
                            this(i,1).crossReferences.(fields{j}) = values(i);
                        end
                    end
                end
                
                this(i,1).promoter35Coordinate = promoter35Coordinate(i);
                this(i,1).promoter35Length     = promoter35Length(i);
                this(i,1).promoter10Coordinate = promoter10Coordinate(i);
                this(i,1).promoter10Length     = promoter10Length(i);
                this(i,1).tssCoordinate        = tssCoordinate(i);
            end
        end
        
        function serializeLinks(this)
            for i = 1:numel(this)
                this(i).genome                                        = this.serializeLinksHelper(this(i).genome);
                
                this(i).genes                                         = this.serializeLinksHelper(this(i).genes);
                this(i).geneCompartments                              = this.serializeLinksHelper(this(i).geneCompartments);
                
                this(i).compartment                                   = this.serializeLinksHelper(this(i).compartment);
                
                this(i).transcriptionFactorProteinMonomers            = this.serializeLinksHelper(this(i).transcriptionFactorProteinMonomers);
                this(i).transcriptionFactorProteinMonomerAffinitys    = this.serializeLinksHelper(this(i).transcriptionFactorProteinMonomerAffinitys);
                this(i).transcriptionFactorProteinMonomerActivitys    = this.serializeLinksHelper(this(i).transcriptionFactorProteinMonomerActivitys);
                this(i).transcriptionFactorProteinMonomerConditions   = this.serializeLinksHelper(this(i).transcriptionFactorProteinMonomerConditions);
                this(i).transcriptionFactorProteinMonomerCompartments = this.serializeLinksHelper(this(i).transcriptionFactorProteinMonomerCompartments);
                
                this(i).transcriptionFactorProteinComplexs            = this.serializeLinksHelper(this(i).transcriptionFactorProteinComplexs);
                this(i).transcriptionFactorProteinComplexAffinitys    = this.serializeLinksHelper(this(i).transcriptionFactorProteinComplexAffinitys);
                this(i).transcriptionFactorProteinComplexActivitys    = this.serializeLinksHelper(this(i).transcriptionFactorProteinComplexActivitys);
                this(i).transcriptionFactorProteinComplexConditions   = this.serializeLinksHelper(this(i).transcriptionFactorProteinComplexConditions);
                this(i).transcriptionFactorProteinComplexCompartments = this.serializeLinksHelper(this(i).transcriptionFactorProteinComplexCompartments);
                
                serializeLinks@edu.stanford.covert.cell.kb.ssRNA(this(i));
            end
        end
        
        function deserializeLinks(this, kb)
            for i = 1:numel(this)
                this(i).genome                                        = this.deserializeLinksHelper(this(i).genome, kb.genome);
                
                this(i).genes                                         = this.deserializeLinksHelper(this(i).genes, kb.genes);
                this(i).geneCompartments                              = this.deserializeLinksHelper(this(i).geneCompartments, kb.compartments);
                
                this(i).compartment                                   = this.deserializeLinksHelper(this(i).compartment, kb.compartments);
                
                this(i).transcriptionFactorProteinMonomers            = this.deserializeLinksHelper(this(i).transcriptionFactorProteinMonomers, kb.proteinMonomers);
                this(i).transcriptionFactorProteinMonomerAffinitys    = this.deserializeLinksHelper(this(i).transcriptionFactorProteinMonomerAffinitys);
                this(i).transcriptionFactorProteinMonomerActivitys    = this.deserializeLinksHelper(this(i).transcriptionFactorProteinMonomerActivitys);
                this(i).transcriptionFactorProteinMonomerConditions   = this.deserializeLinksHelper(this(i).transcriptionFactorProteinMonomerConditions);
                this(i).transcriptionFactorProteinMonomerCompartments = this.deserializeLinksHelper(this(i).transcriptionFactorProteinMonomerCompartments, kb.compartments);
                
                this(i).transcriptionFactorProteinComplexs            = this.deserializeLinksHelper(this(i).transcriptionFactorProteinComplexs, kb.proteinComplexs);
                this(i).transcriptionFactorProteinComplexAffinitys    = this.deserializeLinksHelper(this(i).transcriptionFactorProteinComplexAffinitys);
                this(i).transcriptionFactorProteinComplexActivitys    = this.deserializeLinksHelper(this(i).transcriptionFactorProteinComplexActivitys);
                this(i).transcriptionFactorProteinComplexConditions   = this.deserializeLinksHelper(this(i).transcriptionFactorProteinComplexConditions);
                this(i).transcriptionFactorProteinComplexCompartments = this.deserializeLinksHelper(this(i).transcriptionFactorProteinComplexCompartments, kb.compartments);
                
                deserializeLinks@edu.stanford.covert.cell.kb.ssRNA(this(i), kb);
            end
        end

        function deleteLinks(this)
            for i = 1:numel(this)
                this(i).genome                                        = [];
                
                this(i).genes                                         = [];
                this(i).geneCompartments                              = [];
                
                this(i).compartment                                   = [];
                
                this(i).transcriptionFactorProteinMonomers            = [];
                this(i).transcriptionFactorProteinMonomerAffinitys    = [];
                this(i).transcriptionFactorProteinMonomerActivitys    = [];
                this(i).transcriptionFactorProteinMonomerConditions   = [];
                this(i).transcriptionFactorProteinMonomerCompartments = [];
                
                this(i).transcriptionFactorProteinComplexs            = [];
                this(i).transcriptionFactorProteinComplexAffinitys    = [];
                this(i).transcriptionFactorProteinComplexActivitys    = [];
                this(i).transcriptionFactorProteinComplexConditions   = [];
                this(i).transcriptionFactorProteinComplexCompartments = [];
                
                deleteLinks@edu.stanford.covert.cell.kb.ssRNA(this(i));
            end
        end
        
        function value = get.halfLife(this)
            %retrieve
            if ~isempty(this.halfLife)
                value = this.halfLife;
                return;
            end
            
            %compute
            value = mean([this.genes.halfLife]);
            
            %store
            this.halfLife = value;
        end
        
        function value = get.type(this)
            %retrieve
            if ~isempty(this.type)
                value = this.type;
                return;
            end
            
            %compute
            if ~isempty(this.genes)
                value = this.genes(1).type;
            else
                value = [];
            end
            
            %store
            this.type = value;
        end
        
        function value = get.dnaSequence(this)
            %retrieve
            if ~isempty(this.dnaSequence)
                value = this.dnaSequence;
                return;
            end
            
            %compute
            if this.direction
                value = this.genome.sequence(this.startCoordinate:this.endCoordinate);
            else
                value = seqrcomplement(this.genome.sequence(this.startCoordinate:this.endCoordinate));
            end
            
            %store
            this.dnaSequence = value;
        end
        
        function value = get.sequence(this)
            %retrieve
            if ~isempty(this.sequence)
                value = this.sequence;
                return;
            end
            
            %compute
            value = this.dnaSequence;
            value(this.dnaSequence == 'A') = 'U';
            value(this.dnaSequence == 'C') = 'G';
            value(this.dnaSequence == 'G') = 'C';
            value(this.dnaSequence == 'T') = 'A';
            
            %store
            this.sequence = value;
        end
        
        function value = get.startCoordinate(this)
            %retrieve
            if ~isempty(this.startCoordinate)
                value = this.startCoordinate;
                return;
            end
            
            %compute
            value = min([this.genes.startCoordinate]);
            
            %store
            this.startCoordinate = value;
        end
        
        function value = get.endCoordinate(this)
            %retrieve
            if ~isempty(this.endCoordinate)
                value = this.endCoordinate;
                return;
            end
            
            %compute
            value = max([this.genes.endCoordinate]);
            
            %store
            this.endCoordinate = value;
        end
        
        function value = get.fivePrimeCoordinate(this)
            %retrieve
            if ~isempty(this.fivePrimeCoordinate)
                value = this.fivePrimeCoordinate;
                return;
            end
            
            %compute
            if this.direction
                value = this.startCoordinate;
            else
                value = this.endCoordinate;
            end
            
            %store
            this.fivePrimeCoordinate = value;
        end
        
        function value = get.threePrimeCoordinate(this)
            %retrieve
            if ~isempty(this.threePrimeCoordinate)
                value = this.threePrimeCoordinate;
                return;
            end
            
            %compute
            if this.direction
                value = this.endCoordinate;
            else
                value = this.startCoordinate;
            end
            
            %store
            this.threePrimeCoordinate = value;
        end
        
        function value = get.direction(this)
            %retrieve
            if ~isempty(this.direction)
                value = this.direction;
                return;
            end
            
            %compute
            value = this.genes(1).direction;
            
            %store
            this.direction = value;
        end
        
        function value = get.synthesisRate(this)
            %retrieve
            if ~isempty(this.synthesisRate)
                value = this.synthesisRate;
                return;
            end
            
            %compute
            value = mean(reshape([this.genes.synthesisRate],[],this.numGenes),2)';
            
            %store
            this.synthesisRate = value;
        end
        
        function value = get.expression(this)
            %retrieve
            if ~isempty(this.expression)
                value = this.expression;
                return;
            end
            
            %compute
            value = this.synthesisRate*this.halfLife/log(2);
            
            %store
            this.expression = value;
        end
        
        function value = get.promoter1Sequence(this)
            %retrieve
            if ~isempty(this.promoter1Sequence)
                value = this.promoter1Sequence;
                return;
            end
            
            %compute
            if this.direction
                value = this.genome.sequence(this.fivePrimeCoordinate + this.tssCoordinate + (-1:-1) + 1);
            else
                value = seqcomplement(this.genome.sequence(this.fivePrimeCoordinate - this.tssCoordinate + (1:-1:1) - 1));
            end
            
            %store
            this.promoter1Sequence = value;
        end
        
        function value = get.promoter10Sequence(this)
            %retrieve
            if ~isempty(this.promoter10Sequence)
                value = this.promoter10Sequence;
                return;
            end
            
            %compute
            if this.direction
                value = this.genome.sequence(this.fivePrimeCoordinate + this.promoter10Coordinate + (-this.promoter10Length:-1) + 1);
            else
                value = seqcomplement(this.genome.sequence(this.fivePrimeCoordinate - this.promoter10Coordinate + (this.promoter10Length:-1:1) - 1));
            end
            
            %store
            this.promoter10Sequence = value;
        end
        
        function value = get.promoter35Sequence(this)
            %retrieve
            if ~isempty(this.promoter35Sequence)
                value = this.promoter35Sequence;
                return;
            end
            
            %compute
            if this.direction
                value = this.genome.sequence(this.fivePrimeCoordinate + this.promoter35Coordinate + (-this.promoter35Length:-1) + 1);
            else
                value = seqcomplement(this.genome.sequence(this.fivePrimeCoordinate - this.promoter35Coordinate + (this.promoter35Length:-1:1) - 1));
            end
            
            %store
            this.promoter35Sequence = value;
        end
        
        function value = get.intergenicSequences(this)
            %retrieve
            if ~isempty(this.intergenicSequences)
                value = this.intergenicSequences;
                return;
            end
            
            %compute
            value = cell(0,1);
            for i = 1:numel(this.genes)-1
                if strcmp(this.type,'mRNA')
                    continue;
                end
                if this.genes(i).endCoordinate+1 > this.genes(i+1).startCoordinate-1
                    continue;
                end
                
                if this.direction
                    seq = this.genome.sequence(this.genes(i).endCoordinate+1:this.genes(i+1).startCoordinate-1);
                else
                    seq = seqrcomplement(this.genome.sequence(this.genes(i).endCoordinate+1:this.genes(i+1).startCoordinate-1));
                end
                
                value{i,1} = seq;
                value{i,1}(seq == 'A') = 'U';
                value{i,1}(seq == 'C') = 'G';
                value{i,1}(seq == 'G') = 'C';
                value{i,1}(seq == 'T') = 'A';
            end
            
            %store
            this.intergenicSequences = value;
        end
        
        function value = get.intergenicSequenceLengths(this)
            %retrieve
            if ~isempty(this.intergenicSequenceLengths)
                value = this.intergenicSequenceLengths;
                return;
            end
            
            %compute
            value = zeros(0,1);
            for i = 1:numel(this.genes)-1
                if strcmp(this.type,'mRNA')
                    continue;
                end
                if this.genes(i).endCoordinate+1 > this.genes(i+1).startCoordinate-1
                    continue;
                end
                
                value(i,1) = (this.genes(i+1).startCoordinate-1) - (this.genes(i).endCoordinate+1) + 1;
            end
            
            %store
            this.intergenicSequenceLengths = value;
        end
        
        function value = get.intergenicSequenceBaseCounts(this)
            %retrieve
            if ~isempty(this.intergenicSequenceBaseCounts)
                value = this.intergenicSequenceBaseCounts;
                return;
            end
            
            %compute
            value = zeros(numel(this.intergenicSequences), this.knowledgeBase.numMetabolites);
            
            for i = 1:numel(this.intergenicSequences)
                value(i, this.knowledgeBase.nmpIndexs) = [...
                    sum('A' == this.intergenicSequences{i}) ...
                    sum('C' == this.intergenicSequences{i}) ...
                    sum('G' == this.intergenicSequences{i}) ...
                    sum('U' == this.intergenicSequences{i})];
            end
            
            %store
            this.intergenicSequenceBaseCounts = value;
        end
        
        function value = get.intergenicSequenceDecayReactions(this)
            %retrieve
            if ~isempty(this.intergenicSequenceDecayReactions)
                value = this.intergenicSequenceDecayReactions;
                return;
            end
            
            %compute
            value = this.intergenicSequenceBaseCounts;            
            for i = 1:numel(this.intergenicSequences)
                value(i, this.knowledgeBase.waterIndexs)    = value(i, this.knowledgeBase.waterIndexs)    - max(0, length(this.intergenicSequences{i}) - 1);
                value(i, this.knowledgeBase.hydrogenIndexs) = value(i, this.knowledgeBase.hydrogenIndexs) + max(0, length(this.intergenicSequences{i}) - 1);
            end
            
            %store
            this.intergenicSequenceDecayReactions = value;
        end
        
        function value = get.intergenicSequenceMolecularWeights(this)
            %retrieve
            if ~isempty(this.intergenicSequenceMolecularWeights)
                value = this.intergenicSequenceMolecularWeights;
                return;
            end
            
            %compute
            value = zeros(size(this.intergenicSequences));
            
            for i = 1:numel(this.intergenicSequences)
                value(i) = this.calculateMolecularWeight(...
                    this.intergenicSequences{i}, length(this.intergenicSequences{i}), 'linear', ...
                    this.knowledgeBase.metaboliteMolecularWeights(this.knowledgeBase.nmpIndexs));
            end
            
            %store
            this.intergenicSequenceMolecularWeights = value;
        end
        
        function value = get.intergenicSequenceHalfLives(this)
            %retrieve
            if ~isempty(this.intergenicSequenceHalfLives)
                value = this.intergenicSequenceHalfLives;
                return;
            end
            
            %compute
            value = zeros(size(this.intergenicSequences));
            
            %store
            this.intergenicSequenceHalfLives = value;
        end
        
        function value = get.numGenes(this)
            %retrieve
            if ~isempty(this.numGenes)
                value = this.numGenes;
                return;
            end
            
            %compute
            value = length(this.genes);
            
            %store
            this.numGenes = value;
        end
        
        function value = get.features(this)
            %retrieve
            if ~isempty(this.features)
                value = this.features;
                return;
            end
            
            %compute
            featureStartCoordinates = [this.genome.features.startCoordinate];
            featureEndCoordinates = [this.genome.features.endCoordinate];
            value = this.genome.features(...
                (this.startCoordinate > featureStartCoordinates & this.startCoordinate < featureEndCoordinates) | ...
                (this.endCoordinate   > featureStartCoordinates & this.endCoordinate   < featureEndCoordinates) | ...
                (this.startCoordinate < featureStartCoordinates & this.endCoordinate   > featureEndCoordinates));
            
            %store
            this.features = value;
        end
    end
end