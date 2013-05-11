% Defines a gene
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/7/2009
classdef Gene < edu.stanford.covert.cell.kb.ssRNA

    properties
        genome                      = edu.stanford.covert.cell.kb.Genome.empty(0, 0);
        transcriptionUnits          = edu.stanford.covert.cell.kb.TranscriptionUnit.empty(0, 0);
        proteinMonomers             = edu.stanford.covert.cell.kb.ProteinMonomer.empty(0, 0);
        stableModificationReactions = edu.stanford.covert.cell.kb.Reaction.empty(0, 0);
        aminoAcid                   = edu.stanford.covert.cell.kb.Metabolite.empty(0, 0);
        compartment                 = edu.stanford.covert.cell.kb.Compartment.empty(0, 0);
    end

    properties %(SetAccess = protected)
        halfLife
        sequence

        symbol
        synonyms
        type
        startCodon
        codons
        startCoordinate
        endCoordinate
        direction
        essential
        expression
    end

    %computed properties
    properties %(SetAccess = protected)
        dnaSequence
        fivePrimeCoordinate
        threePrimeCoordinate
        synthesisRate

        matureBaseCount
        matureMolecularWeight
        matureDecayReaction

        aminoacylatedBaseCount
        aminoacylatedMolecularWeight
        aminoacylatedDecayReaction

        features
    end

    methods
        %constructor
        function this = Gene(knowledgeBase, wid, wholeCellModelID, name, ...
                symbol, synonyms, type, startCodon, codons, ...
                startCoordinate, sequenceLength, direction, essential, halfLife, ...
                expression, comments, crossReferences)

            if nargin == 0; return; end;

            this = edu.stanford.covert.cell.kb.Gene.empty(size(wid, 1), 0);
            this(size(wid, 1), 1) = edu.stanford.covert.cell.kb.Gene;
            for i = 1:size(wid, 1)
                this(i, 1).idx = i;
                this(i, 1).knowledgeBase = knowledgeBase;
                this(i, 1).wid = wid(i);
                this(i, 1).wholeCellModelID = wholeCellModelID{i};
                this(i, 1).name = name{i};
                if exist('comments', 'var') && ~isempty(comments); this(i, 1).comments = comments{i}; end;
                if exist('crossReferences', 'var')
                    if size(crossReferences,1) > 1
                        this(i, 1).crossReferences = crossReferences(i);
                    else
                        this(i, 1).crossReferences = struct;
                        fields = fieldnames(crossReferences);
                        for j = 1:size(fields, 1)
                            values = crossReferences.(fields{j});
                            this(i, 1).crossReferences.(fields{j}) = values(i);
                        end
                    end
                end

                this(i, 1).symbol = symbol{i};
                this(i, 1).synonyms = strsplit(';', synonyms{i});
                this(i, 1).type = type{i};
                this(i, 1).startCodon = startCodon(i);
                this(i, 1).codons = strsplit(';', codons{i});
                this(i, 1).startCoordinate = startCoordinate(i);
                this(i, 1).endCoordinate = startCoordinate(i) + sequenceLength(i) - 1;
                this(i, 1).direction = direction(i);
                this(i, 1).essential = essential{i};
                this(i, 1).expression = expression(i, :);
                this(i, 1).halfLife = halfLife(i);
            end
        end

        function serializeLinks(this)
            for i = 1:numel(this)
                this(i).genome                      = this.serializeLinksHelper(this(i).genome);
                this(i).transcriptionUnits          = this.serializeLinksHelper(this(i).transcriptionUnits);
                this(i).proteinMonomers             = this.serializeLinksHelper(this(i).proteinMonomers);
                this(i).stableModificationReactions = this.serializeLinksHelper(this(i).stableModificationReactions);
                this(i).aminoAcid                   = this.serializeLinksHelper(this(i).aminoAcid);
                this(i).compartment                 = this.serializeLinksHelper(this(i).compartment);

                serializeLinks@edu.stanford.covert.cell.kb.ssRNA(this(i));
            end
        end
        
        function deserializeLinks(this, kb)
            for i = 1:numel(this)
                this(i).genome                      = this.deserializeLinksHelper(this(i).genome, kb.genome);
                this(i).transcriptionUnits          = this.deserializeLinksHelper(this(i).transcriptionUnits, kb.transcriptionUnits);
                this(i).proteinMonomers             = this.deserializeLinksHelper(this(i).proteinMonomers, kb.proteinMonomers);
                this(i).stableModificationReactions = this.deserializeLinksHelper(this(i).stableModificationReactions, kb.reactions);
                this(i).aminoAcid                   = this.deserializeLinksHelper(this(i).aminoAcid, kb.metabolites);
                this(i).compartment                 = this.deserializeLinksHelper(this(i).compartment, kb.compartments);
                
                deserializeLinks@edu.stanford.covert.cell.kb.ssRNA(this(i), kb);
            end
        end

        function deleteLinks(this)
            for i = 1:numel(this)
                this(i).genome                      = [];
                this(i).transcriptionUnits          = [];
                this(i).proteinMonomers             = [];
                this(i).stableModificationReactions = [];
                this(i).aminoAcid                   = [];
                this(i).compartment                 = [];

                deleteLinks@edu.stanford.covert.cell.kb.ssRNA(this(i));
            end
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

        function value = get.matureBaseCount(this)
            %retrieve
            if ~isempty(this.matureBaseCount)
                value = this.matureBaseCount;
                return;
            end
            
            %compute
            value = this.baseCount;

            if isempty(this.stableModificationReactions)
                return;
            end

            modifications = zeros(size(value));
            for i = 1:length(this.stableModificationReactions)
                reaction = this.stableModificationReactions(i);

                for j = 1:length(reaction.metabolites)
                    if any(this.knowledgeBase.modifiedNMPIndexs == reaction.metabolites(j).idx) || ...
                            (any(this.knowledgeBase.nmpIndexs == reaction.metabolites(j).idx) && reaction.metaboliteCoefficients(j) < 0)

                        modifications([reaction.metabolites(j).idx]) = ...
                            modifications([reaction.metabolites(j).idx]) + ...
                            reaction.metaboliteCoefficients(j)';
                    end
                end
            end

            value = value + modifications;
            
            %store
            this.matureBaseCount = value;
        end

        function value = get.matureMolecularWeight(this)
            % import classes
            import edu.stanford.covert.util.ConstantUtil;
            
            %retrieve
            if ~isempty(this.matureMolecularWeight)
                value = this.matureMolecularWeight;
                return;
            end
            
            %compute
            molecularWeightHO = ConstantUtil.elements.H + ConstantUtil.elements.O;

            baseCount = this.matureBaseCount;

            metabolicMolecularWeights = this.knowledgeBase.metaboliteMolecularWeights;

            value = baseCount * metabolicMolecularWeights - ...
                molecularWeightHO * max(0, this.sequenceLength - strcmp(this.sequenceTopology, 'linear'));
            
            %store
            this.matureMolecularWeight = value;
        end

        function value = get.matureDecayReaction(this)
            %retrieve
            if ~isempty(this.matureDecayReaction)
                value = this.matureDecayReaction;
                return;
            end
            
            %compute
            value = this.matureBaseCount;
            value(this.knowledgeBase.waterIndexs)    = value(this.knowledgeBase.waterIndexs)    - max(0, this.sequenceLength - strcmp(this.sequenceTopology, 'linear'));
            value(this.knowledgeBase.hydrogenIndexs) = value(this.knowledgeBase.hydrogenIndexs) + max(0, this.sequenceLength - strcmp(this.sequenceTopology, 'linear'));
            
            %store
            this.matureDecayReaction = value;
        end

        function value = get.aminoacylatedBaseCount(this)
            %retrieve
            if ~isempty(this.aminoacylatedBaseCount)
                value = this.aminoacylatedBaseCount;
                return;
            end
            
            %compute
            value = this.matureBaseCount;
            if ~isempty(this.aminoAcid)
                value(this.aminoAcid.idx) = value(this.aminoAcid.idx) + 1;
            end
            
            %store
            this.aminoacylatedBaseCount = value;
        end

        function value = get.aminoacylatedMolecularWeight(this)
            import edu.stanford.covert.util.ConstantUtil;
            
            %retrieve
            if ~isempty(this.aminoacylatedMolecularWeight)
                value = this.aminoacylatedMolecularWeight;
                return;
            end
            
            %compute
            value = this.matureMolecularWeight;
            if ~isempty(this.aminoAcid)
                molecularWeightHO = ConstantUtil.elements.H + ConstantUtil.elements.O;
                value = value + this.aminoAcid.molecularWeight - molecularWeightHO;
            end
            
            %store
            this.aminoacylatedMolecularWeight = value;
        end

        function value = get.aminoacylatedDecayReaction(this)
            %retrieve
            if ~isempty(this.aminoacylatedDecayReaction)
                value = this.aminoacylatedDecayReaction;
                return;
            end
            
            %compute
            value = this.matureDecayReaction;
            if ~isempty(this.aminoAcid)
                value(this.aminoAcid.idx) = value(this.aminoAcid.idx) + 1;
                value(this.knowledgeBase.waterIndexs)    = value(this.knowledgeBase.waterIndexs)    - 1;
                value(this.knowledgeBase.hydrogenIndexs) = value(this.knowledgeBase.hydrogenIndexs) + 1;
            end
            
            %store
            this.aminoacylatedDecayReaction = value;
        end

        function value = get.synthesisRate(this)
            %retrieve
            if ~isempty(this.synthesisRate)
                value = this.synthesisRate;
                return;
            end
            
            %compute
            value = log(2) / this.halfLife * this.expression;
            
            %store
            this.synthesisRate = value;
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