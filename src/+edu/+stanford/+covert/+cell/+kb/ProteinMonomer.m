% Defines a protein monomer
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/7/2009
classdef ProteinMonomer < edu.stanford.covert.cell.kb.Polymer & edu.stanford.covert.cell.kb.Protein
    properties
        gene                             = edu.stanford.covert.cell.kb.Gene.empty(0,0);
        geneCompartments                 = edu.stanford.covert.cell.kb.Compartment.empty(0,0);
        
        prostheticGroups                 = edu.stanford.covert.cell.kb.Metabolite.empty(0,0);
        prostheticGroupCompartments      = edu.stanford.covert.cell.kb.Compartment.empty(0,0);
        prostheticGroupCoefficients      = [];
        
        chaperoneSubstrates              = edu.stanford.covert.cell.kb.ProteinMonomer.empty(0,0);
        chaperoneCompartments            = edu.stanford.covert.cell.kb.Compartment.empty(0,0);
        chaperoneCoefficients            = [];
        
        proteinComplexs                  = edu.stanford.covert.cell.kb.ProteinComplex.empty(0,0);
        
        stableModificationReactions      = edu.stanford.covert.cell.kb.Reaction.empty(0,0);
        
        regulatedTranscriptionUnits      = edu.stanford.covert.cell.kb.TranscriptionUnit.empty(0,0);
        
        compartment                      = edu.stanford.covert.cell.kb.Compartment.empty(0,0);
        
        stimuliRegulators                = edu.stanford.covert.cell.kb.Stimuli.empty(0,0);
        metaboliteRegulators             = edu.stanford.covert.cell.kb.Metabolite.empty(0,0);
        proteinMonomerRegulators         = edu.stanford.covert.cell.kb.ProteinMonomer.empty(0,0);
        proteinComplexRegulators         = edu.stanford.covert.cell.kb.ProteinComplex.empty(0,0);
        
        regulatedProteinMonomers         = edu.stanford.covert.cell.kb.ProteinMonomer.empty(0,0);
        regulatedProteinComplexs         = edu.stanford.covert.cell.kb.ProteinComplex.empty(0,0);
        
        parameters                       = edu.stanford.covert.cell.kb.Parameter.empty(0,0);
    end
    
    properties %(SetAccess = protected)
        empiricalFormula
        smiles
        charge
        pKa
        halfLife
        sequence
        
        topology
        activeSite
        metalBindingSite
        
        nTerminalMethionineCleavage
        signalSequenceType
        signalSequenceLocation
        signalSequenceLength
        
        activationRule
    end

    %computed properties
    properties %(SetAccess = protected)
        molecularWeight
        density
        volume
        pI
        extinctionCoefficient
        absorbanceFactor
        
        negAA
        posAA
        instabilityIndex
        stable
        aliphaticIndex
        gravy
        
        baseCount
        cumulativeBaseCount
        decayReaction
        tRNACount
        nTerminalAA
        
        disulfideBonds
        
        signalSequence
        signalSequenceBaseCount
        signalSequenceMolecularWeight
        signalSequenceHalfLife
        signalSequenceDecayReaction
        
        processedISequence
        processedISequenceLength
        processedISequenceBaseCount
        processedISequenceMolecularWeight
        processedISequenceHalfLife
        processedISequenceDecayReaction
        
        processedIISequence
        processedIISequenceLength
        processedIISequenceBaseCount
        processedIISequenceMolecularWeight
        processedIISequenceHalfLife
        processedIISequenceDecayReaction
        
        foldedSequence
        foldedSequenceLength
        foldedSequenceBaseCount
        foldedSequenceMolecularWeight
        foldedSequenceHalfLife
        foldedSequenceDecayReaction
        
        matureSequence
        matureSequenceLength
        matureSequenceBaseCount
        matureSequenceMolecularWeight
        matureSequenceHalfLife
        matureSequenceDecayReaction
    end

    properties (Constant = true)
        bases = 'ARNDCQEGHILKMFPSTWYV';

        %dipeptide instability weight value (DIWV)
        %1. http://ca.expasy.org/tools/protparam-doc.html
        %2. Kunchur Guruprasad, B.V.Bhasker Reddy and Madhusudan W.Pandit
        %   (1990). Correlation between stability of a protein and its
        %   dipeptide composition: a novel approach for predicting in vivo
        %   stability of a protein from its primary sequence. Protein
        %   Engineering 4(2):155-161. Table 3 copied to
        %   src/+edu/+stanford/+covert/+cell/+kb/DIWV.xlsx.
        DIWV = [
              1.00   1.00   1.00  -7.49 44.94  1.00  1.00   1.00 -7.49  1.00  1.00   1.00  1.00   1.00 20.26  1.00   1.00   1.00  1.00  1.00
              1.00  58.28  13.34   1.00  1.00 20.26  1.00  -7.49 20.26  1.00  1.00   1.00  1.00   1.00 20.26 44.94   1.00  58.28 -6.54  1.00
              1.00   1.00   1.00   1.00 -1.88 -6.54  1.00 -14.03  1.00 44.94  1.00  24.68  1.00 -14.03 -1.88  1.00  -7.49  -9.37  1.00  1.00
              1.00  -6.54   1.00   1.00  1.00  1.00  1.00   1.00  1.00  1.00  1.00  -7.49  1.00  -6.54  1.00 20.26 -14.03   1.00  1.00  1.00
              1.00   1.00   1.00  20.26  1.00 -6.54  1.00   1.00 33.60  1.00 20.26   1.00 33.60   1.00 20.26  1.00  33.60  24.68  1.00 -6.54
              1.00   1.00   1.00  20.26 -6.54 20.26 20.26   1.00  1.00  1.00  1.00   1.00  1.00  -6.54 20.26 44.94   1.00   1.00 -6.54 -6.54
              1.00   1.00   1.00  20.26 44.94 20.26 33.60   1.00 -6.54 20.26  1.00   1.00  1.00   1.00 20.26 20.26   1.00 -14.03  1.00  1.00
             -7.49   1.00  -7.49   1.00  1.00  1.00 -6.54  13.34  1.00 -7.49  1.00  -7.49  1.00   1.00  1.00  1.00  -7.49  13.34 -7.49  1.00
              1.00   1.00  24.68   1.00  1.00  1.00  1.00  -9.37  1.00 44.94  1.00  24.68  1.00  -9.37 -1.88  1.00  -6.54  -1.88 44.94  1.00
              1.00   1.00   1.00   1.00  1.00  1.00 44.94   1.00 13.34  1.00 20.26  -7.49  1.00   1.00 -1.88  1.00   1.00   1.00  1.00 -7.49
              1.00  20.26   1.00   1.00  1.00 33.60  1.00   1.00  1.00  1.00  1.00  -7.49  1.00   1.00 20.26  1.00   1.00  24.68  1.00  1.00
              1.00  33.60   1.00   1.00  1.00 24.68  1.00  -7.49  1.00 -7.49 -7.49   1.00 33.60   1.00 -6.54  1.00   1.00   1.00  1.00 -7.49
             13.34  -6.54   1.00   1.00  1.00 -6.54  1.00   1.00 58.28  1.00  1.00   1.00 -1.88   1.00 44.94 44.94  -1.88   1.00 24.68  1.00
              1.00   1.00   1.00  13.34  1.00  1.00  1.00   1.00  1.00  1.00  1.00 -14.03  1.00   1.00 20.26  1.00   1.00   1.00 33.60  1.00
             20.26  -6.54   1.00  -6.54 -6.54 20.26 18.38   1.00  1.00  1.00  1.00   1.00 -6.54  20.26 20.26 20.26   1.00  -1.88  1.00 20.26
              1.00  20.26   1.00   1.00 33.60 20.26 20.26   1.00  1.00  1.00  1.00   1.00  1.00   1.00 44.94 20.26   1.00   1.00  1.00  1.00
              1.00   1.00 -14.03   1.00  1.00 -6.54 20.26  -7.49  1.00  1.00  1.00   1.00  1.00  13.34  1.00  1.00   1.00 -14.03  1.00  1.00
            -14.03   1.00  13.34   1.00  1.00  1.00  1.00  -9.37 24.68  1.00 13.34   1.00 24.68   1.00  1.00  1.00 -14.03   1.00  1.00 -7.49
             24.68 -15.91   1.00  24.68  1.00  1.00 -6.54  -7.49 13.34  1.00  1.00   1.00 44.94   1.00 13.34  1.00  -7.49  -9.37 13.34  1.00
              1.00   1.00   1.00 -14.03  1.00  1.00  1.00  -7.49  1.00  1.00  1.00  -1.88  1.00   1.00 20.26  1.00  -7.49   1.00 -6.54  1.00];
    end

    methods
        function this = ProteinMonomer(knowledgeBase, wid, wholeCellModelID, name,...
                topology, activeSite, metalBindingSite, dnaFootprint, dnaFootprintBindingStrandedness, dnaFootprintRegionStrandedness, molecularInteraction, ...
                chemicalRegulation, subsystem, generalClassification, proteaseClassification, ...
                transporterClassification, ...
                nTerminalMethionineCleavage, signalSequenceType, signalSequenceLocation, signalSequenceLength, ...
                activationRule, ...
                comments, crossReferences)

            if nargin == 0; return; end;

            this = edu.stanford.covert.cell.kb.ProteinMonomer.empty(size(wid, 1), 0);
            this(size(wid, 1), 1) = edu.stanford.covert.cell.kb.ProteinMonomer;
            for i = 1:size(wid, 1)
                this(i, 1).idx = i;
                this(i, 1).knowledgeBase = knowledgeBase;
                this(i, 1).wid = wid(i);
                this(i, 1).wholeCellModelID = wholeCellModelID{i};
                this(i, 1).name = name{i};
                if exist('comments', 'var') && ~isempty(comments); this(i, 1).comments = comments{i}; end;
                if exist('crossReferences', 'var')
                    if size(crossReferences, 1) > 1
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

                this(i, 1).topology = topology{i};
                this(i, 1).activeSite = activeSite{i};
                this(i, 1).metalBindingSite = metalBindingSite{i};
                this(i, 1).dnaFootprint = dnaFootprint(i);
                this(i, 1).dnaFootprintBindingStrandedness = dnaFootprintBindingStrandedness{i};
                this(i, 1).dnaFootprintRegionStrandedness = dnaFootprintRegionStrandedness{i};
                this(i, 1).molecularInteraction = molecularInteraction{i};
                this(i, 1).chemicalRegulation = chemicalRegulation{i};
                this(i, 1).subsystem = subsystem{i};
                this(i, 1).generalClassification = generalClassification{i};
                this(i, 1).proteaseClassification = proteaseClassification;
                this(i, 1).transporterClassification = transporterClassification{i};
                this(i, 1).nTerminalMethionineCleavage = nTerminalMethionineCleavage(i);
                this(i, 1).signalSequenceType = signalSequenceType{i};
                this(i, 1).signalSequenceLocation = signalSequenceLocation{i};
                this(i, 1).signalSequenceLength = signalSequenceLength(i);
                this(i, 1).activationRule = activationRule{i};
            end
        end

        function serializeLinks(this)
            for i = 1:numel(this)
                this(i).gene                             = this.serializeLinksHelper(this(i).gene);
                this(i).geneCompartments                 = this.serializeLinksHelper(this(i).geneCompartments);

                this(i).prostheticGroups                 = this.serializeLinksHelper(this(i).prostheticGroups);
                this(i).prostheticGroupCompartments      = this.serializeLinksHelper(this(i).prostheticGroupCompartments);
                this(i).prostheticGroupCoefficients      = this.serializeLinksHelper(this(i).prostheticGroupCoefficients);

                this(i).chaperoneSubstrates              = this.serializeLinksHelper(this(i).chaperoneSubstrates);
                this(i).chaperoneCompartments            = this.serializeLinksHelper(this(i).chaperoneCompartments);
                this(i).chaperoneCoefficients            = this.serializeLinksHelper(this(i).chaperoneCoefficients);

                this(i).proteinComplexs                  = this.serializeLinksHelper(this(i).proteinComplexs);

                this(i).stableModificationReactions      = this.serializeLinksHelper(this(i).stableModificationReactions);

                this(i).regulatedTranscriptionUnits      = this.serializeLinksHelper(this(i).regulatedTranscriptionUnits);

                this(i).compartment                      = this.serializeLinksHelper(this(i).compartment);

                this(i).stimuliRegulators                = this.serializeLinksHelper(this(i).stimuliRegulators);
                this(i).metaboliteRegulators             = this.serializeLinksHelper(this(i).metaboliteRegulators);
                this(i).proteinMonomerRegulators         = this.serializeLinksHelper(this(i).proteinMonomerRegulators);
                this(i).proteinComplexRegulators         = this.serializeLinksHelper(this(i).proteinComplexRegulators);
                
                this(i).regulatedProteinMonomers         = this.serializeLinksHelper(this(i).regulatedProteinMonomers);
                this(i).regulatedProteinComplexs         = this.serializeLinksHelper(this(i).regulatedProteinComplexs);
                
                this(i).parameters                       = this.serializeLinksHelper(this(i).parameters);
                                
                serializeLinks@edu.stanford.covert.cell.kb.Polymer(this(i));
                serializeLinks@edu.stanford.covert.cell.kb.Protein(this(i));
            end
        end
        
        function deserializeLinks(this, kb)
            for i = 1:numel(this)
                this(i).gene                             = this.deserializeLinksHelper(this(i).gene, kb.genes);
                this(i).geneCompartments                 = this.deserializeLinksHelper(this(i).geneCompartments, kb.compartments);
                
                this(i).prostheticGroups                 = this.deserializeLinksHelper(this(i).prostheticGroups, kb.metabolites);
                this(i).prostheticGroupCompartments      = this.deserializeLinksHelper(this(i).prostheticGroupCompartments, kb.compartments);
                this(i).prostheticGroupCoefficients      = this.deserializeLinksHelper(this(i).prostheticGroupCoefficients);
                
                this(i).chaperoneSubstrates              = this.deserializeLinksHelper(this(i).chaperoneSubstrates, kb.proteinMonomers);
                this(i).chaperoneCompartments            = this.deserializeLinksHelper(this(i).chaperoneCompartments, kb.compartments);
                this(i).chaperoneCoefficients            = this.deserializeLinksHelper(this(i).chaperoneCoefficients);
                
                this(i).proteinComplexs                  = this.deserializeLinksHelper(this(i).proteinComplexs, kb.proteinComplexs);
                
                this(i).stableModificationReactions      = this.deserializeLinksHelper(this(i).stableModificationReactions, kb.reactions);
                
                this(i).regulatedTranscriptionUnits      = this.deserializeLinksHelper(this(i).regulatedTranscriptionUnits, kb.transcriptionUnits);
                
                this(i).compartment                      = this.deserializeLinksHelper(this(i).compartment, kb.compartments);
                
                this(i).stimuliRegulators                = this.deserializeLinksHelper(this(i).stimuliRegulators, kb.stimulis);
                this(i).metaboliteRegulators             = this.deserializeLinksHelper(this(i).metaboliteRegulators, kb.metabolites);
                this(i).proteinMonomerRegulators         = this.deserializeLinksHelper(this(i).proteinMonomerRegulators, kb.proteinMonomers);
                this(i).proteinComplexRegulators         = this.deserializeLinksHelper(this(i).proteinComplexRegulators, kb.proteinComplexs);
                
                this(i).regulatedProteinMonomers         = this.deserializeLinksHelper(this(i).regulatedProteinMonomers, kb.proteinMonomers);
                this(i).regulatedProteinComplexs         = this.deserializeLinksHelper(this(i).regulatedProteinComplexs, kb.proteinComplexs);
                
                this(i).parameters                       = this.deserializeLinksHelper(this(i).parameters, kb.parameters);
                
                deserializeLinks@edu.stanford.covert.cell.kb.Polymer(this(i), kb);
                deserializeLinks@edu.stanford.covert.cell.kb.Protein(this(i), kb);
            end
        end

        function deleteLinks(this)
            for i = 1:numel(this)
                this(i).gene                             = [];
                this(i).geneCompartments                 = [];

                this(i).prostheticGroups                 = [];
                this(i).prostheticGroupCompartments      = [];
                this(i).prostheticGroupCoefficients      = [];

                this(i).chaperoneSubstrates              = [];
                this(i).chaperoneCompartments            = [];
                this(i).chaperoneCoefficients            = [];

                this(i).proteinComplexs                  = [];

                this(i).stableModificationReactions      = [];

                this(i).regulatedTranscriptionUnits      = [];

                this(i).compartment                      = [];

                this(i).stimuliRegulators                = [];
                this(i).metaboliteRegulators             = [];
                this(i).proteinMonomerRegulators         = [];
                this(i).proteinComplexRegulators         = [];

                this(i).regulatedProteinMonomers         = [];
                this(i).regulatedProteinComplexs         = [];

                this(i).parameters                       = [];

                deleteLinks@edu.stanford.covert.cell.kb.Polymer(this(i));
                deleteLinks@edu.stanford.covert.cell.kb.Protein(this(i));
            end
        end

        function value = get.empiricalFormula(~)
            throw(MException('ProteinMonomer:error property is not defined'));
        end

        function value = get.smiles(~)
            throw(MException('ProteinMonomer:error', 'property is not defined'));
        end

        function value = get.charge(~)
            throw(MException('ProteinMonomer:error', 'property is not defined'));
        end

        function value = get.pKa(~)
            throw(MException('ProteinMonomer:error', 'property is not defined'));
        end

        function value = get.sequence(this)
            %retrieve
            if ~isempty(this.sequence)
                value = this.sequence;
                return;
            end
            
            %calculate
            value = nt2aa(this.gene.dnaSequence,'GeneticCode',this.knowledgeBase.translationTable);
            stopCodonIdx = strfind(value, '*');
            if ~isempty(stopCodonIdx)
                value = value(1:stopCodonIdx-1);
            end
            
            %store
            this.sequence = value;
        end

        function value = get.signalSequence(this)
            %retrieve
            if ~isempty(this.signalSequence)
                value = this.signalSequence;
                return;
            end
            
            %calculate
            value = this.sequence(1:this.signalSequenceLength);
            
            %store
            this.signalSequence = value;
        end

        function value = get.processedISequence(this)
            %retrieve
            if ~isempty(this.processedISequence)
                value = this.processedISequence;
                return;
            end
            
            %calculate
            value = this.sequence(this.nTerminalMethionineCleavage+1:end);
            
            %store
            this.processedISequence = value;
        end

        function value = get.processedIISequence(this)
            %retrieve
            if ~isempty(this.processedIISequence)
                value = this.processedIISequence;
                return;
            end
            
            %calculate
            if this.signalSequenceLength>0
                value = this.sequence(this.signalSequenceLength+1:end);
            else
                value = this.sequence(this.nTerminalMethionineCleavage+1:end);
            end
            
            %store
            this.processedIISequence = value;
        end

        function value = get.foldedSequence(this)
            %retrieve
            if ~isempty(this.foldedSequence)
                value = this.foldedSequence;
                return;
            end
            
            %calculate
            value = this.processedIISequence;
            
            %store
            this.foldedSequence = value;
        end

        function value = get.matureSequence(this)
            %retrieve
            if ~isempty(this.matureSequence)
                value = this.matureSequence;
                return;
            end
            
            %calculate
            value = this.foldedSequence;

            if isempty(this.stableModificationReactions); return; end;

            ligations = zeros(1, this.knowledgeBase.numMetabolites);
            for i = 1:length(this.stableModificationReactions)
                reaction = this.stableModificationReactions(i);
                if ~strcmp(reaction.type,'ligation') || isempty(reaction.metabolites); continue; end;
                ligations([reaction.metabolites.idx]) = ...
                    ligations([reaction.metabolites.idx])+...
                    reaction.metaboliteCoefficients';
            end

            for i = 1:length(this.bases)
                value = [value repmat(this.bases(i), 1, -ligations(this.knowledgeBase.aminoAcidIndexs(i)))]; %#ok<AGROW>
            end
            
            %store
            this.matureSequence = value;
        end

        function value = get.processedISequenceLength(this)
            %retrieve
            if ~isempty(this.processedISequenceLength)
                value = this.processedISequenceLength;
                return;
            end
            
            %calculate
            value = length(this.processedISequence);
            
            %store
            this.processedISequenceLength = value;
        end

        function value = get.processedIISequenceLength(this)
            %retrieve
            if ~isempty(this.processedIISequenceLength)
                value = this.processedIISequenceLength;
                return;
            end
            
            %calculate
            value = length(this.processedIISequence);
            
            %store
            this.processedIISequenceLength = value;
        end

        function value = get.foldedSequenceLength(this)
            %retrieve
            if ~isempty(this.foldedSequenceLength)
                value = this.foldedSequenceLength;
                return;
            end
            
            %calculate
            value = length(this.foldedSequence);
            
            %store
            this.foldedSequenceLength = value;
        end

        function value = get.matureSequenceLength(this)
            %retrieve
            if ~isempty(this.matureSequenceLength)
                value = this.matureSequenceLength;
                return;
            end
            
            %calculate
            value = length(this.matureSequence);
            
            %store
            this.matureSequenceLength = value;
        end

        function value = get.baseCount(this)
            %retrieve
            if ~isempty(this.baseCount)
                value = this.baseCount;
                return;
            end
            
            %calculate
            value = this.computeBaseCount(this.sequence, this.knowledgeBase.numMetabolites, this.knowledgeBase.aminoAcidIndexs, true);
            
            %store
            this.baseCount = value;
        end

        function value = get.signalSequenceBaseCount(this)
            %retrieve
            if ~isempty(this.signalSequenceBaseCount)
                value = this.signalSequenceBaseCount;
                return;
            end
            
            %calculate
            value = this.computeBaseCount(this.signalSequence, this.knowledgeBase.numMetabolites, this.knowledgeBase.aminoAcidIndexs, false);
            
            %store
            this.signalSequenceBaseCount = value;
        end

        function value = get.processedISequenceBaseCount(this)
            %retrieve
            if ~isempty(this.processedISequenceBaseCount)
                value = this.processedISequenceBaseCount;
                return;
            end
            
            %calculate
            value = this.computeBaseCount(this.processedISequence, this.knowledgeBase.numMetabolites, this.knowledgeBase.aminoAcidIndexs, false);
            
            %store
            this.processedISequenceBaseCount = value;
        end

        function value = get.processedIISequenceBaseCount(this)
            %retrieve
            if ~isempty(this.processedIISequenceBaseCount)
                value = this.processedIISequenceBaseCount;
                return;
            end
            
            %calculate
            value = this.computeBaseCount(this.processedIISequence, this.knowledgeBase.numMetabolites, this.knowledgeBase.aminoAcidIndexs, false);

            if strcmp(this.signalSequenceType,'lipoprotein')
                value(this.knowledgeBase.cysteineIndexs) = ...
                    value(this.knowledgeBase.cysteineIndexs)-1;
                value(this.knowledgeBase.diacylglycerolCysteineIndexs) = ...
                    value(this.knowledgeBase.diacylglycerolCysteineIndexs)+1;
            end
            
            %store
            this.processedIISequenceBaseCount = value;
        end

        function value = get.foldedSequenceBaseCount(this)
            %retrieve
            if ~isempty(this.foldedSequenceBaseCount)
                value = this.foldedSequenceBaseCount;
                return;
            end
            
            %calculate
            value = this.processedIISequenceBaseCount;

            if isempty(this.prostheticGroups); return; end;

            value([this.prostheticGroups.idx]) = ...
                value([this.prostheticGroups.idx]) + ...
                max(1, this.prostheticGroupCoefficients)';
            
            %store
            this.foldedSequenceBaseCount = value;
        end

        function value = get.matureSequenceBaseCount(this)
            %retrieve
            if ~isempty(this.matureSequenceBaseCount)
                value = this.matureSequenceBaseCount;
                return;
            end
            
            %calculate
            value = this.foldedSequenceBaseCount;
            
            if ~isempty(this.stableModificationReactions)                
                adductions = zeros(size(value));
                ligations  = zeros(size(value));
                for i = 1:length(this.stableModificationReactions)
                    reaction = this.stableModificationReactions(i);
                    if isempty(reaction.metabolites); continue; end;
                    
                    switch reaction.type
                        case 'adduction'
                            adductions([reaction.metabolites.idx]) = ...
                                adductions([reaction.metabolites.idx]) + ...
                                reaction.metaboliteCoefficients';
                        case 'ligation'
                            ligations([reaction.metabolites.idx]) = ...
                                ligations([reaction.metabolites.idx]) + ...
                                reaction.metaboliteCoefficients';
                    end
                end
                
                value([this.knowledgeBase.aminoAcidIndexs; this.knowledgeBase.modifiedAminoAcidIndexs]) = ...
                    value([this.knowledgeBase.aminoAcidIndexs; this.knowledgeBase.modifiedAminoAcidIndexs]) + ...
                    adductions([this.knowledgeBase.aminoAcidIndexs; this.knowledgeBase.modifiedAminoAcidIndexs]);
                
                value(this.knowledgeBase.aminoAcidIndexs) = ...
                    value(this.knowledgeBase.aminoAcidIndexs) + ...
                    -ligations(this.knowledgeBase.aminoAcidIndexs);
            end
            
            %store
            this.matureSequenceBaseCount = value;
        end

        function value = get.cumulativeBaseCount(this)
            %retrieve
            if ~isempty(this.cumulativeBaseCount)
                value = this.cumulativeBaseCount;
                return;
            end
            
            %calculate
            aminoAcidIndexs = this.knowledgeBase.aminoAcidIndexs;
            value = zeros(length(this.sequence), this.knowledgeBase.numMetabolites);
            for i = 1:length(this.bases)
                value(:,aminoAcidIndexs(i)) = this.sequence == this.bases(i);
            end
            value = cumsum(value,2);
            
            %store
            this.cumulativeBaseCount = value;
        end

        function value = get.tRNACount(~)
            throw(MException('ProteinMonomer:error','property is not defined'));
        end

        function value = get.molecularWeight(this)
            %retrieve
            if ~isempty(this.molecularWeight)
                value = this.molecularWeight;
                return;
            end
            
            %calculate
            value = this.computeMolecularWeight(this.baseCount, this.sequenceLength, this.knowledgeBase.metaboliteMolecularWeights);
            
            %store
            this.molecularWeight = value;
        end

        function value = get.signalSequenceMolecularWeight(this)
            %retrieve
            if ~isempty(this.signalSequenceMolecularWeight)
                value = this.signalSequenceMolecularWeight;
                return;
            end
            
            %calculate
            value = this.computeMolecularWeight(this.signalSequenceBaseCount, this.signalSequenceLength, this.knowledgeBase.metaboliteMolecularWeights);
            
            %store
            this.signalSequenceMolecularWeight = value;
        end

        function value = get.processedISequenceMolecularWeight(this)
            %retrieve
            if ~isempty(this.processedISequenceMolecularWeight)
                value = this.processedISequenceMolecularWeight;
                return;
            end
            
            %calculate
            value = this.computeMolecularWeight(this.processedISequenceBaseCount, this.processedISequenceLength, this.knowledgeBase.metaboliteMolecularWeights);
            
            %store
            this.processedISequenceMolecularWeight = value;
        end

        function value = get.processedIISequenceMolecularWeight(this)
            %retrieve
            if ~isempty(this.processedIISequenceMolecularWeight)
                value = this.processedIISequenceMolecularWeight;
                return;
            end
            
            %calculate
            value = this.computeMolecularWeight(this.processedIISequenceBaseCount, this.processedIISequenceLength, this.knowledgeBase.metaboliteMolecularWeights);
            
            %store
            this.processedIISequenceMolecularWeight = value;
        end

        function value = get.foldedSequenceMolecularWeight(this)
            %retrieve
            if ~isempty(this.foldedSequenceMolecularWeight)
                value = this.foldedSequenceMolecularWeight;
                return;
            end
            
            %calculate
            value = this.processedIISequenceMolecularWeight;

            if ~isempty(this.prostheticGroups)
                value = value + [this.prostheticGroups.molecularWeight] * ...
                    max(1, this.prostheticGroupCoefficients);
            end
            
            %store
            this.foldedSequenceMolecularWeight = value;
        end

        function value = get.matureSequenceMolecularWeight(this)
            % import classes
            import edu.stanford.covert.util.ConstantUtil;
            
            %retrieve
            if ~isempty(this.matureSequenceMolecularWeight)
                value = this.matureSequenceMolecularWeight;
                return;
            end
            
            %calculate
            value = this.foldedSequenceMolecularWeight;

            if ~isempty(this.stableModificationReactions)
                adductions = zeros(this.knowledgeBase.numMetabolites, 1);
                ligations  = zeros(this.knowledgeBase.numMetabolites, 1);
                for i = 1:length(this.stableModificationReactions)
                    reaction = this.stableModificationReactions(i);
                    if isempty(reaction.metabolites); continue; end;

                    switch reaction.type
                        case 'adduction'
                            adductions([reaction.metabolites.idx]) = ...
                                adductions([reaction.metabolites.idx]) + ...
                                reaction.metaboliteCoefficients;
                        case 'ligation'
                            ligations([reaction.metabolites.idx]) = ...
                                ligations([reaction.metabolites.idx]) + ...
                                reaction.metaboliteCoefficients;
                    end
                end

                molecularWeightH2O = 2 * ConstantUtil.elements.H + ConstantUtil.elements.O;

                value = value ...
                    + this.knowledgeBase.metaboliteMolecularWeights([this.knowledgeBase.aminoAcidIndexs; this.knowledgeBase.modifiedAminoAcidIndexs])'...
                    * adductions([this.knowledgeBase.aminoAcidIndexs;this.knowledgeBase.modifiedAminoAcidIndexs]) ...
                    - (this.knowledgeBase.metaboliteMolecularWeights(this.knowledgeBase.aminoAcidIndexs) - molecularWeightH2O)'...
                    * ligations(this.knowledgeBase.aminoAcidIndexs);
            end
            
            %store
            this.matureSequenceMolecularWeight = value;
        end

        function value = get.decayReaction(this)
            %retrieve
            if ~isempty(this.decayReaction)
                value = this.decayReaction;
                return;
            end
            
            %calculate
            value = this.computeDecayReaction(this.baseCount, this.sequenceLength, this.knowledgeBase.waterIndexs);
            
            %store
            this.decayReaction = value;
        end

        function value = get.signalSequenceDecayReaction(this)
            %retrieve
            if ~isempty(this.signalSequenceDecayReaction)
                value = this.signalSequenceDecayReaction;
                return;
            end
            
            %calculate
            value = this.computeDecayReaction(this.signalSequenceBaseCount, this.signalSequenceLength, this.knowledgeBase.waterIndexs);
            
            %store
            this.signalSequenceDecayReaction = value;
        end

        function value = get.processedISequenceDecayReaction(this)
            %retrieve
            if ~isempty(this.processedISequenceDecayReaction)
                value = this.processedISequenceDecayReaction;
                return;
            end
            
            %calculate
            value = this.computeDecayReaction(this.processedISequenceBaseCount, this.processedISequenceLength, this.knowledgeBase.waterIndexs);
            
            %store
            this.processedISequenceDecayReaction = value;
        end

        function value = get.processedIISequenceDecayReaction(this)
            %retrieve
            if ~isempty(this.processedIISequenceDecayReaction)
                value = this.processedIISequenceDecayReaction;
                return;
            end
            
            %calculate
            value = this.computeDecayReaction(this.processedIISequenceBaseCount, this.processedIISequenceLength, this.knowledgeBase.waterIndexs);
            
            %store
            this.processedIISequenceDecayReaction = value;
        end

        function value = get.foldedSequenceDecayReaction(this)
            %retrieve
            if ~isempty(this.foldedSequenceDecayReaction)
                value = this.foldedSequenceDecayReaction;
                return;
            end
            
            %calculate
            value = this.computeDecayReaction(this.foldedSequenceBaseCount, this.foldedSequenceLength, this.knowledgeBase.waterIndexs);
            
            %store
            this.foldedSequenceDecayReaction = value;
        end

        function value = get.matureSequenceDecayReaction(this)
            %retrieve
            if ~isempty(this.matureSequenceDecayReaction)
                value = this.matureSequenceDecayReaction;
                return;
            end
            
            %calculate
            value = this.computeDecayReaction(this.matureSequenceBaseCount, this.matureSequenceLength, this.knowledgeBase.waterIndexs);
            
            %store
            this.matureSequenceDecayReaction = value;
        end

        function value = get.halfLife(this)
            %retrieve
            if ~isempty(this.halfLife)
                value = this.halfLife;
                return;
            end
            
            %calculate
            value = this.computeHalfLife(this.nTerminalAA);
            
            %store
            this.halfLife = value;
        end
        
        function value = get.signalSequenceHalfLife(this)
            %retrieve
            if ~isempty(this.signalSequenceHalfLife)
                value = this.signalSequenceHalfLife;
                return;
            end
            
            %calculate
            value = 0;
            
            %store
            this.signalSequenceHalfLife = value;
        end
        
        function value = get.processedISequenceHalfLife(this)
            %retrieve
            if ~isempty(this.processedISequenceHalfLife)
                value = this.processedISequenceHalfLife;
                return;
            end
            
            %calculate
            value = this.computeHalfLife(this.processedISequence(1));
            
            %store
            this.processedISequenceHalfLife = value;
        end
        
        function value = get.processedIISequenceHalfLife(this)
            %retrieve
            if ~isempty(this.processedIISequenceHalfLife)
                value = this.processedIISequenceHalfLife;
                return;
            end
            
            %calculate
            value = this.computeHalfLife(this.processedIISequence(1));
            
            %store
            this.processedIISequenceHalfLife = value;
        end
        
        function value = get.foldedSequenceHalfLife(this)
            %retrieve
            if ~isempty(this.foldedSequenceHalfLife)
                value = this.foldedSequenceHalfLife;
                return;
            end
            
            %calculate
            value = this.computeHalfLife(this.foldedSequence(1));
            
            %store
            this.foldedSequenceHalfLife = value;
        end
        
        function value = get.matureSequenceHalfLife(this)
            %retrieve
            if ~isempty(this.matureSequenceHalfLife)
                value = this.matureSequenceHalfLife;
                return;
            end
            
            %calculate
            value = this.computeHalfLife(this.matureSequence(1));
            
            %store
            this.matureSequenceHalfLife = value;
        end

        %http://www3.interscience.wiley.com/journal/121602274/abstract?CRETRY=1&SRETRY=0
        function value = get.density(this)
            %retrieve
            if ~isempty(this.density)
                value = this.density;
                return;
            end
            
            %calculate
            value = 1.35;
            
            %store
            this.density = value;
        end

        function value = get.volume(this)
            %retrieve
            if ~isempty(this.volume)
                value = this.volume;
                return;
            end
            
            %calculate
            value = this.molecularWeight/this.density;
            
            %store
            this.volume = value;
        end

        %http://isoelectric.ovh.org/files/practise-isoelectric-point.html
        function value = get.pI(this)
            %retrieve
            if ~isempty(this.pI)
                value = this.pI;
                return;
            end
            
            %calculate
            numAsp = sum('D' == this.sequence);
            numGlu = sum('E' == this.sequence);
            numCys = sum('C' == this.sequence);
            numTyr = sum('Y' == this.sequence);
            numHis = sum('H' == this.sequence);
            numLys = sum('K' == this.sequence);
            numArg = sum('R' == this.sequence);

            pH = 6.5;             %starting point pI = 6.5 - theoretically it should be 7, but average protein pI is 6.5 so we increase the probability
            pHprev = 0.0;         %of finding the solution
            pHnext = 14.0;        %0-14 is possible pH range
            E = 0.01;             %epsilon means precision [pI = pH \pm E]

            %the infinite loop
            while true

                % we are using pK values form Wikipedia as they give quite good approximation
                % if you want you can change it
                QN1 = -1/(1+10^(3.65-pH));       %C-terminal charge
                QN2 = -numAsp/(1+10^(3.9-pH));   %D charge
                QN3 = -numGlu/(1+10^(4.07-pH));  %E charge
                QN4 = -numCys/(1+10^(8.18-pH));  %C charge
                QN5 = -numTyr/(1+10^(10.46-pH)); %Y charge
                QP1 = numHis/(1+10^(pH-6.04));   %H charge
                QP2 = 1/(1+10^(pH-8.2));         %NH2 charge
                QP3 = numLys/(1+10^(pH-10.54));  %K charge
                QP4 = numArg/(1+10^(pH-12.48));  %R charge

                NQ = QN1+QN2+QN3+QN4+QN5+QP1+QP2+QP3+QP4; %net charge in given pH


                if pH >= 14.0
                    throw(MException('ProteinMonomer:error','pH higher than 14'));
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%   BISECTION   %%%%%%%%%%%%%%%%%%%%%%%%

                %we are out of range, thus the new pH value must be smaller
                if NQ<0
                    temp = pH;
                    pH = pH-((pH-pHprev)/2);
                    pHnext = temp;

                    %we used to small pH value, so we have to increase it
                else
                    temp = pH;
                    pH = pH + ((pHnext-pH)/2);
                    pHprev = temp;
                end

                %terminal condition, finding isoelectric point with given precision
                if (pH-pHprev<E) && (pHnext-pH<E)
                    break;
                end
            end

            value = pH;
            
            %store
            this.pI = value;
        end

        %http://ca.expasy.org/tools/protparam-doc.html
        function value = get.extinctionCoefficient(this)
            %retrieve
            if ~isempty(this.extinctionCoefficient)
                value = this.extinctionCoefficient;
                return;
            end
            
            %calculate
            value = ...
                sum('W' == this.sequence) * 5500 + ...
                sum('Y' == this.sequence) * 1490 + ...
                sum('C' == this.sequence) * 125;
            
            %store
            this.extinctionCoefficient = value;
        end

        %http://ca.expasy.org/tools/protparam-doc.html
        function value = get.absorbanceFactor(this)
            %retrieve
            if ~isempty(this.absorbanceFactor)
                value = this.absorbanceFactor;
                return;
            end
            
            %calculate
            value = this.extinctionCoefficient/this.molecularWeight;
            
            %store
            this.absorbanceFactor = value;
        end

        %http://en.wikipedia.org/wiki/Amino_acid
        function value = get.negAA(this)
            %retrieve
            if ~isempty(this.negAA)
                value = this.negAA;
                return;
            end
            
            %calculate
            value = sum('D' == this.sequence) + sum('E' == this.sequence);
            
            %store
            this.negAA = value;
        end

        %http://en.wikipedia.org/wiki/Amino_acid
        function value = get.posAA(this)
            %retrieve
            if ~isempty(this.posAA)
                value = this.posAA;
                return;
            end
            
            %calculate
            value = ...
                sum('R' == this.sequence) + ...
                sum('H' == this.sequence) + ...
                sum('K' == this.sequence);
            
            %store
            this.posAA = value;
        end

        function value = get.nTerminalAA(this)
            %retrieve
            if ~isempty(this.nTerminalAA)
                value = this.nTerminalAA;
                return;
            end
            
            %calculate
            value = this.sequence(1);
            
            %store
            this.nTerminalAA = value;
        end

        function value = get.disulfideBonds(~)
            throw(MException('ProteinMonomer:error','property is not defined'));
        end

        %Instability index (II)
        %http://ca.expasy.org/tools/protparam-doc.html
        function value = get.instabilityIndex(this)
            %retrieve
            if ~isempty(this.instabilityIndex)
                value = this.instabilityIndex;
                return;
            end
            
            %calculate
            value = 0;
            idx2 = aa2int(this.sequence(1));
            for i = 1:length(this.sequence)-1
                idx1 = idx2;
                idx2 = aa2int(this.sequence(i+1));
                if idx1>size(this.DIWV,1) || idx2>size(this.DIWV,2); continue; end;
                value = value+this.DIWV(idx1,idx2);
            end
            value = 10/length(this.sequence)*value;
            
            %store
            this.instabilityIndex = value;
        end

        %http://ca.expasy.org/tools/protparam-doc.html
        function value = get.stable(this)
            %retrieve
            if ~isempty(this.stable)
                value = this.stable;
                return;
            end
            
            %calculate
            value = (this.instabilityIndex<40);
            
            %store
            this.sequence = value;
        end

        %http://ca.expasy.org/tools/protparam-doc.html
        function value = get.aliphaticIndex(this)
            %retrieve
            if ~isempty(this.sequence)
                value = this.sequence;
                return;
            end
            
            %calculate
            value = 100*(...
                1.0*sum('A' == this.sequence)+...
                2.9*sum('V' == this.sequence)+...
                3.9*sum('I' == this.sequence)+...
                3.9*sum('L' == this.sequence)) / ...
                this.sequenceLength;
            
            %store
            this.stable = value;
        end

        %GRAVY (Grand Average of Hydropathy)
        %http://ca.expasy.org/tools/protparam-doc.html
        function value = get.gravy(this)
            %retrieve
            if ~isempty(this.gravy)
                value = this.gravy;
                return;
            end
            
            %calculate
            value = ...
                (1.8*sum('A' == this.sequence)+...
                -4.5*sum('R' == this.sequence)+...
                -3.5*sum('N' == this.sequence)+...
                -3.5*sum('D' == this.sequence)+...
                 2.5*sum('C' == this.sequence)+...
                -3.5*sum('Q' == this.sequence)+...
                -3.5*sum('E' == this.sequence)+...
                -0.4*sum('G' == this.sequence)+...
                -3.2*sum('H' == this.sequence)+...
                 4.5*sum('I' == this.sequence)+...
                 3.8*sum('L' == this.sequence)+...
                -3.9*sum('K' == this.sequence)+...
                 1.9*sum('M' == this.sequence)+...
                 2.8*sum('F' == this.sequence)+...
                -1.6*sum('P' == this.sequence)+...
                -0.8*sum('S' == this.sequence)+...
                -0.7*sum('T' == this.sequence)+...
                -0.9*sum('W' == this.sequence)+...
                -1.3*sum('Y' == this.sequence)+...
                 4.2*sum('V' == this.sequence)) /...
                 this.sequenceLength;
             
             %store
            this.gravy = value;
        end               
    end   

    methods (Static = true)
        function value = computeBaseCount(sequence, numMetabolites, aminoAcidIndexs, nTerminalFormylMethionine)
            import edu.stanford.covert.cell.kb.ProteinMonomer;

            value = zeros(1,numMetabolites);
            for i = 1:length(ProteinMonomer.bases)
                value(aminoAcidIndexs(i)) = sum(ProteinMonomer.bases(i) == sequence);
            end

            if nTerminalFormylMethionine
                value(aminoAcidIndexs(13)) = value(aminoAcidIndexs(13)) - 1;
                value(aminoAcidIndexs(21)) = value(aminoAcidIndexs(21)) + 1;
            end
        end

        function value = computeDecayReaction(baseCount, sequenceLength, waterIndexs)
            value = baseCount;
            value(waterIndexs) = ...
                + value(waterIndexs) ...
                - max(0, sequenceLength - 1);
        end

        function value = computeMolecularWeight(baseCount, sequenceLength, metaboliteMolecularWeights)
            % import classes
            import edu.stanford.covert.util.ConstantUtil;

            molecularWeightH2O = 2 * ConstantUtil.elements.H + ConstantUtil.elements.O;

            if ~any(baseCount)
                value = 0;
                return;
            end

            value = baseCount * metaboliteMolecularWeights ...
                - molecularWeightH2O * max(0, sequenceLength - 1);
        end

        %N-end rule
        %http://ca.expasy.org/tools/protparam-doc.html
        function value = computeHalfLife(nTerminalAA) %#ok<INUSD>
%             switch nTerminalAA
%                 case {'R','L','K','F','W','Y'}, value = 2*60;
%                 case {'P'}, value = Inf;
%                 otherwise, value = 10*60*60;
%             end
            value = 20 * 60 * 60;
        end
    end
end