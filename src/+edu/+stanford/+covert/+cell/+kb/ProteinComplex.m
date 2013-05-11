% Defines a protein complex
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/7/2009
classdef ProteinComplex < edu.stanford.covert.cell.kb.PhysicalObject & edu.stanford.covert.cell.kb.Protein
    properties
        prostheticGroups            = edu.stanford.covert.cell.kb.Metabolite.empty(0, 0);
        prostheticGroupCompartments = edu.stanford.covert.cell.kb.Compartment.empty(0, 0);
        prostheticGroupCoefficients = [];

        chaperoneSubstrates         = edu.stanford.covert.cell.kb.ProteinMonomer.empty(0, 0);
        chaperoneCompartments       = edu.stanford.covert.cell.kb.Compartment.empty(0, 0);
        chaperoneCoefficients       = [];
        
        metabolites                 = edu.stanford.covert.cell.kb.Metabolite.empty(0, 0);
        metaboliteCompartments      = edu.stanford.covert.cell.kb.Compartment.empty(0, 0);
        metaboliteCoefficients      = [];

        rnas                        = edu.stanford.covert.cell.kb.Gene.empty(0, 0);
        rnaCompartments             = edu.stanford.covert.cell.kb.Compartment.empty(0, 0);
        rnaCoefficients             = [];

        proteinMonomers             = edu.stanford.covert.cell.kb.ProteinMonomer.empty(0, 0);
        proteinMonomerCompartments  = edu.stanford.covert.cell.kb.Compartment.empty(0, 0);
        proteinMonomerCoefficients  = [];

        proteinComplexs             = edu.stanford.covert.cell.kb.ProteinComplex.empty(0, 0);
        proteinComplexCompartments  = edu.stanford.covert.cell.kb.Compartment.empty(0, 0);
        proteinComplexCoefficients  = [];

        stableModificationReactions = edu.stanford.covert.cell.kb.Reaction.empty(0, 0);

        regulatedTranscriptionUnits = edu.stanford.covert.cell.kb.TranscriptionUnit.empty(0, 0);

        compartment                 = edu.stanford.covert.cell.kb.Compartment.empty(0, 0);

        stimuliRegulators           = edu.stanford.covert.cell.kb.Stimuli.empty(0, 0);
        metaboliteRegulators        = edu.stanford.covert.cell.kb.Metabolite.empty(0, 0);
        proteinMonomerRegulators    = edu.stanford.covert.cell.kb.ProteinMonomer.empty(0, 0);
        proteinComplexRegulators    = edu.stanford.covert.cell.kb.ProteinComplex.empty(0, 0);

        regulatedProteinMonomers    = edu.stanford.covert.cell.kb.ProteinMonomer.empty(0, 0);
        regulatedProteinComplexs    = edu.stanford.covert.cell.kb.ProteinComplex.empty(0, 0);

        parameters                  = edu.stanford.covert.cell.kb.Parameter.empty(0, 0);

        complexFormationProcess      = edu.stanford.covert.cell.kb.Process.empty(0, 0); %process where complex is formed
    end

    properties %(SetAccess = protected)
        empiricalFormula
        smiles
        charge
        pKa
        halfLife
        activationRule
        disulfideBondMonomers = cell(0, 1);
        disulfideBondCysteines = zeros(0, 2);
    end

    %computed properties
    properties %(SetAccess = protected)
        rnaComposition
        proteinMonomerComposition

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

        matureBaseCount
        matureMolecularWeight
        matureHalfLife
        matureDecayReaction

        numDistinctSubunits
        numSubunits

        numNucleicAcids
        numAminoAcids
        numRNASubunits
        numProteinMonomerSubunits
    end

    methods
        function this = ProteinComplex(knowledgeBase, wid, wholeCellModelID, name,...
                dnaFootprint, dnaFootprintBindingStrandedness, dnaFootprintRegionStrandedness, molecularInteraction, chemicalRegulation, ...
                subsystem, generalClassification, proteaseClassification, ...
                transporterClassification, ...
                activationRule,...
                disulfideBonds,...
                comments, crossReferences)

            if nargin == 0; return; end;

            this = edu.stanford.covert.cell.kb.ProteinComplex.empty(size(wid, 1), 0);
            this(size(wid, 1), 1) = edu.stanford.covert.cell.kb.ProteinComplex;
            for i = 1:size(wid, 1)
                this(i, 1).idx = i;
                this(i, 1).knowledgeBase = knowledgeBase;
                this(i, 1).wid = wid(i);
                this(i, 1).wholeCellModelID = wholeCellModelID{i};
                this(i, 1).name = name{i};
                if exist('comments','var') && ~isempty(comments); this(i,1).comments = comments{i}; end;
                if exist('crossReferences','var')
                    if size(crossReferences, 1) > 1
                        this(i, 1).crossReferences = crossReferences(i);
                    else
                        this(i, 1).crossReferences = struct;
                        fields = fieldnames(crossReferences);
                        for j = 1:size(fields, 1)
                            values = crossReferences.(fields{j});
                            this(i,1).crossReferences.(fields{j}) = values(i);
                        end
                    end
                end

                this(i, 1).dnaFootprint = dnaFootprint(i);
                this(i, 1).dnaFootprintBindingStrandedness = dnaFootprintBindingStrandedness{i};
                this(i, 1).dnaFootprintRegionStrandedness = dnaFootprintRegionStrandedness{i};
                this(i, 1).molecularInteraction = molecularInteraction{i};
                this(i, 1).chemicalRegulation = chemicalRegulation{i};
                this(i, 1).subsystem = subsystem{i};
                this(i, 1).generalClassification = generalClassification{i};
                this(i, 1).proteaseClassification = proteaseClassification{i};
                this(i, 1).transporterClassification = transporterClassification{i};
                this(i, 1).activationRule = activationRule{i};
                
                if disulfideBonds{i}
                    disulfideBondsArr = strsplit(';', disulfideBonds{i});
                    
                    this(i, 1).disulfideBondMonomers = cell(length(disulfideBondsArr), 1);
                    this(i, 1).disulfideBondCysteines = zeros(length(disulfideBondsArr), 2);
                    for j = 1:length(disulfideBondsArr)
                        tmp = strsplit(': ', disulfideBondsArr{j});
                        this(i, 1).disulfideBondMonomers{j, 1} = tmp{1};
                        tmp = strsplit('-', tmp{2});
                        this(i, 1).disulfideBondCysteines(j, :) = [str2double(tmp{1}(2:end)) str2double(tmp{2}(2:end))];
                    end
                end
            end
        end

        function serializeLinks(this)
            for i = 1:numel(this)
                this(i).prostheticGroups            = this.serializeLinksHelper(this(i).prostheticGroups);
                this(i).prostheticGroupCompartments = this.serializeLinksHelper(this(i).prostheticGroupCompartments);
                this(i).prostheticGroupCoefficients = this.serializeLinksHelper(this(i).prostheticGroupCoefficients);

                this(i).chaperoneSubstrates         = this.serializeLinksHelper(this(i).chaperoneSubstrates);
                this(i).chaperoneCompartments       = this.serializeLinksHelper(this(i).chaperoneCompartments);
                this(i).chaperoneCoefficients       = this.serializeLinksHelper(this(i).chaperoneCoefficients);

                this(i).metabolites                 = this.serializeLinksHelper(this(i).metabolites);
                this(i).metaboliteCompartments      = this.serializeLinksHelper(this(i).metaboliteCompartments);
                this(i).metaboliteCoefficients      = this.serializeLinksHelper(this(i).metaboliteCoefficients);
               
                this(i).rnas                        = this.serializeLinksHelper(this(i).rnas);
                this(i).rnaCompartments             = this.serializeLinksHelper(this(i).rnaCompartments);
                this(i).rnaCoefficients             = this.serializeLinksHelper(this(i).rnaCoefficients);

                this(i).proteinMonomers             = this.serializeLinksHelper(this(i).proteinMonomers);
                this(i).proteinMonomerCompartments  = this.serializeLinksHelper(this(i).proteinMonomerCompartments);
                this(i).proteinMonomerCoefficients  = this.serializeLinksHelper(this(i).proteinMonomerCoefficients);

                this(i).proteinComplexs             = this.serializeLinksHelper(this(i).proteinComplexs);
                this(i).proteinComplexCompartments  = this.serializeLinksHelper(this(i).proteinComplexCompartments);
                this(i).proteinComplexCoefficients  = this.serializeLinksHelper(this(i).proteinComplexCoefficients);

                this(i).stableModificationReactions = this.serializeLinksHelper(this(i).stableModificationReactions);

                this(i).regulatedTranscriptionUnits = this.serializeLinksHelper(this(i).regulatedTranscriptionUnits);

                this(i).compartment                 = this.serializeLinksHelper(this(i).compartment);

                this(i).stimuliRegulators           = this.serializeLinksHelper(this(i).stimuliRegulators);
                this(i).metaboliteRegulators        = this.serializeLinksHelper(this(i).metaboliteRegulators);
                this(i).proteinMonomerRegulators    = this.serializeLinksHelper(this(i).proteinMonomerRegulators);
                this(i).proteinComplexRegulators    = this.serializeLinksHelper(this(i).proteinComplexRegulators);

                this(i).regulatedProteinMonomers    = this.serializeLinksHelper(this(i).regulatedProteinMonomers);
                this(i).regulatedProteinComplexs    = this.serializeLinksHelper(this(i).regulatedProteinComplexs);

                this(i).parameters                  = this.serializeLinksHelper(this(i).parameters);
                
                this(i).complexFormationProcess     = this.serializeLinksHelper(this(i).complexFormationProcess);

                serializeLinks@edu.stanford.covert.cell.kb.PhysicalObject(this(i));
                serializeLinks@edu.stanford.covert.cell.kb.Protein(this(i));
            end
        end
        
        function deserializeLinks(this, kb)
            for i = 1:numel(this)
                this(i).prostheticGroups            = this.deserializeLinksHelper(this(i).prostheticGroups, kb.metabolites);
                this(i).prostheticGroupCompartments = this.deserializeLinksHelper(this(i).prostheticGroupCompartments, kb.compartments);
                this(i).prostheticGroupCoefficients = this.deserializeLinksHelper(this(i).prostheticGroupCoefficients);

                this(i).chaperoneSubstrates         = this.deserializeLinksHelper(this(i).chaperoneSubstrates, kb.proteinMonomers);
                this(i).chaperoneCompartments       = this.deserializeLinksHelper(this(i).chaperoneCompartments, kb.compartments);
                this(i).chaperoneCoefficients       = this.deserializeLinksHelper(this(i).chaperoneCoefficients);

                this(i).metabolites                 = this.deserializeLinksHelper(this(i).metabolites, kb.metabolites);
                this(i).metaboliteCompartments      = this.deserializeLinksHelper(this(i).metaboliteCompartments, kb.compartments);
                this(i).metaboliteCoefficients      = this.deserializeLinksHelper(this(i).metaboliteCoefficients);
               
                this(i).rnas                        = this.deserializeLinksHelper(this(i).rnas, kb.genes);
                this(i).rnaCompartments             = this.deserializeLinksHelper(this(i).rnaCompartments, kb.compartments);
                this(i).rnaCoefficients             = this.deserializeLinksHelper(this(i).rnaCoefficients);

                this(i).proteinMonomers             = this.deserializeLinksHelper(this(i).proteinMonomers, kb.proteinMonomers);
                this(i).proteinMonomerCompartments  = this.deserializeLinksHelper(this(i).proteinMonomerCompartments, kb.compartments);
                this(i).proteinMonomerCoefficients  = this.deserializeLinksHelper(this(i).proteinMonomerCoefficients);

                this(i).proteinComplexs             = this.deserializeLinksHelper(this(i).proteinComplexs, kb.proteinComplexs);
                this(i).proteinComplexCompartments  = this.deserializeLinksHelper(this(i).proteinComplexCompartments, kb.compartments);
                this(i).proteinComplexCoefficients  = this.deserializeLinksHelper(this(i).proteinComplexCoefficients);

                this(i).stableModificationReactions = this.deserializeLinksHelper(this(i).stableModificationReactions, kb.reactions);

                this(i).regulatedTranscriptionUnits = this.deserializeLinksHelper(this(i).regulatedTranscriptionUnits, kb.transcriptionUnits);

                this(i).compartment                 = this.deserializeLinksHelper(this(i).compartment, kb.compartments);

                this(i).stimuliRegulators           = this.deserializeLinksHelper(this(i).stimuliRegulators, kb.stimulis);
                this(i).metaboliteRegulators        = this.deserializeLinksHelper(this(i).metaboliteRegulators, kb.metabolites);
                this(i).proteinMonomerRegulators    = this.deserializeLinksHelper(this(i).proteinMonomerRegulators, kb.proteinMonomers);
                this(i).proteinComplexRegulators    = this.deserializeLinksHelper(this(i).proteinComplexRegulators, kb.proteinComplexs);

                this(i).regulatedProteinMonomers    = this.deserializeLinksHelper(this(i).regulatedProteinMonomers, kb.proteinMonomers);
                this(i).regulatedProteinComplexs    = this.deserializeLinksHelper(this(i).regulatedProteinComplexs, kb.proteinComplexs);

                this(i).parameters                  = this.deserializeLinksHelper(this(i).parameters, kb.parameters);
                
                this(i).complexFormationProcess     = this.deserializeLinksHelper(this(i).complexFormationProcess, kb.processes);
                
                deserializeLinks@edu.stanford.covert.cell.kb.PhysicalObject(this(i), kb);
                deserializeLinks@edu.stanford.covert.cell.kb.Protein(this(i), kb);
            end
        end

        function deleteLinks(this)
            for i = 1:numel(this)
                this(i).prostheticGroups            = [];
                this(i).prostheticGroupCompartments = [];
                this(i).prostheticGroupCoefficients = [];

                this(i).chaperoneSubstrates         = [];
                this(i).chaperoneCompartments       = [];
                this(i).chaperoneCoefficients       = [];

                this(i).metabolites                 = [];
                this(i).metaboliteCompartments      = [];
                this(i).metaboliteCoefficients      = [];
               
                this(i).rnas                        = [];
                this(i).rnaCompartments             = [];
                this(i).rnaCoefficients             = [];

                this(i).proteinMonomers             = [];
                this(i).proteinMonomerCompartments  = [];
                this(i).proteinMonomerCoefficients  = [];

                this(i).proteinComplexs             = [];
                this(i).proteinComplexCompartments  = [];
                this(i).proteinComplexCoefficients  = [];

                this(i).stableModificationReactions = [];

                this(i).regulatedTranscriptionUnits = [];

                this(i).compartment                 = [];

                this(i).stimuliRegulators           = [];
                this(i).metaboliteRegulators        = [];
                this(i).proteinMonomerRegulators    = [];
                this(i).proteinComplexRegulators    = [];

                this(i).regulatedProteinMonomers    = [];
                this(i).regulatedProteinComplexs    = [];

                this(i).parameters                  = [];

                deleteLinks@edu.stanford.covert.cell.kb.PhysicalObject(this(i));
                deleteLinks@edu.stanford.covert.cell.kb.Protein(this(i));
            end
        end

        function value = get.empiricalFormula(this)
            throw(MException('ProteinComplex:error', 'property is not defined'));
        end

        function value = get.smiles(this)
            throw(MException('ProteinComplex:error', 'property is not defined'));
        end

        function value = get.charge(this)
            throw(MException('ProteinComplex:error', 'property is not defined'));
        end

        function value = get.pKa(this)
            throw(MException('ProteinComplex:error', 'property is not defined'));
        end

        function value = get.halfLife(this)
            %retrieve
            if ~isempty(this.halfLife)
                value = this.halfLife;
                return;
            end
            
            %calculate
            %rnaHalfLives = zeros(1,length(this.rnas));
            %for i = 1:length(this.rnas)
            %    rnaHalfLives(i) = this.rnas(i).halfLife;
            %end

            monomerHalfLives = zeros(1,length(this.proteinMonomers));
            for i = 1:length(this.proteinMonomers)
                monomerHalfLives(i) = this.proteinMonomers(i).matureSequenceHalfLife;
            end

            complexHalfLives = zeros(1,length(this.proteinComplexs));
            for i = 1:length(this.proteinComplexs)
                complexHalfLives(i) = this.proteinComplexs(i).matureHalfLife;
            end

            halfLives = [monomerHalfLives complexHalfLives];
            coefficients = [this.proteinMonomerCoefficients; this.proteinComplexCoefficients];

            if ~isempty(coefficients)
                value = halfLives * coefficients / sum(coefficients);
            else
                value = Inf;
            end
            
            %store
            this.halfLife = value;
        end

        function value = get.matureHalfLife(this)
            %retrieve
            if ~isempty(this.matureHalfLife)
                value = this.matureHalfLife;
                return;
            end
            
            %calculate
            value = this.halfLife;
            
            %store
            this.matureHalfLife = value;
        end

        %leave on composition tree
        function value = get.rnaComposition(this)
            %retrieve
            if ~isempty(this.rnaComposition)
                value = this.rnaComposition;
                return;
            end
            
            %calculate
            value = this.rnaComposition_Helper();
            
            %store
            this.rnaComposition = value;
        end

        %helper method because MATLAB won't let you call another objects
        %get.rnaComposition method from within get.rnaComposition
        function value = rnaComposition_Helper(this)
            value = zeros(1, this.knowledgeBase.numGenes, this.knowledgeBase.numCompartments);

            if ~isempty(this.rnas)
                value(sub2ind(size(value),...
                    ones(1, numel(this.rnas)),...
                    double([this.rnas.idx]),...
                    double([this.rnaCompartments.idx]))) = ...
                    this.rnaCoefficients;
            end

            if ~isempty(this.proteinComplexs)
                value = value + ...
                    sum(permute(reshape([this.proteinComplexs.rnaComposition],this.knowledgeBase.numGenes, [], this.knowledgeBase.numCompartments),[2 1 3]) .* ...
                    repmat(this.proteinComplexCoefficients, [1 this.knowledgeBase.numGenes this.knowledgeBase.numCompartments]),1);
            end
        end

        %leave on composition tree
        function value = get.proteinMonomerComposition(this)
            %retrieve
            if ~isempty(this.proteinMonomerComposition)
                value = this.proteinMonomerComposition;
                return;
            end
            
            %calculate
            value = this.proteinMonomerComposition_Helper();
            
            %store
            this.proteinMonomerComposition = value;
        end

        %helper method because MATLAB won't let you call another objects
        %get.proteinMonomerComposition method from within
        %get.proteinMonomerComposition
        function value = proteinMonomerComposition_Helper(this)
            value = zeros(1, this.knowledgeBase.numProteinMonomers, this.knowledgeBase.numCompartments);

            if ~isempty(this.proteinMonomers)
                value(sub2ind(size(value),...
                    ones(1, numel(this.proteinMonomers)),...
                    double([this.proteinMonomers.idx]),...
                    double([this.proteinMonomerCompartments.idx]))) = ...
                    this.proteinMonomerCoefficients;
            end

            if ~isempty(this.proteinComplexs)
                value = value + ...
                    sum(permute(reshape([this.proteinComplexs.proteinMonomerComposition],this.knowledgeBase.numProteinMonomers, [], this.knowledgeBase.numCompartments),[2 1 3]) .* ...
                    repmat(this.proteinComplexCoefficients, [1 this.knowledgeBase.numProteinMonomers this.knowledgeBase.numCompartments]),1);
            end
        end

        function value = get.molecularWeight(this)
            %retrieve
            if ~isempty(this.molecularWeight)
                value = this.molecularWeight;
                return;
            end
            
            %calculate
            metaboliteMolecularWeights = [this.metabolites.molecularWeight];
            
            rnaMolecularWeights = zeros(1, length(this.rnas));
            for i = 1:length(this.rnas)
                rnaMolecularWeights(i) = this.rnas(i).matureMolecularWeight;
            end

            monomerMolecularWeights = zeros(1, length(this.proteinMonomers));
            for i = 1:length(this.proteinMonomers)
                monomerMolecularWeights(i) = this.proteinMonomers(i).matureSequenceMolecularWeight;
            end

            complexMolecularWeights = zeros(1, length(this.proteinComplexs));
            for i = 1:length(this.proteinComplexs)
                complexMolecularWeights(i) = this.proteinComplexs(i).matureMolecularWeight;
            end

            molecularWeights = [metaboliteMolecularWeights rnaMolecularWeights monomerMolecularWeights complexMolecularWeights];
            coefficients = [this.metaboliteCoefficients; this.rnaCoefficients; this.proteinMonomerCoefficients; this.proteinComplexCoefficients];

            if ~isempty(coefficients)
                value = molecularWeights * coefficients;
            else
                value = 0;
            end
            
            %store
            this.molecularWeight = value;
        end

        function value = get.matureMolecularWeight(this)
            %retrieve
            if ~isempty(this.matureMolecularWeight)
                value = this.matureMolecularWeight;
                return;
            end
            
            %calculate
            value = this.molecularWeight;
            if ~isempty(this.prostheticGroups)
                value = value + [this.prostheticGroups.molecularWeight] * ...
                    max(1, this.prostheticGroupCoefficients');
            end
            
            %store
            this.matureMolecularWeight = value;
        end

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
            value = this.molecularWeight / this.density;
            
            %store
            this.volume = value;
        end

        function value = get.pI(this)
            %retrieve
            if ~isempty(this.pI)
                value = this.pI;
                return;
            end
            
            %calculate
            value = this.meanSubunitPropertyValue('pI', true);
            
            %store
            this.pI = value;
        end

        function value = get.extinctionCoefficient(this)
            %retrieve
            if ~isempty(this.extinctionCoefficient)
                value = this.extinctionCoefficient;
                return;
            end
            
            %calculate
            value = this.meanSubunitPropertyValue('extinctionCoefficient', true);
            
            %store
            this.extinctionCoefficient = value;
        end

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

        function value = get.negAA(this)
            %retrieve
            if ~isempty(this.negAA)
                value = this.negAA;
                return;
            end
            
            %calculate
            value = this.meanSubunitPropertyValue('negAA', false);
            
            %store
            this.negAA = value;
        end

        function value = get.posAA(this)
            %retrieve
            if ~isempty(this.posAA)
                value = this.posAA;
                return;
            end
            
            %calculate
            value = this.meanSubunitPropertyValue('posAA', false);
            
            %store
            this.posAA = value;
        end

        %Instability index (II)
        function value = get.instabilityIndex(this)
            %retrieve
            if ~isempty(this.instabilityIndex)
                value = this.instabilityIndex;
                return;
            end
            
            %calculate
            value = this.meanSubunitPropertyValue('instabilityIndex', true);
            
            %store
            this.instabilityIndex = value;
        end

        function value = get.stable(this)
            %retrieve
            if ~isempty(this.stable)
                value = this.stable;
                return;
            end
            
            %calculate
            value = (this.instabilityIndex<40);
            
            %store
            this.stable = value;
        end

        %Aliphatic index
        function value = get.aliphaticIndex(this)
            %retrieve
            if ~isempty(this.aliphaticIndex)
                value = this.aliphaticIndex;
                return;
            end
            
            %calculate
            value = this.meanSubunitPropertyValue('aliphaticIndex', true);
            
            %store
            this.aliphaticIndex = value;
        end

        %GRAVY (Grand Average of Hydropathy)
        function value = get.gravy(this)
            %retrieve
            if ~isempty(this.gravy)
                value = this.gravy;
                return;
            end
            
            %calculate
            value = this.meanSubunitPropertyValue('gravy', true);
            
            %store
            this.gravy = value;
        end

        %base count (numbers of amino acids)
        function value = get.baseCount(this)   
            %retrieve
            if ~isempty(this.baseCount)
                value = this.baseCount;
                return;
            end
            
            %calculate
            rnaBaseCounts = zeros(length(this.rnas), this.knowledgeBase.numMetabolites);
            for i = 1:length(this.rnas)
                rnaBaseCounts(i, :) = this.rnas(i).matureBaseCount;
            end

            monomerBaseCounts = zeros(length(this.proteinMonomers), this.knowledgeBase.numMetabolites);
            for i = 1:length(this.proteinMonomers)
                monomerBaseCounts(i, :) = this.proteinMonomers(i).matureSequenceBaseCount;
            end

            complexBaseCounts = zeros(length(this.proteinComplexs), this.knowledgeBase.numMetabolites);
            for i = 1:length(this.proteinComplexs)
                complexBaseCounts(i, :) = this.proteinComplexs(i).matureBaseCount;
            end

            baseCounts = [rnaBaseCounts; monomerBaseCounts; complexBaseCounts];
            coefficients = [this.rnaCoefficients; this.proteinMonomerCoefficients; this.proteinComplexCoefficients];

            if ~isempty(coefficients)
                value = coefficients' * baseCounts;
            else
                value = zeros(1, this.knowledgeBase.numMetabolites);
            end
            
            if ~isempty(this.metabolites)
                value([this.metabolites.idx]) = ...
                    value([this.metabolites.idx]) + ...
                    this.metaboliteCoefficients';
            end
            
            %store
            this.baseCount = value;
        end

        function value = get.matureBaseCount(this)
            %retrieve
            if ~isempty(this.matureBaseCount)
                value = this.matureBaseCount;
                return;
            end
            
            %calculate
            value = this.baseCount;
            if ~isempty(this.prostheticGroups)
                value([this.prostheticGroups.idx]) = ...
                    value([this.prostheticGroups.idx]) + ...
                    max(1, this.prostheticGroupCoefficients);
            end
            
            %store
            this.matureBaseCount = value;
        end

        %base count (numbers of amino acids)
        function value = get.cumulativeBaseCount(this)
            throw(MException('ProteinComplex:error','property is not defined'));
        end

        %codon count
        function value = get.tRNACount(this)
            %retrieve
            if ~isempty(this.tRNACount)
                value = this.tRNACount;
                return;
            end
            
            %calculate
            monomerTRNACounts = zeros(length(this.proteinMonomers),this.knowledgeBase.numMetabolites);
            for i = 1:length(this.proteinMonomers)
                switch this.proteinMonomers(i).gene.type
                    case 'mRNA', monomerTRNACounts(i,:) = this.proteinMonomers(i).matureSequenceBaseCount;
                end
            end

            complexTRNACounts = zeros(length(this.proteinComplexs),this.knowledgeBase.numMetabolites);
            for i = 1:length(this.proteinComplexs)
                complexTRNACounts(i,:) = this.proteinComplexs.matureBaseCount;
            end

            TRNACounts = [monomerTRNACounts;complexTRNACounts];
            coefficients = [this.proteinMonomerCoefficients;this.proteinComplexCoefficients];

            if ~isempty(coefficients)
                value = coefficients'*TRNACounts;
            else
                value = zeros(1,this.knowledgeBase.numMetabolites);
            end
            
            %store
            this.tRNACount = value;
        end

        function value = get.decayReaction(this)
            %retrieve
            if ~isempty(this.decayReaction)
                value = this.decayReaction;
                return;
            end
            
            %calculate
            value = this.decayReaction_Helper();
            
            %store
            this.decayReaction = value;
        end
        
        function value = decayReaction_Helper(this)
            value = zeros(1, this.knowledgeBase.numMetabolites);
            
            for i = 1:length(this.proteinComplexs)
                value = value + this.proteinComplexCoefficients(i) * this.proteinComplexs(i).decayReaction;
            end
            
            if ~isempty(this.metabolites)
                value([this.metabolites.idx]) =   ...
                    + value([this.metabolites.idx]) ...
                    + this.metaboliteCoefficients';
            end
        end

        function value = get.matureDecayReaction(this)
            %retrieve
            if ~isempty(this.matureDecayReaction)
                value = this.matureDecayReaction;
                return;
            end
            
            %calculate
            value = this.decayReaction;
            if ~isempty(this.prostheticGroups)
                value([this.prostheticGroups.idx]) = ...
                    + value([this.prostheticGroups.idx]) ...
                    + max(1, this.prostheticGroupCoefficients);
            end
            
            %store
            this.matureDecayReaction = value;
        end

        %number of distinct subunits
        function value = get.numDistinctSubunits(this)
            %retrieve
            if ~isempty(this.numDistinctSubunits)
                value = this.numDistinctSubunits;
                return;
            end
            
            %calculate
            value = length(this.rnas)+length(this.proteinMonomers)+length(this.proteinComplexs);
            
            %store
            this.numDistinctSubunits = value;
        end

        %number of subunits
        function value = get.numSubunits(this)
            %retrieve
            if ~isempty(this.numSubunits)
                value = this.numSubunits;
                return;
            end
            
            %calculate
            value = sum([this.proteinMonomerCoefficients;this.proteinComplexCoefficients]);
            if ~isempty(value)
                value = 0;
            end
            
            %store
            this.numSubunits = value;
        end

        function value = get.numNucleicAcids(this)
            %retrieve
            if ~isempty(this.numNucleicAcids)
                value = this.numNucleicAcids;
                return;
            end            
            
            %calculate
           value = this.computeNumNucleicAcids();
           
           %store
            this.numNucleicAcids = value;
        end

        function value = computeNumNucleicAcids(this)
            value = 0;

            for i = 1:length(this.rnas)
                value = value + this.rnas(i).sequenceLength * this.rnaCoefficients(i);
            end

            for i = 1:length(this.proteinComplexs)
                value = value + this.proteinComplexs(i).computeNumNucleicAcids() * this.proteinComplexCoefficients(i);
            end
        end

        function value = get.numAminoAcids(this)
            %retrieve
            if ~isempty(this.numAminoAcids)
                value = this.numAminoAcids;
                return;
            end
            
            %calculate
            value = this.computeNumAminoAcids();
            
            %store
            this.numAminoAcids = value;
        end

        function value = computeNumAminoAcids(this)
            value = 0;

            for i = 1:length(this.proteinMonomers)
                value = value + this.proteinMonomers(i).matureSequenceLength * this.proteinMonomerCoefficients(i);
            end

            for i = 1:length(this.proteinComplexs)
                value = value + this.proteinComplexs(i).computeNumAminoAcids() * this.proteinComplexCoefficients(i);
            end
        end

        function value = get.numRNASubunits(this)
            %retrieve
            if ~isempty(this.numRNASubunits)
                value = this.numRNASubunits;
                return;
            end
            
            %calculate
            value = this.computeNumRNASubunits();
            
            %store
            this.numRNASubunits = value;
        end

        function value = computeNumRNASubunits(this)
            value = sum(this.rnaCoefficients);

            for i = 1:length(this.proteinComplexs)
                value = value + this.proteinComplexs(i).computeNumRNASubunits() * this.proteinComplexCoefficients(i);
            end
        end

        function value = get.numProteinMonomerSubunits(this)
            %retrieve
            if ~isempty(this.numProteinMonomerSubunits)
                value = this.numProteinMonomerSubunits;
                return;
            end
            
            %calculate
            value = this.computeNumProteinMonomerSubunits();
            
            %store
            this.numProteinMonomerSubunits = value;
        end

        function value = computeNumProteinMonomerSubunits(this)
            value = sum(this.proteinMonomerCoefficients);

            for i = 1:length(this.proteinComplexs)
                value = value + this.proteinComplexs(i).computeNumProteinMonomerSubunits() * this.proteinComplexCoefficients(i);
            end
        end

        %average/sum property value over subunits
        function value = meanSubunitPropertyValue(this, property, average)
            values = [];

            if ~isempty(this.rnas)
                values = [values this.rnas.(property)];
            end

            if ~isempty(this.proteinMonomers)
                values = [values this.proteinMonomers.(property)];
            end

            if ~isempty(this.proteinComplexs)
                values = [values this.proteinComplexs.(property)];
            end

            coefficients = [this.rnaCoefficients;this.proteinMonomerCoefficients;this.proteinComplexCoefficients];
            if ~isempty(coefficients)
                if average
                    value = values*coefficients/sum(coefficients);
                else
                    value = values*coefficients;
                end
            else
                if average
                    value = NaN;
                else
                    value = 0;
                end
            end
        end
    end
end