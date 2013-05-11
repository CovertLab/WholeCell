% Defines a metabolite
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/7/2009
classdef Metabolite < edu.stanford.covert.cell.kb.PhysicalObject
    properties
        reactions           = edu.stanford.covert.cell.kb.Reaction.empty(0, 0);

        coenzymeReactions   = edu.stanford.covert.cell.kb.Reaction.empty(0, 0);

        biomassCompartments = edu.stanford.covert.cell.kb.Compartment.empty(0, 0);
        biomassCoefficients = [];

        mediaCompartments   = edu.stanford.covert.cell.kb.Compartment.empty(0, 0);
        mediaConcentrations = []; %(mM)
        mediaInitialTimes   = [];
        mediaFinalTimes     = [];

        regulatedProteinMonomers = edu.stanford.covert.cell.kb.ProteinMonomer.empty(0, 0);
        regulatedProteinComplexs = edu.stanford.covert.cell.kb.ProteinComplex.empty(0, 0);
    end

    properties %(SetAccess = protected)
        traditionalName
        iupacName
        category
        subcategory
        empiricalFormula
        smiles
        charge
        pKa
        logP
        logD
        halfLife
        exchangeLowerBound
        exchangeUpperBound
    end

    %computed properties
    properties %(SetAccess = protected)
        molecularWeight
        density
        volume
        pI
        extinctionCoefficient
        absorbanceFactor
    end

    methods
        function this = Metabolite(knowledgeBase, wid,wholeCellModelID, name, ...
                traditionalName, iupacName, ...
                category, subcategory, ...
                empiricalFormula, smiles, charge, hydrophobic, pKa, pI, logP, logD, volume, molecularWeight, ...
                exchangeLowerBound, exchangeUpperBound, ...
                comments, crossReferences)

            if nargin == 0; return; end;

            this = edu.stanford.covert.cell.kb.Metabolite.empty(size(wid,1),0);
            this(size(wid,1),1) = edu.stanford.covert.cell.kb.Metabolite;
            for i = 1:size(wid,1)
                this(i, 1).idx = i;
                this(i, 1).knowledgeBase = knowledgeBase;
                this(i, 1).wid = wid(i);
                this(i, 1).wholeCellModelID = wholeCellModelID{i};
                this(i, 1).name = name{i};
                if exist('comments','var') && ~isempty(comments); this(i, 1).comments = comments{i}; end;
                if exist('crossReferences','var')
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
                this(i,1).traditionalName = traditionalName{i};
                this(i,1).iupacName = iupacName{i};
                this(i,1).category = category{i};
                this(i,1).subcategory = subcategory{i};
                this(i,1).empiricalFormula = this.parseEmpiricalFormula(empiricalFormula{i});
                this(i,1).smiles = smiles{i};
                this(i,1).charge = charge(i);
                this(i,1).hydrophobic = hydrophobic(i);
                this(i,1).pKa = edu.stanford.covert.util.parseDoubles(', ', pKa{i});
                %this(i,1).pI = pI(i);
                this(i,1).logP = logP(i);
                this(i,1).logD = logD(i);
                %this(i,1).volume = volume(i);
                this(i,1).exchangeLowerBound = exchangeLowerBound(i);
                this(i,1).exchangeUpperBound = exchangeUpperBound(i);
            end
        end

        function serializeLinks(this)
            for i = 1:numel(this)
                this(i).reactions           = this.serializeLinksHelper(this(i).reactions);

                this(i).coenzymeReactions   = this.serializeLinksHelper(this(i).coenzymeReactions);

                this(i).biomassCompartments = this.serializeLinksHelper(this(i).biomassCompartments);
                this(i).biomassCoefficients = this.serializeLinksHelper(this(i).biomassCoefficients);

                this(i).mediaCompartments   = this.serializeLinksHelper(this(i).mediaCompartments);
                this(i).mediaConcentrations = this.serializeLinksHelper(this(i).mediaConcentrations);
                this(i).mediaInitialTimes   = this.serializeLinksHelper(this(i).mediaInitialTimes);
                this(i).mediaFinalTimes     = this.serializeLinksHelper(this(i).mediaFinalTimes);

                this(i).regulatedProteinMonomers = this.serializeLinksHelper(this(i).regulatedProteinMonomers);
                this(i).regulatedProteinComplexs = this.serializeLinksHelper(this(i).regulatedProteinComplexs);

                serializeLinks@edu.stanford.covert.cell.kb.PhysicalObject(this(i));
            end
        end
        
        function deserializeLinks(this, kb)
            for i = 1:numel(this)
                this(i).reactions           = this.deserializeLinksHelper(this(i).reactions, kb.reactions);

                this(i).coenzymeReactions   = this.deserializeLinksHelper(this(i).coenzymeReactions, kb.reactions);

                this(i).biomassCompartments = this.deserializeLinksHelper(this(i).biomassCompartments, kb.compartments);
                this(i).biomassCoefficients = this.deserializeLinksHelper(this(i).biomassCoefficients);

                this(i).mediaCompartments   = this.deserializeLinksHelper(this(i).mediaCompartments, kb.compartments);
                this(i).mediaConcentrations = this.deserializeLinksHelper(this(i).mediaConcentrations);
                this(i).mediaInitialTimes   = this.deserializeLinksHelper(this(i).mediaInitialTimes);
                this(i).mediaFinalTimes     = this.deserializeLinksHelper(this(i).mediaFinalTimes);

                this(i).regulatedProteinMonomers = this.deserializeLinksHelper(this(i).regulatedProteinMonomers, kb.proteinMonomers);
                this(i).regulatedProteinComplexs = this.deserializeLinksHelper(this(i).regulatedProteinComplexs, kb.proteinComplexs);
                
                deserializeLinks@edu.stanford.covert.cell.kb.PhysicalObject(this(i), kb);
            end
        end

        function deleteLinks(this)
            for i = 1:numel(this)
                this(i).reactions           = [];

                this(i).coenzymeReactions   = [];

                this(i).biomassCompartments = [];
                this(i).biomassCoefficients = [];

                this(i).mediaCompartments   = [];
                this(i).mediaConcentrations = [];
                this(i).mediaInitialTimes   = [];
                this(i).mediaFinalTimes     = [];

                this(i).regulatedProteinMonomers = [];
                this(i).regulatedProteinComplexs = [];

                deleteLinks@edu.stanford.covert.cell.kb.PhysicalObject(this(i));
            end
        end

        function value = get.halfLife(this)
            throw(MException('Metabolite:error', 'property is not defined'));
        end

        function value = get.molecularWeight(this)
            %retrieve
            if ~isempty(this.molecularWeight)
                value = this.molecularWeight;
                return;
            end
                
            %compute
            fields = fieldnames(this.empiricalFormula);
            value = 0;
            for i = 1:length(fields)
                if ~isfield(edu.stanford.covert.util.ConstantUtil.elements,fields{i})
                    value = 0;
                    return;
                end
                value = value+...
                    this.empiricalFormula.(fields{i})*edu.stanford.covert.util.ConstantUtil.elements.(fields{i});
            end
            
            %store
            this.molecularWeight = value;
        end

        function value = get.density(this)
            throw(MException('Metabolite:error', 'property is not defined'));
        end

        function value = get.volume(this)
            throw(MException('Metabolite:error', 'property is not defined'));
        end

        function value = get.pI(this)
            value = mean(this.pKa);
        end

        function value = get.extinctionCoefficient(this)
            throw(MException('Metabolite:error', 'property is not defined'));
        end

        function value = get.absorbanceFactor(this)
            throw(MException('Metabolite:error', 'property is not defined'));
        end
    end
    
    methods (Static = true)
        function structure = parseEmpiricalFormula(string)
            if ~regexp(string, '^(([A-Z][a-z]*)(\d+))+$');
                throw(MException('Metabolite:invalidEmpiricalFormula', 'Empirical formula must match pattern ''^(([A-Z][a-z]*)(\d+))+$''.'))
            end
            structure = struct;
            tokens = regexp(string, '([A-Z][a-z]*)(\d+)', 'tokens');
            for i = 1:length(tokens)
                tokens{i}{1} = tokens{i}{1};
                tokens{i}{2} = str2double(tokens{i}{2});
                if isfield(structure,tokens{i}{1})
                    structure.(tokens{i}{1}) = structure.(tokens{i}{1}) + tokens{i}{2};
                else
                    structure.(tokens{i}{1}) = tokens{i}{2};
                end
            end
        end
    end
end