% Defines a compartment
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 11/9/2009
classdef Compartment < edu.stanford.covert.cell.kb.PhysicalObject
    properties
        genes              = edu.stanford.covert.cell.kb.Gene.empty(0, 0);
        transcriptionUnits = edu.stanford.covert.cell.kb.TranscriptionUnit.empty(0, 0);
        proteinMonomers    = edu.stanford.covert.cell.kb.ProteinMonomer.empty(0, 0);
        proteinComplexs    = edu.stanford.covert.cell.kb.ProteinComplex.empty(0, 0);
    end

    properties %(SetAccess = protected)
        empiricalFormula
        smiles
        charge
        pKa
        halfLife
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
        function this = Compartment(knowledgeBase, wid, wholeCellModelID, name, ...
                comments, crossReferences)

            if nargin == 0; return; end;

            this = edu.stanford.covert.cell.kb.Compartment.empty(size(wid, 1), 0);
            this(size(wid, 1), 1) = edu.stanford.covert.cell.kb.Compartment;
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
            end
        end

        function serializeLinks(this)
            for i = 1:numel(this)
                this(i).genes              = this.serializeLinksHelper(this(i).genes);
                this(i).transcriptionUnits = this.serializeLinksHelper(this(i).transcriptionUnits);
                this(i).proteinMonomers    = this.serializeLinksHelper(this(i).proteinMonomers);
                this(i).proteinComplexs    = this.serializeLinksHelper(this(i).proteinComplexs);

                serializeLinks@edu.stanford.covert.cell.kb.PhysicalObject(this(i));
            end
        end        
        
        function deserializeLinks(this, kb)
            for i = 1:numel(this)
                this(i).genes              = this.deserializeLinksHelper(this(i).genes, kb.genes);
                this(i).transcriptionUnits = this.deserializeLinksHelper(this(i).transcriptionUnits, kb.transcriptionUnits);
                this(i).proteinMonomers    = this.deserializeLinksHelper(this(i).proteinMonomers, kb.proteinMonomers);
                this(i).proteinComplexs    = this.deserializeLinksHelper(this(i).proteinComplexs, kb.proteinComplexs);
                
                deserializeLinks@edu.stanford.covert.cell.kb.PhysicalObject(this(i), kb);
            end
        end

		function deleteLinks(this)
            for i = 1:numel(this)
                this(i).genes              = [];
                this(i).transcriptionUnits = [];
                this(i).proteinMonomers    = [];
                this(i).proteinComplexs    = [];

                deleteLinks@edu.stanford.covert.cell.kb.PhysicalObject(this(i));
            end
        end

        function value = get.halfLife(this)
            throw(MException('Compartment:error', 'property is not defined'));
        end

        function value = get.molecularWeight(this)
            throw(MException('Compartment:error', 'property is not defined'));
        end

        function value = get.density(this)
            throw(MException('Compartment:error', 'property is not defined'));
        end

        function value = get.volume(this)
            throw(MException('Compartment:error', 'property is not defined'));
        end

        function value = get.pI(this)
            throw(MException('Compartment:error', 'property is not defined'));
        end

        function value = get.extinctionCoefficient(this)
            throw(MException('Compartment:error', 'property is not defined'));
        end

        function value = get.absorbanceFactor(this)
            throw(MException('Compartment:error', 'property is not defined'));
        end
    end
end