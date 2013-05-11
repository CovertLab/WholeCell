% Defines a physical object (chemical, dsDNA, protein etc.). Base class for
% - Metabolite
% - Polymer
%   - NucleicAcid
%   - ProteinMonomer
% - ProteinComplex
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/7/2009
classdef PhysicalObject < edu.stanford.covert.cell.kb.KnowledgeBaseObject
    properties
        compartments = edu.stanford.covert.cell.kb.Compartment.empty(0, 0);
    end

    properties (Abstract = true) %SetAccess = protected
        empiricalFormula        %
        smiles                  %
        charge                  %
        pKa                     %
        halfLife                %s
    end

    %computed properties
    properties (Abstract = true) %SetAccess = protected
        molecularWeight         %g/mol
        density                 %g/cm^3
        volume                  %cm^3/mol
        pI                      %
        extinctionCoefficient   %at 260 nm
        absorbanceFactor        %mmol L^-1 at 260nm
    end

    %computed properties
    properties %(SetAccess = protected)
        halfLifeTimeConstant    %s^-1
    end

    properties %(SetAccess = protected)
        hydrophobic = false
    end

    methods
        function this = PhysicalObject(knowledgeBase, wid, wholeCellModelID, name, ...
                comments, crossReferences)

            if nargin == 0; return; end;

            this = edu.stanford.covert.cell.kb.PhysicalObject.empty(size(wid, 1), 0);
            this(size(wid, 1), 1) = edu.stanford.covert.cell.kb.PhysicalObject;
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
                        for j = 1:size(fields,1)
                            values = crossReferences.(fields{j});
                            this(i, 1).crossReferences.(fields{j}) = values(i);
                        end
                    end
                end
            end
        end
        
        function serializeLinks(this)
            for i = 1:numel(this)
                this(i).compartments = this.serializeLinksHelper(this(i).compartments);
                
                serializeLinks@edu.stanford.covert.cell.kb.KnowledgeBaseObject(this(i));
            end
        end
        
        function deserializeLinks(this, kb)
            for i = 1:numel(this)
                this(i).compartments = this.deserializeLinksHelper(this(i).compartments, kb.compartments);
                
                deserializeLinks@edu.stanford.covert.cell.kb.KnowledgeBaseObject(this(i), kb);
            end
        end

        function deleteLinks(this)
            for i = 1:numel(this)
                this(i).compartments = [];
                
                deleteLinks@edu.stanford.covert.cell.kb.KnowledgeBaseObject(this(i));
            end
        end
        
        function value = get.halfLifeTimeConstant(this)
            %retrieve
            if isempty(this.halfLifeTimeConstant)
                value = this.halfLifeTimeConstant;
                return;
            end
            
            %compute
            value = ln(2) / this.halfLife;
            
            %store
            this.halfLifeTimeConstant = value;
        end
    end
end