% Defines a stimulus
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 11/30/2009
classdef Stimuli < edu.stanford.covert.cell.kb.KnowledgeBaseObject
    properties
        compartments    = edu.stanford.covert.cell.kb.Compartment.empty(0, 0);
        values          = [];
        initialTimes    = [];
        finalTimes      = [];

        regulatedProteinMonomers = edu.stanford.covert.cell.kb.ProteinMonomer.empty(0, 0);
        regulatedProteinComplexs = edu.stanford.covert.cell.kb.ProteinComplex.empty(0, 0);
    end

    methods
        function this = Stimuli(knowledgeBase, wid, wholeCellModelID, name,...
                comments, crossReferences)

            if nargin == 0; return; end;

            this = edu.stanford.covert.cell.kb.Stimuli.empty(size(wid, 1), 0);
            if size(wid,1) == 0; return; end;

            this(size(wid, 1), 1) = edu.stanford.covert.cell.kb.Stimuli;
            for i = 1:size(wid, 1)
                this(i, 1).idx = i;
                this(i, 1).knowledgeBase = knowledgeBase;
                this(i, 1).wid = wid(i);
                this(i, 1).wholeCellModelID = wholeCellModelID{i};
                this(i, 1).name = name{i};
                if exist('comments', 'var') && ~isempty(comments); this(i,1).comments = comments{i}; end;
                if exist('crossReferences', 'var')
                    if size(crossReferences, 1) > 1
                        this(i, 1).crossReferences = crossReferences(i);
                    else
                        this(i, 1).crossReferences = struct;
                        fields = fieldnames(crossReferences);
                        for j = 1:size(fields,1)
                            vals = crossReferences.(fields{j});
                            this(i,1).crossReferences.(fields{j}) = vals(i);
                        end
                    end
                end
            end
        end

        function serializeLinks(this)
            for i = 1:numel(this)
                this(i).compartments             = this.serializeLinksHelper(this(i).compartments);
                this(i).values                   = this.serializeLinksHelper(this(i).values);
                this(i).initialTimes             = this.serializeLinksHelper(this(i).initialTimes);
                this(i).finalTimes               = this.serializeLinksHelper(this(i).finalTimes);

                this(i).regulatedProteinMonomers = this.serializeLinksHelper(this(i).regulatedProteinMonomers);
                this(i).regulatedProteinComplexs = this.serializeLinksHelper(this(i).regulatedProteinComplexs);

                serializeLinks@edu.stanford.covert.cell.kb.KnowledgeBaseObject(this(i));
            end
        end
        
        function deserializeLinks(this, kb)
            for i = 1:numel(this)
                this(i).compartments             = this.deserializeLinksHelper(this(i).compartments, kb.compartments);
                this(i).values                   = this.deserializeLinksHelper(this(i).values);
                this(i).initialTimes             = this.deserializeLinksHelper(this(i).initialTimes);
                this(i).finalTimes               = this.deserializeLinksHelper(this(i).finalTimes);

                this(i).regulatedProteinMonomers = this.deserializeLinksHelper(this(i).regulatedProteinMonomers, kb.proteinMonomers);
                this(i).regulatedProteinComplexs = this.deserializeLinksHelper(this(i).regulatedProteinComplexs, kb.proteinComplexs);
                
                deserializeLinks@edu.stanford.covert.cell.kb.KnowledgeBaseObject(this(i), kb);
            end
        end

        function deleteLinks(this)
            for i = 1:numel(this)
                this(i).compartments             = [];
                this(i).values                   = [];
                this(i).initialTimes             = [];
                this(i).finalTimes               = [];

                this(i).regulatedProteinMonomers = [];
                this(i).regulatedProteinComplexs = [];

                deleteLinks@edu.stanford.covert.cell.kb.KnowledgeBaseObject(this(i));
            end
        end
    end
end