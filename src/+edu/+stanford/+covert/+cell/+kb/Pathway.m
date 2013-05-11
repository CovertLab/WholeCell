% Defines a pathway
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/7/2009
classdef Pathway < edu.stanford.covert.cell.kb.KnowledgeBaseObject
    properties
        reactions = edu.stanford.covert.cell.kb.Reaction.empty(0, 0);
    end

    methods
        function this = Pathway(knowledgeBase, wid, wholeCellModelID, name, ...
                comments, crossReferences)

            if nargin == 0; return; end;

            this = edu.stanford.covert.cell.kb.Pathway.empty(size(wid, 1), 0);
            this(size(wid, 1), 1) = edu.stanford.covert.cell.kb.Pathway;
            for i = 1:size(wid, 1)
                this(i, 1).idx = i;
                this(i, 1).knowledgeBase = knowledgeBase;
                this(i, 1).wid = wid(i);
                this(i, 1).wholeCellModelID = wholeCellModelID{i};
                this(i, 1).name = name{i};
                if exist('comments', 'var') && ~isempty(comments); this(i, 1).comments = comments{i}; end;
                if exist('crossReferences', 'var')
                    if size(crossReferences, 1) > 1
                        this(i,1).crossReferences = crossReferences(i);
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
                this(i).reactions = this.serializeLinksHelper(this(i).reactions);

                serializeLinks@edu.stanford.covert.cell.kb.KnowledgeBaseObject(this(i));
            end
        end
        
        function deserializeLinks(this, kb)
            for i = 1:numel(this)
                this(i).reactions = this.deserializeLinksHelper(this(i).reactions, kb.reactions);
                
                deserializeLinks@edu.stanford.covert.cell.kb.KnowledgeBaseObject(this(i), kb);
            end
        end

        function deleteLinks(this)
            for i = 1:numel(this)
                this(i).reactions = [];

                deleteLinks@edu.stanford.covert.cell.kb.KnowledgeBaseObject(this(i));
            end
        end
    end
end