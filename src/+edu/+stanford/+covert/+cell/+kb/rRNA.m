% Defines a rRNA polymer
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/7/2009
classdef rRNA < edu.stanford.covert.cell.kb.ssRNA & edu.stanford.covert.cell.kb.Enzyme
    methods
        function this = rRNA(knowledgeBase, wid, wholeCellModelID, name, ...
                sequence, ...
                molecularInteraction, chemicalRegulation, subsystem, ...
                generalClassification, proteaseClassification, ...
                transporterClassification, ...
                comments, crossReferences)

            if nargin == 0; return; end;

            this = edu.stanford.covert.cell.kb.rRNA.empty(size(wid,1),0);
            this(size(wid,1),1) = edu.stanford.covert.cell.kb.ssRNA;
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

                this(i,1).sequence = sequence(i);
                this(i,1).molecularInteraction = molecularInteraction{i};
                this(i,1).chemicalRegulation = chemicalRegulation{i};
                this(i,1).subsystem = subsystem{i};
                this(i,1).generalClassification = generalClassification{i};
                this(i,1).proteaseClassification = proteaseClassification;
                this(i,1).transporterClassification = transporterClassification{i};
            end
        end
        
        function serializeLinks(this)
            for i = 1:numel(this)
                this(i).serializeLinks@edu.stanford.covert.cell.kb.ssRNA();
                this(i).serializeLinks@edu.stanford.covert.cell.kb.Enzyme();
            end
        end
        
        function deserializeLinks(this, kb)
            for i = 1:numel(this)
                deserializeLinks@edu.stanford.covert.cell.kb.ssRNA(this(i), kb);
                deserializeLinks@edu.stanford.covert.cell.kb.Enzyme(this(i), kb);
            end
        end

        function deleteLinks(this)
            for i = 1:numel(this)
                this(i).deleteLinks@edu.stanford.covert.cell.kb.ssRNA();
                this(i).deleteLinks@edu.stanford.covert.cell.kb.Enzyme();
            end
        end
    end
end