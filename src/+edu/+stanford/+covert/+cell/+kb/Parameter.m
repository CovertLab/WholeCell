% Defines a parameter
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 11/17/2009
classdef Parameter < edu.stanford.covert.cell.kb.KnowledgeBaseObject
    properties
        process          = edu.stanford.covert.cell.kb.Process.empty(0, 0);
        state           = edu.stanford.covert.cell.kb.State.empty(0, 0);
        reactions       = edu.stanford.covert.cell.kb.Reaction.empty(0, 0);
        proteinMonomers = edu.stanford.covert.cell.kb.ProteinMonomer.empty(0, 0);
        proteinComplexs = edu.stanford.covert.cell.kb.ProteinComplex.empty(0, 0);
    end

    properties %(SetAccess = protected)
        index
        defaultValue
        units
        experimentallyConstrained
    end

    methods
        function this = Parameter(knowledgeBase, wid, wholeCellModelID, name, ...
                index, defaultValue, units, experimentallyConstrained, ...
                comments, crossReferences)

            if nargin == 0; return; end;

            this = edu.stanford.covert.cell.kb.Parameter.empty(size(wid, 1), 0);
            this(size(wid, 1), 1) = edu.stanford.covert.cell.kb.Parameter;
            for i = 1:size(wid, 1)
                this(i, 1).idx = i;
                this(i, 1).knowledgeBase = knowledgeBase;
                this(i, 1).wid = wid(i);
                this(i, 1).wholeCellModelID = wholeCellModelID{i};
                this(i, 1).name = name{i};
                this(i, 1).index = index{i};
                this(i, 1).defaultValue = this.parseDefaultValue(defaultValue{i});
                this(i, 1).units = units{i};
                this(i, 1).experimentallyConstrained = experimentallyConstrained(i);
                if exist('comments','var') && ~isempty(comments); this(i,1).comments = comments{i}; end;
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
                this(i).process         = this.serializeLinksHelper(this(i).process);
                this(i).state           = this.serializeLinksHelper(this(i).state);
                this(i).reactions       = this.serializeLinksHelper(this(i).reactions);
                this(i).proteinMonomers = this.serializeLinksHelper(this(i).proteinMonomers);
                this(i).proteinComplexs = this.serializeLinksHelper(this(i).proteinComplexs);
                serializeLinks@edu.stanford.covert.cell.kb.KnowledgeBaseObject(this(i))
            end
        end
        
        function deserializeLinks(this, kb)
            for i = 1:numel(this)
                this(i).process         = this.deserializeLinksHelper(this(i).process, kb.processes);
                this(i).state           = this.deserializeLinksHelper(this(i).state, kb.states);
                this(i).reactions       = this.deserializeLinksHelper(this(i).reactions, kb.reactions);
                this(i).proteinMonomers = this.deserializeLinksHelper(this(i).proteinMonomers, kb.proteinMonomers);
                this(i).proteinComplexs = this.deserializeLinksHelper(this(i).proteinComplexs, kb.proteinComplexs);                
                deserializeLinks@edu.stanford.covert.cell.kb.KnowledgeBaseObject(this(i), kb);
            end
        end

        function deleteLinks(this)
            for i = 1:numel(this)
                this(i).process         = [];
                this(i).state           = [];
                this(i).reactions       = [];
                this(i).proteinMonomers = [];
                this(i).proteinComplexs = [];
                deleteLinks@edu.stanford.covert.cell.kb.KnowledgeBaseObject(this(i))
            end
        end
    end

    methods (Static = true)
        function parsedValue = parseDefaultValue(rawValue)
            if isempty(rawValue)
                parsedValue = {};
                return;
            end

            if iscell(rawValue) && numel(rawValue) == 1
                rawValue = rawValue{1};
            end

            rawValue = strsplit(';', rawValue)';
            parsedValue = zeros(size(rawValue));
            for i = 1:length(rawValue)
                parsedValue(i) = str2double(rawValue{i});
            end
            if any(isnan(parsedValue))
                parsedValue = rawValue;
            end
        end
    end
end