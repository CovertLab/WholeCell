% Base class for knowledge base objects
% - KnowledgeBase
% - PhysicalObject
% - Enzyme
% - Reaction
% - Pathway
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 11/26/2012
classdef KnowledgeBaseObject < handle
    properties
        idx = uint32(0);                 %index in array

        wid = uint32(0);                 %unique id in database
        wholeCellModelID                 %unique id in model
		insertDate
		modifiedDate
		insertUser
		modifiedUser
        name
        comments
        crossReferences = struct;

        knowledgeBase = edu.stanford.covert.cell.kb.KnowledgeBase.empty(0, 0);
        references    = edu.stanford.covert.cell.kb.Reference.empty(0, 0);
    end

    methods
        function this = KnowledgeBaseObject(knowledgeBase, wid, wholeCellModelID, name, comments, crossReferences)
            if nargin == 0; return; end;
            this.knowledgeBase = knowledgeBase;
            this.wid = wid;
            this.wholeCellModelID = wholeCellModelID;

            if nargin < 4; return; end;
            this.name = name;
            this.comments = comments;
            this.crossReferences = crossReferences;
        end
        
        function calcIndices(this)
            for i = 1:numel(this)
                this(i).idx = i;
            end
        end

        function serializeLinks(this)
            for i = 1:numel(this)
                this(i).knowledgeBase = this.serializeLinksHelper(this(i).knowledgeBase);
                this(i).references    = this.serializeLinksHelper(this(i).references);
            end
        end
        
        function deserializeLinks(this, kb)
            for i = 1:numel(this)
                this(i).knowledgeBase = this.deserializeLinksHelper(this(i).knowledgeBase, kb);
                this(i).references = this.deserializeLinksHelper(this(i).references, kb.references);
            end
        end
        
        function classAndIndexs = serializeLinksHelper(~, vals)
            if isa(vals, 'edu.stanford.covert.cell.kb.KnowledgeBaseObject')
                classAndIndexs = {class(vals) [vals.idx]};
            else
                classAndIndexs = vals;
            end
        end
        
        function vals = deserializeLinksHelper(~, classAndIndexs, all_objs)
            if nargin >= 3 && iscell(classAndIndexs) && isequal(size(classAndIndexs), [1 2]) && isequal(classAndIndexs{1}(1:min(28, end)), 'edu.stanford.covert.cell.kb.')
                if ~iscell(all_objs)
                    all_objs = {all_objs};
                end
                
                valid = false;
                for i = 1:numel(all_objs)
                    if isequal(class(all_objs{i}), classAndIndexs{1})
                        vals = all_objs{i}(classAndIndexs{2});
                        valid = true;
                        break;
                    end
                end
                
                if ~valid
                    throw(MException('KnowledgeBaseObject:deserializeLinksHelper:invalidClass', 'invalid class %s', classAndIndexs{1}));
                end
            else
                vals = classAndIndexs;
            end
        end

        function deleteLinks(this)
            for i = 1:length(this)
                this(i).knowledgeBase = [];
                this(i).references    = [];
            end
        end

        function set.idx(this, value)
            this.idx = uint32(value);
        end

        function set.wid(this, wid)
            this.wid = uint32(wid);
        end
    end
    
    methods
        function invalidate(this)
            metaClass = metaclass(this);
            for i = 1:numel(metaClass.Properties)
                if ~isempty(metaClass.Properties{i}.GetMethod) && isempty(metaClass.Properties{i}.SetMethod) && isequal(metaClass.Properties{i}.SetAccess, 'public')
                    for j = 1:numel(this)
                        this(j).(metaClass.Properties{i}.Name) = [];
                    end
                end
            end
        end
    end
end