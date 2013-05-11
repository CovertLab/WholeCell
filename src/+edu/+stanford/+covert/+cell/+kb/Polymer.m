% Defines a polymer (dna, rna, protein, etc.). Base class for
% - NucleicAcid
%   - dsDNA
%   - dsRNA
%   - ssDNA
%   - ssRNA
% - ProteinMonomer
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/7/2009
classdef Polymer < edu.stanford.covert.cell.kb.PhysicalObject
    properties %(SetAccess = protected)
        sequenceTopology = 'linear';
    end
    
    properties (Abstract = true) %SetAccess = protected
        sequence
    end
    
    %computed properties
    properties (Abstract = true) %SetAccess = protected
        baseCount
        cumulativeBaseCount
        decayReaction
    end
    
    properties %(SetAccess = protected)
        sequenceLength
    end
    
    methods
        function this = Polymer(knowledgeBase, wid, wholeCellModelID, name, ...
                comments, crossReferences)
            
            if nargin == 0; return; end;
            
            this = edu.stanford.covert.cell.kb.Polymer.empty(size(wid,1),0);
            this(size(wid,1),1) = edu.stanford.covert.cell.kb.Polymer;
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
            end
        end
        
        function value = get.sequenceLength(this)
            %retrieve
            if ~isempty(this.sequenceLength)
                value = this.sequenceLength;
                return;
            end
            
            %compute
            value = length(this.sequence);
            
            %store
            this.sequenceLength = value;
        end
    end
end