% Defines a genome
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/7/2009
classdef Genome < edu.stanford.covert.cell.kb.dsDNA

    properties
        genes              = edu.stanford.covert.cell.kb.Gene.empty(0, 0);
        transcriptionUnits = edu.stanford.covert.cell.kb.TranscriptionUnit.empty(0, 0);
        features           = edu.stanford.covert.cell.kb.GenomeFeature.empty(0, 0);
    end

    properties %(SetAccess = protected)
        halfLife
        sequence
    end

    methods
        %constructor
        function this = Genome(knowledgeBase, wid, wholeCellModelID, name, ...
                sequenceTopology, sequence, ...
                comments, crossReferences)

            if nargin == 0; return; end;

            this = edu.stanford.covert.cell.kb.Genome.empty(size(wid, 1), 0);
            this(size(wid, 1), 1) = edu.stanford.covert.cell.kb.Genome;
            for i = 1:size(wid, 1)
                this(i, 1).idx = i;
                this(i, 1).knowledgeBase = knowledgeBase;
                this(i, 1).wid = wid(i);
                if ~isempty(wholeCellModelID); this(i, 1).wholeCellModelID = wholeCellModelID{i}; end;
                if ~isempty(name); this(i, 1).name = name{i}; end;
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

                this(i, 1).sequenceTopology = sequenceTopology{i};
                this(i, 1).sequence = upper(sequence{i});
            end
        end

        function serializeLinks(this)
            for i = 1:numel(this)
                this(i).genes              = this.serializeLinksHelper(this(i).genes);
                this(i).transcriptionUnits = this.serializeLinksHelper(this(i).transcriptionUnits);
                this(i).features           = this.serializeLinksHelper(this(i).features);

                serializeLinks@edu.stanford.covert.cell.kb.dsDNA(this(i));
            end
        end
        
        function deserializeLinks(this, kb)
            for i = 1:numel(this)
                this(i).genes              = this.deserializeLinksHelper(this(i).genes, kb.genes);
                this(i).transcriptionUnits = this.deserializeLinksHelper(this(i).transcriptionUnits, kb.transcriptionUnits);
                this(i).features           = this.deserializeLinksHelper(this(i).features, kb.genomeFeatures);
                
                deserializeLinks@edu.stanford.covert.cell.kb.dsDNA(this(i), kb);
            end
        end

        function deleteLinks(this)
            for i = 1:numel(this)
                this(i).genes              = [];
                this(i).transcriptionUnits = [];
                this(i).features           = [];

                deleteLinks@edu.stanford.covert.cell.kb.dsDNA(this(i));
            end
        end

        function value = get.halfLife(this)
            throw(MException('Genome:error', 'property is not defined'));
        end
    end
end