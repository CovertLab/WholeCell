% Defines a reference
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/7/2009
classdef Reference < edu.stanford.covert.cell.kb.KnowledgeBaseObject
    properties
        type
        pmid
        isbn
        authors
        editors
        year
        title
        publication
        volume
        issue
        pages
        publisher
        url
        citations

        processes          = edu.stanford.covert.cell.kb.Process.empty(0,0);
        parameters         = edu.stanford.covert.cell.kb.Parameter.empty(0,0);
        compartments       = edu.stanford.covert.cell.kb.Compartment.empty(0,0);
        metabolites        = edu.stanford.covert.cell.kb.Metabolite.empty(0,0);
        genes              = edu.stanford.covert.cell.kb.Gene.empty(0,0);
        transcriptionUnits = edu.stanford.covert.cell.kb.ProteinComplex.empty(0,0);
        genomeFeatures     = edu.stanford.covert.cell.kb.GenomeFeature.empty(0,0);
        proteinMonomers    = edu.stanford.covert.cell.kb.TranscriptionUnit.empty(0,0);
        proteinComplexs    = edu.stanford.covert.cell.kb.ProteinMonomer.empty(0,0);
        reactions          = edu.stanford.covert.cell.kb.Reaction.empty(0,0);
        pathways           = edu.stanford.covert.cell.kb.Pathway.empty(0,0);
        stimulis           = edu.stanford.covert.cell.kb.Stimuli.empty(0,0);
        notes              = edu.stanford.covert.cell.kb.Note.empty(0,0);
    end
   
    methods
        function this = Reference(knowledgeBase, wid, wholeCellModelID, name,...
                type, pmid, isbn, ...
                authors, editors, year, title, publication, volume, issue, pages, ...
                publisher, url, citations, ...
                comments, crossReferences)

            if nargin == 0; return; end;

            this = edu.stanford.covert.cell.kb.Reference.empty(size(wid,1),0);
            this(size(wid,1),1) = edu.stanford.covert.cell.kb.Reference;
            for i = 1:size(wid,1)
                this(i,1).idx = i;
                this(i,1).knowledgeBase = knowledgeBase;
                this(i,1).wid = wid(i);
                this(i,1).wholeCellModelID = wholeCellModelID{i};
                if ~isempty(name); this(i,1).name = name{i}; end;
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

				this(i,1).type = type{i};
                this(i,1).pmid = pmid(i);
                this(i,1).isbn = isbn(i);
                this(i,1).authors = authors{i};
                this(i,1).editors = editors{i};
                this(i,1).year = year(i);
                this(i,1).title = title{i};
                this(i,1).publication = publication{i};
                this(i,1).volume = volume(i);
                this(i,1).issue = issue(i);
                this(i,1).pages = pages{i};
                this(i,1).publisher = publisher{i};
                this(i,1).url = url{i};
                this(i,1).citations = citations{i};
            end
        end

        function serializeLinks(this)
            for i = 1:numel(this)
                this(i).processes          = this.serializeLinksHelper(this(i).processes);
                this(i).parameters         = this.serializeLinksHelper(this(i).parameters);
                this(i).compartments       = this.serializeLinksHelper(this(i).compartments);
                this(i).metabolites        = this.serializeLinksHelper(this(i).metabolites);
                this(i).genes              = this.serializeLinksHelper(this(i).genes);
                this(i).transcriptionUnits = this.serializeLinksHelper(this(i).transcriptionUnits);
                this(i).genomeFeatures     = this.serializeLinksHelper(this(i).genomeFeatures);
                this(i).proteinMonomers    = this.serializeLinksHelper(this(i).proteinMonomers);
                this(i).proteinComplexs    = this.serializeLinksHelper(this(i).proteinComplexs);
                this(i).reactions          = this.serializeLinksHelper(this(i).reactions);
                this(i).pathways           = this.serializeLinksHelper(this(i).pathways);
                this(i).stimulis           = this.serializeLinksHelper(this(i).stimulis);
                this(i).notes              = this.serializeLinksHelper(this(i).notes);

                serializeLinks@edu.stanford.covert.cell.kb.KnowledgeBaseObject(this(i));
            end
        end
        
        function deserializeLinks(this, kb)
            for i = 1:numel(this)
                this(i).processes          = this.deserializeLinksHelper(this(i).processes, kb.processes);
                this(i).parameters         = this.deserializeLinksHelper(this(i).parameters, kb.parameters);
                this(i).compartments       = this.deserializeLinksHelper(this(i).compartments, kb.compartments);
                this(i).metabolites        = this.deserializeLinksHelper(this(i).metabolites, kb.metabolites);
                this(i).genes              = this.deserializeLinksHelper(this(i).genes, kb.genes);
                this(i).transcriptionUnits = this.deserializeLinksHelper(this(i).transcriptionUnits, kb.transcriptionUnits);
                this(i).genomeFeatures     = this.deserializeLinksHelper(this(i).genomeFeatures, kb.genomeFeatures);
                this(i).proteinMonomers    = this.deserializeLinksHelper(this(i).proteinMonomers, kb.proteinMonomers);
                this(i).proteinComplexs    = this.deserializeLinksHelper(this(i).proteinComplexs, kb.proteinComplexs);
                this(i).reactions          = this.deserializeLinksHelper(this(i).reactions, kb.reactions);
                this(i).pathways           = this.deserializeLinksHelper(this(i).pathways, kb.pathways);
                this(i).stimulis           = this.deserializeLinksHelper(this(i).stimulis, kb.stimulis);
                this(i).notes              = this.deserializeLinksHelper(this(i).notes, kb.notes);
                
                deserializeLinks@edu.stanford.covert.cell.kb.KnowledgeBaseObject(this(i), kb);
            end
        end

        function deleteLinks(this)
            for i = 1:numel(this)
                this(i).processes            = [];
                this(i).parameters         = [];
                this(i).compartments       = [];
                this(i).metabolites        = [];
                this(i).genes              = [];
                this(i).transcriptionUnits = [];
                this(i).genomeFeatures     = [];
                this(i).proteinMonomers    = [];
                this(i).proteinComplexs    = [];
                this(i).reactions          = [];
                this(i).pathways           = [];
                this(i).stimulis           = [];
                this(i).notes              = [];

                deleteLinks@edu.stanford.covert.cell.kb.KnowledgeBaseObject(this(i));
            end
        end
    end
end