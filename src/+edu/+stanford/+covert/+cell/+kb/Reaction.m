% Defines a reaction
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/7/2009
classdef Reaction < edu.stanford.covert.cell.kb.KnowledgeBaseObject
    properties
        process                             = edu.stanford.covert.cell.kb.Process.empty(0,0);
        state                              = edu.stanford.covert.cell.kb.State.empty(0,0);

        stimulis                           = edu.stanford.covert.cell.kb.Stimuli.empty(0,0);
        stimuliCompartments                = edu.stanford.covert.cell.kb.Compartment.empty(0,0);
        stimuliCoefficients                = [];

        metabolites                        = edu.stanford.covert.cell.kb.Metabolite.empty(0,0);
        metaboliteCompartments             = edu.stanford.covert.cell.kb.Compartment.empty(0,0);
        metaboliteCoefficients             = [];

        rnas                               = edu.stanford.covert.cell.kb.TranscriptionUnit.empty(0,0);
        rnaCompartments                    = edu.stanford.covert.cell.kb.Compartment.empty(0,0);
        rnaCoefficients                    = [];

        proteinMonomers                    = edu.stanford.covert.cell.kb.ProteinMonomer.empty(0,0);
        proteinMonomerCompartments         = edu.stanford.covert.cell.kb.Compartment.empty(0,0);
        proteinMonomerCoefficients         = [];

        proteinComplexs                    = edu.stanford.covert.cell.kb.ProteinComplex.empty(0,0);
        proteinComplexCompartments         = edu.stanford.covert.cell.kb.Compartment.empty(0,0);
        proteinComplexCoefficients         = [];

        enzymes                            = edu.stanford.covert.cell.kb.Enzyme.empty(0,0);
        enzymeCompartments                 = edu.stanford.covert.cell.kb.Compartment.empty(0,0);

        coenzymes                          = edu.stanford.covert.cell.kb.Metabolite.empty(0,0);
        coenzymeCompartments               = edu.stanford.covert.cell.kb.Compartment.empty(0,0);
        coenzymeCoefficients               = [];

        stableModifications                = [];
        stableModificationCompartments     = edu.stanford.covert.cell.kb.Compartment.empty(0,0);
        stableModificationPositions        = [];

        parameters                         = edu.stanford.covert.cell.kb.Parameter.empty(0,0);

        pathways                           = edu.stanford.covert.cell.kb.Pathway.empty(0,0);
    end

    properties %(SetAccess = protected)
		type
        ecNumber
        spontaneous
        direction
        deltaG
        keq
        rateLawForward
        kmForward
        vmaxForward
        vmaxUnitsForward
        rateLawBackward
        kmBackward
        vmaxBackward
        vmaxUnitsBackward
        optimalpH
        optimalTemp
        activators
        inhibitors
        lowerBound
        upperBound
        boundUnits
    end

    %computed properties
    properties %(SetAccess = protected)
        enzymeMolecularWeight
    end

    methods
        function this = Reaction(knowledgeBase, wid, wholeCellModelID, name, ...
                type, ecNumber, spontaneous, direction, deltaG, keq, ...
                rateLawForward, kmForward, vmaxForward, vmaxUnitsForward, ...
                rateLawBackward, kmBackward, vmaxBackward, vmaxUnitsBackward,  ...
                optimalpH, optimalTemp, ...
                activators, inhibitors, ...
                lowerBound, upperBound, boundUnits, ...
                comments, crossReferences)

            if nargin == 0; return; end;

            this = edu.stanford.covert.cell.kb.Reaction.empty(size(wid,1),0);
            this(size(wid,1),1) = edu.stanford.covert.cell.kb.Reaction;
            for i=1:size(wid,1)
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
                        for j=1:size(fields,1)
                            values = crossReferences.(fields{j});
                            this(i,1).crossReferences.(fields{j}) = values(i);
                        end
                    end
                end

                this(i,1).type = type{i};
                this(i,1).ecNumber = ecNumber{i};
                this(i,1).spontaneous = spontaneous{i};
                this(i,1).direction = direction{i};
                this(i,1).deltaG = deltaG(i);
                this(i,1).keq = keq(i);
                this(i,1).rateLawForward = rateLawForward{i};
                this(i,1).kmForward = kmForward(i);
                this(i,1).vmaxForward = vmaxForward(i);
                this(i,1).vmaxUnitsForward = vmaxUnitsForward{i};
                this(i,1).rateLawBackward = rateLawBackward{i};
                this(i,1).kmBackward = kmBackward(i);
                this(i,1).vmaxBackward = vmaxBackward(i);
                this(i,1).vmaxUnitsBackward = vmaxUnitsBackward{i};
                this(i,1).optimalpH = optimalpH(i);
                this(i,1).optimalTemp = optimalTemp(i);
                this(i,1).activators = activators{i};
                this(i,1).inhibitors = inhibitors{i};
                this(i,1).lowerBound = lowerBound(i);
                this(i,1).upperBound = upperBound(i);
                this(i,1).boundUnits = boundUnits{i};
            end
        end

        function serializeLinks(this)
            for i = 1:numel(this)
                this(i).process                            = this.serializeLinksHelper(this(i).process);
                this(i).state                              = this.serializeLinksHelper(this(i).state);

                this(i).stimulis                           = this.serializeLinksHelper(this(i).stimulis);
                this(i).stimuliCompartments                = this.serializeLinksHelper(this(i).stimuliCompartments);
                this(i).stimuliCoefficients                = this.serializeLinksHelper(this(i).stimuliCoefficients);

                this(i).metabolites                        = this.serializeLinksHelper(this(i).metabolites);
                this(i).metaboliteCompartments             = this.serializeLinksHelper(this(i).metaboliteCompartments);
                this(i).metaboliteCoefficients             = this.serializeLinksHelper(this(i).metaboliteCoefficients);

                this(i).rnas                               = this.serializeLinksHelper(this(i).rnas);
                this(i).rnaCompartments                    = this.serializeLinksHelper(this(i).rnaCompartments);
                this(i).rnaCoefficients                    = this.serializeLinksHelper(this(i).rnaCoefficients);

                this(i).proteinMonomers                    = this.serializeLinksHelper(this(i).proteinMonomers);
                this(i).proteinMonomerCompartments         = this.serializeLinksHelper(this(i).proteinMonomerCompartments);
                this(i).proteinMonomerCoefficients         = this.serializeLinksHelper(this(i).proteinMonomerCoefficients);

                this(i).proteinComplexs                    = this.serializeLinksHelper(this(i).proteinComplexs);
                this(i).proteinComplexCompartments         = this.serializeLinksHelper(this(i).proteinComplexCompartments);
                this(i).proteinComplexCoefficients         = this.serializeLinksHelper(this(i).proteinComplexCoefficients);

                this(i).enzymes                            = this.serializeLinksHelper(this(i).enzymes);
                this(i).enzymeCompartments                 = this.serializeLinksHelper(this(i).enzymeCompartments);

                this(i).coenzymes                          = this.serializeLinksHelper(this(i).coenzymes);
                this(i).coenzymeCompartments               = this.serializeLinksHelper(this(i).coenzymeCompartments);
                this(i).coenzymeCoefficients               = this.serializeLinksHelper(this(i).coenzymeCoefficients);

                this(i).stableModifications                = this.serializeLinksHelper(this(i).stableModifications);
                this(i).stableModificationCompartments     = this.serializeLinksHelper(this(i).stableModificationCompartments);
                this(i).stableModificationPositions        = this.serializeLinksHelper(this(i).stableModificationPositions);

                this(i).parameters                         = this.serializeLinksHelper(this(i).parameters);

                this(i).pathways                           = this.serializeLinksHelper(this(i).pathways);

                serializeLinks@edu.stanford.covert.cell.kb.KnowledgeBaseObject(this(i));
            end
        end
        
        function deserializeLinks(this, kb)
            for i = 1:numel(this)
                this(i).process                            = this.deserializeLinksHelper(this(i).process, kb.processes);
                this(i).state                              = this.deserializeLinksHelper(this(i).state, kb.states);

                this(i).stimulis                           = this.deserializeLinksHelper(this(i).stimulis, kb.stimulis);
                this(i).stimuliCompartments                = this.deserializeLinksHelper(this(i).stimuliCompartments, kb.compartments);
                this(i).stimuliCoefficients                = this.deserializeLinksHelper(this(i).stimuliCoefficients);

                this(i).metabolites                        = this.deserializeLinksHelper(this(i).metabolites, kb.metabolites);
                this(i).metaboliteCompartments             = this.deserializeLinksHelper(this(i).metaboliteCompartments, kb.compartments);
                this(i).metaboliteCoefficients             = this.deserializeLinksHelper(this(i).metaboliteCoefficients);

                this(i).rnas                               = this.deserializeLinksHelper(this(i).rnas, kb.transcriptionUnits);
                this(i).rnaCompartments                    = this.deserializeLinksHelper(this(i).rnaCompartments, kb.compartments);
                this(i).rnaCoefficients                    = this.deserializeLinksHelper(this(i).rnaCoefficients);

                this(i).proteinMonomers                    = this.deserializeLinksHelper(this(i).proteinMonomers, kb.proteinMonomers);
                this(i).proteinMonomerCompartments         = this.deserializeLinksHelper(this(i).proteinMonomerCompartments, kb.compartments);
                this(i).proteinMonomerCoefficients         = this.deserializeLinksHelper(this(i).proteinMonomerCoefficients);

                this(i).proteinComplexs                    = this.deserializeLinksHelper(this(i).proteinComplexs, kb.proteinComplexs);
                this(i).proteinComplexCompartments         = this.deserializeLinksHelper(this(i).proteinComplexCompartments, kb.compartments);
                this(i).proteinComplexCoefficients         = this.deserializeLinksHelper(this(i).proteinComplexCoefficients);

                this(i).enzymes                            = this.deserializeLinksHelper(this(i).enzymes, {kb.proteinMonomers; kb.proteinComplexs});
                this(i).enzymeCompartments                 = this.deserializeLinksHelper(this(i).enzymeCompartments, kb.compartments);

                this(i).coenzymes                          = this.deserializeLinksHelper(this(i).coenzymes, kb.metabolites);
                this(i).coenzymeCompartments               = this.deserializeLinksHelper(this(i).coenzymeCompartments, kb.compartments);
                this(i).coenzymeCoefficients               = this.deserializeLinksHelper(this(i).coenzymeCoefficients);

                this(i).stableModifications                = this.deserializeLinksHelper(this(i).stableModifications, {kb.genes; kb.proteinMonomers});
                this(i).stableModificationCompartments     = this.deserializeLinksHelper(this(i).stableModificationCompartments, kb.compartments);
                this(i).stableModificationPositions        = this.deserializeLinksHelper(this(i).stableModificationPositions);

                this(i).parameters                         = this.deserializeLinksHelper(this(i).parameters, kb.parameters);

                this(i).pathways                           = this.deserializeLinksHelper(this(i).pathways, kb.pathways);
                
                deserializeLinks@edu.stanford.covert.cell.kb.KnowledgeBaseObject(this(i), kb);
            end
        end

        function deleteLinks(this)
            for i = 1:numel(this)
                this(i).process                             = [];
                this(i).state                              = [];

                this(i).stimulis                           = [];
                this(i).stimuliCompartments                = [];
                this(i).stimuliCoefficients                = [];

                this(i).metabolites                        = [];
                this(i).metaboliteCompartments             = [];
                this(i).metaboliteCoefficients             = [];

                this(i).rnas                               = [];
                this(i).rnaCompartments                    = [];
                this(i).rnaCoefficients                    = [];

                this(i).proteinMonomers                    = [];
                this(i).proteinMonomerCompartments         = [];
                this(i).proteinMonomerCoefficients         = [];

                this(i).proteinComplexs                    = [];
                this(i).proteinComplexCompartments         = [];
                this(i).proteinComplexCoefficients         = [];

                this(i).enzymes                            = [];
                this(i).enzymeCompartments                 = [];

                this(i).coenzymes                          = [];
                this(i).coenzymeCompartments               = [];
                this(i).coenzymeCoefficients               = [];

                this(i).stableModifications                = [];
                this(i).stableModificationCompartments     = [];
                this(i).stableModificationPositions        = [];

                this(i).parameters                         = [];

                this(i).pathways                           = [];

                deleteLinks@edu.stanford.covert.cell.kb.KnowledgeBaseObject(this(i));
            end
        end

        function value = get.enzymeMolecularWeight(this)
            %retrieve
            if ~isempty(this.enzymeMolecularWeight)
                value = this.enzymeMolecularWeight;
                return;
            end
            
            %compute
            value = this.enzymes.molecularWeight;
            
            %store
            this.enzymeMolecularWeight = value;
        end
    end
end