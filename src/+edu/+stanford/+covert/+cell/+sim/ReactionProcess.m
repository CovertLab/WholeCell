% Builds on base process class. Provides additional support for processes
% represented as reactions in database.
% - reaction Stoichiometry Matrix
% - reaction Catalysis Matrix
% - reaction Modification Matrix
% - enzyme bounds
% - reaction bounds
%
% Used to implement the following processes
% - metabolism (flux-balance analysis)
% - RNA, protein modification
% - terminal organelle assembly
% - tRNA aminoacylation
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/5/2010
classdef ReactionProcess < edu.stanford.covert.cell.sim.Process
    %Whole cell model IDs of the process's stimuli, subsrates, enzymes and
    %reactions. Set by initializeConstants from the reactions associated
    %with the process in the knowledge base
    properties
        stimuliWholeCellModelIDs   = {};  %stimuli Whole Cell model IDs
        substrateWholeCellModelIDs = {};  %substrate Whole Cell model IDs
        enzymeWholeCellModelIDs    = {};  %enzyme Whole Cell model IDs
        reactionWholeCellModelIDs  = {};  %reaction Whole Cell model IDs
    end

    %constants set by initialize constants method from the knowledge base.
    %these constants describe the reactions implemented by the process
    properties        
        reactionNames               %reaction names
        reactionTypes               %reaction types
        reactionStoichiometryMatrix %stoichiometries of metabolic reactions (substrates X reactions X compartments)
        reactionCatalysisMatrix     %reactions and enzyme (reactions X monomers and complexes X compartments)
        reactionModificationMatrix  %reactions X (RNAs, protein monomers) X compartments
        reactionCoenzymeMatrix      %reactions X metabolites X compartments
        enzymeBounds                %maximal flux bounds per enzyme for enzyme catalyzed reactions where kinetic data is available (reactions/enzyme/s)
        reactionBounds              %maximal flux bounds for exchange reactions, (reactions/(gram dry biomass)/s)
    end
    
    %constructor
    methods
        %sets process meta data
        function this = ReactionProcess(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.Process(wholeCellModelID, name);
        end
    end

    %communication between process/simulation
    methods
        %Ues knowledge base to initialize process's constants including:
        %- whole cell model IDs, names, and types of reactions associated with process
        %- whole cell model IDs of stimuli, substrates, and enzymes
        %  associated with reactions
        %- mappings between process's reactions, stimuli, substrates, and
        %  enzymes and the simulation's state
        %- adjacency matrices describing
        %  - substrate stoichiometry of reactions
        %  - substrate coenzymes of reactions
        %  - substrates modified by each reaction
        %  - enzymes which catalyze each reaction
        %- quantitative properties of reactions:
        %  - reaction bounds
        %  - enzyme kinetic rates
        function initializeConstants(this, knowledgeBase, simulation, options)
            %options
            if ~exist('options','var') || ~isstruct(options)
                options = struct(...
                    'retainSubstrateCompartments', false, ...
                    'retainEnzymeCompartments', false, ...
                    'retainModificationCompartments', false);
            else
                if ~isfield(options, 'retainSubstrateCompartments')
                    options.retainSubstrateCompartments = false;
                end
                if ~isfield(options, 'retainEnzymeCompartments')
                    options.retainEnzymeCompartments = false;
                end
                if ~isfield(options, 'retainModificationCompartments')
                    options.retainModificationCompartments = false;
                end
            end

            %reactions
            process = findobj(knowledgeBase.processes, 'wholeCellModelID', this.wholeCellModelID);
            reactions = process.reactions;
            reactionGlobalIndexs = [reactions.idx]';
            
            this.reactionWholeCellModelIDs        = {reactions.wholeCellModelID}';            
            this.reactionNames                    = {reactions.name}';
            this.reactionTypes                    = {reactions.type}';
            
            this.reactionCatalysisMatrix          = knowledgeBase.reactionCatalysisMatrix(reactionGlobalIndexs, :, :);
            this.reactionModificationMatrix       = knowledgeBase.reactionModificationMatrix(reactionGlobalIndexs, :, :);
            this.reactionCoenzymeMatrix           = knowledgeBase.reactionCoenzymeMatrix(reactionGlobalIndexs, :, :);
            this.enzymeBounds                     = knowledgeBase.enzymeBounds(reactionGlobalIndexs, :);
            this.reactionBounds                   = knowledgeBase.reactionBounds(reactionGlobalIndexs, :);

            reactionStimuliStoichiometryMatrix    = knowledgeBase.reactionStimuliStoichiometryMatrix(:, reactionGlobalIndexs, :);
            reactionMetaboliteStoichiometryMatrix = knowledgeBase.reactionMetaboliteStoichiometryMatrix(:, reactionGlobalIndexs, :);
            reactionRNAStoichiometryMatrix        = knowledgeBase.reactionRNAStoichiometryMatrix(:, reactionGlobalIndexs, :);
            reactionMonomerStoichiometryMatrix    = knowledgeBase.reactionProteinMonomerStoichiometryMatrix(:, reactionGlobalIndexs, :);
            reactionComplexStoichiometryMatrix    = knowledgeBase.reactionProteinComplexStoichiometryMatrix(:, reactionGlobalIndexs, :);

            %substrates
            stimuliIndexs    = find(sum(sum(abs(reactionStimuliStoichiometryMatrix), 3), 2));
            metaboliteIndexs = find(sum(sum(abs(reactionMetaboliteStoichiometryMatrix), 3), 2) | sum(sum(abs(this.reactionCoenzymeMatrix), 3), 1)');
            rnaIndexs        = find(sum(sum(abs(reactionRNAStoichiometryMatrix), 3), 2));
            monomerIndexs    = find(sum(sum(abs(reactionMonomerStoichiometryMatrix), 3), 2));
            complexIndexs    = find(sum(sum(abs(reactionComplexStoichiometryMatrix), 3), 2));

            this.substrateWholeCellModelIDs = [
                this.stimulus.wholeCellModelIDs(stimuliIndexs);
                this.metabolite.wholeCellModelIDs(metaboliteIndexs);
                this.rna.wholeCellModelIDs(this.rna.matureIndexs(rnaIndexs));
                this.monomer.wholeCellModelIDs(this.monomer.matureIndexs(monomerIndexs));
                this.complex.wholeCellModelIDs(this.complex.matureIndexs(complexIndexs))];
            this.reactionStoichiometryMatrix = [
                reactionStimuliStoichiometryMatrix(stimuliIndexs, :, :);
                reactionMetaboliteStoichiometryMatrix(metaboliteIndexs, :, :);
                reactionRNAStoichiometryMatrix(rnaIndexs, :, :);
                reactionMonomerStoichiometryMatrix(monomerIndexs, :, :);
                reactionComplexStoichiometryMatrix(complexIndexs, :, :)];
            this.reactionCoenzymeMatrix = [
                zeros(size(this.reactionCoenzymeMatrix, 1), length(stimuliIndexs),         size(this.reactionCoenzymeMatrix, 3)) ...
                this.reactionCoenzymeMatrix(:, metaboliteIndexs, :)  ...
                zeros(size(this.reactionCoenzymeMatrix, 1), length(rnaIndexs),             size(this.reactionCoenzymeMatrix, 3)) ...
                zeros(size(this.reactionCoenzymeMatrix, 1), length(monomerIndexs),         size(this.reactionCoenzymeMatrix, 3)) ...
                zeros(size(this.reactionCoenzymeMatrix, 1), length(complexIndexs),         size(this.reactionCoenzymeMatrix, 3))];

            %enzymes
            enzymeIndexs  = find(sum(sum(this.reactionCatalysisMatrix, 3), 1))';
            this.reactionCatalysisMatrix = this.reactionCatalysisMatrix(:, enzymeIndexs, :);
            enzymeRNAIndexs     = [];
            enzymeMonomerIndexs = enzymeIndexs(enzymeIndexs <= length(this.monomer.matureIndexs));
            enzymeComplexIndexs = enzymeIndexs(enzymeIndexs >  length(this.monomer.matureIndexs)) - length(this.monomer.matureIndexs);

            this.enzymeWholeCellModelIDs = { ...
                this.rna.wholeCellModelIDs{this.rna.matureIndexs(enzymeRNAIndexs)} ...
                this.monomer.wholeCellModelIDs{this.monomer.matureIndexs(enzymeMonomerIndexs)} ...
                this.complex.wholeCellModelIDs{this.complex.matureIndexs(enzymeComplexIndexs)}}';

            %super class method
            initializeConstants@edu.stanford.covert.cell.sim.Process(this, knowledgeBase, simulation, options);

            %substrates
            if ~options.retainSubstrateCompartments
                numSubstrates = size(this.reactionStoichiometryMatrix, 1);
                numReactions  = size(this.reactionStoichiometryMatrix, 2);

                substrateCompartmentIndexs = [
                    this.substrateStimulusCompartmentIndexs;
                    this.substrateMetaboliteCompartmentIndexs;
                    this.substrateRNACompartmentIndexs;
                    this.substrateMonomerCompartmentIndexs;
                    this.substrateComplexCompartmentIndexs];

                this.reactionStoichiometryMatrix =  this.reactionStoichiometryMatrix(sub2ind(...
                    size(this.reactionStoichiometryMatrix),...
                    repmat((1:numSubstrates)', [1 numReactions]),...
                    repmat(1:numReactions, [numSubstrates 1]),...
                    repmat(substrateCompartmentIndexs, [1 numReactions])));

                this.reactionCoenzymeMatrix =  this.reactionCoenzymeMatrix(sub2ind(...
                    size(this.reactionCoenzymeMatrix),...
                    repmat((1:numReactions)', [1 numSubstrates]),...
                    repmat(1:numSubstrates, [numReactions 1]),...
                    repmat(substrateCompartmentIndexs', [numReactions 1])));
            end

            %enzymes
            if ~options.retainEnzymeCompartments
                numReactions = size(this.reactionCatalysisMatrix, 1);
                numEnzymes   = size(this.reactionCatalysisMatrix, 2);
                enzymeCompartments = zeros(1, numEnzymes);
                enzymeCompartments(this.enzymeMonomerLocalIndexs) = this.enzymeMonomerCompartmentIndexs;
                enzymeCompartments(this.enzymeComplexLocalIndexs) = this.enzymeComplexCompartmentIndexs;
                this.reactionCatalysisMatrix = this.reactionCatalysisMatrix(sub2ind(...
                    size(this.reactionCatalysisMatrix),...
                    repmat((1:numReactions)', [1 numEnzymes]),...
                    repmat(1:numEnzymes, [numReactions 1]),...
                    repmat(enzymeCompartments, [numReactions 1])));
            end

            %modifications
            if ~options.retainModificationCompartments
                numReactions = size(this.reactionModificationMatrix, 1);
                numGenes = size(this.reactionModificationMatrix, 2);
                geneCompartments = zeros(1, numGenes);
                geneCompartments(this.gene.mRNAIndexs) = this.monomer.compartments(this.monomer.matureIndexs);
                geneCompartments(setdiff(1:numGenes, this.gene.mRNAIndexs)) = this.compartment.cytosolIndexs;

                this.reactionModificationMatrix = this.reactionModificationMatrix(sub2ind(...
                    size(this.reactionModificationMatrix),...
                    repmat((1:numReactions)', [1 numGenes]),...
                    repmat(1:numGenes, [numReactions 1]),...
                    repmat(geneCompartments, [numReactions 1])));
            end
        end
    end

    %helper methods of initialize constants
    methods
        %Computes indices of process's reactions within simulation
        function value = reactionIndexs(this, wholeCellModelIDs)
            value = this.componentIndexs(wholeCellModelIDs, 'reaction');
        end
    end

    %get/set methods of annotation properties
    methods
        %Annotates as fixed constants:
        function value = computeFixedConstantsNames(this)
            value = [this.computeFixedConstantsNames@edu.stanford.covert.cell.sim.Process();
                'reactionNames';
                'reactionTypes';
                'reactionStoichiometryMatrix';
                'reactionCatalysisMatrix';
                'reactionModificationMatrix';
                'reactionCoenzymeMatrix';
                'enzymeBounds';
                'reactionBounds';
                ];
        end
    end

    methods (Access = protected)
        function initializeConstants_overrideReactions(this, simulation, wholeCellModelIDs)
            %map old substrates onto new
            [tfs, idxs] = ismember(wholeCellModelIDs, this.reactionWholeCellModelIDs);
            assert(all(tfs));

            %update properties
            this.reactionWholeCellModelIDs   = this.reactionWholeCellModelIDs(idxs);
            this.reactionNames               = this.reactionNames(idxs);
            this.reactionTypes               = this.reactionTypes(idxs);
            this.reactionStoichiometryMatrix = this.reactionStoichiometryMatrix(:, idxs);
            this.reactionCatalysisMatrix     = this.reactionCatalysisMatrix(idxs, :);
            this.reactionModificationMatrix  = this.reactionModificationMatrix(idxs, :);
            this.reactionCoenzymeMatrix      = this.reactionCoenzymeMatrix(idxs, :);
            this.reactionBounds              = this.reactionBounds(idxs, :);
            this.enzymeBounds                = this.enzymeBounds(idxs, :);

            %trim substrates
            idxs = any(this.reactionStoichiometryMatrix, 2) | any(this.reactionCoenzymeMatrix, 1)';
            this.initializeConstants_overrideSubstrates(this.substrateWholeCellModelIDs(idxs));

            %trim enzymes
            idxs = any(this.reactionCatalysisMatrix, 1);
            this.initializeConstants_overrideEnzymes(this.enzymeWholeCellModelIDs(idxs));
        end

        function initializeConstants_overrideSubstrates(this, wholeCellModelIDs)           
            %map old substrates onto new
            [~, idxs1, idxs2] = intersect(this.substrateWholeCellModelIDs, wholeCellModelIDs);

            %map old reactionStoichiometryMatrix, reactionCoenzymeMatrix
            %onto new
            reactionStoichiometryMatrix = zeros(length(wholeCellModelIDs), size(this.reactionStoichiometryMatrix, 2), size(this.reactionStoichiometryMatrix, 3));
            reactionCoenzymeMatrix      = zeros(size(this.reactionCoenzymeMatrix, 1), length(wholeCellModelIDs), size(this.reactionCoenzymeMatrix, 3));

            reactionStoichiometryMatrix(idxs2, :, :) = this.reactionStoichiometryMatrix(idxs1, :, :);
            reactionCoenzymeMatrix(:, idxs2, :) = this.reactionCoenzymeMatrix(:, idxs1, :);

            this.reactionStoichiometryMatrix = reactionStoichiometryMatrix;
            this.reactionCoenzymeMatrix      = reactionCoenzymeMatrix;

            %super class method
            this.initializeConstants_overrideSubstrates@edu.stanford.covert.cell.sim.Process(wholeCellModelIDs);
        end

        function initializeConstants_overrideEnzymes(this, wholeCellModelIDs)
            %sort
            wholeCellModelIDs = sort(wholeCellModelIDs);

            %map old enzymes onto new
            [~, idxs1, idxs2] = intersect(this.enzymeWholeCellModelIDs, wholeCellModelIDs);

            %map old reactionCatalysisMatrix onto new
            reactionCatalysisMatrix = zeros(size(this.reactionCatalysisMatrix, 1), length(wholeCellModelIDs), size(this.reactionCatalysisMatrix, 3));
            reactionCatalysisMatrix(:, idxs2, :) = this.reactionCatalysisMatrix(:, idxs1, :);
            this.reactionCatalysisMatrix = reactionCatalysisMatrix;

            %super class method
            this.initializeConstants_overrideEnzymes@edu.stanford.covert.cell.sim.Process(wholeCellModelIDs);
        end
    end
end