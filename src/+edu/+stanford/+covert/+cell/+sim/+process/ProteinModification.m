%Protein Modification
%
% @wholeCellModelID Process_ProteinModification
% @name             Protein Modification
% @description
%   Biology
%   ==================
%   This process simulates protein modifications including
%   - Adductions:
%      - Ser/Thr/Tyr phosphorylation   (MG_109_DIMER)
%      - lipoate ligation              (MG_270_MONOMER)
%   - Ligations:
%      - C-terminal glutamate ligation (MG_012_MONOMER)
%
%   Many protein require covalent modifications to function properly. These
%   modifications are enzymatic catalyzed.
%
%   Knowledge Base
%   ==================
%   In this model of M. genitalium, there are currently 16 proteins that require
%   phosphorylation at anywhere from one to 11 sites each. Only one protein
%   requires lipoate ligation (pdhC), and only one requires glutamate ligation
%   (rplF). All other proteins require no modification and are guaranteed to
%   progress through this phase of protein maturation in a single time step.
%
%   The knowledge base representation of the stoichiometry of protein
%   monomer modification reactions includes both the unmodified amino acid
%   of the unmodified protein monomer on the left-hand-side and the
%   modified amino acid of the resulting modified protein monomer on the
%   right-hand-side.
%
%   Representation
%   ==================
%   The counts of unmodified and fully modified proteins are represented by the
%   unmodifiedMonomers and modifiedMonomers properties. Intermediate modified
%   states (proteins which have some, but not all of their requisite
%   modifications) are not represented here. The molecular weights of the
%   unmodified and fully modified protein mononomers are computed by the
%   knowledge base protein monomer class.
%
%   This process uses the reactionModificationMatrix, reactionCatalysisMatrix,
%   reactionStoichiometryMatrix, and enzymeBounds properties to represent the
%   modifications required to mature each protein. These properties are
%   initialized from the knowledge base by initializeConstants.
%   reactionModificationMatrix represents the protein monomer modified by each
%   reaction. reactionCatalysisMatrix represents the enzyme which catalyzes each
%   reaction. enzymeBounds represents the kcat of the enzyme for each reaction.
%   reactionStoichiometryMatrix represents the free metabolites reactants and
%   products of each reaction. Note reactionStoichiometryMatrix used in this
%   reaction is different from that of the superclass. The
%   reactionStoichiometryMatrix used here does not include either the unmodified
%   amino acid of the unmodified protein monomer on the left-hand-side of
%   reactions or the modified amino acid of the resulting modified protein
%   monomer on the right-hand-side of reactions; these amino acids are
%   represented within the unmodified and modified protein monomers. For this
%   process to mature a protein monomer, each of the reactions which modifies
%   that protein monomer must proceed.
%
%   Initialization
%   ==================
%   All protein monomers are initialized to the mature state. This is
%   accomplished by the simulation class initializeState method.
%
%   Simulation
%   ==================
%   While(true)
%     1. Calculate numbers of protein monomers that can be modified based on
%        substrate, enzyme, and unmodified protein monomer availability and
%        kinetics.
%     2. Randomly select protein monomer to modify weighted by limits
%        calculated in step (1).
%     3. Update substrates, enzymes, unmodified protein monomers
%     4. Repeat until insufficient resources to further modify protein
%        monomers
%   End
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/4/2010
classdef  ProteinModification < edu.stanford.covert.cell.sim.ReactionProcess
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
            'monomerCompartments';
            };
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {   %names of simulation state properties redundant with timecourses in this or other processes or the simulation
            'unmodifiedMonomers';
            'modifiedMonomers'};
    end
    
    %IDs, names, and local indices
    properties
        %substrateWholeCellModelIDs = {'ATP';'ADP';'H2O';...};
        substrateIndexs_aminoAcids          %indices of amino acids within substrates
        substrateIndexs_modifiedAminoAcids  %indices of modified amino acids within substrates
        substrateIndexs_atp                 %index of ATP within substrates
        substrateIndexs_adp                 %index of ADP within substrates
        substrateIndexs_amp                 %index of AMP within substrates
        substrateIndexs_hydrogen            %index of H within substrates
        substrateIndexs_phosphate           %index of Pi within substrates
        substrateIndexs_glutamate           %index of glutamate within substrates
        substrateIndexs_lipoylAmp           %index of lipoyl AMP within substrates
        substrateIndexs_lipoylLys           %index of lipoyl Lysine within substrates
        
        enzymeIndexs_serineThreonineKinase  %index of serine/threonine kinase (MG_109_DIMER) within enzymes
        enzymeIndexs_lipoylTransferase      %index of lipoyl transferase (MG_270_MONOMER) within enzymes
        enzymeIndexs_glutamateLigase        %index of glutamate ligase (MG_012_MONOMER) within enzymes
        
        reactionIndexs_adduction            %indices of adduction reactions within reactionWholeCellModelIDs
        reactionIndexs_ligation             %indices of ligation reactions within reactionWholeCellModelIDs
        
        unmodifiedMonomerWholeCellModelIDs  %whole cell model ids of unmodified protein monomers
        modifiedMonomerWholeCellModelIDs    %whole cell model ids of modified protein monomers
        
        monomerIndexs_modified              %indices of monomers that are modified within unmodifiedMonomerWholeCellModelIDs
        monomerIndexs_notmodified           %indices of monomers that are not modified within unmodifiedMonomerWholeCellModelIDs
        
        speciesIndexs_enzymes               %indexs of enzyme species within speciesReactantByproductMatrix, speciesReactantMatrix
    end
    
    %fixed biological constants
    properties
        monomerCompartments               %final compartments of protein monomers (when mature after translocation and incorporation into special compartments such as the terminal organelle)
        
        speciesReactantByproductMatrix    %stoichiometry of susbtrates, enzymes, RNA in modifications necessary to mature each modifying RNA: [RNA modifications] X [subtrates; enzymes; unmodified RNAs]
        speciesReactantMatrix             %reactant stoichiometry of susbtrates, enzymes, RNA in modifications necessary to mature each modifying RNA: [RNA modifications] X [subtrates; enzymes; unmodified RNAs]
    end
    
    %local state
    properties
        unmodifiedMonomers %counts of unmodified protein monomers (monomers X 1 X time)
        modifiedMonomers   %counts of modified protein monomers   (monomers X 1 X time)
    end
    
    %constructor
    methods
        function this = ProteinModification(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.ReactionProcess(wholeCellModelID, name);
        end
    end
    
    %communication between process/simulation
    methods
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)
            this.initializeConstants@edu.stanford.covert.cell.sim.ReactionProcess(...
                knowledgeBase, simulation, varargin{:});
            
            %% metabolite indices
            this.substrateIndexs_aminoAcids         = this.substrateIndexs({'GLU'; 'LYS'; 'SER'; 'THR'; 'TYR';});
            this.substrateIndexs_modifiedAminoAcids = this.substrateIndexs({'LIPOYLLYS'; 'pSER'; 'pTHR'; 'pTYR'});
            this.substrateIndexs_atp                = this.substrateIndexs({'ATP'});
            this.substrateIndexs_adp                = this.substrateIndexs({'ADP'});
            this.substrateIndexs_amp                = this.substrateIndexs({'AMP'});
            this.substrateIndexs_hydrogen           = this.substrateIndexs({'H'});
            this.substrateIndexs_phosphate          = this.substrateIndexs({'PI'});
            this.substrateIndexs_glutamate          = this.substrateIndexs({'GLU'});
            this.substrateIndexs_lipoylAmp          = this.substrateIndexs({'LIPOYLAMP'});
            this.substrateIndexs_lipoylLys          = this.substrateIndexs({'LIPOYLLYS'});
            
            %% enzyme indices
            this.enzymeIndexs_serineThreonineKinase  = this.enzymeIndexs({'MG_109_DIMER'});
            this.enzymeIndexs_lipoylTransferase      = this.enzymeIndexs({'MG_270_MONOMER'});
            this.enzymeIndexs_glutamateLigase        = this.enzymeIndexs({'MG_012_MONOMER'});
            
            %% protein monomers
            this.unmodifiedMonomerWholeCellModelIDs = this.monomer.wholeCellModelIDs(this.monomer.foldedIndexs);
            this.modifiedMonomerWholeCellModelIDs   = this.monomer.wholeCellModelIDs(this.monomer.matureIndexs);
            
            %% protein monomer modification reactions
            %indices
            this.reactionIndexs_adduction = find(strcmp(this.reactionTypes, 'adduction'));
            this.reactionIndexs_ligation  = find(strcmp(this.reactionTypes, 'ligation'));
            
            %Note in constrast to the knowledge base,
            %reactionStoichiometryMatrix doesn't include the unmodified
            %amino acid of the unmodified protein monomer on the
            %left-hand-side or the modified amino acid of the resulting
            %modified protein monomer.
            this.reactionStoichiometryMatrix(this.substrateIndexs_aminoAcids,         this.reactionIndexs_adduction) = 0;
            this.reactionStoichiometryMatrix(this.substrateIndexs_modifiedAminoAcids, this.reactionIndexs_adduction) = 0;
            
            %In constrast to other processes where reactionModificationMatrix
            %maps genes onto reactions, here reactionModificationMatrix
            %maps protein monomers onto reactions. Also un contrast to
            %other processes where entries in reactionModificationMatrix
            %contain the positions of the modification in the modified
            %macromolecule, here entries of reactionModificationMatrix are simply 0 or 1
            this.reactionModificationMatrix = this.reactionModificationMatrix(:, this.gene.mRNAIndexs)>0;
            
            this.monomerIndexs_modified    = find( any(this.reactionModificationMatrix, 1));
            this.monomerIndexs_notmodified = find(~any(this.reactionModificationMatrix, 1));
            
            this.initializeSpeciesNetwork();
            
            %compartments of protein monomers as seen by this process
            %  Note: proteins undergo modification before inclusion into
            %  the terminal organelle. Consequently, this process modifies
            %  proteins in the membrane and cytosol rather than in the
            %  terminal organelle membrane or cytosol.
            this.monomerCompartments = knowledgeBase.proteinMonomerCompartments;
            this.monomerCompartments(this.monomerCompartments == this.compartment.terminalOrganelleCytosolIndexs)  = this.compartment.cytosolIndexs;
            this.monomerCompartments(this.monomerCompartments == this.compartment.terminalOrganelleMembraneIndexs) = this.compartment.membraneIndexs;
        end
        
        %initialize reaction stoichiomety to be used in evolveState:
        %  [Monomer modifications] X [subtrates; enzymes; unmodified monomers]
        function initializeSpeciesNetwork(this)
            numMetabolites = size(this.reactionStoichiometryMatrix, 1);
            numEnzymes     = size(this.reactionCatalysisMatrix, 2);
            numProteins    = numel(this.monomerIndexs_modified);
            
            this.speciesReactantByproductMatrix = [this.reactionModificationMatrix(:, this.monomerIndexs_modified)' * [...
                -this.reactionStoichiometryMatrix' ...
                this.reactionCatalysisMatrix ./ repmat(this.enzymeBounds(:, 2) * this.stepSizeSec, 1, numEnzymes)] ...
                eye(numProteins)];
            this.speciesReactantMatrix = [this.reactionModificationMatrix(:, this.monomerIndexs_modified)' * max(0, [...
                -this.reactionStoichiometryMatrix' ...
                this.reactionCatalysisMatrix ./ repmat(this.enzymeBounds(:, 2) * this.stepSizeSec, 1, numEnzymes)]) ...
                eye(numProteins)];
            
            this.speciesIndexs_enzymes = (1:numEnzymes)' + numMetabolites;
        end
        
        %retrieve state from simulation
        function copyFromState(this)
            this.copyFromState@edu.stanford.covert.cell.sim.ReactionProcess();
            
            %compartments
            mComps = this.monomer.compartments(this.monomer.foldedIndexs);
            mComps(mComps == this.compartment.terminalOrganelleCytosolIndexs)  = this.compartment.cytosolIndexs;
            mComps(mComps == this.compartment.terminalOrganelleMembraneIndexs) = this.compartment.membraneIndexs;
            
            % indices
            numTimePoints = size(this.unmodifiedMonomers, 3);
            if numTimePoints == 1
                unmodifiedMonomerIndexs = sub2ind(size(this.monomer.counts), ...
                    this.monomer.foldedIndexs, mComps);
                modifiedMonomerIndexs = sub2ind(size(this.monomer.counts), ...
                    this.monomer.matureIndexs, mComps);
            else
                unmodifiedMonomerIndexs = sub2ind(size(this.monomer.counts), ...
                    repmat(this.monomer.foldedIndexs, [1 1 numTimePoints]),...
                    repmat(mComps, [1 1 numTimePoints]), ...
                    permute(repmat((1:numTimePoints)', [1 size(this.unmodifiedMonomers, 1) 1]), [2 3 1]));
                modifiedMonomerIndexs = sub2ind(size(this.monomer.counts), ...
                    repmat(this.monomer.matureIndexs, [1 1 numTimePoints]),...
                    repmat(mComps, [1 1 numTimePoints]), ...
                    permute(repmat((1:numTimePoints)', [1 size(this.unmodifiedMonomers, 1) 1]), [2 3 1]));
            end
            
            %monomers
            this.unmodifiedMonomers = this.monomer.counts(unmodifiedMonomerIndexs);
            this.modifiedMonomers   = this.monomer.counts(modifiedMonomerIndexs);
        end
        
        %send state to simulation
        function copyToState(this)
            this.copyToState@edu.stanford.covert.cell.sim.ReactionProcess();
            
            %compartments
            mComps = this.monomer.compartments(this.monomer.foldedIndexs);
            mComps(mComps == this.compartment.terminalOrganelleCytosolIndexs)  = this.compartment.cytosolIndexs;
            mComps(mComps == this.compartment.terminalOrganelleMembraneIndexs) = this.compartment.membraneIndexs;
            
            % indices
            numTimePoints = size(this.unmodifiedMonomers, 3);
            if numTimePoints == 1
                unmodifiedMonomerIndexs = sub2ind(size(this.monomer.counts), ...
                    this.monomer.foldedIndexs, mComps);
                modifiedMonomerIndexs = sub2ind(size(this.monomer.counts), ...
                    this.monomer.matureIndexs, mComps);
            else
                unmodifiedMonomerIndexs = sub2ind(size(this.monomer.counts), ...
                    repmat(this.monomer.foldedIndexs, [1 1 numTimePoints]),...
                    repmat(mComps, [1 1 numTimePoints]),...
                    permute(repmat((1:numTimePoints)', [1 size(this.unmodifiedMonomers, 1) 1]), [2 3 1]));
                modifiedMonomerIndexs = sub2ind(size(this.monomer.counts), ...
                    repmat(this.monomer.matureIndexs, [1 1 numTimePoints]),...
                    repmat(mComps, [1 1 numTimePoints]),...
                    permute(repmat((1:numTimePoints)', [1 size(this.unmodifiedMonomers, 1) 1]), [2 3 1]));
            end
            
            %monomers
            this.monomer.counts(unmodifiedMonomerIndexs) = this.unmodifiedMonomers;
            this.monomer.counts(modifiedMonomerIndexs)   = this.modifiedMonomers;
        end
    end
    
    %allocate memory for state
    methods
        function allocateMemoryForState(this, numTimePoints)
            this.allocateMemoryForState@edu.stanford.covert.cell.sim.ReactionProcess(numTimePoints);
            
            this.unmodifiedMonomers = zeros(length(this.unmodifiedMonomerWholeCellModelIDs), 1, numTimePoints);
            this.modifiedMonomers   = zeros(length(this.modifiedMonomerWholeCellModelIDs),   1, numTimePoints);
        end
    end
    
    methods
        %Calculate
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, ~, states)
            %substrate and byproducts
            nRxn = this.reactionModificationMatrix * states.monomerProductions;
            bmProd = max(0, -this.reactionStoichiometryMatrix) * nRxn;
            byProd = max(0,  this.reactionStoichiometryMatrix) * nRxn;
            
            %enzymes
            nRxn = this.reactionModificationMatrix * states.monomerProductions0;
            minEnzExp = 2 * this.reactionCatalysisMatrix' * (nRxn ./ this.enzymeBounds(:, 2));
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
        end
        
        %initialization: monomers initialized to mature/bound/inactivated state
        %by simulation initializeState method
        function initializeState(~)
        end
        
        %resource requirements
        function result = calcResourceRequirements_Current(this)
            result = max(0, -this.reactionStoichiometryMatrix) * min(...
                ceil(this.reactionCatalysisMatrix * this.enzymes * this.stepSizeSec), ...
                this.reactionModificationMatrix * this.unmodifiedMonomers);
        end
        
        %simulation
        function evolveState(this)
            %% update protein monomers that don't require modification
            this.modifiedMonomers(this.monomerIndexs_notmodified) = ...
                this.modifiedMonomers(this.monomerIndexs_notmodified) + ...
                this.unmodifiedMonomers(this.monomerIndexs_notmodified);
            this.unmodifiedMonomers(this.monomerIndexs_notmodified) = 0;
            
            %% stop early if no proteins require modification
            if ~any(this.unmodifiedMonomers)
                return;
            end
            
            %% simulate protein monomer modification
            numProteins = numel(this.monomerIndexs_modified);
            
            species = [
                this.substrates;
                this.enzymes;
                this.unmodifiedMonomers(this.monomerIndexs_modified)]';
            
            reactionFluxes = zeros(size(this.speciesReactantByproductMatrix,1),1);
            anyFlux = false;
            
            limits = species(ones(numProteins, 1), :) ./ this.speciesReactantMatrix;
            reactionLimits = min(limits, [], 2)';
            reactionLimits(isinf(reactionLimits) | isnan(reactionLimits) | reactionLimits < 0) = 0;
            isReactionInactive = reactionLimits <= 0;
            
            while true
                %compute maximum number of each protein monomer species that can be modified
                limits = species(ones(numProteins, 1), :) ./ max(0, this.speciesReactantByproductMatrix);
                reactionLimits = min(...
                    this.randStream.stochasticRound(min(limits(:, this.speciesIndexs_enzymes), [], 2)), ...
                    min(limits(:, [1:this.speciesIndexs_enzymes(1)-1 this.speciesIndexs_enzymes(end)+1:end]), [], 2))';
                reactionLimits(isReactionInactive | isinf(reactionLimits) | isnan(reactionLimits) | reactionLimits < 1) = 0;
                
                %stop if no more substrates can be modified
                if ~any(reactionLimits)
                    break;
                end
                anyFlux = true;
                
                %pick reaction
                selectedReaction = this.randStream.randsample(numel(reactionLimits), 1, true, reactionLimits);
                reactionFluxes(selectedReaction) = reactionFluxes(selectedReaction) + 1;
                
                %update metabolites, enzymes, unmodified substrates
                species = species - this.speciesReactantByproductMatrix(selectedReaction,:);
            end
            
            %stop if no reactions can proceed
            if ~anyFlux
                return;
            end
            
            % update substrates
            this.substrates = this.substrates + ...
                this.reactionStoichiometryMatrix * this.reactionModificationMatrix(:,this.monomerIndexs_modified) * reactionFluxes;
            
            %update protein monomers that require modification
            this.unmodifiedMonomers(this.monomerIndexs_modified) = this.unmodifiedMonomers(this.monomerIndexs_modified) - reactionFluxes;
            this.modifiedMonomers(this.monomerIndexs_modified)   = this.modifiedMonomers(this.monomerIndexs_modified)   + reactionFluxes;
        end
    end
    
    %get methods of dependent local state
    methods
        function value = getDryWeight(this)
            if size(this.unmodifiedMonomers, 3) == 1
                value = this.getDryWeight@edu.stanford.covert.cell.sim.ReactionProcess() + (...
                    this.monomer.molecularWeights(this.monomer.foldedIndexs)' * this.unmodifiedMonomers + ...
                    this.monomer.molecularWeights(this.monomer.matureIndexs)' * this.modifiedMonomers) ...
                    / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            else
                value = this.getDryWeight@edu.stanford.covert.cell.sim.ReactionProcess() + (...
                    permute(this.monomer.molecularWeights(this.monomer.foldedIndexs)' * permute(this.unmodifiedMonomers, [1 3 2]), [1 3 2]) + ...
                    permute(this.monomer.molecularWeights(this.monomer.matureIndexs)' * permute(this.modifiedMonomers,   [1 3 2]), [1 3 2])) ...
                    / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            end
        end
    end
end
