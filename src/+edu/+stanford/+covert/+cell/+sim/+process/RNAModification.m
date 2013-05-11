%RNA Modification
%
% @wholeCellModelID Process_RNAModification
% @name             RNA Modification
% @description
%   Biology
%   ==================
%   This process simulates r/tRNA modifications including
%   - rRNA methylation                                       MG_252_DIMER, MG_346_DIMER, MG_380_MONOMER, MG_463_MONOMER
%   - rRNA pseudouridation                                   MG_209_MONOMER, MG_370_MONOMER
%   - tRNA methylation                                       MG_347_DIMER, MG_445_DIMER
%   - tRNA pseudouridation                                   MG_182_DIMER
%   - tRNA lysidine synthetase                               MG_084_TETRAMER
%   - tRNA sulfur transfer                                   MG_295_MONOMER, MG_372_MONOMER
%   - tRNA uridine 5-carboxymethylaminomethyl modification   MG_008_379_TETRAMER
%
%   These modifications are believe to help RNAs fold properly and achieve their
%   catalytically active structure. The modifications are also believed to
%   improve RNA stability. In addition some tRNA modifications, namely the
%   modifications near the wobble position are believed to enhance codon
%   recognition. These modifications are enzymatically catalyzed by 13 proteins.
%
%   Knowledge Base
%   ==================
%   As of 8/11/2010 the M. genitalium knowledge base contained 91 RNA
%   modification reactions involving 38 RNAs and 13 enzymes. Each of the 38
%   RNAs is involved in 1-7 reactions. All other RNAs require no
%   modification and are guaranteed to progress through this phase of RNA
%   maturation in a single time step.
%
%   The knowledge base representation of the stoichiometry of RNA
%   modification reactions includes both the unmodified nucleic acid of the
%   unmodified RNA on the left-hand-side and the modified nucleic acid of
%   the resulting modified RNA on the right-hand-side.
%
%   Representation
%   ==================
%   The counts of unmodified and fully modified RNAs are represented by the
%   unmodifiedRNAs and modifiedRNAs properties. Intermediate modified
%   states (RNAs which have some, but not all of their requisite
%   modifications) are not represented here. The molecular weights of the
%   unmodified and fully modified RNAs are computed by the knowledge base RNA
%   classes.
%
%   This process uses the reactionModificationMatrix, reactionCatalysisMatrix,
%   reactionStoichiometryMatrix, and enzymeBounds properties to represent the
%   modifications required to mature each RNA. These properties are
%   initialized from the knowledge base by initializeConstants.
%   reactionModificationMatrix represents the RNA modified by each reaction.
%   reactionCatalysisMatrix represents the enzyme which catalyzes each reaction.
%   enzymeBounds represents the kcat of the enzyme for each reaction.
%   reactionStoichiometryMatrix represents the free metabolites reactants and
%   products of each reaction. Note reactionStoichiometryMatrix used in this
%   reaction is different from that of the superclass. The
%   reactionStoichiometryMatrix used here does not include either the unmodified
%   nucleic acid of the unmodified RNA on the left-hand-side of reactions or the
%   modified nucleic acid of the resulting modified RNA on the right-hand-side
%   of reactions; these nucleic acids are represented within the unmodified and
%   modified RNAs. For this process to mature a RNA, each of the reactions which
%   modifies that RNA must proceed.
%
%   Initialization
%   ==================
%   All RNAs are initialized to the mature state. This is accomplished by the
%   simulation class initializeState method.
%
%   Simulation
%   ==================
%   While(true)
%     1. Calculate numbers of RNAs that can be modified based on
%        substrate, enzyme, and unmodified RNA availability and kinetics.
%     2. Randomly select RNA to modify weighted by limits calculated in
%        step (1).
%     3. Update substrates, enzymes, unmodified RNAs
%     4. Repeat until insufficient resources to further modify RNAs
%   End
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/10/2010
classdef RNAModification < edu.stanford.covert.cell.sim.ReactionProcess

    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
            'speciesReactantByproductMatrix'
            'speciesReactantMatrix'
			};		
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {   %names of simulation state properties redundant with timecourses in this or other processes or the simulation
            'unmodifiedRNAs';
            'modifiedRNAs'};
    end

    %IDs, names, and local indices
    properties
        substrateIndexs_amet                             %index   within substrates of amet
        substrateIndexs_fthf5                            %index   within substrates of 5-fthf
        substrateIndexs_cys                              %index   within substrates of cysteine
        substrateIndexs_gly                              %index   within substrates of glycine
        substrateIndexs_lys                              %index   within substrates of lysine
        substrateIndexs_atp                              %index   within substrates of ATP
        substrateIndexs_water                            %index   within substrates of water
        substrateIndexs_hydrogen                         %index   within substrates of hydrogen
        substrateIndexs_nmp                              %indices within substrates of NMPs
        substrateIndexs_modifiedNMP                      %indices within substrates of modified NMPs

        enzymeIndexs_rRNA16SDimethyladenosineTransferase %index within enzymes of   MG_463_MONOMER       rRNA 16S dimethyladenosine transferase
        enzymeIndexs_rRNA16SMethyltransferaseGidB        %index within enzymes of   MG_380_MONOMER       rRNA 16S methyltransferase GidB
        enzymeIndexs_rRNA23SMethyltransferaseI           %index within enzymes of   MG_252_DIMER         rRNA 23S methyltransferase I
        enzymeIndexs_rRNA23SMethyltransferaseII          %index within enzymes of   MG_346_DIMER         rRNA 23S methyltransferase II
        enzymeIndexs_rRNA23SPseudouridineSynthaseI       %index within enzymes of   MG_209_MONOMER       rRNA 23S pseudouridine synthase I
        enzymeIndexs_rRNA23SPseudouridineSynthaseII      %index within enzymes of   MG_370_MONOMER       rRNA 23S pseudouridine synthase II
        enzymeIndexs_tRNAGuanineN1Methyltransferase      %index within enzymes of   MG_445_DIMER         tRNA guanine-N1-methyltransferase
        enzymeIndexs_tRNAGuanineN7Methyltransferase      %index within enzymes of   MG_347_DIMER         tRNA guanine-N7-methyltransferase
        enzymeIndexs_tRNALysidineSynthetase              %index within enzymes of   MG_084_TETRAMER      tRNA lysidine synthetase
        enzymeIndexs_tRNAPseudouridineSynthase           %index within enzymes of   MG_182_DIMER         tRNA pseudouridine synthase
        enzymeIndexs_tRNAUracil2Sulfurtransferase        %index within enzymes of   MG_295_MONOMER       tRNA uracil-2-sulfurtransferase
        enzymeIndexs_tRNAUracil4Sulfurtransferase        %index within enzymes of   MG_372_DIMER         tRNA uracil-4-sulfurtransferase
        enzymeIndexs_tRNAUracil5Carboxymethylaminomethyl %index within enzymes of   MG_008_379_TETRAMER  tRNA uracil-5-carboxymethylaminomethyl

        unmodifiedRNAWholeCellModelIDs                   %whole cell model ids of unmodified RNAs
        modifiedRNAWholeCellModelIDs                     %whole cell model ids of modified RNAs
        
        rnaIndexs_modified    %indices of RNAs that are modified within unmodifiedRNAWholeCellModelIDs
        rnaIndexs_unmodified  %indices of RNAs that are not modified within unmodifiedRNAWholeCellModelIDs
        rnaMask_modified      %boolean indicating RNAs that are modified
        rnaMask_unmodified    %boolean indicating RNAs that are not modified
        
        speciesIndexs_enzymes %indexs of enzyme species within speciesReactantByproductMatrix, speciesReactantMatrix
    end

    %fixed biological constants
    properties
        speciesReactantByproductMatrix %stoichiometry of susbtrates, enzymes, RNA in modifications necessary to mature each modifying RNA: [RNA modifications] X [subtrates; enzymes; unmodified RNAs]
        speciesReactantMatrix          %reactant stoichiometry of susbtrates, enzymes, RNA in modifications necessary to mature each modifying RNA: [RNA modifications] X [subtrates; enzymes; unmodified RNAs]
    end

    %global state (copied locally for convenience)
    properties
        unmodifiedRNAs %counts of unmodified RNAs (RNAs X 1 X time)
        modifiedRNAs   %counts of modified RNAs   (RNAs X 1 X time)
    end

    %constructor
    methods
        function this = RNAModification(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.ReactionProcess(wholeCellModelID, name);
        end
    end

    %communication between process/simulation
    methods
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)
            this.initializeConstants@edu.stanford.covert.cell.sim.ReactionProcess(knowledgeBase, simulation, varargin{:});

            %metabolites indices
            this.substrateIndexs_amet        = this.substrateIndexs({'AMET'});
            this.substrateIndexs_fthf5       = this.substrateIndexs({'FTHF5'});
            this.substrateIndexs_cys         = this.substrateIndexs({'CYS'});
            this.substrateIndexs_gly         = this.substrateIndexs({'GLY'});
            this.substrateIndexs_lys         = this.substrateIndexs({'LYS'});
            this.substrateIndexs_atp         = this.substrateIndexs({'ATP'});
            this.substrateIndexs_water       = this.substrateIndexs({'H2O'});
            this.substrateIndexs_hydrogen    = this.substrateIndexs({'H'});
            this.substrateIndexs_nmp         = this.substrateIndexs({'AMP'; 'CMP'; 'GMP'; 'UMP'});
            this.substrateIndexs_modifiedNMP = this.substrateMetaboliteLocalIndexs(...
                strcmp({knowledgeBase.metabolites(this.substrateMetaboliteGlobalIndexs).category}, ...
                'modified ribonucleotide monophosphate'));

            %enzyme indices
            this.enzymeIndexs_rRNA16SDimethyladenosineTransferase = this.enzymeIndexs({'MG_463_MONOMER'});      %dimethyladenosine transferase
            this.enzymeIndexs_rRNA16SMethyltransferaseGidB        = this.enzymeIndexs({'MG_380_MONOMER'});      %methyltransferase GidB
            this.enzymeIndexs_rRNA23SMethyltransferaseI           = this.enzymeIndexs({'MG_252_DIMER'});        %23S rRNA methyltransferase; G2251
            this.enzymeIndexs_rRNA23SMethyltransferaseII          = this.enzymeIndexs({'MG_346_DIMER'});        %23S rRNA methyltransferase; U2552
            this.enzymeIndexs_rRNA23SPseudouridineSynthaseI       = this.enzymeIndexs({'MG_209_MONOMER'});      %23S rRNA pseudouridine synthase; U955, U2504, U2580
            this.enzymeIndexs_rRNA23SPseudouridineSynthaseII      = this.enzymeIndexs({'MG_370_MONOMER'});      %23S rRNA pseudouridine synthase; U1911, U1915, U1917
            this.enzymeIndexs_tRNAGuanineN1Methyltransferase      = this.enzymeIndexs({'MG_445_DIMER'});        %tRNA (guanine-N1)-methyltransferase
            this.enzymeIndexs_tRNAGuanineN7Methyltransferase      = this.enzymeIndexs({'MG_347_DIMER'});        %tRNA (guanine-N(7)-)-methyltransferase
            this.enzymeIndexs_tRNALysidineSynthetase              = this.enzymeIndexs({'MG_084_TETRAMER'});     %tRNA(Ile)-lysidine synthetase
            this.enzymeIndexs_tRNAPseudouridineSynthase           = this.enzymeIndexs({'MG_182_DIMER'});        %tRNA pseudouridine synthase A
            this.enzymeIndexs_tRNAUracil2Sulfurtransferase        = this.enzymeIndexs({'MG_295_MONOMER'});      %tRNA U34 sulfurtransferase
            this.enzymeIndexs_tRNAUracil4Sulfurtransferase        = this.enzymeIndexs({'MG_372_DIMER'});        %thiamine biosynthesis/tRNA modification protein ThiI
            this.enzymeIndexs_tRNAUracil5Carboxymethylaminomethyl = this.enzymeIndexs({'MG_008_379_TETRAMER'}); %tRNA uridine 5-carboxymethylaminomethyl modification enzyme

            %RNAs
            this.unmodifiedRNAWholeCellModelIDs = this.rna.wholeCellModelIDs(this.rna.processedIndexs);
            this.modifiedRNAWholeCellModelIDs   = this.rna.wholeCellModelIDs(this.rna.matureIndexs);

            %RNA modification reactions
            %- Note in constrast to the knowledge base,
            %  reactionStoichiometryMatrix doesn't include the unmodified
            %  nucleic acid of the unmodified RNA on the left-hand-side or
            %  the modified nucleic acid of the resulting modified RNA.
            %- In constrast to other processes where
            %  reactionModificationMatrix maps genes onto reactions, here
            %  reactionModificationMatrix maps RNAs onto reactions
            %- In contrast to other processes where entries in
            %  reactionModificationMatrix contain the positions of the
            %  modification in the modified macromolecule, here entries of
            %  reactionModificationMatrix are simply 0 or 1
            this.reactionStoichiometryMatrix(this.substrateIndexs_nmp, :) = ...
                max(this.reactionStoichiometryMatrix(this.substrateIndexs_nmp, :), 0);
            this.reactionStoichiometryMatrix(this.substrateIndexs_modifiedNMP, :) = 0;
            this.reactionModificationMatrix = (this.reactionModificationMatrix > 0) * ...
                this.rna.matureRNAGeneComposition;
            
            this.rnaMask_modified = any(this.reactionModificationMatrix, 1);
            this.rnaMask_unmodified = ~any(this.reactionModificationMatrix, 1);
            
            this.rnaIndexs_modified = find(this.rnaMask_modified);
            this.rnaIndexs_unmodified = find(this.rnaMask_unmodified);
            
            this.initializeSpeciesNetwork();
        end
        
        %initialize reaction stoichiomety to be used in evolveState:
        %  [RNA modifications] X [subtrates; enzymes; unmodified RNAs]
        function initializeSpeciesNetwork(this)            
            numMetabolites = size(this.reactionStoichiometryMatrix, 1);
            numEnzymes     = size(this.reactionCatalysisMatrix, 2);
            numRNAs        = sum(this.rnaMask_modified);
            
            this.speciesIndexs_enzymes = (1:numEnzymes)' + numMetabolites;
            
            this.speciesReactantByproductMatrix = [this.reactionModificationMatrix(:, this.rnaMask_modified)' * [...
                -this.reactionStoichiometryMatrix'...
                this.reactionCatalysisMatrix ./ repmat(this.enzymeBounds(:,2) * this.stepSizeSec, 1, numEnzymes)]...
                eye(numRNAs)];
            this.speciesReactantMatrix = [this.reactionModificationMatrix(:, this.rnaMask_modified)' * max(0, [...
                -this.reactionStoichiometryMatrix'...
                this.reactionCatalysisMatrix ./ repmat(this.enzymeBounds(:, 2) * this.stepSizeSec, 1, numEnzymes)])...
                eye(numRNAs)];
        end

        %retrieve state from simulation
        function copyFromState(this)
            this.copyFromState@edu.stanford.covert.cell.sim.ReactionProcess();

            this.unmodifiedRNAs = this.rna.counts(this.rna.processedIndexs, this.compartment.cytosolIndexs, :);
            this.modifiedRNAs   = this.rna.counts(this.rna.matureIndexs,    this.compartment.cytosolIndexs, :);
        end

        %send state to simulation
        function copyToState(this)
            this.copyToState@edu.stanford.covert.cell.sim.ReactionProcess();

            this.rna.counts(this.rna.processedIndexs, this.compartment.cytosolIndexs, :) = this.unmodifiedRNAs;
            this.rna.counts(this.rna.matureIndexs,    this.compartment.cytosolIndexs, :) = this.modifiedRNAs;
        end
    end

    %allocate memory for state
    methods
        function allocateMemoryForState(this, numTimePoints)
            this.allocateMemoryForState@edu.stanford.covert.cell.sim.ReactionProcess(numTimePoints);

            this.unmodifiedRNAs = zeros(length(this.unmodifiedRNAWholeCellModelIDs), 1, numTimePoints);
            this.modifiedRNAs   = zeros(length(this.modifiedRNAWholeCellModelIDs),   1, numTimePoints);
        end
    end

    %model
    methods
        %Calculate
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, ~, states)
            %substrate and byproducts
            nRxn = this.reactionModificationMatrix * states.rnaProductions;

            bmProd = max(0, -this.reactionStoichiometryMatrix) * nRxn;
            byProd = max(0,  this.reactionStoichiometryMatrix) * nRxn;
            
            %enzymes
            minEnzExp = 2 * this.reactionCatalysisMatrix' * ...
                ((this.reactionModificationMatrix * states.rnaProductions0) ./ this.enzymeBounds(:, 2));
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
        end      

        %initialization: RNAs initialized to mature/aminoacylated state by
        %simulation initializeState method
        function initializeState(~)
        end

        %resource requirements
        function result = calcResourceRequirements_Current(this)
            result = ...
                max(0, -this.reactionStoichiometryMatrix) * min(...
                ceil(this.reactionCatalysisMatrix * this.enzymes * this.stepSizeSec), ...
                (this.reactionModificationMatrix * this.unmodifiedRNAs));
        end

        %simulation
        function evolveState(this)
            %% update RNAs that don't require modification
            this.modifiedRNAs(this.rnaMask_unmodified) = ...
                this.modifiedRNAs(this.rnaMask_unmodified) + ...
                this.unmodifiedRNAs(this.rnaMask_unmodified);
            this.unmodifiedRNAs(this.rnaMask_unmodified) = 0;
            
            %% stop early if no RNAs need to be modified
            if ~any(this.unmodifiedRNAs)
                return;
            end

            %% simulate RNA modification
            numRNAs = numel(this.rnaIndexs_modified);
            
            species = [
                this.substrates;
                this.enzymes;
                this.unmodifiedRNAs(this.rnaMask_modified)]';
            
            reactionFluxes = zeros(size(this.speciesReactantByproductMatrix, 1), 1);
            anyFlux = false;
            
            limits = species(ones(numRNAs, 1), :) ./ this.speciesReactantMatrix;
            limits(:, this.substrateIndexs_water) = NaN;
            limits(:, this.substrateIndexs_hydrogen) = NaN;
            reactionLimits = min(limits, [], 2)';
            reactionLimits(isinf(reactionLimits) | isnan(reactionLimits) | reactionLimits < 0) = 0;
            isReactionInactive = reactionLimits <= 0;
            
            while true
                %compute maximum number of each RNA species that can be modified
                limits = species(ones(numRNAs, 1), :) ./ max(0, this.speciesReactantByproductMatrix);
                limits(:, this.substrateIndexs_water) = NaN;
                limits(:, this.substrateIndexs_hydrogen) = NaN;
                reactionLimits = min(...
                    this.randStream.stochasticRound(min(limits(:, this.speciesIndexs_enzymes), [], 2)), ...
                    min(limits(:, [1:this.speciesIndexs_enzymes(1)-1 this.speciesIndexs_enzymes(end)+1:end]), [], 2))';
                reactionLimits(isReactionInactive | isinf(reactionLimits) | isnan(reactionLimits) | reactionLimits < 1) = 0;

                %stop if no more substrates can be modified
                if ~any(reactionLimits); break; end;
                anyFlux = true;
                
                %pick reaction
                selectedReaction = this.randStream.randsample(numel(reactionLimits), 1, true, reactionLimits);
                reactionFluxes(selectedReaction) = ...
                    reactionFluxes(selectedReaction) + 1;

                %update metabolites, enzymes, unmodified RNAs
                species = species - this.speciesReactantByproductMatrix(selectedReaction, :);
            end

            %stop if no reactions can proceed
            if ~anyFlux
                return;
            end
            
            % update substrates
            this.substrates = this.substrates + ...
                this.reactionStoichiometryMatrix * this.reactionModificationMatrix(:, this.rnaMask_modified) * reactionFluxes;

            %update RNAs that require modification
            this.unmodifiedRNAs(this.rnaMask_modified) = this.unmodifiedRNAs(this.rnaMask_modified) - reactionFluxes;
            this.modifiedRNAs(this.rnaMask_modified)   = this.modifiedRNAs(this.rnaMask_modified)   + reactionFluxes;
        end
    end

    %get methods of dependent local state
    methods
        function value = getDryWeight(this)
            if size(this.unmodifiedRNAs, 3) == 1
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + (...
                    this.rna.molecularWeights(this.rna.processedIndexs)' * this.unmodifiedRNAs + ...
                    this.rna.molecularWeights(this.rna.matureIndexs)'    * this.modifiedRNAs) ...
                    / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            else
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + (...
                    permute(this.rna.molecularWeights(this.rna.processedIndexs)' * permute(this.unmodifiedRNAs, [1 3 2]),[1 3 2]) + ...
                    permute(this.rna.molecularWeights(this.rna.matureIndexs)'    * permute(this.modifiedRNAs,   [1 3 2]),[1 3 2])) ...
                    / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            end
        end
    end
end
