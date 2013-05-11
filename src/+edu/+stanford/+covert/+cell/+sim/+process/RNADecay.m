%RNA Decay
%
% @wholeCellModelID Process_RNADecay
% @name             RNA Decay
% @description
%   Biology
%   ===============
%   In presence of ribonucleases such as ribonuclease R (MG_104_MONOMER) RNAs
%   have relatively short half lives compared to that of other macromolecules
%   (eg. protein, DNA) and the M. genitalium cell cycle length. The relatively
%   short half lives of RNAs enables the small M. genitalium with its very small
%   pool of RNAs and particularly mRNAs to sample a broader range of
%   configurations of the RNA pool over a shorter period that would be possible
%   with longer half lifes. This helps the cell more finely tune the expression
%   of proteins, more efficiently execute cell-cycle dependent events, and
%   respond to the external environment. This enhanced fitness due to short RNA
%   half lifes comes at a large energetic cost however.
%
%   In addition to ribonucleases, aminoacylated RNAs require peptidyl tRNA
%   hydrolase (MG_083_MONOMER) to release their conjugated amino acids.
%
%   This process decays all species of RNA, and at all maturation states
%   including aminoacylated states.
%
%   Knowledge Base
%   ===============
%   The knowledge base contains experimentally measured half lifes of many RNA
%   species measured largely in E. coli and mapped to M. genitalium by homology.
%   These half lifes are refined, by simulation.fitConstants to make them
%   consistent with other experimental data used to fit the model. Prior to
%   fitting missing half lifes are imputed either as the average of that of all
%   measured RNA species.
%
%      Type   Avg Half Life (m)
%      ====   =================
%      mRNA   4.5 +/- 2.0
%      rRNA   150
%      sRNA   89
%      tRNA   45
%
%   Representation
%   ===============
%   The substrates, enzymes, and RNAs properties represent the counts of
%   metabolites, ribonuclease R and peptidyl tRNA hydrolase enzymes, and RNAs.
%   This process contains no intermediate representation of RNA degradation; RNA
%   degradation is treated as an all-or-nothing event that either proceeds to
%   complete with a time step or doesn't progress at all.
%
%   decayRates represents the decay rate of each RNA species in seconds.
%   decayRates is informed by experimentally measured RNA half lifes organized
%   in the knowledge base, and fit by simulation.fitConstants. decayReactions
%   represents the metabolites required to decay each RNA species, and the
%   metabolic byproducts of the decay of each RNA species. decayReactions is
%   computed by the knowledge RNA classes based on the sequence, processing, and
%   modifications of each RNA species.
%
%   Initialization
%   ===============
%   All RNAs are initialized to the mature state. This is accomplished by the
%   simulation class initializeState method.
%
%   Simulation
%   ===============
%   This process models RNA decay as an enyzme-dependent poisson process with
%   rate parameter:
%     lambda = RNAs .* decayRates * stepSizeSec
%
%   Algorithm
%   +++++++++++++++
%   1. Stochastically select RNAs to decay based on poission distribution with
%      lambda = RNAs .* decayRates * stepSizeSec
%   2. (Ignore limits to decay posed by availability of metabolite reactants
%      since the only reactant is water, and water is abundantly available)
%   3. Limit RNA decay by available enzyme activity
%      a. All RNAs require ribonuclease R to decay
%      b. Additionally, only decay aminoacylated tRNAs up to the limit of
%         available peptidyl tRNA hydrolase activity.
%   4. Update counts of RNAs
%   5. Update counts of metabolic byproducts of RNA decay
%
% Author: Markus Covert, mcovert@stanford.edu
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/30/2010
classdef RNADecay < edu.stanford.covert.cell.sim.Process

    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
		    'peptidylTRNAHydrolaseSpecificRate';
            'ribonucleaseRFragmentLength';
            'decayReactions';
            };			
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {   %names of simulation state properties redundant with timecourses in this or other processes or the simulation
            'RNAs'};
    end

    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs = {};   %whole cell model IDs of stimuli

        substrateWholeCellModelIDs = {}; %whole cell model IDs of substrates
        substrateIndexs_hydrogen         %index within substrates of hydrogen
        substrateIndexs_water            %index within substrates of water
        substrateIndexs_methionine       %index within substrates of methionine
        substrateIndexs_fmethionine      %index within substrates of formylmethionine
        substrateIndexs_glutamate        %index within substrates of glutamate
        substrateIndexs_glutamine        %index within substrates of glutamine
        substrateIndexs_formate          %index within substrates of formate
        substrateIndexs_ammonia          %index within substrates of ammonia
        substrateIndexs_aminoAcids       %index within substrates of amino acids
        substrateIndexs_nmps             %index within substrates of NMPs

        enzymeWholeCellModelIDs = {      %enzyme whole cell model ids
            'MG_104_MONOMER';            %ribonuclease R
            'MG_083_MONOMER'};           %peptidyl-tRNA hydrolase
        enzymeIndexs_ribonucleaseR         = 1; %index within enzymes of ribonuclease R
        enzymeIndexs_peptidylTRNAHydrolase = 2; %index within enzymes of peptidyl-tRNA hydrolase
        
        matureTRNAIndexs                 %indices of mature tRNAs within RNAs
        matureTMRNAIndexs                %indices of mature tmRNAs within RNAs
    end
    
    %fixed biological constants
    properties
        peptidylTRNAHydrolaseSpecificRate   %0.700 [PUB_0026]
        ribonucleaseRFragmentLength         %5 [PUB_0039]
        decayReactions                      %adjacency matrix -- RNAs X (reactants and products of RNA decay)
    end

    %global state (copied locally for convenience)
    properties
        RNAs                                %counts of RNAs
    end
            
    %global state (referenced locally for convenience)
    properties
        transcripts                         %New Transcripts state class
    end

    %constructor
    methods
        function this = RNADecay(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.Process(wholeCellModelID, name);
        end
    end

    %communication between process/simulation
    methods
        %set references to state objects
        function storeObjectReferences(this, simulation)
            this.storeObjectReferences@edu.stanford.covert.cell.sim.Process(simulation);
            
            this.transcripts = simulation.state('Transcript');
            this.states = [this.states; {this.transcripts}];
        end
        
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)
            s = this.rna;
            g = this.gene;
            
            %include all metabolites involved in RNA decay
            decayReactions_nascentRNA       = knowledgeBase.nascentRNAs.decayReactions;
            decayReactions_processedRNA     = knowledgeBase.processedRNAs.decayReactions;
            decayReactions_intergenicRNA    = knowledgeBase.intergenicRNAs.decayReactions;
            decayReactions_matureRNA        = knowledgeBase.matureRNAs.decayReactions;
            decayReactions_aminoacylatedRNA = knowledgeBase.aminoacylatedRNAs.decayReactions;  
            
            this.substrateWholeCellModelIDs = unique([simulation.state('Metabolite').wholeCellModelIDs(...
                any(decayReactions_nascentRNA,       1) | ...
                any(decayReactions_processedRNA,     1) | ...
                any(decayReactions_intergenicRNA,    1) | ...
                any(decayReactions_matureRNA,        1) | ...
                any(decayReactions_aminoacylatedRNA, 1));
                'H';'H2O';'NH3';'FOR';
                'ALA';'ARG';'ASN';'ASP';'CYS';'GLN';'GLU';'GLY';'HIS';'ILE';'LEU';'LYS';'MET';'PHE';'PRO';'SER';'THR';'TRP';'TYR';'VAL';'FMET']);

            %super class method
            this.initializeConstants@edu.stanford.covert.cell.sim.Process(...
                knowledgeBase, simulation, varargin{:});
            
            %substrate indices
            this.substrateIndexs_hydrogen    = this.substrateIndexs({'H'});
            this.substrateIndexs_water       = this.substrateIndexs({'H2O'});
            this.substrateIndexs_methionine  = this.substrateIndexs({'MET'});
            this.substrateIndexs_fmethionine = this.substrateIndexs({'FMET'});
            this.substrateIndexs_glutamate   = this.substrateIndexs({'GLU'});
            this.substrateIndexs_glutamine   = this.substrateIndexs({'GLN'});
            this.substrateIndexs_ammonia     = this.substrateIndexs({'NH3'});
            this.substrateIndexs_formate     = this.substrateIndexs({'FOR'});
            this.substrateIndexs_aminoAcids  = this.substrateIndexs({'ALA';'ARG';'ASN';'ASP';'CYS';'GLN';'GLU';'GLY';'HIS';'ILE';'LEU';'LYS';'MET';'PHE';'PRO';'SER';'THR';'TRP';'TYR';'VAL';'FMET'});
            this.substrateIndexs_nmps        = this.substrateIndexs({'AMP';'CMP';'GMP';'UMP'});
            
            this.matureTRNAIndexs = find(any(s.matureRNAGeneComposition(g.tRNAIndexs, :), 1))';
            this.matureTMRNAIndexs = s.getIndexs('MG_0004');
            
            this.decayReactions = zeros(numel(this.rna.wholeCellModelIDs), numel(this.substrateWholeCellModelIDs));
            this.decayReactions(s.nascentIndexs, :)       = decayReactions_nascentRNA(:,       this.substrateMetaboliteGlobalIndexs);
            this.decayReactions(s.processedIndexs, :)     = decayReactions_processedRNA(:,     this.substrateMetaboliteGlobalIndexs);
            this.decayReactions(s.intergenicIndexs, :)    = decayReactions_intergenicRNA(:,    this.substrateMetaboliteGlobalIndexs);
            this.decayReactions(s.matureIndexs, :)        = decayReactions_matureRNA(:,        this.substrateMetaboliteGlobalIndexs);
            this.decayReactions(s.boundIndexs, :)         = decayReactions_matureRNA(:,        this.substrateMetaboliteGlobalIndexs);
            this.decayReactions(s.misfoldedIndexs, :)     = decayReactions_matureRNA(:,        this.substrateMetaboliteGlobalIndexs);
            this.decayReactions(s.damagedIndexs, :)       = decayReactions_matureRNA(:,        this.substrateMetaboliteGlobalIndexs);
            this.decayReactions(s.aminoacylatedIndexs, :) = decayReactions_aminoacylatedRNA(:, this.substrateMetaboliteGlobalIndexs);           
        end

        %retrieve state from simulation
        function copyFromState(this)
            this.copyFromState@edu.stanford.covert.cell.sim.Process();

            this.RNAs = this.rna.counts(:, this.compartment.cytosolIndexs, :);
        end

        %send state to simulation
        function copyToState(this)
            this.copyToState@edu.stanford.covert.cell.sim.Process();

            this.rna.counts(:, this.compartment.cytosolIndexs, :) = this.RNAs;
        end
    end
    
    %memory alloction for state
    methods
        function allocateMemoryForState(this, numTimePoints)
            this.allocateMemoryForState@edu.stanford.covert.cell.sim.Process(numTimePoints);
            
            this.RNAs = zeros(size(this.rna.counts, 1), 1, numTimePoints);
        end
    end

    %model
    methods
        %Calculate
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, ~, states)
            import edu.stanford.covert.util.ComputationUtil;
            invMat = this.rna.intergenicRNAMatrix * ...
                ComputationUtil.invertCompositionMatrix(this.rna.nascentRNAMatureRNAComposition);
            
            %% substrate and byproducts
            %data
            matureRNAReactions = this.decayReactions(this.rna.aminoacylatedIndexs, :); %same as mature for (m,r)RNA and for non-aminoacylated sRNA
            intergenicRNAReactions = this.decayReactions(this.rna.intergenicIndexs, :);
            intergenicRNADecays = invMat * states.rnaProductions;

            %RNA decay
            bmProd = ...
                + max(0, -matureRNAReactions)'     * states.rnaDecays ...
                + max(0, -intergenicRNAReactions)' * intergenicRNADecays;
            byProd = ...
                + max(0,  matureRNAReactions)'     * states.rnaDecays ...
                + max(0,  intergenicRNAReactions)' * intergenicRNADecays;
            
            %% enzymes
            
            %data
            rnaDecays = states.rnaDecays0;
            aminoacylatedRNADecays = rnaDecays([this.matureTRNAIndexs; this.matureTMRNAIndexs]);
            intergenicRNADecays = invMat * states.rnaProductions0;

            %RNA decay
            minEnzExp = zeros(size(this.enzymeWholeCellModelIDs));
            minEnzExp(this.enzymeIndexs_ribonucleaseR) = ...
                2 * (sum(rnaDecays) + sum(intergenicRNADecays));
            minEnzExp(this.enzymeIndexs_peptidylTRNAHydrolase) = ...
                2 * sum(aminoacylatedRNADecays) / this.peptidylTRNAHydrolaseSpecificRate;
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
        end

        %initialization: RNAs intialized to mature/aminoacylated state by
        %simulation intializeState method
        function initializeState(~)
        end

        %resource requirements
        function result = calcResourceRequirements_Current(this)
            if this.enzymes(this.enzymeIndexs_ribonucleaseR) == 0
                result = zeros(size(this.substrates));
                return;
            end
            
            result = max(0, -this.decayReactions' * (min(1, this.rna.decayRates) .* this.RNAs));
            if ~isempty(this.transcripts.abortedTranscripts)
                result(this.substrateIndexs_water) = ...
                    + result(this.substrateIndexs_water) ...
                    + sum(this.transcripts.abortedTranscripts(:, 2) - 1);
            end
        end

        %simulation
        function evolveState(this)
            % import classes
            import edu.stanford.covert.cell.kb.ssRNA;
            
            % numbers of enzymes
            ribonucleaseR = this.enzymes(this.enzymeIndexs_ribonucleaseR);
            peptidylTRNAHydrolase = this.randStream.stochasticRound(...
                this.enzymes(this.enzymeIndexs_peptidylTRNAHydrolase) ...
                * this.peptidylTRNAHydrolaseSpecificRate ...
                * this.stepSizeSec);
            
            %Ribonuclease R required for decay, terminate early if no ribonuclease R
            if ribonucleaseR == 0
                return;
            end
            
            %% decay all aborted transcripts
            abortedSeqs = this.transcripts.abortedSequences;
            abortedTfs = false(size(abortedSeqs));
            for i = 1:numel(abortedSeqs)
                substrateCost = ssRNA.computeDecayReaction(ssRNA.computeBaseCount(...
                    abortedSeqs{i}, numel(this.substrates), this.substrateIndexs_nmps), ...
                    numel(abortedSeqs{i}), 'linear', ...
                    this.substrateIndexs_water, this.substrateIndexs_hydrogen)';
                if any(this.substrates < -substrateCost)
                    break;
                end
                abortedTfs(i) = true;
                this.substrates = this.substrates + substrateCost;
            end
            this.transcripts.abortedTranscripts = this.transcripts.abortedTranscripts(~abortedTfs, :);
            
            %% Stochastically decay free RNA as poisson process
            decayingRNAs = min(this.randStream.random('poisson', ...
                this.RNAs .* min(1e6, this.rna.decayRates * this.stepSizeSec)), ...
                this.RNAs);
            if ~any(decayingRNAs)
                return;
            end
            
            %Require peptidyl tRNA hydrolase to decay aminoacylated tRNAs
            tmp = decayingRNAs(this.rna.aminoacylatedIndexs);
            tmp2 = zeros(size(tmp));
            while any(tmp)
                if peptidylTRNAHydrolase <= 0
                    break;
                end
                idx = this.randStream.randsample(numel(tmp), 1, true, tmp);
                tmp(idx) = tmp(idx) - 1;
                tmp2(idx) = tmp2(idx) + 1;
                peptidylTRNAHydrolase = peptidylTRNAHydrolase - 1;
            end
            decayingRNAs(this.rna.aminoacylatedIndexs) = tmp2;
            
            %require substrates (water) to decay RNAs
            tmp = decayingRNAs;
            tmp2 = zeros(size(tmp));
            water = this.substrates(this.substrateIndexs_water);
            waterReqs = max(0, -this.decayReactions(:, this.substrateIndexs_water));
            while any(tmp)
                idx = this.randStream.randsample(numel(tmp), 1, true, tmp);
                if water < waterReqs(idx);
                    break;
                end
                water = water - waterReqs(idx);
                tmp(idx) = tmp(idx) - 1;
                tmp2(idx) = tmp2(idx) + 1;
            end
            decayingRNAs = tmp2;
            
            %update numbers of RNAs
            this.RNAs = this.RNAs - decayingRNAs;
            
            %update metabolites
            %- water, hydrogen for hydrolysis
            %- nucleotide monophosphate salvage
            this.substrates = this.substrates + this.decayReactions' * decayingRNAs;
            
            this.rna.counts(:, this.compartment.cytosolIndexs) = this.RNAs;
            if any(any(this.rna.updateExternalState(-decayingRNAs, true)))
                throw(MException('RNADecay:error', 'All RNAs should have been degraded'));
            end
            this.RNAs = this.rna.counts(:, this.compartment.cytosolIndexs);
        end
    end

    %get methods of dependent local state
    methods
        function value = getDryWeight(this)
            if size(this.RNAs, 3) == 1
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + ...
                    this.rna.molecularWeights' * this.RNAs / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            else
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + ...
                    permute(this.rna.molecularWeights' * permute(this.RNAs,[1 3 2]),[1 3 2]) / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            end
        end
    end
end
