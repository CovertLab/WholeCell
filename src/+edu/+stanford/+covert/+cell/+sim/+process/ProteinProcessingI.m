%Protein Processing I
%
% @wholeCellModelID Process_ProteinProcessingI
% @name             Protein Processing I
% @description
%   Biology
%   ===========
%   Following translation, nascent peptides are deformylated, cleaved,
%   translocated, folded, and modified. First, peptide deformylase (MG_106)
%   deformylates the N-terminal formylmethionine of each nascent peptide.
%   Second, methionine aminopeptidase (MG_172) cleaves the N-terminal methionine
%   of 35 peptides. Second (see protein translocation process), 117 integral
%   membrane, lipoproteins, and extracellular proteins bind the SecA translocase
%   and are translocated into and through the plasma membrane via the
%   SecYEGDF-YidC pore. Third (see protein processing II process), diacylglyceryl
%   is transferred to the C-terminal cysteine of the signal sequence of each
%   lipoprotein by diacylglyceryl transferase, and the signal sequence of each
%   lipoprotein is cleaved by signal peptidase II. Next (see protein
%   folding process), 85 peptides bind inorganic ions at particular sites
%   and 64 peptides fold with the assistance of the chaperones and chaperonins
%   DnaJ, DnaK, GroEL, GroES, and GrpE. All peptides require trigger factor
%   to properly fold. Finally (see protein modification process), 20 peptide
%   species are modified at 63 sites by 3 enzymes – serine/threonine protein
%   kinase, lipoate ligase, and alpha glutamate ligase.
%
%   This process simulates
%   - N-terminal peptide deformylation, and
%   - N-terminal amino acid cleavage.
%
%   Knowledge Base
%   ===========
%   Every M. genitalium protein requires deformylation, and roughly 7%
%   require N-terminal amino acid cleavage. N-terminal methionine cleavages were
%   reconstructed by mapping N-terminal methionine cleavages observed in
%   Shewanella oneidensis MR-1 [PUB_0280] onto homologous M. genitalium genes.
%   The N-terminal methionine cleavage state of each protein monomer is
%   organized in the knowledge base, and loaded into the
%   nascentMonomerNTerminalMethionineCleavages property of this process by the
%   initializeConstants method.
%
%   Representation
%   ===========
%   substrates, enzymes, unprocessedMonomers, and processedMonomers represent
%   the counts of metabolites, the deformylase and methionine aminopeptidase,
%   and nascent and deformylated, cleaved protein monomers. The compartment
%   dimension of each of these properties has length 1. That is the process only
%   accesses the counts of these objects in the relevant compartments. This
%   process doesn't represent any additional intermediate processed states. This
%   process treats N-terminal deformylation and methionine cleavage as an
%   all-or-nothing event that either proceeds to complete within a single time
%   step, or does not occur at all.
%
%   nascentMonomerNTerminalMethionineCleavages is a boolean which represents
%   whether or not the N-terminal methionine of each protein monomer must be
%   cleaved.
%
%   Initialization
%   ===============
%   All protein monomers are initialized to the mature state. This is
%   implemented by the simulation class initializeState method.
%
%   Simulation
%   ===========
%   1. Compute the maximum number of peptides that can be processed based on the
%      availability of the two enzymes.
%   2. Randomly select peptides to be processed, weighted by the counts of each
%      protein species.
%   3. Update the counts of unprocessed and processed proteins monomers.
%      Decrement the counts of available enzyme activity.
%   4. Compute the maximum number of peptides which don't require cleavage that
%      can be processed based on the availability of peptide deformylase
%      activity.
%   5. Randomly select peptides to be deformylated, weighted by the counts of
%      each protein species.
%   6. Update the counts of unprocessed and processed proteins monomers.
%      Decrement the counts of available deformylase activity.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/30/2010
classdef ProteinProcessingI < edu.stanford.covert.cell.sim.Process
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
		    'deformylaseSpecificRate';
            'methionineAminoPeptidaseSpecificRate';
            'nascentMonomerNTerminalMethionineCleavages';            
		    };
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {   %names of simulation state properties redundant with timecourses in this or other processes or the simulation
            'unprocessedMonomers';
            'processedMonomers'};
    end

    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs = {};  %whole cell model IDs of stimuli

        substrateWholeCellModelIDs = {  %whole cell model IDs of substrates
            'H2O';
            'H';
            'MET';
            'FOR'};
        substrateIndexs_water      = 1; %index within substrates of water
        substrateIndexs_hydrogen   = 2; %index within substrates of hydrogen
        substrateIndexs_methionine = 3; %index within substrates of methionine
        substrateIndexs_formate    = 4; %index within substrates of formate

        enzymeWholeCellModelIDs={       %enzyme whole cell model ids
            'MG_106_DIMER';             %peptide deformylase
            'MG_172_MONOMER'};          %methionine aminopeptidase, type I
        enzymeIndexs_deformylase              = 1; %index within enzymes of peptide deformylase
        enzymeIndexs_methionineAminoPeptidase = 2; %index within enzymes of methionine aminopeptidase, type I

        unprocessedMonomerWholeCellModelIDs %whole cell model IDs of unprocessed protein monomers
        processedMonomerWholeCellModelIDs   %whole cell model IDs of processed protein monomers
    end
    
    %fixed biological constants
    properties
        deformylaseSpecificRate                     %number of reactions per second (38) [PUB_0021]
        methionineAminoPeptidaseSpecificRate        %number of reactions per second (6) [PUB_0026]
        nascentMonomerNTerminalMethionineCleavages  %whether N-terminal amino acid needs to be cleaved        
    end

    %global state (stored locally for convenience)
    properties
        unprocessedMonomers %counts of unprocessed protein monomers
        processedMonomers   %counts of processed protein monomers
    end

    %constructor
    methods
        function this = ProteinProcessingI(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.Process(wholeCellModelID, name);
        end
    end

    %communication between process/simulation
    methods
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)
            this.initializeConstants@edu.stanford.covert.cell.sim.Process(knowledgeBase, simulation, varargin{:});

            this.unprocessedMonomerWholeCellModelIDs = this.monomer.wholeCellModelIDs(this.monomer.nascentIndexs);
            this.processedMonomerWholeCellModelIDs   = this.monomer.wholeCellModelIDs(this.monomer.processedIIndexs);

            this.nascentMonomerNTerminalMethionineCleavages = knowledgeBase.proteinMonomerNTerminalMethionineCleavages;
        end

        %retrieve state from simulation
        function copyFromState(this)
            this.copyFromState@edu.stanford.covert.cell.sim.Process();

            this.unprocessedMonomers = this.monomer.counts(this.monomer.nascentIndexs,    this.compartment.cytosolIndexs, :);
            this.processedMonomers   = this.monomer.counts(this.monomer.processedIIndexs, this.compartment.cytosolIndexs, :);
        end

        %send state to simulation
        function copyToState(this)
            this.copyToState@edu.stanford.covert.cell.sim.Process();

            this.monomer.counts(this.monomer.nascentIndexs,    this.compartment.cytosolIndexs, :) = this.unprocessedMonomers;
            this.monomer.counts(this.monomer.processedIIndexs, this.compartment.cytosolIndexs, :) = this.processedMonomers;
        end
    end

    %allocate memory for state
    methods
        function allocateMemoryForState(this, numTimePoints)
            this.allocateMemoryForState@edu.stanford.covert.cell.sim.Process(numTimePoints);

            this.unprocessedMonomers = zeros(length(this.unprocessedMonomerWholeCellModelIDs), 1, numTimePoints);
            this.processedMonomers   = zeros(length(this.processedMonomerWholeCellModelIDs),   1, numTimePoints);
        end
    end

    %model
    methods
        %Calculate
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, ~, states)
            %% initialize
            bmProd = zeros(size(this.substrateWholeCellModelIDs));
            byProd = zeros(size(this.substrateWholeCellModelIDs));
            minEnzExp = zeros(size(this.enzymeWholeCellModelIDs));
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
            
            %% substrate and byproducts
            %deformylation
            %formyl-L-methionyl peptide + H2O <==> methionyl peptide + formate + H+
            bmProd(this.substrateIndexs_water)      = bmProd(this.substrateIndexs_water)    + sum(states.monomerProductions);
            byProd(this.substrateIndexs_formate)    = byProd(this.substrateIndexs_formate)  + sum(states.monomerProductions);
            byProd(this.substrateIndexs_hydrogen)   = byProd(this.substrateIndexs_hydrogen) + sum(states.monomerProductions);
            
            %N-terminal methionine cleavage
            bmProd(this.substrateIndexs_water)      = bmProd(this.substrateIndexs_water)      + states.monomerProductions' * this.nascentMonomerNTerminalMethionineCleavages;
            byProd(this.substrateIndexs_methionine) = byProd(this.substrateIndexs_methionine) + states.monomerProductions' * this.nascentMonomerNTerminalMethionineCleavages;
            
            %% enzymes
            %deformylation
            minEnzExp(this.enzymeIndexs_deformylase) = 2 * sum(states.monomerProductions0) / this.deformylaseSpecificRate;
            
            %N-terminal methionine cleavage
            minEnzExp(this.enzymeIndexs_methionineAminoPeptidase) = ...
                2 * (states.monomerProductions0' * this.nascentMonomerNTerminalMethionineCleavages) / ...
                this.methionineAminoPeptidaseSpecificRate;
        end

        %initialization: monomers initialized to mature/bound/inactivated state
        %by simulation initializeState method
        function initializeState(~)
        end

        %resource requirements
        function result = calcResourceRequirements_Current(this)
            %Numbers of enzymes
            deformylaseLimit = ceil(this.enzymes(this.enzymeIndexs_deformylase) * this.deformylaseSpecificRate * this.stepSizeSec);
            cleavageLimit = ceil(this.enzymes(this.enzymeIndexs_methionineAminoPeptidase) * ...
                this.methionineAminoPeptidaseSpecificRate * this.stepSizeSec);
            
            %initialize
            result = zeros(size(this.substrates));
            
            %upper bound on metabolites
            result(this.substrateIndexs_water) = min(...
                + min(deformylaseLimit, sum(this.unprocessedMonomers)) ...                                          %deformylation: formyl-L-methionyl peptide + H2O <==> methionyl peptide + formate + H+
                + min(cleavageLimit, this.nascentMonomerNTerminalMethionineCleavages' * this.unprocessedMonomers)); %N-terminal methionine cleavage
        end

        %simulation
        %- Deformylation: formyl-L-methionyl peptide + H2O <==> methionyl peptide + formate + H+
        %- N-terminal methionine cleavage
        function evolveState(this)
            %Numbers of enzymes
            deformylaseLimit = this.enzymes(this.enzymeIndexs_deformylase) * this.deformylaseSpecificRate * this.stepSizeSec;
            cleavageLimit = this.enzymes(this.enzymeIndexs_methionineAminoPeptidase) * ...
                this.methionineAminoPeptidaseSpecificRate * this.stepSizeSec;
            
            %Simulate deformylation and cleavage
            if cleavageLimit > 0 && any(this.unprocessedMonomers(this.nascentMonomerNTerminalMethionineCleavages))
                %compute maximum of transformations that can occur, limited by
                %- unprocessed protein monomers
                %- enzymes (kinetics, availability)
                transformations = this.unprocessedMonomers;
                transformations = transformations * ...
                    min(1, deformylaseLimit / sum(transformations));
                transformations(this.nascentMonomerNTerminalMethionineCleavages) = ...
                    transformations(this.nascentMonomerNTerminalMethionineCleavages) * ...
                    min(1, cleavageLimit / sum(transformations(this.nascentMonomerNTerminalMethionineCleavages)));
                transformations = this.randStream.stochasticRound(transformations);
                
                %water availability
                if sum(transformations) + sum(transformations(this.nascentMonomerNTerminalMethionineCleavages)) > ...
                       this.substrates(this.substrateIndexs_water)
                   water_cleavage = this.randStream.stochasticRound(this.substrates(this.substrateIndexs_water) * ...
                       2 * sum(transformations(this.nascentMonomerNTerminalMethionineCleavages)) / ...
                       (sum(transformations) + sum(transformations(this.nascentMonomerNTerminalMethionineCleavages))));
                   water_deformylations = this.substrates(this.substrateIndexs_water) - water_cleavage;
                   if any(transformations(this.nascentMonomerNTerminalMethionineCleavages))
                       transformations(this.nascentMonomerNTerminalMethionineCleavages) = ...
                           min(transformations(this.nascentMonomerNTerminalMethionineCleavages), ...
                           this.randStream.mnrnd(floor(water_cleavage / 2), ...
                           transformations(this.nascentMonomerNTerminalMethionineCleavages) / ...
                           sum(transformations(this.nascentMonomerNTerminalMethionineCleavages)))');
                   end
                   if any(transformations(~this.nascentMonomerNTerminalMethionineCleavages))
                       transformations(~this.nascentMonomerNTerminalMethionineCleavages) = ...
                           min(transformations(~this.nascentMonomerNTerminalMethionineCleavages), ...
                           this.randStream.mnrnd(water_deformylations, ...
                           transformations(~this.nascentMonomerNTerminalMethionineCleavages) / ...
                           sum(transformations(~this.nascentMonomerNTerminalMethionineCleavages)))');
                   end
                end
                
                %update numbers of monomers
                this.processedMonomers   = this.processedMonomers   + transformations;
                this.unprocessedMonomers = this.unprocessedMonomers - transformations;
                
                %update enzyme availability
                deformylaseLimit = max(0, deformylaseLimit - sum(transformations));
                
                %update metabolites
                this.substrates([this.substrateIndexs_water; this.substrateIndexs_formate; this.substrateIndexs_hydrogen]) = ...
                    this.substrates([this.substrateIndexs_water; this.substrateIndexs_formate; this.substrateIndexs_hydrogen]) + ...
                    [-1; 1; 1] * sum(transformations);
                
                this.substrates([this.substrateIndexs_water; this.substrateIndexs_methionine]) = ...
                    this.substrates([this.substrateIndexs_water; this.substrateIndexs_methionine]) + ...
                    [-1; 1] * sum(transformations(this.nascentMonomerNTerminalMethionineCleavages));
            end

            %% Continue to simulate deformylation only

            %compute maximum of transformations that can occur, limited by
            %- unprocessed protein monomers
            %- enzymes (kinetics, availability)
            transformations = this.unprocessedMonomers;
            transformations(this.nascentMonomerNTerminalMethionineCleavages) = 0;
            if deformylaseLimit > 0 && any(transformations)
                transformations = transformations * ...
                    min(1, deformylaseLimit / sum(transformations));

                transformations = this.randStream.stochasticRound(transformations);
                if sum(transformations) > this.substrates(this.substrateIndexs_water)
                    transformations = min(transformations, this.randStream.mnrnd(this.substrates(this.substrateIndexs_water), transformations / sum(transformations))');
                end

                %update numbers of monomers
                this.processedMonomers   = this.processedMonomers   + transformations;
                this.unprocessedMonomers = this.unprocessedMonomers - transformations;

                %update metabolites
                this.substrates([this.substrateIndexs_water; this.substrateIndexs_formate; this.substrateIndexs_hydrogen]) = ...
                    this.substrates([this.substrateIndexs_water; this.substrateIndexs_formate; this.substrateIndexs_hydrogen]) + ...
                    [-1; 1; 1] * sum(transformations);
            end
        end
    end

    %get methods of dependent local state
    methods
        function value = getDryWeight(this)
            if size(this.unprocessedMonomers, 3) == 1
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + (...
                    this.monomer.molecularWeights(this.monomer.nascentIndexs)'    * this.unprocessedMonomers + ...
                    this.monomer.molecularWeights(this.monomer.processedIIndexs)' * this.processedMonomers) ...
                    / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            else
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + (...
                    permute(this.monomer.molecularWeights(this.monomer.nascentIndexs)'    * permute(this.unprocessedMonomers,[1 3 2]),[1 3 2]) + ...
                    permute(this.monomer.molecularWeights(this.monomer.processedIIndexs)' * permute(this.processedMonomers,  [1 3 2]),[1 3 2])) ...
                    / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            end
        end
    end
end
