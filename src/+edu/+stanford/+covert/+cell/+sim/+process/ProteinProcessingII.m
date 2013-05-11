%Protein Processing II
%
% @wholeCellModelID Process_ProteinProcessingII
% @name             Protein Processing II
% @description
%   Biology
%   ===============
%   Proteins are produced in the cytoplasm, and integral membrane proteins,
%   lipoproteins, and extracellular proteins must be translocated into and
%   through the cell membrane to reach their intended destination. To be
%   recognized by the translocation machinery lipoprotein and secreted proteins
%   contain type II N-terminal signal sequences, which are typically 10-15 amino
%   acids long and positively charged. Following translocation lipoproteins must
%   be anchored to the outter leaflet of the cell membrane by the addition of
%   diacylglyceryl by diacylglyceryl transferase (Lgt, MG_086) to what will
%   become the C-terminal cysteine, and their signal sequence must be cleaved by
%   signal peptidase II (LspA, MG_210) at lipoboxes (L[ASI][GA]C). M. genitalium
%   does not contain an apolipoprotein transacylase [PUB_0655, PUB_0656,
%   PUB_0657]. The signal peptides of secreted proteins is similarly cleaved.
%
%   This process simulates lipoprotein and secreted protein maturation:
%   - diacylglyceryl transfer – Lgt (MG_086) catalyzes transfer of
%     diacylglycerol group to sulfhydryl group of lipobox cysteine
%     [PUB_0654].
%   - signal peptide cleavage – LspA (MG_210) cleaves lipoprotein at
%     lipobox cysteine [PUB_0654].
%
%   Knowledge Base
%   ===============
%   The localization of each protein, and the signal peptide type and length of
%   lipoproteins and secreted proteins was compiled from several sources (see
%   protein translocation process). The information was organized in the
%   knowledge base, and is the type of each signal peptide is encoded in the
%   this process's lipoproteinMonomerIndexs, secretedMonomerIndexs, and
%   unprocessedMonomerIndexs properties by initializeConstants.
%
%   Representation
%   ===============
%   The substrates and enzymes represents the counts of available metabolites
%   and diacylglyceryl transferase and signal peptidase enzymes.
%   unprocessedMonomers, processedMonomers, and signalSequenceMonomers represent
%   the counts of unanchored, uncleaved protein monomers; anchored, cleaved
%   protein monomers; and the separted, free signal sequences. This process does
%   not represent any intermediate states in the anchoring and cleavage of lipo-
%   and secreted proteins. This process consideres anchoring and cleavage to be
%   an all-or-nothing event.
%
%   The lipoproteinMonomerIndexs and secretedMonomerIndexs properties represent
%   the indices of lipo-and secreted proteins and their released signal
%   sequences within unprocessedMonomers, processedMonomers, and
%   signalSequenceMonomers. unprocessedMonomerIndexs indicating the indices of
%   non-lipo-, non-secreted proteins and their released signal sequences (not
%   used in simulation; only allocated for convenience and parallelism) within
%   unprocessedMonomers, processedMonomers, and signalSequenceMonomers.
%
%   Initialization
%   ===============
%   All protein monomers are initialized to the mature state. This is
%   accomplished by the simulation class initializeState method.
%
%   Simulation
%   ===============
%   1. Compute the maximum number of peptides that can be processed based on the
%      availability of metabolites, and of the two enzymes.
%   2. Randomly select peptides to be processed, weighted by the counts of each
%      protein species.
%   3. Update the counts of unprocessed and processed proteins monomers.
%      Decrement the counts of available metabolites and enzyme activity.
%   4. Compute the maximum number of peptides which don't require anchoring that
%      can be processed based on the availability of signal peptidase II
%      activity.
%   5. Randomly select peptides to be cleaved, weighted by the counts of
%      each protein species.
%   6. Update the counts of unprocessed and processed proteins monomers.
%      Decrement the counts of available signal peptidase II activity.
%   7. Transition proteins which don't requiring anchoring or cleavage (eg.
%      cytosolic and integral membrane proteins). That is set processedMonomers
%      equal to its sum with unprocessedMonomers for these monomers, and set
%      unprocessedMonomers to zero for these monomers.
%
%   References
%   ===============
%   1. Chambaud I, Wróblewski H, Blanchard A (1999). Interactions between
%      mycoplasma lipoproteins and the host immune system. Trends
%      Microbiol. 7(12): 493-9. [PUB_0654].
%   2. Chambaud I, Heilig R, Ferris S, Barbe V, Samson D, Galisson F,
%      Moszer I, Dybvig K, Wroblewski H, Viari A, Rocha EP, Blanchard A
%      (2001). The complete genome sequence of the murine respiratory
%      pathogen. Mycoplasma pulmonis. Nucleic Acids Res. 29(10):2145-53.
%      [PUB_0655]
%   3. Muhlradt PF, Kiess M, Meyer H, Süssmuth R, Jung G (1997). Isolation,
%      structure elucidation, and synthesis of a macrophage stimulatory
%      lipopeptide from Mycoplasma fermentans acting at picomolar
%      concentration. J Exp Med. 185(11):1951-8. [PUB_0656].
%   4. Piec G, Mirkovitch J, Palacio S, Muhlradt PF, Felix R (1999). Effect
%      of MALP-2, a lipopeptide from Mycoplasma fermentans, on bone
%      resorption in vitro. Infect Immun. 67(12):6281-5. [PUB_0657].
%   5. Sankaran K, Wu HC (1994). Lipid modification of bacterial
%      prolipoprotein. Transfer of diacylglyceryl moiety from
%      phosphatidylglycerol. J Biol Chem. 269(31):19701-6. [PUB_0266]
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/9/2010
classdef  ProteinProcessingII < edu.stanford.covert.cell.sim.Process
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
		    'lipoproteinSignalPeptidaseSpecificRate';
            'lipoproteinDiacylglycerylTransferaseSpecificRate';
		    };		
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {   %names of simulation state properties redundant with timecourses in this or other processes or the simulation
            'unprocessedMonomers';
            'processedMonomers';
            'signalSequenceMonomers'};
    end

    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs = {}; %whole cell model IDs of stimuli

        substrateWholeCellModelIDs = { %whole cell model IDs of substrates
            'H2O';
            'diacylglycerolCys';
            'PG160';
            'SNGLYP';
            'H'};
        substrateIndexs_water             = 1; %index within substrates of water
        substrateIndexs_diacylglycerolCys = 2; %index within substrates of diacylglycerol cysteine
        substrateIndexs_PG160             = 3; %index within substrates of PG160
        substrateIndexs_SNGLYP            = 4; %index within substrates of sn-glycerol-1-phosphate
        substrateIndexs_hydrogen          = 5; %index within substrates of hydrogen

        enzymeWholeCellModelIDs = { %enzyme whole cell model ids
            'MG_210_MONOMER'        %prolipoprotein signal peptidase, signal peptidase II
            'MG_086_MONOMER'        %prolipoprotein diacylglyceryl transferase
            };
        enzymeIndexs_signalPeptidase           = 1; %index within enzymes of prolipoprotein signal peptidase, signal peptidase II
        enzymeIndexs_diacylglycerylTransferase = 2; %index within enzymes of prolipoprotein diacylglyceryl transferase

        monomerCompartments                    %compartment indices (within compartment.wholeCellModelIDs) of monomers

        lipoproteinMonomerIndexs               %indices within unprocessedMonomers of lipoproteins
        secretedMonomerIndexs                  %indices within unprocessedMonomers of secreted proteins
        unprocessedMonomerIndexs               %indices within unprocessedMonomers of non-lipoproteins, non-secreted proteins

        unprocessedMonomerWholeCellModelIDs    %whole cell model IDs of unprocessed monomers
        processedMonomerWholeCellModelIDs      %whole cell model IDs of processed monomers
        signalSequenceMonomerWholeCellModelIDs %whole cell model IDs of signal sequences
    end
    
    %fixed biological constants
    properties
        lipoproteinSignalPeptidaseSpecificRate             %number of reactions per second (11) [PUB_0008]
        lipoproteinDiacylglycerylTransferaseSpecificRate   %number of reactions per second (0.0165) [PUB_0157, PUB_0158]
    end

    %global state (stored locally for convenience)
    properties
        unprocessedMonomers     %counts of unprocessed protein monomers [unprocessed X 1]
        processedMonomers       %counts of processed protein monomers [processed X 1]
        signalSequenceMonomers  %counts of cleaved signal sequences [signal sequence X 1]
    end

    %constructor
    methods
        function this = ProteinProcessingII(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.Process(wholeCellModelID, name);
        end
    end

    %communication between process/simulation
    methods
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)
            this.initializeConstants@edu.stanford.covert.cell.sim.Process(...
                knowledgeBase, simulation, varargin{:});

            this.unprocessedMonomerWholeCellModelIDs    = this.monomer.wholeCellModelIDs(this.monomer.processedIIndexs);
            this.processedMonomerWholeCellModelIDs      = this.monomer.wholeCellModelIDs(this.monomer.processedIIIndexs);
            this.signalSequenceMonomerWholeCellModelIDs = this.monomer.wholeCellModelIDs(this.monomer.signalSequenceIndexs);

            this.monomerCompartments = double(knowledgeBase.proteinMonomerCompartments);
            this.monomerCompartments(this.monomerCompartments == this.compartment.terminalOrganelleCytosolIndexs)  = this.compartment.cytosolIndexs;
            this.monomerCompartments(this.monomerCompartments == this.compartment.terminalOrganelleMembraneIndexs) = this.compartment.membraneIndexs;

            this.lipoproteinMonomerIndexs = find(strcmp(...
                {knowledgeBase.proteinMonomers.signalSequenceType},'lipoprotein'))';
            this.secretedMonomerIndexs = find(strcmp(...
                {knowledgeBase.proteinMonomers.signalSequenceType},'secretory'))';
            this.unprocessedMonomerIndexs = setdiff(...
                (1:knowledgeBase.numProteinMonomers)',...
                [this.lipoproteinMonomerIndexs;this.secretedMonomerIndexs]);
        end

        %retrieve state from simulation
        function copyFromState(this)
            this.copyFromState@edu.stanford.covert.cell.sim.Process();

            numTimePoints = size(this.unprocessedMonomers, 3);
            
            if numTimePoints == 1
                this.unprocessedMonomers = this.monomer.counts(sub2ind(...
                    size(this.monomer.counts), ...
                    this.monomer.processedIIndexs, ...
                    this.monomerCompartments));
                this.processedMonomers = this.monomer.counts(sub2ind(...
                    size(this.monomer.counts), ...
                    this.monomer.processedIIIndexs, ...
                    this.monomerCompartments));
                this.signalSequenceMonomers = this.monomer.counts(...
                    this.monomer.signalSequenceIndexs, ...
                    this.compartment.cytosolIndexs);
            else
                this.unprocessedMonomers = this.monomer.counts(sub2ind(...
                    size(this.monomer.counts), ...
                    repmat(this.monomer.processedIIndexs, [1 1 numTimePoints]), ...
                    repmat(this.monomerCompartments,[1 1 numTimePoints]),...
                    permute(repmat((1:numTimePoints)',[1 size(this.processedMonomers,1) 1]),[2 3 1])));
                this.processedMonomers = this.monomer.counts(sub2ind(...
                    size(this.monomer.counts), ...
                    repmat(this.monomer.processedIIIndexs, [1 1 numTimePoints]), ...
                    repmat(this.monomerCompartments,[1 1 numTimePoints]),...
                    permute(repmat((1:numTimePoints)',[1 size(this.processedMonomers,1) 1]),[2 3 1])));
                this.signalSequenceMonomers = this.monomer.counts(...
                    this.monomer.signalSequenceIndexs, ...
                    this.compartment.cytosolIndexs, :);
            end
        end

        %send state to simulation
        function copyToState(this)
            this.copyToState@edu.stanford.covert.cell.sim.Process();

            numTimePoints = size(this.unprocessedMonomers, 3);

            if numTimePoints == 1
                this.monomer.counts(sub2ind(...
                    size(this.monomer.counts), ...
                    this.monomer.processedIIndexs, ...
                    this.monomerCompartments)) = ...
                    this.unprocessedMonomers;
                this.monomer.counts(sub2ind(...
                    size(this.monomer.counts), ...
                    this.monomer.processedIIIndexs, ...
                    this.monomerCompartments)) = ...
                    this.processedMonomers;
                this.monomer.counts(this.monomer.signalSequenceIndexs, ...
                    this.compartment.cytosolIndexs) = ...
                    this.signalSequenceMonomers;
            else
                this.monomer.counts(sub2ind(...
                    size(this.monomer.counts), ...
                    repmat(this.monomer.processedIIndexs, [1 1 numTimePoints]), ...
                    repmat(this.monomerCompartments,[1 1 numTimePoints]),...
                    permute(repmat((1:numTimePoints)',[1 size(this.unprocessedMonomers,1) 1]),[2 3 1]))) = ...
                    this.unprocessedMonomers;
                this.monomer.counts(sub2ind(...
                    size(this.monomer.counts), ...
                    repmat(this.monomer.processedIIIndexs, [1 1 numTimePoints]), ...
                    repmat(this.monomerCompartments,[1 1 numTimePoints]),...
                    permute(repmat((1:numTimePoints)',[1 size(this.processedMonomers,1) 1]),[2 3 1]))) = ...
                    this.processedMonomers;
                this.monomer.counts(this.monomer.signalSequenceIndexs, ...
                    this.compartment.cytosolIndexs, :) = ...
                    this.signalSequenceMonomers;
            end
        end
    end

    %allocate memory for state
    methods
        function allocateMemoryForState(this, numTimePoints)
            this.allocateMemoryForState@edu.stanford.covert.cell.sim.Process(numTimePoints);

            numCompartments = this.compartment.count;
            this.unprocessedMonomers    = zeros(length(this.unprocessedMonomerWholeCellModelIDs),    numCompartments, numTimePoints);
            this.processedMonomers      = zeros(length(this.processedMonomerWholeCellModelIDs),      numCompartments, numTimePoints);
            this.signalSequenceMonomers = zeros(length(this.signalSequenceMonomerWholeCellModelIDs), 1, numTimePoints);
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
            cleavages = sum(states.monomerProductions([this.lipoproteinMonomerIndexs; this.secretedMonomerIndexs]));
            transfers = sum(states.monomerProductions(this.lipoproteinMonomerIndexs));

            %signal peptide cleavage
            bmProd(this.substrateIndexs_water) = bmProd(this.substrateIndexs_water) + cleavages;

            %diacylglyceryl transfer
            %CYS[c] + PG160[m] ==> diacylglycerolCys[c] + SNGLYP[c] + H[c]
            bmProd(this.substrateIndexs_PG160)    = bmProd(this.substrateIndexs_PG160)    + transfers;
            byProd(this.substrateIndexs_SNGLYP)   = byProd(this.substrateIndexs_SNGLYP)   + transfers;
            byProd(this.substrateIndexs_hydrogen) = byProd(this.substrateIndexs_hydrogen) + transfers;
            
            %% enzymes
            minEnzExp(this.enzymeIndexs_signalPeptidase) = ...
                2 * sum(states.monomerProductions0([this.lipoproteinMonomerIndexs; this.secretedMonomerIndexs])) / ...
                this.lipoproteinSignalPeptidaseSpecificRate;
            minEnzExp(this.enzymeIndexs_diacylglycerylTransferase) = ...
                2 * sum(states.monomerProductions0(this.lipoproteinMonomerIndexs)) / ...
                this.lipoproteinDiacylglycerylTransferaseSpecificRate;
        end
        
        %initialization: monomers initialized to mature/bound/inactivated state
        %by simulation initializeState method
        function initializeState(~)
        end

        %resource requirements
        function result = calcResourceRequirements_Current(this)
            %Enzyme availability
            peptidaseLimit = ceil(this.enzymes(this.enzymeIndexs_signalPeptidase) * ...
                this.lipoproteinSignalPeptidaseSpecificRate * ...
                this.stepSizeSec);
            transferaseLimit = ceil(this.enzymes(this.enzymeIndexs_diacylglycerylTransferase) * ...
                this.lipoproteinDiacylglycerylTransferaseSpecificRate * ...
                this.stepSizeSec);
            
            %initialize
            result = zeros(size(this.substrates));
            
            %signal peptide cleavage
            result(this.substrateIndexs_water) = ...
                min(peptidaseLimit, sum(this.unprocessedMonomers([this.lipoproteinMonomerIndexs; this.secretedMonomerIndexs])));
            
            %diacylglyceryl transfer
            %CYS[c] + PG160[m] ==> diacylglycerolCys[c] + SNGLYP[c] + H[c]
            result(this.substrateIndexs_PG160) = ...
                min(transferaseLimit, sum(this.unprocessedMonomers(this.lipoproteinMonomerIndexs)));
        end

        %simulation
        function evolveState(this)            
            %pass along monomers that don't need to be processed
            this.processedMonomers(this.unprocessedMonomerIndexs) = ...
                this.processedMonomers(this.unprocessedMonomerIndexs) + ...
                this.unprocessedMonomers(this.unprocessedMonomerIndexs);
            this.unprocessedMonomers(this.unprocessedMonomerIndexs) = 0;
            
            %stop early if no proteins need to be processed
            if ~any(this.unprocessedMonomers([this.lipoproteinMonomerIndexs; this.secretedMonomerIndexs]))
                return;
            end
            
            %Indices of transformed monomers
            peptidaseIndexs = [this.lipoproteinMonomerIndexs; this.secretedMonomerIndexs];
            transferaseIndexs = this.lipoproteinMonomerIndexs;
            
            %Enzyme availability
            peptidaseLimit = this.enzymes(this.enzymeIndexs_signalPeptidase) * ...
                this.lipoproteinSignalPeptidaseSpecificRate * ...
                this.stepSizeSec;
            
            %% Simulate signal peptide cleavage and diacylglyceryl transfer
            % - lipoproteins or secreted
            if any(this.unprocessedMonomers(transferaseIndexs))
                %Enzyme availability
                transferaseLimit = this.enzymes(this.enzymeIndexs_diacylglycerylTransferase) * ...
                    this.lipoproteinDiacylglycerylTransferaseSpecificRate * ...
                    this.stepSizeSec;
                
                %compute maximum of transformations that can occur, limited by
                %- unprocessed protein monomers
                %- enzymes (kinetics, availability)
                %- metabolites (PG160, water)
                transformations = this.unprocessedMonomers;
                transformations(peptidaseIndexs) = transformations(peptidaseIndexs) * ...
                    min([1, peptidaseLimit / sum(transformations(peptidaseIndexs))]);
                transformations(transferaseIndexs) = transformations(transferaseIndexs) * ...
                    min([1, transferaseLimit / sum(transformations(transferaseIndexs))]);
                transformations = this.randStream.stochasticRound(transformations);
                if sum(transformations(peptidaseIndexs)) > this.substrates(this.substrateIndexs_water)
                   transformations(peptidaseIndexs) = min(transformations(peptidaseIndexs), ...
                       this.randStream.mnrnd(this.substrates(this.substrateIndexs_water), ...
                       transformations(peptidaseIndexs) / sum(transformations(peptidaseIndexs)))');
                end
                if sum(transformations(transferaseIndexs)) > this.substrates(this.substrateIndexs_PG160)
                    transformations(transferaseIndexs) = min(transformations(transferaseIndexs), ...
                        this.randStream.mnrnd(this.substrates(this.substrateIndexs_PG160), ...
                        transformations(transferaseIndexs) / sum(transformations(transferaseIndexs)))');
                end
                
                %update numbers of monomers
                this.signalSequenceMonomers(peptidaseIndexs) = ...
                    this.signalSequenceMonomers(peptidaseIndexs) + ...
                    transformations(peptidaseIndexs);
                this.processedMonomers   = this.processedMonomers   + transformations;
                this.unprocessedMonomers = this.unprocessedMonomers - transformations;
                
                %update enzyme availability
                peptidaseLimit = peptidaseLimit - sum(transformations(peptidaseIndexs));
                
                %update metabolites
                this.substrates(this.substrateIndexs_water) = ...
                    this.substrates(this.substrateIndexs_water) - ...
                    sum(transformations(peptidaseIndexs));
                this.substrates([this.substrateIndexs_PG160; this.substrateIndexs_SNGLYP; this.substrateIndexs_hydrogen]) = ...
                    this.substrates([this.substrateIndexs_PG160; this.substrateIndexs_SNGLYP; this.substrateIndexs_hydrogen]) + ...
                    [-1;1;1] * sum(transformations(transferaseIndexs));
            end
            
            %% Continue to simulate signal peptide cleavage - extracellular proteins only
            transformations = this.unprocessedMonomers;
            transformations(peptidaseIndexs) = transformations(peptidaseIndexs) * ...
                min([1, peptidaseLimit / sum(transformations(peptidaseIndexs))]);
            transformations(transferaseIndexs) = 0;
            if any(transformations)
                %compute maximum of transformations that can occur, limited by
                %- unprocessed protein monomers
                %- enzymes (kinetics, availability)
                %- metabolites (PG160)                
                transformations = this.randStream.stochasticRound(transformations);
                if sum(transformations(peptidaseIndexs)) > this.substrates(this.substrateIndexs_water)
                   transformations(peptidaseIndexs) = min(transformations(peptidaseIndexs), ...
                       this.randStream.mnrnd(this.substrates(this.substrateIndexs_water), ...
                       transformations(peptidaseIndexs) / sum(transformations(peptidaseIndexs)))');
                end
                
                %update numbers of monomers
                this.signalSequenceMonomers(peptidaseIndexs) = ...
                    this.signalSequenceMonomers(peptidaseIndexs) + ...
                    transformations(peptidaseIndexs);
                this.processedMonomers   = this.processedMonomers   + transformations;
                this.unprocessedMonomers = this.unprocessedMonomers - transformations;
                
                %update metabolites
                this.substrates(this.substrateIndexs_water) = ...
                    this.substrates(this.substrateIndexs_water) - ...
                    sum(transformations(peptidaseIndexs));
            end
        end
    end

    %get methods of dependent local state
    methods
        function value = getDryWeight(this)
            if size(this.unprocessedMonomers, 3) == 1
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + (...
                    this.monomer.molecularWeights(this.monomer.processedIIndexs)'     * this.unprocessedMonomers + ...
                    this.monomer.molecularWeights(this.monomer.processedIIIndexs)'    * this.processedMonomers   + ...
                    this.monomer.molecularWeights(this.monomer.signalSequenceIndexs)' * this.signalSequenceMonomers) ...
                    / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            else
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + (...
                    permute(this.monomer.molecularWeights(this.monomer.processedIIndexs)'     * permute(this.unprocessedMonomers,   [1 3 2]), [1 3 2]) + ...
                    permute(this.monomer.molecularWeights(this.monomer.processedIIIndexs)'    * permute(this.processedMonomers,     [1 3 2]), [1 3 2]) + ...
                    permute(this.monomer.molecularWeights(this.monomer.signalSequenceIndexs)' * permute(this.signalSequenceMonomers,[1 3 2]), [1 3 2])) ...
                    / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            end
        end
    end
end
