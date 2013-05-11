%Transcription
%
% @wholeCellModelID Process_Transcription
% @name             Transcription
% @description
%   Biology
%   ===============
%   Transcription is the first step in the synthesis of functional gene
%   products where RNA polymerase and several accessory enzymes translate
%   transcription units, or regions of the DNA containing 1 or more genes, into
%   RNA molecules. Following transcription the RNA molecules follow one of
%   two pathways:
%   - mRNAs are used as templates for translation (see translation process)
%   - r/s/tRNAs which are transcribed as molecules containing multiple genes are
%     cleavaged into their individual genes (see RNA processing process), are
%     modified at several bases (see RNA modification process) to improve their
%     stability and enhance their catalytic activity, and finally act as
%     ribozymes (r/sRNAs) or as adaptors between mRNAs and the amino acids they
%     code for (tRNAs).
%
%   Transcription begins with the recruitment of RNA polymerase to a promoter
%   with the help of the sigma initiation factor and possiblity transcription
%   factors. Next elongation factors are recruited, RNA begins to be
%   polymerized, and sigma factor is released. Finally the RNA polymerase
%   reaches  a terminator at the the end of the transcription unit, and
%   with the help of termination factors releases the polymerized RNA and
%   dissociates from the DNA. Termination in E. coli occurs via either a
%   Rho-dependent (50%) or Rho-independent mechanism (50%). Rho-dependent
%   termination is catalyzed by the hexameric ATP-dependent helicase Rho.
%   Rho is not essential in B. subtilis [PUB_0234]. Rho-independent
%   termination occur via the intrinsic properties RNA which disrupt RNA
%   polymerase-DNA binding. Terminator hairpins are not predicted in M.
%   genitalium. In E. coli termination is incompetition with
%   antitermination. However antitermination has not been reported in any
%   mcyoplasma [PUB_0182].
%
%   As soon as the RNA begins to polymerized, even prior to termination, the
%   mRNA transcripts may be bound by ribosomes and polymerized. For simplicity,
%   our model doesn't represent this phenomenon.
%
%   Knowledge Base
%   ===============
%   The transcription unit structure was compiled from several sources:
%   - Primary reports of cotranscribed genes [PUB_0176, PUB_0182, PUB_0186,
%     PUB_0188, PUB_0188, PUB_0244, PUB_0247, PUB_0248, PUB_0249, PUB_0251]
%   - OperonDB database of cotranscribed genes [PUB_0250]
%   - Conservation of gene order across multiple species
%   - Related function of adjacent genes
%   - Expression levels as measured by microarrays of adjacent genes [PUB_0569]
%   - Strandedness of adjacent genes
%   - Weiner, Hermann, and Browning computational model of promoters and
%     transcription unit start sites [PUB_0411]
%
%   The transcription unit structure is organized in the knowledge base, and is
%   loaded into this process by the initializeConstants method.
%
%   The knowledge base also contains the measured expression and half-lives of
%   many transcripts. These values are loaded by initializeConstants and fit by
%   simulation.fitConstants to be consistent with other processes.
%
%   Representation
%   ===============
%   substrates, enzymes, boundEnzymes, and RNAs represent the counts of
%   metabolites, free transcription enzymes, transcription enzymes bound to RNAs
%   and RNA polymerases, and nascent RNAs.
%
%   rnaPolymerases.states represents the current state / pseudostate (actively
%   transcribing, specifically bound, non-specifically bound, free,
%   non-existent) of each RNA polymerase, where each state is indicated by the
%   enumeration:
%   - rnaPolymeraseActivelyTranscribingValue
%   - rnaPolymeraseSpecificallyBoundValue
%   - rnaPolymeraseNonSpecificallyBoundValue
%   - rnaPolymeraseFreeValue
%   - rnaPolymeraseNotExistValue (state exists as a way to account for memory
%     allocated for future RNA polymerases)
%   For actively transcribing polymerases rnaPolymerases.states also represents
%   the position of the polymerase along the transcription unit.
%
%   That is entries of rnaPolymerases.states with the following values
%   corresponding to these states:
%     >= RNAPolymerases.activelyTranscribingValue: RNA polymerases position on genome actively transcribing
%     == RNAPolymerases.specificallyBoundValue:    RNA polymerase specifically bound
%     == RNAPolymerases.nonSpecificallyBoundValue: RNA polymerase non-specifically bound
%     == RNAPolymerases.freeValue:                 RNA polymerase free
%     == RNAPolymerases.notExistValue:             RNA polymerase doesn't exist
%
%   transcripts.boundTranscriptionUnits represents the particular transcription
%   unit to which each actively transcribing and specifically bound polymerase
%   is bound.
%
%   rnaPolymeraseStateExpectations represents the expected occupancies of the
%   RNA polymerase states.
%
%   transcriptionUnitBindingProbabilities represents the relative affinity of
%   RNA polymerases for the promoters of each transcription unit.
%   transcriptionFactorBindingProbFoldChange represents the fold change affect of
%   transcription factors on the relative affinities of the RNA polymerse for
%   the promoters. RNA polymerases are assigned to
%   transcription units weighted by the product of
%   transcriptionUnitBindingProbabilities and
%   transcriptionFactorBindingProbFoldChange.
%
%   Initialization
%   ===============
%   All RNAs are initialized to the mature state. This is implemented by the
%   simulation class initializeState method.
%
%   RNA polymerases are initialized to their steady state:
%   - Each RNA polymerase is randomly assigned (with replacement) to one of the
%     actively transcribing, specifically bound, non-specifically bound, or free
%     states weighted by the expected occupancy of each state
%     (rnaPolymeraseStateExpectations)
%   - Actively transcribing and specifically bound polymerases randomly assigned
%     to transcription units weighted by their transcription rates
%     (transcriptionUnitBindingProbabilities)
%   - Each transcription unit to which an actively transcribing polymerase has
%     been assigned is divided into 1 segment for each polymerase
%   - Actively transcribing polymerases randomly assigned to positions within
%     the assigned segment of their assigned transcription unit (positions near
%     the segment border are not allowed to prevent polymerases from being too
%     close to each other) with uniform probably.
%
%   Simulation
%   ===============
%   Evolves the state of RNA polymerase using a markov chain model with four
%   states:
%   - actively translating
%   - specifically bound
%   - non-specifically bound
%   - free
%
%   Transition probabilities are designed to maintain the occupancy of each
%   state within a narrow window around their expected values. Transition
%   probability are determined by four logistic control functions. These can
%   be tuned with the constants
%   - rnaPolymeraseStateExpectations
%
%   RNA polymerase are created in the free state.
%
%   Actively transcribing state:
%   1. Release sigma factor if after first second of elongation
%   2. Elongate transcript according to nucleic acid limits (substrates)
%      if elongation factors are available
%   3. If transcription complete and termination factor available
%      - release transcript
%      - transition RNA polymerase to free state
%      - increment gene expression
%      Otherwise remain in active state
%
%   Specifically bound State:
%   - Can transition to active, specifically bound, non-specifically bound,
%     or free states
%   - Transition into state only if a free sigma factor is available
%     1. Decrement number of free sigma factors
%     2. Pick a transcription unit (tu) to bind to according to
%     Expression transcription unit i~prob(ribosome releases tu i|ribosome active)
%                      =prob(ribosome within RNA polymerase elongation rate bases of length of tu i|ribosome active)
%                      =prob(ribosome within RNA polymerase elongation rate bases of length of tu i|ribosome active, bound to tu i)*prob(ribosome bound to tu i|ribosome active)
%                      =[(RNA polymerase transcription rate)/(length of tu i)] * [(length of tu i)*prob(binding tu i|binding)]
%     prob(binding tu i | binding)~expression tu i
%
%   Non-specifically bound state:
%   - Can transition to specifically bound, non-specifically bound, or free
%     states
%
%   Free state:
%   - Can transition to specifically bound, non-specifically bound, or free
%     states
%
%   Algorithm
%   +++++++++++++++
%   1. Randomly transition RNA polymerases among activlely transcribing,
%      specifically bound, non-specifically, bound, and free states weighted by
%      state transition probabilities. Update rnaPolymerases.states.
%   2. Randomly assign RNA polymerases entering the specifically bound to
%      specific transcription units weighted by the product of
%      transcriptionUnitBindingProbabilities and
%      transcriptionFactorBindingProbFoldChange. Update
%      transcripts.boundTranscriptionUnits.
%   3. Assign RNA polymerase entering the actively transcribing state sigma
%      factors. Update enzymes and boundEnzymes.
%   4. Simulate RNA polymerization by actively transcribing RNA polymerases with
%      the aid of elongation factors. Allocate available nucleic acids among the
%      actively transcribing RNA polymerases. Release sigma factors from RNA
%      polymerases that started at the beginning of the transcription unit and
%      progressed. Update rnaPolymerases.states. Update substrates. Update enzymes
%      and boundEnzymes.
%   5. If termination factors are available dissociate RNA polymerases which
%      have reached the terminus of the transcription they're bound to, and
%      release RNAs. Update rnaPolymerases.states and
%      transcripts.boundTranscriptionUnits. Increment RNAs.
%
%   References
%   ===========
%   1. McClure, W. R. 1985. Mechanism and control of transcription
%      initiation in prokaryotes. Annu. Rev. Biochem. 54:171-204.
%      [PUB_0775]
%   2. Ciampi MS (2006). Rho-dependent terminators and transcription
%      termination. Microbiology. 152(9):2515-28 [PUB_0233].
%   3. Nudler E, Gottesman ME (2002). Transcription termination and
%      anti-termination in E. coli. Genes Cells. 7(8):755-68. [PUB_0662]
%   4. Washio T, Sasayama J, Tomita M (1998). Analysis of complete genomes
%      suggests that many prokaryotes do not rely on hairpin formation in
%      transcription termination. Nucleic Acids Res. 26(23):5456-63
%      [PUB_0234]
%   5. Peterson JD, Umayam LA, Dickinson T, Hickey EK, White O (2001). The
%      Comprehensive Microbial Resource. Nucleic Acids Res. 29(1):123-5.
%      [PUB_0182]
%   6. Shepherd N, Dennis P, Bremer H (2001). Cytoplasmic RNA Polymerase in
%      Escherichia coli. J Bacteriol. 183(8): 2527-34. [PUB_0784]
%   7. Klumpp S, Hwa T (2008). Growth-rate-dependent partitioning of RNA
%      polymerases in bacteria. Proc Natl Acad Sci U S A. 105(21):
%      20245-50. [PUB_0785]
%   8. Grigorova IL, Phleger NJ, Mutalik VK, Gross CA (2006). Insights into
%      transcriptional regulation and sigma competition from an equilibrium
%      model of RNA polymerase binding to DNA. Proc Natl Acad Sci U S A.
%      103(14): 5332-7. [PUB_0786]
%
% Author: Markus Covert, mcovert@stanford.edu
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/4/2010
%
%TODO: require Mg2+ or Mn2+ as cofactor for transcription
classdef Transcription < edu.stanford.covert.cell.sim.Process & edu.stanford.covert.cell.sim.ChromosomeProcessAspect
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
            'rnaPolymeraseElongationRate';
            'transcriptionUnitBaseCounts';
            'transcriptionUnitBindingProbabilities';
            };
        fittedConstantNames__      = {   %names of fitted constant properties
            'transcriptionUnitBindingProbabilities'
            };
        localStateNames__          = {   %names of simulation state properties redundant with timecourses in this or other processes or the simulation
            'RNAs'};
    end
    
    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs = {}; %whole cell model IDs of stimuli
        
        substrateWholeCellModelIDs = { %whole cell model IDs of substrates
            'ATP'; 'CTP'; 'GTP'; 'UTP';
            'AMP'; 'CMP'; 'GMP'; 'UMP';
            'ADP'; 'PPI'; 'H2O'; 'H'};
        substrateIndexs_ntp         = (1:4)'; %indices within substrates of NTPs
        substrateIndexs_nmp         = (5:8)'; %indices within substrates of NMPs
        substrateIndexs_atp         = 1;      %index within substrates of ATP
        substrateIndexs_adp         = 9;      %index within substrates of ADP
        substrateIndexs_diphosphate = 10;     %index within substrates of diphosphate
        substrateIndexs_water       = 11;     %index within substrates of water
        substrateIndexs_hydrogen    = 12;     %index within substrates of hydrogen
        
        enzymeWholeCellModelIDs = {       %enzyme whole cell model ids
            'MG_249_MONOMER';             %RNA polymerase sigma factor RpoD
            'MG_282_MONOMER';             %transcription elongation factor GreA
            'MG_141_MONOMER';             %transcription termination factor NusA
            'MG_027_MONOMER';             %transcription termination/antitermination protein NusB
            'RNA_POLYMERASE';             %DNA-directed RNA polymerase
            'RNA_POLYMERASE_HOLOENZYME'}; %DNA-directed RNA polymerase holoenzyme
        enzymeIndexs_transcriptionFactors    = (1:4)'; %indices within enzymes of transcription factors
        enzymeIndexs_sigmaFactor             = 1;      %index within enzymes of sigma factor RpoD
        enzymeIndexs_elongationFactor        = 2;      %index within enzymes of elongation factor GreA
        enzymeIndexs_terminationFactor       = 3;      %index within enzymes of termination factor NusA
        enzymeIndexs_antiterminationFactor   = 4;      %index within enzymes of termination/antitermination protein NusB
        enzymeIndexs_rnaPolymerase           = 5;      %index within enzymes of DNA-directed RNA polymerase
        enzymeIndexs_rnaPolymeraseHoloenzyme = 6;      %index within enzymes of DNA-directed RNA polymerase holoenzyme
        
        complexIndexs_DnaA_ATP                         %indices within protein complexes of DnaA-ATP
        transcriptionUnitIndexs_DnaAR12345Boxes        %indices of transcription units containing the functional DnaA boxes R1-5
    end
    
    %fixed biological constants
    properties
        rnaPolymeraseElongationRate                   %RNA polymerase elongation rate (50 nucleotides per second per RNA polymerase) [PUB_0562, PUB_0563]
        transcriptionUnitBaseCounts                   %transcription unit base counts
        transcriptionUnitBindingProbabilities         %transcription unit binding probabilities
        stateTransitionProbabilities                  %transition probabilities among RNA polymerase states
    end
    
    %global state, stored locally for convenience
    properties
        RNAs                                          %copy number of transcription units
    end
    
    %global state (referenced locally for convenience)
    properties
        rnaPolymerases   %RNA Polymerase state class
        transcripts      %New Transcripts state class
    end
    
    %constructor
    methods
        function this = Transcription(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.Process(wholeCellModelID, name);
        end
    end
    
    %communication between process/simulation
    methods
        %set references to state objects
        function storeObjectReferences(this, simulation)
            this.storeObjectReferences@edu.stanford.covert.cell.sim.Process(simulation);
            this.storeObjectReferences@edu.stanford.covert.cell.sim.ChromosomeProcessAspect(simulation);
            
            this.rnaPolymerases = simulation.state('RNAPolymerase');
            this.transcripts = simulation.state('Transcript');
            
            this.states = [this.states; {this.rnaPolymerases; this.transcripts}];
        end
        
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)
            %call super class method
            this.initializeConstants@edu.stanford.covert.cell.sim.Process(knowledgeBase, simulation, varargin{:});
            this.initializeConstants@edu.stanford.covert.cell.sim.ChromosomeProcessAspect(knowledgeBase, simulation, varargin{:});
            
            this.transcriptionUnitBaseCounts = reshape([knowledgeBase.transcriptionUnits.baseCount],[],length(knowledgeBase.transcriptionUnits))';
            this.transcriptionUnitBaseCounts = this.transcriptionUnitBaseCounts(:, this.substrateMetaboliteGlobalIndexs);
            
            this.transcriptionUnitBindingProbabilities  = zeros(size(this.transcripts.transcriptionUnitLengths));
            
            this.complexIndexs_DnaA_ATP = find(ismember({knowledgeBase.proteinComplexs.wholeCellModelID}', {
                'MG_469_1MER_ATP'      %DnaA-ATP 1mer
                'MG_469_2MER_1ATP_ADP' %DnaA 2mer-(1)ATP-(1)ADP
                'MG_469_2MER_ATP'      %DnaA-ATP 2mer
                'MG_469_3MER_2ATP_ADP' %DnaA 3mer-(2)ATP-(1)ADP
                'MG_469_3MER_ATP'      %DnaA-ATP 3mer
                'MG_469_4MER_3ATP_ADP' %DnaA 4mer-(3)ATP-(1)ADP
                'MG_469_4MER_ATP'      %DnaA-ATP 4mer
                'MG_469_5MER_4ATP_ADP' %DnaA 5mer-(4)ATP-(1)ADP
                'MG_469_5MER_ATP'      %DnaA-ATP 5mer
                'MG_469_6MER_5ATP_ADP' %DnaA 6mer-(5)ATP-(1)ADP
                'MG_469_6MER_ATP'      %DnaA-ATP 6mer
                'MG_469_7MER_6ATP_ADP' %DnaA 7mer-(6)ATP-(1)ADP
                'MG_469_7MER_ATP'}));  %DnaA-ATP 7mer;
            [~, idxs] = ismember({
                'Functional box R1'
                'Functional box R2'
                'Functional box R3'
                'Functional box R4'
                'Functional box R5'}, {knowledgeBase.genomeFeatures.name});
            tus = [knowledgeBase.genomeFeatures(idxs).transcriptionUnits];
            this.transcriptionUnitIndexs_DnaAR12345Boxes = unique([tus.idx])';
        end
        
        %retrieve state from simulation
        function copyFromState(this)
            this.copyFromState@edu.stanford.covert.cell.sim.Process();
            
            this.RNAs = this.rna.counts(this.rna.nascentIndexs, this.compartment.cytosolIndexs, :);
        end
        
        %send state to simulation
        function copyToState(this)
            this.copyToState@edu.stanford.covert.cell.sim.Process();
            
            this.rna.counts(this.rna.nascentIndexs, this.compartment.cytosolIndexs, :) = this.RNAs;
        end
    end
    
    %allocate memory for state
    methods
        function allocateMemoryForState(this, numTimePoints)
            this.allocateMemoryForState@edu.stanford.covert.cell.sim.Process(numTimePoints);
            
            numTranscriptionUnits = length(this.transcripts.transcriptionUnitLengths);
            
            this.RNAs = zeros(numTranscriptionUnits, 1, numTimePoints);
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
            %Contributions, by process, to RNA metabolism:
            %                 NTPs + H2O --(transcription)--> nascent RNA + PPi + H
            %          nascent RNA + H2O --(processing)-----> processed RNA + NMPs
            %processed RNA + metabolites --(modification)---> modified RNA + metabolites
            %               modified RNA --(decay)----------> NMPs + modified NMPs
            %               NMPs + 2 ATP --(charging)-------> NTPs + 2 ADP
            %        modified NMPs + H2O --(catabolism)-----> modified bases + Pi
            %
            %Here we only account for the contributions of transcription.
            invMat = edu.stanford.covert.util.ComputationUtil.invertCompositionMatrix(this.rna.nascentRNAMatureRNAComposition);
            
            transcribedNTPs = this.transcriptionUnitBaseCounts(:, this.substrateIndexs_nmp)' * invMat * states.rnaProductions;
            bmProd(this.substrateIndexs_ntp)         = bmProd(this.substrateIndexs_ntp)         + transcribedNTPs;
            byProd(this.substrateIndexs_diphosphate) = byProd(this.substrateIndexs_diphosphate) + sum(transcribedNTPs);
            bmProd(this.substrateIndexs_water)       = bmProd(this.substrateIndexs_water)       + sum(states.rnaProductions);
            byProd(this.substrateIndexs_hydrogen)    = byProd(this.substrateIndexs_hydrogen)    + sum(states.rnaProductions);
            
            %% enzymes
            %RNA polymerase
            fractionInitiating = this.rnaPolymeraseElongationRate * sum(invMat * states.rnaProductions0) / ... %fraction of active RNA polymerases that are initiating
                (this.transcripts.transcriptionUnitLengths' * invMat * states.rnaProductions0);
            minEnzExp(this.enzymeIndexs_rnaPolymerase) = 4.5867 * ...
                (this.transcripts.transcriptionUnitLengths' * invMat * states.rnaProductions0) / ...
                this.rnaPolymeraseElongationRate / ...
                this.rnaPolymerases.stateExpectations(this.rnaPolymerases.activelyTranscribingIndex) / ...
                (1 - fractionInitiating);
            
            %sigma, elongation, termination factors
            minEnzExp(this.enzymeIndexs_sigmaFactor) = minEnzExp(this.enzymeIndexs_rnaPolymerase) * ...
                (this.rnaPolymerases.stateExpectations(this.rnaPolymerases.activelyTranscribingIndex) * fractionInitiating + ...
                this.rnaPolymerases.stateExpectations(this.rnaPolymerases.specificallyBoundIndex));
            minEnzExp(this.enzymeIndexs_elongationFactor)  = 2;
            minEnzExp(this.enzymeIndexs_terminationFactor) = 2;
        end
        
        %Initialize RNA polymerase, bound transcription units, transcription factor states
        %Assumptions:
        %- Metabolic cost of the transcription of these nucleic acids is negligible
        %- Probability that so many RNA polymerases bind a given transcription unit
        %  that they can't pack along the transcription unit without steric
        %  interactions is negligible. Thus we don't try to handle this case
        %  separately
        function initializeState(this)
            rnaPols = this.rnaPolymerases;
            trnspts = this.transcripts;
            
            %% states
            trnsptStrnds = 2 - trnspts.transcriptionUnitDirections;
            trnsptLens = trnspts.transcriptionUnitLengths;
            trnsptDirs = 2 * trnspts.transcriptionUnitDirections - 1;
            trnspt5Coords = trnspts.transcriptionUnitFivePrimeCoordinates;
            trnsptStarts = ...
                trnspt5Coords - (1 - trnspts.transcriptionUnitDirections) .* (trnspts.transcriptionUnitLengths-1);
            trnsptStarts(trnsptDirs == 1) = ...
                trnsptStarts(trnsptDirs == 1) - ...
                this.enzymeDNAFootprints5Prime(this.enzymeIndexs_rnaPolymerase);
            trnsptStarts(trnsptDirs == -1) = ...
                trnsptStarts(trnsptDirs == -1) + ...
                this.enzymeDNAFootprints5Prime(this.enzymeIndexs_rnaPolymerase);
            
            %% enzymes
            freTFs = this.enzymes(this.enzymeIndexs_transcriptionFactors);
            bndTFs = this.boundEnzymes(this.enzymeIndexs_transcriptionFactors);
            frePls = this.enzymes(this.enzymeIndexs_rnaPolymerase);
            bndPls = this.boundEnzymes(this.enzymeIndexs_rnaPolymerase);
            freHls = this.enzymes(this.enzymeIndexs_rnaPolymeraseHoloenzyme);
            bndHls = this.boundEnzymes(this.enzymeIndexs_rnaPolymeraseHoloenzyme);
            pls = frePls + bndPls + freHls + bndHls;
            
            %% allocate
            rnaPols.states = repmat(rnaPols.notExistValue, [2 * pls, 1, 1]);
            rnaPols.positionStrands = zeros([2 * pls, 2, 1]);
            trnspts.boundTranscriptionUnits = repmat(trnspts.nullTranscriptValue, [2 * pls, 1, 1]);
            trnspts.boundTranscriptProgress = repmat(trnspts.nullTranscriptValue, [2 * pls, 1, 1]);
            trnspts.boundTranscriptChromosome = repmat(trnspts.nullTranscriptValue, [2 * pls, 1, 1]);
            
            %% initialize
            pState = rnaPols.stateExpectations;
            probs = this.computeRNAPolymeraseTUBindingProbabilities();
            pBinds = trnspts.transcriptionUnitLengths .* probs(:, 1);
            pSpBinds = pBinds;
            
            for i = 1:pls
                %check whether conditions necessary for actively transcribing state are met
                if frePls == 0 || ~any(pBinds) || freTFs(this.enzymeIndexs_elongationFactor) == 0
                    pState(rnaPols.activelyTranscribingIndex) = 0;
                end
                
                %check whether conditions necessary for specifically bound state are met
                if frePls == 0 || ~any(pSpBinds) || freTFs(this.enzymeIndexs_sigmaFactor) == 0
                    pState(rnaPols.specificallyBoundIndex) = 0;
                end
                
                %check whether conditions necessary for non-specifically bound state are met
                if frePls == 0
                    pState(rnaPols.nonSpecificallyBoundIndex) = 0;
                end
                
                if ~any(pState)
                    pState(rnaPols.freeIndex) = 1;
                end
                
                %randomly select state
                state = this.randStream.randsample(4, 1, true, pState);
                
                %Actively transcribing, specifically bound: randomly select
                %transcription unit
                switch state
                    case rnaPols.activelyTranscribingIndex
                        tmp_pBinds = pBinds;
                        posStrnds = zeros(0, 2);
                        while any(tmp_pBinds)
                            iTU = this.randStream.randsample(numel(pBinds), 1, true, tmp_pBinds);
                            tmp_pBinds(iTU) = 0;
                            [~, posStrnds] = this.bindProteinToChromosomeStochastically(...
                                this.enzymeIndexs_rnaPolymerase, 1, [trnsptStarts(iTU)+1 trnsptStrnds(iTU)], trnsptLens(iTU)-2);
                            if ~isempty(posStrnds)
                                break;
                            end
                        end
                        if isempty(posStrnds)
                            throw(MException('Transcription:error', 'unable to bind chromosome'));
                        end
                        
                        posStrnds( isodd(posStrnds(:, 2)), 1) = posStrnds( isodd(posStrnds(:, 2)), 1) + this.enzymeDNAFootprints5Prime(this.enzymeIndexs_rnaPolymerase);
                        posStrnds(iseven(posStrnds(:, 2)), 1) = posStrnds(iseven(posStrnds(:, 2)), 1) + this.enzymeDNAFootprints3Prime(this.enzymeIndexs_rnaPolymerase);
                        if trnspts.transcriptionUnitDirections(iTU)
                            state = posStrnds(:, 1) - trnspt5Coords(iTU) + 1;
                        else
                            state = trnspt5Coords(iTU) - posStrnds(:, 1) + 1;
                        end
                        
                        rnaPols.states(i) = state;
                        rnaPols.positionStrands(i, :) = posStrnds;
                        trnspts.boundTranscriptProgress(i) = state;
                        trnspts.boundTranscriptionUnits(i) = iTU;
                        trnspts.boundTranscriptChromosome(i) = 1; %initialize to first chromosome
                        
                        frePls = frePls - 1;
                        bndPls = bndPls + 1;
                        
                        %freTFs(this.enzymeIndexs_elongationFactor) = freTFs(this.enzymeIndexs_elongationFactor) - 1;
                        %bndTFs(this.enzymeIndexs_elongationFactor) = bndTFs(this.enzymeIndexs_elongationFactor) + 1;
                    case rnaPols.specificallyBoundIndex
                        while any(pSpBinds)
                            iTU = this.randStream.randsample(numel(pSpBinds), 1, true, pSpBinds);
                            pSpBinds(iTU) = 0;
                            [~, ~, posStrnds] = this.bindProteinToChromosome([trnspt5Coords(iTU)-trnsptDirs(iTU) trnsptStrnds(iTU)], this.enzymeIndexs_rnaPolymeraseHoloenzyme, 1);
                            if ~isempty(posStrnds)
                                break;
                            end
                        end
                        if isempty(posStrnds)
                            throw(MException('Transcription:error', 'unable to bind chromosome'));
                        end
                        
                        rnaPols.states(i) = rnaPols.specificallyBoundValue;
                        rnaPols.positionStrands(i, :) = posStrnds;
                        trnspts.boundTranscriptProgress(i) = trnspts.nullTranscriptValue;
                        trnspts.boundTranscriptionUnits(i) = iTU;
                        trnspts.boundTranscriptChromosome(i) = 1; %initialize to first chromosome
                        
                        bndHls = bndHls + 1;
                        frePls = frePls - 1;
                        freTFs(this.enzymeIndexs_sigmaFactor) = freTFs(this.enzymeIndexs_sigmaFactor) - 1;
                    case rnaPols.nonSpecificallyBoundIndex
                        [~, posStrnds] = this.bindProteinToChromosomeStochastically(this.enzymeIndexs_rnaPolymerase, 1);
                        if isempty(posStrnds)
                            throw(MException('Transcription:error', 'unable to bind chromosome'));
                        end
                        posStrnds( isodd(posStrnds(:, 2)), 1) = posStrnds( isodd(posStrnds(:, 2)), 1) + this.enzymeDNAFootprints5Prime(this.enzymeIndexs_rnaPolymerase);
                        posStrnds(iseven(posStrnds(:, 2)), 1) = posStrnds(iseven(posStrnds(:, 2)), 1) + this.enzymeDNAFootprints3Prime(this.enzymeIndexs_rnaPolymerase);
                        
                        rnaPols.states(i) = rnaPols.nonSpecificallyBoundValue;
                        rnaPols.positionStrands(i, :) = posStrnds;
                        trnspts.boundTranscriptionUnits(i) = trnspts.nullTranscriptValue;
                        trnspts.boundTranscriptProgress(i) = trnspts.nullTranscriptValue;
                        trnspts.boundTranscriptChromosome(i) = trnspts.nullTranscriptValue;
                        
                        frePls = frePls - 1;
                        bndPls = bndPls + 1;
                    case rnaPols.freeIndex
                        rnaPols.states(i) = rnaPols.freeValue;
                        rnaPols.positionStrands(i, :) = 0;
                        trnspts.boundTranscriptionUnits(i) = trnspts.nullTranscriptValue;
                        trnspts.boundTranscriptProgress(i) = trnspts.nullTranscriptValue;
                        trnspts.boundTranscriptChromosome(i) = trnspts.nullTranscriptValue;
                end
            end
            
            %% store states of enzymes
            this.enzymes(this.enzymeIndexs_transcriptionFactors)      = freTFs;
            this.boundEnzymes(this.enzymeIndexs_transcriptionFactors) = bndTFs;
            this.enzymes(this.enzymeIndexs_rnaPolymerase)             = frePls;
            this.boundEnzymes(this.enzymeIndexs_rnaPolymerase)        = bndPls;
            this.enzymes(this.enzymeIndexs_rnaPolymeraseHoloenzyme)   = freHls;
            this.boundEnzymes(this.enzymeIndexs_rnaPolymeraseHoloenzyme) = bndHls;
            
            %% decrement counts of RNAs
            this.copyToState();
            
            comp = this.compartment;
            rna = this.rna;
            
            matureRnaWt = max(0, sum(rna.dryWeight) - trnspts.dryWeight);
            initRnaCnts = rna.counts;
            while sum(rna.dryWeight) > matureRnaWt
                idx = this.randStream.randsample(size(rna.counts, 1), 1, false, rna.counts(:, comp.cytosolIndexs));
                rna.counts(idx, comp.cytosolIndexs) = rna.counts(idx, comp.cytosolIndexs) - 1;
            end
            rna.counts(:, comp.cytosolIndexs) = ...
                + rna.counts(:, comp.cytosolIndexs) ...
                + rna.updateExternalState(-(initRnaCnts(:, comp.cytosolIndexs) - rna.counts(:, comp.cytosolIndexs)), true);
            
            this.copyFromState();
        end
        
        %resource requirements
        function result = calcResourceRequirements_Current(this)
            result = zeros(size(this.substrates));
            result(this.substrateIndexs_ntp) = 2 * ...
                sum(this.metabolite.nmpComposition, 2) * this.rnaPolymerases.nActive * ...
                this.rnaPolymeraseElongationRate;
            result(this.substrateIndexs_water) = this.rnaPolymerases.nActive;
        end
        
        %simulation
        function evolveState(this)
            %% define and allocate variables
            
            %states
            c = this.chromosome;
            rnaPols = this.rnaPolymerases;
            trnspts = this.transcripts;
            
            %numbers of enzymes
            nFreTfs  = this.enzymes(this.enzymeIndexs_transcriptionFactors);
            nBndTFs  = this.boundEnzymes(this.enzymeIndexs_transcriptionFactors);
            nFrePols = this.enzymes(this.enzymeIndexs_rnaPolymerase);
            nBndPols = this.boundEnzymes(this.enzymeIndexs_rnaPolymerase);
            nFreHols = this.enzymes(this.enzymeIndexs_rnaPolymeraseHoloenzyme);
            nBndHols = this.boundEnzymes(this.enzymeIndexs_rnaPolymeraseHoloenzyme);
            nTotPols = nFrePols + nBndPols + nFreHols + nBndHols;
            
            %properties of RNA polymerases
            pStTrns = this.stateTransitionProbabilities; %probabilities of state transitions
            
            %transcription unit properties
            tuDirs = 2 * trnspts.transcriptionUnitDirections - 1;
            tuStrnds = c.transcriptionUnitStrands;
            tuLens = trnspts.transcriptionUnitLengths;
            tu5Coords = trnspts.transcriptionUnitFivePrimeCoordinates;
            tuSeqs = trnspts.transcriptionUnitSequences;  %sequences
            
            %% Update states of RNA polymerases
            
            %allocate space to store states of new RNA polymerases
            if nTotPols > size(rnaPols.states, 1)
                rnaPols.states = [
                    rnaPols.states;
                    rnaPols.notExistValue(ones(2 * nTotPols - size(rnaPols.states, 1), 1), 1)];
                rnaPols.positionStrands = [
                    rnaPols.positionStrands;
                    zeros(2 * nTotPols - size(rnaPols.positionStrands, 1), 2)];
                trnspts.boundTranscriptionUnits = [
                    trnspts.boundTranscriptionUnits;
                    trnspts.nullTranscriptValue(ones(2 * nTotPols - size(trnspts.boundTranscriptChromosome, 1), 1), 1)];
                trnspts.boundTranscriptProgress = [
                    trnspts.boundTranscriptProgress;
                    trnspts.nullTranscriptValue(ones(2 * nTotPols - size(trnspts.boundTranscriptChromosome, 1), 1), 1)];
                trnspts.boundTranscriptChromosome = [
                    trnspts.boundTranscriptChromosome;
                    trnspts.nullTranscriptValue(ones(2 * nTotPols - size(trnspts.boundTranscriptChromosome, 1), 1), 1)];
            end
            
            %indices of RNA polymerases in various states
            actPols = find(rnaPols.states >= rnaPols.activelyTranscribingValue);
            sbPols = find(rnaPols.states == rnaPols.specificallyBoundValue);
            freeNsbPols = find(...
                rnaPols.states == rnaPols.nonSpecificallyBoundValue | ...
                rnaPols.states == rnaPols.freeValue);
            sbPols = sbPols(this.randStream.randperm(numel(sbPols)));
            freeNsbPols = freeNsbPols(this.randStream.randperm(numel(freeNsbPols)));
            
            %number of new RNA polymerases
            nNewPols = nTotPols - (numel(actPols) + numel(sbPols) + numel(freeNsbPols));
            
            %dissociate any free RNA polymerase holoenzymes (created say by
            %their release from chromosome by another protein) into RNA
            %polymerase core and sigma factor
            nFrePols = nFrePols + nFreHols;
            nFreTfs(this.enzymeIndexs_sigmaFactor) = nFreTfs(this.enzymeIndexs_sigmaFactor) + nFreHols;
            nFreHols = 0;
            
            %initiating specifically bound polymerases
            nInitMax = this.randStream.stochasticRound(numel(sbPols) * pStTrns(rnaPols.activelyTranscribingIndex, rnaPols.specificallyBoundIndex)); %#ok<*PROP>
            if nInitMax > 0
                releasedProteins = this.releaseProteinFromSites(rnaPols.positionStrands(sbPols(1:nInitMax), :), false, this.enzymeIndexs_rnaPolymeraseHoloenzyme, true, true);
                if ~isequal(releasedProteins(this.enzymeIndexs_rnaPolymeraseHoloenzyme), nInitMax)
                    throw(MException('Transcription:error', 'Unable to release protein'));
                end
                
                %try to initiate RNA polymerases
                positions = tu5Coords(trnspts.boundTranscriptionUnits(sbPols(1:nInitMax)));
                tfs = this.bindProteinToChromosome([positions rnaPols.positionStrands(sbPols(1:nInitMax), 2)], ...
                    this.enzymeIndexs_rnaPolymeraseHoloenzyme);
                rnaPols.states(sbPols(tfs)) = rnaPols.activelyTranscribingValue;
                trnspts.boundTranscriptProgress(sbPols(tfs)) = trnspts.newTranscriptValue;
                rnaPols.positionStrands(sbPols(tfs), 1) = positions(tfs, 1);
                
                %rebind RNA polymerases that couldn't move forward
                if ~all(this.bindProteinToChromosome(rnaPols.positionStrands(sbPols(~tfs), :), ...
                        this.enzymeIndexs_rnaPolymeraseHoloenzyme))
                    throw(MException('Transcription:error', 'Unable to unbind protein'));
                end
            end
            
            %specifically bound polymerases becoming non-specifically bound / free
            nUnsb = max(0, numel(sbPols) - nInitMax - this.randStream.stochasticRound(numel(sbPols) * ...
                pStTrns(rnaPols.specificallyBoundIndex, rnaPols.specificallyBoundIndex)));
            if nUnsb > 0
                idxs = sbPols(end-nUnsb+1:end);
                posStrnds = rnaPols.positionStrands(idxs, :);
                releasedProteins = this.releaseProteinFromSites(posStrnds, false, this.enzymeIndexs_rnaPolymeraseHoloenzyme, true, true);
                if ~isequal(releasedProteins(this.enzymeIndexs_rnaPolymeraseHoloenzyme), nUnsb)
                    throw(MException('Transcription:error', 'Unable to unbind protein'));
                end
                
                rnaPols.states(idxs) = rnaPols.freeValue;
                rnaPols.positionStrands(idxs, :) = 0;
                trnspts.boundTranscriptionUnits(idxs) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptProgress(idxs) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptChromosome(idxs) = trnspts.nullTranscriptValue;
                
                nBndHols = nBndHols - nUnsb;                
                nFrePols = nFrePols + nUnsb;
                nFreTfs(this.enzymeIndexs_sigmaFactor) = nFreTfs(this.enzymeIndexs_sigmaFactor) + nUnsb;
            end
            
            %free polymerases/non-specifically bound becoming specifically bound
            nFreeNsbPols = this.randStream.stochasticRound(min([
                numel(freeNsbPols) * pStTrns(rnaPols.specificallyBoundIndex, rnaPols.freeIndex)
                nFreTfs(this.enzymeIndexs_sigmaFactor)]));
            if nFreeNsbPols > 0
                %unbind from non-specifically bound site
                idxs = freeNsbPols(1:nFreeNsbPols);
                nNsbPols = sum(rnaPols.states(idxs) == rnaPols.nonSpecificallyBoundValue);
                releasedProteins = this.releaseProteinFromSites(rnaPols.positionStrands(idxs, :), false, this.enzymeIndexs_rnaPolymerase, true, true);
                if ~isequal(releasedProteins(this.enzymeIndexs_rnaPolymerase), nNsbPols)
                    throw(MException('Transcription:error', 'Unable to release protein'));
                end
                
                rnaPols.states(idxs) = rnaPols.freeValue;
                rnaPols.positionStrands(idxs, :) = 0;
                trnspts.boundTranscriptionUnits(idxs) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptProgress(idxs) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptChromosome(idxs) = trnspts.nullTranscriptValue;
                
                nFrePols = nFrePols + nNsbPols;
                nBndPols = nBndPols - nNsbPols;
                
                %bind to new specifically bound site
                pBnds = this.computeRNAPolymeraseTUBindingProbabilities();
                
                if nnz(c.polymerizedRegions) == 2
                    [~, iTU, posStrnds, nFreeNsbPols] = this.bindProteinToChromosome( ...
                        [tu5Coords-tuDirs tuStrnds], ...
                        this.enzymeIndexs_rnaPolymeraseHoloenzyme, nFreeNsbPols, pBnds(:, 1), ...
                        true, true, 1, false, []);
                    iChr = 1;
                else
                    [~, iTU, posStrnds, nFreeNsbPols] = this.bindProteinToChromosome( ...
                        [tu5Coords-tuDirs tuStrnds
                        tu5Coords-tuDirs tuStrnds + 2], ...
                        this.enzymeIndexs_rnaPolymeraseHoloenzyme, nFreeNsbPols, pBnds(:), ...
                        true, true, 1, false, []);
                    iTU = mod(iTU - 1, numel(tuLens)) + 1;
                    iChr = ceil(posStrnds(:, 2) / 2);
                end
                
                idxs = freeNsbPols(1:nFreeNsbPols);
                rnaPols.states(idxs) = rnaPols.specificallyBoundValue;
                rnaPols.positionStrands(idxs, :) = posStrnds;
                trnspts.boundTranscriptionUnits(idxs) = iTU;
                trnspts.boundTranscriptProgress(idxs) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptChromosome(idxs) = iChr;
                
                nBndHols = nBndHols + nFreeNsbPols;
                nFrePols = nFrePols - nFreeNsbPols;
                nFreTfs(this.enzymeIndexs_sigmaFactor) = nFreTfs(this.enzymeIndexs_sigmaFactor) - nFreeNsbPols;
            end
            
            %free polymerases/non-specifically bound becoming free/non-specifically bound
            nsbPols = find(rnaPols.states == rnaPols.nonSpecificallyBoundValue);
            if ~isempty(nsbPols)
                %set non-specifically bound to free
                releasedProteins = this.releaseProteinFromSites(rnaPols.positionStrands(nsbPols, :), false, this.enzymeIndexs_rnaPolymerase, true, true);
                if ~isequal(releasedProteins(this.enzymeIndexs_rnaPolymerase), numel(nsbPols))
                    throw(MException('Transcription:error', 'Unable to release protein'));
                end
                
                rnaPols.states(nsbPols) = rnaPols.freeValue;
                rnaPols.positionStrands(nsbPols, :) = 0;
                trnspts.boundTranscriptionUnits(nsbPols) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptProgress(nsbPols) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptChromosome(nsbPols) = trnspts.nullTranscriptValue;
                
                nFrePols = nFrePols + numel(nsbPols);
                nBndPols = nBndPols - numel(nsbPols);
            end
            
            freePols = find(rnaPols.states == rnaPols.freeValue);
            nNsbPols = this.randStream.stochasticRound(min([...
                nFrePols
                numel(freePols) * pStTrns(rnaPols.nonSpecificallyBoundIndex, rnaPols.freeIndex) / (pStTrns(rnaPols.nonSpecificallyBoundIndex, rnaPols.freeIndex) + pStTrns(rnaPols.freeIndex, rnaPols.freeIndex))]));
            if nNsbPols > 0
                %set some to non-specifically bound state
                [nNsbPols posStrnds] = this.bindProteinToChromosomeStochastically(this.enzymeIndexs_rnaPolymerase, nNsbPols);
                posStrnds( isodd(posStrnds(:, 2)), 1) = posStrnds( isodd(posStrnds(:, 2)), 1) + this.enzymeDNAFootprints5Prime(this.enzymeIndexs_rnaPolymerase);
                posStrnds(iseven(posStrnds(:, 2)), 1) = posStrnds(iseven(posStrnds(:, 2)), 1) + this.enzymeDNAFootprints3Prime(this.enzymeIndexs_rnaPolymerase);
                
                rnaPols.states(freePols(1:nNsbPols)) = rnaPols.nonSpecificallyBoundValue;
                rnaPols.positionStrands(freePols(1:nNsbPols), :) = posStrnds;
                trnspts.boundTranscriptionUnits(freePols(1:nNsbPols)) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptProgress(freePols(1:nNsbPols)) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptChromosome(freePols(1:nNsbPols)) = trnspts.nullTranscriptValue;
                
                nFrePols = nFrePols - nNsbPols;
                nBndPols = nBndPols + nNsbPols;
            end
            
            %set newly created RNA polymerases to free state
            if nNewPols > 0
                newPols = find(rnaPols.states == rnaPols.notExistValue, nNewPols, 'first');
                rnaPols.states(newPols) = rnaPols.freeValue;
                rnaPols.positionStrands(newPols, :) = 0;
                trnspts.boundTranscriptionUnits(newPols) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptProgress(newPols) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptChromosome(newPols) = trnspts.nullTranscriptValue;
            end
            
            %% Transcribe active sequences
            %compute progress of transcription and cost
            usedNTPs = zeros(4, 1);
            if ~isempty(actPols) && nFreTfs(this.enzymeIndexs_elongationFactor) && this.rnaPolymeraseElongationRate > 0
                posStrnds = rnaPols.positionStrands(actPols, :);
                iTUs = trnspts.boundTranscriptionUnits(actPols);
                
                %temporarily release RNA polymerase from old position strands
                releasedProteins = ...
                    this.releaseProteinFromSites(posStrnds, false, this.enzymeIndexs_rnaPolymerase, true, true) + ...
                    this.releaseProteinFromSites(posStrnds, false, this.enzymeIndexs_rnaPolymeraseHoloenzyme, true, true);
                if ~isequal(releasedProteins(this.enzymeIndexs_rnaPolymerase) + releasedProteins(this.enzymeIndexs_rnaPolymeraseHoloenzyme), numel(actPols))
                    throw(MException('Transcription:error', 'Unable to unbind protein'));
                end
                
                %extract active sequences (part of sequence that should be transcribed)
                actSeqs = cell(size(actPols)); %active sequences (those can be transcribed)
                for idx = 1:numel(actPols)
                    i = actPols(idx);
                    left = rnaPols.states(i);
                    actSeqs{idx} = tuSeqs{trnspts.boundTranscriptionUnits(i)}(left:end);
                end
                actSeqs = char(actSeqs);
                actSeqs = actSeqs(:, 1:min(end, this.rnaPolymeraseElongationRate));
                
                %elongation limits by RNA polymerase kinetic rate and transcription unit extent
                elngMax = min(this.rnaPolymeraseElongationRate, ...
                    tuLens(trnspts.boundTranscriptionUnits(actPols)) - rnaPols.states(actPols) + 1);
                                
                %elongation limits by RNA-polymerase / RNA-polymerase interactions
                for strnd = 1 : 4
                    idxs = find(posStrnds(:, 2) == strnd);
                    if numel(idxs) <= 1
                        continue;
                    end
                    [~, order] = sort(posStrnds(idxs, 1));
                    
                    if isodd(strnd)
                        elngMax(idxs(order)) = min(...
                            elngMax(idxs(order)), ...
                            diff([posStrnds(idxs(order), 1); posStrnds(idxs(order(1)), 1) + c.sequenceLen]) - ...
                            this.enzymeDNAFootprints(this.enzymeIndexs_rnaPolymerase));
                    else
                        elngMax(idxs(order)) = min(...
                            elngMax(idxs(order)), ...
                            diff([posStrnds(idxs(order(end)), 1) - c.sequenceLen; posStrnds(idxs(order), 1)]) - ...
                            this.enzymeDNAFootprints(this.enzymeIndexs_rnaPolymerase));
                    end
                end
                
                %elongation limits by RNA-polymerase DNA interactions (Note: independent of other RNA polymerases)
                [~, ~, ~, tmp_elngMax] = c.isRegionAccessible(...
                    posStrnds, (elngMax + 1) .* tuDirs(iTUs), ...
                    [], this.enzymeGlobalIndexs(this.enzymeIndexs_rnaPolymerase), true, [], true, false);
                if ~all(tmp_elngMax)
                    throw(MException('Transcription:error', 'Unable to bind proteins'));
                end
                elngMax = min(abs(tmp_elngMax) - 1, elngMax);
                
                %elongation limits by NTPs
                for i = 1:size(actSeqs, 1)
                    actSeqs(i, elngMax(i) + 1:end) = ' ';
                end
                [elngProg, this.substrates(this.substrateIndexs_ntp), usedNTPs] = ...
                    edu.stanford.covert.cell.sim.util.polymerize(...
                    actSeqs, this.substrates(this.substrateIndexs_ntp), 'ACGU', ' ', 0, 0, this.randStream);
                
                %release sigma factor after initiation
                nInits = sum(...
                    rnaPols.states(actPols) == rnaPols.activelyTranscribingValue & ...
                    elngProg > 0);
                nBndHols = nBndHols - nInits;
                nFreTfs(this.enzymeIndexs_sigmaFactor) = nFreTfs(this.enzymeIndexs_sigmaFactor) + nInits;
                nBndPols = nBndPols + nInits;
                
                %rebind RNA polymerase to new positions
                rnaPols.states(actPols) = rnaPols.states(actPols) + elngProg;
                rnaPols.positionStrands(actPols, 1) = posStrnds(:, 1) + elngProg .* tuDirs(iTUs);
                trnspts.boundTranscriptProgress(actPols) = rnaPols.states(actPols);
                
                initPols = rnaPols.states(actPols) == rnaPols.activelyTranscribingValue;
                if ...
                        ~all(this.bindProteinToChromosome(rnaPols.positionStrands(actPols(~initPols), :), this.enzymeIndexs_rnaPolymerase)) || ...
                        ~all(this.bindProteinToChromosome(rnaPols.positionStrands(actPols( initPols), :), this.enzymeIndexs_rnaPolymeraseHoloenzyme))
                    throw(MException('Transcription:error', 'Unable to bind proteins'));
                end
            end
            
            %release polymerases which have completed transcription, and increment
            %expression of the corresponding transcription unit
            if ~isempty(actPols) && nFreTfs(this.enzymeIndexs_elongationFactor) && nFreTfs(this.enzymeIndexs_terminationFactor)
                trmPols = actPols(rnaPols.states(actPols) > ...
                    tuLens(trnspts.boundTranscriptionUnits(actPols)));
                nTrmPols = numel(trmPols);
                if nTrmPols > this.substrates(this.substrateIndexs_water)
                   order = this.randStream.randperm(nTrmPols);
                   nTrmPols = this.substrates(this.substrateIndexs_water);
                   trmPols = trmPols(order(1:nTrmPols));
                end
                
                if nTrmPols > 0
                    if isscalar(trmPols)
                        this.RNAs(trnspts.boundTranscriptionUnits(trmPols)) = ...
                            this.RNAs(trnspts.boundTranscriptionUnits(trmPols)) + 1;
                    else
                        this.RNAs = this.RNAs + ...
                            histc(trnspts.boundTranscriptionUnits(trmPols, 1), 1:numel(tuLens));
                    end
                    
                    %release proteins from chromsome
                    releasedProteins = this.releaseProteinFromSites(rnaPols.positionStrands(trmPols, :), false, this.enzymeIndexs_rnaPolymerase, true, true);
                    if ~isequal(releasedProteins(this.enzymeIndexs_rnaPolymerase), numel(trmPols))
                        throw(MException('Transcription:error', 'Unable to unbind protein'));
                    end
                    
                    nFrePols = nFrePols + nTrmPols;
                    nBndPols = nBndPols - nTrmPols;
                    
                    this.substrates(this.substrateIndexs_water)    = this.substrates(this.substrateIndexs_water)    - nTrmPols;
                    this.substrates(this.substrateIndexs_hydrogen) = this.substrates(this.substrateIndexs_hydrogen) + nTrmPols;
                    
                    rnaPols.states(trmPols) = rnaPols.freeValue;
                    rnaPols.positionStrands(trmPols, :) = 0;
                    trnspts.boundTranscriptionUnits(trmPols) = trnspts.nullTranscriptValue;
                    trnspts.boundTranscriptProgress(trmPols) = trnspts.nullTranscriptValue;
                    trnspts.boundTranscriptChromosome(trmPols) = trnspts.nullTranscriptValue;
                end
            end
            
            %% diphosphate released by polymerization
            this.substrates(this.substrateIndexs_diphosphate) = ...
                this.substrates(this.substrateIndexs_diphosphate) + ...
                sum(usedNTPs);
            
            %% store enzyme state
            this.enzymes(this.enzymeIndexs_transcriptionFactors) = nFreTfs;
            this.boundEnzymes(this.enzymeIndexs_transcriptionFactors) = nBndTFs;
            this.enzymes(this.enzymeIndexs_rnaPolymerase) = nFrePols;
            this.boundEnzymes(this.enzymeIndexs_rnaPolymerase) = nBndPols;
            this.enzymes(this.enzymeIndexs_rnaPolymeraseHoloenzyme) = nFreHols;
            this.boundEnzymes(this.enzymeIndexs_rnaPolymeraseHoloenzyme) = nBndHols;
        end
    end
    
    %model helper functions
    methods    
        function pBinds = computeRNAPolymeraseTUBindingProbabilities(this, protectDnaABoxes)
            c = this.chromosome;
            r = this.rnaPolymerases;
            
            %relative probability RNA polymerase binds each transcription
            %unit
            pBinds = this.transcriptionUnitBindingProbabilities(:, ones(1, 2)) .* ...
                r.transcriptionFactorBindingProbFoldChange .* ...
                r.supercoilingBindingProbFoldChange;
            
            %if functional DnaA boxes R1-5 are occupied by DnaA-ATP,
            %don't allow RNA polymerase to bind since it would knockoff
            %DnaA-ATP
            if nargin == 1 || protectDnaABoxes
                [posStrnds, complexs] = find(...
                    c.complexBoundSites(c.transcriptionUnitStartCoordinates(this.transcriptionUnitIndexs_DnaAR12345Boxes) + ...
                    (1:c.transcriptionUnitLengths(this.transcriptionUnitIndexs_DnaAR12345Boxes))' - 1, :));
                if any(ismembc(complexs, this.complexIndexs_DnaA_ATP) & posStrnds(:, 2) <= 2)
                    pBinds(this.transcriptionUnitIndexs_DnaAR12345Boxes, 1) = 0;
                end
                if any(ismembc(complexs, this.complexIndexs_DnaA_ATP) & posStrnds(:, 2) > 2)
                    pBinds(this.transcriptionUnitIndexs_DnaAR12345Boxes, 2) = 0;
                end
            end
        end
        
        function pStProbs = calcStateTransitionProbabilities(this, nPols, ntpProdRate, tuBindingProbs, method)
            import edu.stanford.covert.util.ComputationUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            rna = this.rna;
            rnaPols = this.rnaPolymerases;
            
            pStProbs = zeros(4, 4);
            pStProbs(rnaPols.nonSpecificallyBoundIndex, rnaPols.activelyTranscribingIndex) = ...
                ntpProdRate / (nPols * rnaPols.stateExpectations(rnaPols.activelyTranscribingIndex)) / ...
                (tuBindingProbs' * rna.lengths(rna.nascentIndexs));
            if nargin < 5 || isequal(method, 'handTuned')
                pStProbs(rnaPols.activelyTranscribingIndex, rnaPols.specificallyBoundIndex) = 1;
                pStProbs(rnaPols.nonSpecificallyBoundIndex, rnaPols.specificallyBoundIndex) = 0;
                pStProbs(rnaPols.specificallyBoundIndex, rnaPols.nonSpecificallyBoundIndex) = 0.10;
            else
                pStProbs(rnaPols.activelyTranscribingIndex, rnaPols.specificallyBoundIndex) = ...
                    pStProbs(rnaPols.nonSpecificallyBoundIndex, rnaPols.activelyTranscribingIndex) * ...
                    rnaPols.stateExpectations(rnaPols.activelyTranscribingIndex) / ...
                    rnaPols.stateExpectations(rnaPols.specificallyBoundIndex);
                pStProbs(rnaPols.specificallyBoundIndex, rnaPols.nonSpecificallyBoundIndex) = ...
                    pStProbs(rnaPols.nonSpecificallyBoundIndex, rnaPols.activelyTranscribingIndex) * ...
                    rnaPols.stateExpectations(rnaPols.activelyTranscribingIndex) / ...
                    rnaPols.stateExpectations(rnaPols.nonSpecificallyBoundIndex);
            end
            
            pStProbs(rnaPols.activelyTranscribingIndex, rnaPols.activelyTranscribingIndex) = ...
                1 - pStProbs(rnaPols.nonSpecificallyBoundIndex, rnaPols.activelyTranscribingIndex);
            pStProbs(rnaPols.specificallyBoundIndex, rnaPols.specificallyBoundIndex) = ...
                + 1 ...
                - pStProbs(rnaPols.activelyTranscribingIndex, rnaPols.specificallyBoundIndex)...
                - pStProbs(rnaPols.nonSpecificallyBoundIndex, rnaPols.specificallyBoundIndex);
            pStProbs(rnaPols.nonSpecificallyBoundIndex, rnaPols.nonSpecificallyBoundIndex) = ...
                1 - pStProbs(rnaPols.specificallyBoundIndex, rnaPols.nonSpecificallyBoundIndex);
            pStProbs(:, rnaPols.freeIndex) = pStProbs(:, rnaPols.nonSpecificallyBoundIndex);
        end
    end
    
    %get methods of dependent local state
    methods
        function value = getDryWeight(this)
            if size(this.RNAs, 3) == 1
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + ...
                    this.rna.molecularWeights(this.rna.nascentIndexs)' * this.RNAs ...
                    / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            else
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + ...
                    permute(this.rna.molecularWeights(this.rna.nascentIndexs)' * permute(this.RNAs,[1 3 2]),[1 3 2]) ...
                    / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            end
        end
    end
end
