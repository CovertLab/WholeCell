%ReplicationInitiation
%
% @wholeCellModelID Process_ReplicationInitiation
% @name             ReplicationInitiation
% @description
%   Biology
%   =======================================
%   Chromosomal replication begins with the formation of large DnaA-ATP
%   polymers, totaling approximately 30 DnaA molecues, at several sites denoted
%   R1-5 near the OriC. This process simulates the binding and unbinding of
%   DnaA-ATP and DnaA-ADP to these and 2000 additional sites throughout the
%   chromosome throughout the cell cycle. Although binding occurs throughout the
%   cell cycle, due to the cell's limited amount of DnaA, the titration affects
%   of the additional 2000 sites, and the cooperativity of DnaA polymerization
%   at the OriC sites, DnaA complexation formation at the OriC only occurs
%   approximately 2/3 through the cell cycle, providing robust control of
%   replication initiation.
%
%   DnaA Boxes (Mycoplasma genitalium)
%   =======================================
%   All the DnaA box positions based on the M. genitalium motif described in
%   Cordova 2002. 
%   - 9mer sites (high affinity) are the exact matches of the motif (and reverse
%     complement). 
%   - 8mer sites (medium affinity) are matches of the motif (and reverse
%     complement) with 1 incorrect base. 
%
%   We assume that in the oriC, 5 boxes are present: one 9mer, three
%   8mers, and one 7mer, to mimic E. coli's R1-R4, R5. These boxes should reside
%   between MG_469 and MG_470 (bases: 578581-579224). There are 9 (8-9mer) boxes
%   in this region, but we only recognize 4, so we ignore the boxes at positions
%   578837, 578855, 578881, 578966, and 579139.
%
%   R5 is a 7mer, so it is a very weak binder of DnaA. Essentially it is only
%   bound by cooperativity given the presence of the initiator complex. Since we
%   do not know its exact mechanism/purpose we say it just binds after the
%   complex is formed, and its binding triggers initiation.
%
%   Knowledge Base
%   =======================================
%   The DnaA boxes are represented in the knowledge base as genomic features and
%   loaded into this class by the initializeConstants method. The knowledge base
%   also contains the footprint sizes of the DnaA complexes; these are used by
%   the chromosomes object to determine whether DnaA complexes can bind to
%   specific chromosomal regions.
%
%   Representation
%   =======================================
%   The substrates, enzymes, and boundEnzymes properties represent the counts of
%   free metabolites, free DnaA, and DnaA bound to the chromosome. The
%   complexBoundSites property of the chromosomes object represent the specific
%   chromosomal locations of bound DnaA. The ATP/ADP bound and polymerization
%   status of each bound DnaA molecule is indicated by the specific identity of
%   the bound DnaA complex (DnaA-ATP 1-7 mer; DnaA-ADP + DnaA-ATP 0-6 mer).   
%
%   Initialization
%   =======================================
%   The process is initialized to a steady state by the initializeState method.
%   The steady state (amounts of free, 8mer/9mer bound DnaA-AxP) is found using
%   non-linear constrained optimization where we try to identify a state which
%   is a stable point and which maximizes the amount of 9-mer bound ATP. In the
%   initializeState method we make the simplifications that there is no free
%   DnaA (all DnaA is ATP or ADP bound) and that there are no DnaA polymers at
%   the functional R1-4 OriC boxes.
%
%   Simulation
%   =======================================
%   This process follows the general ideas in Atlas et al. 2008. The process
%   consists of several subfunctions executed in a deterministic order:
%   - Activate free DnaA to DnaA-ATP (activateFreeDnaA)
%     Deterministically form DnaA-ATP complexes upto the limit of available DnaA
%     monomers and ATP. The kinetics of DnaA activation are not known, and are
%     not modeled.
%   - Dissociate free DnaA-ATP polymers into monomers and hydrolyze ATPs (inactivateFreeDnaAATP)
%   - polymerized DnaA-ATP (polymerizeDnaAATP)
%     If chromosomes are supercoiled, stochastically polymerize R1-4 DnaA boxes
%     (which are bound by DnaA-ATP monomers/polymers (of up to length 6)) by 1
%     additional DnaA-ATP molecule at rate
%       kbATP * numFreeDnaAATP / V * C
%     where C is a cooperativity constant which depends on the polymerization
%     status of the other R1-4 boxes
%   - polymerized DnaA-ADP (polymerizeDnaAADP)
%     Similar to DnaA-ATP polymerization, but with slower kinetic rate, kbADP
%   - Bind DnaA-ATP (bindDnaAATP)
%     Stochastically bind DnaA-ATP to free DnaA boxes at rate
%       kbATP * numFreeDnaAATP / V
%   - Bind DnaA-ADP (bindDnaAADP)
%     Similar to DnaA-ATP binding, but with slower kinetic rate, kbADP
%   - Stochastically release bound DnaA-ATP with uniform probability (releaseDnaAAxP)
%     Stochastically release bound DnaA-ATP monomers, and stochastically
%     depolymerize R1-4 boxes (except those which have polymer lengths equal to
%     the minimum of that over the R1-4 boxes) at rate kd1ATP.
%   - Stochastically release bound DnaA-ADP with uniform probability (releaseDnaAAxP)
%     Stochastically release bound DnaA-ADP monomers, and stochastically
%     depolymerize R1-4 boxes at rate kd1ADP.
%   - Reactivate free DnaA from free DnaA-ADPs (reactivateFreeDnaAADP)
%     Deterministically reactivate free DnaA-ATP from free DnaA-ATP a rate
%        numFreeDnaAADP * (k_Regen * membraneConc) /
%                       (K_Regen_P4 + membraneConc)
%
%   Replication-dependent bound DnaA-ATP inactivation is modeled differently
%   here than by Atlas et al 2008. Atlas et, 2008 included a global term for the
%   effect of active beta-clamps on bound DnaA-ATP inactivation. Because this
%   model is evaluated as part of a larger model and in particular the exact
%   position of active beta-clamps are known, we are able to model the local
%   affects of beta-clamps on bound DnaA-ATP, which is to release the bound
%   protein from DNA. However, because we cannot distinguish free DnaA-ATP from
%   DnaA-ATP released by beta-clamps we only model the release of these proteins
%   from the DNA, and not their hydrolysis to DnaA-ADP.
%
%   References
%   =======================================
%   1. Atlas, J.C., Nikolaev, E.V., Browning, S.T., Shuler, M.L. (2008).
%      Incorporating genome-wide DNA sequence information into a dynamic
%      whole-cell model of E. coli: application to DNA replication. Systems
%      Biology, IET 2: 369-382.
%   2. Browning, S.T., Castellanos, M., Shuler, M.L. (2004). Robust control of
%      Initiation of prokaryotic chromosome replication: essential considerations
%      for a minimal cell. Biotechnology and Bioengineering 88: 575-584.
%      All rate constants are from Browning (2004).
%   3. Cordova, C.M.M., Lartigue, C., Sirand-Pugnet, P., Renaudin, J., Cunha,
%      R.A.F., Blanchard, A. (2002). Identification of the origin of replication
%      of the Mycoplasma pulmonis chromosome and its use in oriC replicative
%      plasmids. Journal of Bacteriology 184: 5426-5435.
%   4. Margulies, C., Kaguni, J.M. (1996). Ordered and sequential binding of DnaA
%      protein to oriC, the chromosomal origin of escherichia coli. Journal of
%      biological chemistry 271: 17035-17040.
%
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/10/2010
classdef ReplicationInitiation < edu.stanford.covert.cell.sim.Process & edu.stanford.covert.cell.sim.ChromosomeProcessAspect
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
		    'dnaABoxStartPositions';
            'siteCooperativity';
            'stateCooperativity';
            'kb1ATP';
            'kb2ATP';
            'kd1ATP';
            'kb1ADP';
            'kb2ADP';
            'kd1ADP';
            'k_Regen';
            'K_Regen_P4';
            'k_inact';
		    };
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {}; %names of simulation state properties redundant with timecourses in this or other processes or the simulation
    end

    %enumerations
    properties (Constant)
        dnaABoxStatus_NotExist     = -1; %site on the second chromosome that has not yet been replicated
        dnaABoxStatus_NotBound     =  0; %site not bound by DnaA-ATP or DnaA-ADP
        dnaABoxStatus_DnaAATPBound =  1; %site bound by DnaA-ATP
        dnaABoxStatus_DnaAADPBound =  2; %site bound by DnaA-ADP
    end

    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs = {}; %whole cell model IDs of stimuli
        
        substrateWholeCellModelIDs = { %whole cell model IDs of substrates
            'ATP';'ADP';'PI';'H2O';'H'};
        substrateIndexs_atp        = 1;
        substrateIndexs_adp        = 2;
        substrateIndexs_phosphate  = 3;
        substrateIndexs_water      = 4;
        substrateIndexs_hydrogen   = 5;
        
        enzymeWholeCellModelIDs = {    %whole cell model IDs of enzymes
            'MG_469_1MER_ADP'      %DnaA-ADP 1mer
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
            'MG_469_7MER_ATP'      %DnaA-ATP 7mer
            'MG_469_MONOMER'       %chromosomal replication initiator protein DnaA
            };
        enzymeIndexs_DnaA             = 15;
        enzymeIndexs_DnaA_1mer_ADP    = 1;
        enzymeIndexs_DnaA_1mer_ATP    = 2;
        enzymeIndexs_DnaA_Nmer_ATP    = (2:2:14)';
        enzymeIndexs_DnaA_Nmer_ADP    = (1:2:13)';
        enzymeIndexs_DnaA_polymer_ATP = (4:2:14)';
        enzymeIndexs_DnaA_polymer_ADP = (3:2:13)';
        
        dnaABoxIndexs_7mer   %indices of 7mer binding sites within dnaABindingSites
        dnaABoxIndexs_8mer   %indices of 8mer binding sites within dnaABindingSites
        dnaABoxIndexs_9mer   %indices of 9mer binding sites within dnaABindingSites
        dnaABoxIndexs_R12345 %indices of R1-5 functional binding sites within dnaABindingSites
        dnaABoxIndexs_R1234  %indices of R1-4 functional binding sites within dnaABindingSites
        dnaABoxIndexs_R5     %indices of R5 functional binding sites within dnaABindingSites
    end

    %fixed biological constants
    properties
        dnaABoxStartPositions%positions of all DnaA binding sites on chromosome
        
        siteCooperativity    %factor by which DnaAATP to oriC site binding probability increases when other sites are bound (70)
        stateCooperativity   %factor by which DnaAATP to oriC site binding probability increases when x*4 sized DnaA complex has formed at oriC (20)
        
        kb1ATP               %rate for DnaA-ATP binding high affinity DnaA boxes (25 nM/h) --ASSUME THESE ARE ACTUALLY 1/(nM*h) %Browning (2004)[PUB_0448]
        kb2ATP               %rate for DnaA-ATP binding medium affinity DnaA boxes (0.6 nM/h) --ASSUME THESE ARE ACTUALLY 1/(nM*h) %Browning (2004)[PUB_0448]
        kd1ATP               %rate for DnaA-ATP dissociating from dna (20 1/h) %Browning (2004)[PUB_0448]
        kb1ADP               %rate for DnaA-ADP binding high affinity DnaA boxes (2.5 nM/h) --ASSUME THESE ARE ACTUALLY 1/(nM*h) %Browning (2004)[PUB_0448]
        kb2ADP               %rate for DnaA-ADP binding medium affinity DnaA boxes (0.61 nM/h) --ASSUME THESE ARE ACTUALLY 1/(nM*h) %Browning (2004)[PUB_0448]
        kd1ADP               %rate for DnaA-ADP dissociating from dna (20 1/h) %Browning (2004)[PUB_0448]
        k_Regen              %rate for DnaA-ADP to DnaA-ATP regeneration (2.3026 1/h)
        K_Regen_P4           %rate for DnaA-ADP to DnaA-ATP regeneration catalyzed by membrane lipids (0.018 g/L) (Atlas 2008 xml file)[PUB_0447]
        k_inact              %rate for DnaA-ATP to DnaA-ADP inactivation (4.24e14 1/s) (Browning 2004) %units unknown, likely 1/g [PUB_0448]
        
        dnaARelease_remainingDnaAIndexs %local indices of DnaA-AxP complexes remaining after bound DnaA-AxP complex release
        dnaARelease_remainingDnaAMatrix %DnaA-AxP complexes remaining after bound DnaA-AxP complex release [remaining X bound]
        dnaARelease_releasedDnaAIndexs  %local indices of DnaA-AxP complexes released by DnaA-AxP complex release
        dnaARelease_releasedDnaAMatrix  %DnaA-AxP complexes released by bound DnaA-AxP complex release [released X bound]
    end
    
    %state references
    properties
        mass
    end

    methods
        %constructor
        function this = ReplicationInitiation(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.Process(wholeCellModelID, name);
        end
    end

    %communication between process/simulation
    methods
        function storeObjectReferences(this, simulation)
            this.storeObjectReferences@edu.stanford.covert.cell.sim.Process(simulation);
            this.storeObjectReferences@edu.stanford.covert.cell.sim.ChromosomeProcessAspect(simulation);
            
            this.mass = simulation.state('Mass');
            this.states = [this.states; {this.mass}];
        end
        
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)
            this.initializeConstants@edu.stanford.covert.cell.sim.Process(...
                knowledgeBase, simulation, varargin{:});
            this.initializeConstants@edu.stanford.covert.cell.sim.ChromosomeProcessAspect(...
                knowledgeBase, simulation, varargin{:});
            
            %DnaA binding sites
            dnaABoxes = findobj(knowledgeBase.genomeFeatures, 'type', 'DnaA box');
            dnaABoxSubtypes = {dnaABoxes.subtype}';
            dnaABoxNames = {dnaABoxes.name}';
            
            this.dnaABoxStartPositions = ...
                ceil([dnaABoxes.startCoordinate]' + [dnaABoxes.sequenceLength]'/2 - this.enzymeDNAFootprints(this.enzymeIndexs_DnaA_1mer_ATP)/2);
            
            this.dnaABoxIndexs_7mer = find(strcmp(dnaABoxSubtypes, '7mer'));
            this.dnaABoxIndexs_8mer = find(strcmp(dnaABoxSubtypes, '8mer'));
            this.dnaABoxIndexs_9mer = find(strcmp(dnaABoxSubtypes, '9mer'));
            
            [~, this.dnaABoxIndexs_R12345] = ismember({
                'Functional box R1';
                'Functional box R2';
                'Functional box R3';
                'Functional box R4';
                'Functional box R5';
                }, dnaABoxNames);
            this.dnaABoxIndexs_R1234 = this.dnaABoxIndexs_R12345(1:4);
            this.dnaABoxIndexs_R5 = this.dnaABoxIndexs_R12345(5);
            
            %DnaA-AxP release reactions
            this.calcDnaAxpReleaseReactions();
        end
        
        function sampleDnaABoxes(this, nBoxes)
            validateattributes(nBoxes, {'numeric'}, {...
                'integer', ...
                '>=', numel(this.dnaABoxIndexs_R12345), ...
                '<=', numel(this.dnaABoxStartPositions)});
            
            idxs = sort([this.randStream.randomlySelectNRows(...
                setdiff((1:numel(this.dnaABoxStartPositions))', this.dnaABoxIndexs_R12345), ...
                nBoxes - numel(this.dnaABoxIndexs_R12345))
                this.dnaABoxIndexs_R12345]);
            
            this.dnaABoxStartPositions(idxs);
            
            [tfs, tmpIdxs] = ismember(this.dnaABoxIndexs_7mer, idxs);
            this.dnaABoxIndexs_7mer = tmpIdxs(tfs);            
            [tfs, tmpIdxs] = ismember(this.dnaABoxIndexs_8mer, idxs);
            this.dnaABoxIndexs_8mer = tmpIdxs(tfs);            
            [tfs, tmpIdxs] = ismember(this.dnaABoxIndexs_9mer, idxs);
            this.dnaABoxIndexs_9mer = tmpIdxs(tfs);
            [tfs, tmpIdxs] = ismember(this.dnaABoxIndexs_R12345, idxs);
            this.dnaABoxIndexs_R12345 = tmpIdxs(tfs);
            [tfs, tmpIdxs] = ismember(this.dnaABoxIndexs_R1234, idxs);
            this.dnaABoxIndexs_R1234 = tmpIdxs(tfs);
            [tfs, tmpIdxs] = ismember(this.dnaABoxIndexs_R5, idxs);
            this.dnaABoxIndexs_R5 = tmpIdxs(tfs);
        end
        
        function calcDnaAxpReleaseReactions(this)
            this.dnaARelease_remainingDnaAIndexs = zeros(size(this.enzymeComplexGlobalIndexs));
            this.dnaARelease_remainingDnaAIndexs(this.enzymeIndexs_DnaA_1mer_ATP) = 0;
            this.dnaARelease_remainingDnaAIndexs(this.enzymeIndexs_DnaA_1mer_ADP) = 0;
            this.dnaARelease_remainingDnaAIndexs(this.enzymeIndexs_DnaA_Nmer_ATP(2:end)) = this.enzymeIndexs_DnaA_Nmer_ATP(1:end-1);
            this.dnaARelease_remainingDnaAIndexs(this.enzymeIndexs_DnaA_Nmer_ADP(2:end)) = this.enzymeIndexs_DnaA_Nmer_ATP(1:end-1);
            
            this.dnaARelease_remainingDnaAMatrix = zeros(numel(this.enzymeWholeCellModelIDs)); %remaining bound X original bound
            this.dnaARelease_remainingDnaAMatrix(sub2ind(size(this.dnaARelease_remainingDnaAMatrix), this.enzymeIndexs_DnaA_Nmer_ATP(1:end-1), this.enzymeIndexs_DnaA_Nmer_ATP(2:end))) = 1;
            this.dnaARelease_remainingDnaAMatrix(sub2ind(size(this.dnaARelease_remainingDnaAMatrix), this.enzymeIndexs_DnaA_Nmer_ATP(1:end-1), this.enzymeIndexs_DnaA_Nmer_ADP(2:end))) = 1;
            
            this.dnaARelease_releasedDnaAIndexs = zeros(size(this.enzymeComplexGlobalIndexs));
            this.dnaARelease_releasedDnaAIndexs(this.enzymeIndexs_DnaA_1mer_ATP) = this.enzymeIndexs_DnaA_1mer_ATP;
            this.dnaARelease_releasedDnaAIndexs(this.enzymeIndexs_DnaA_1mer_ADP) = this.enzymeIndexs_DnaA_1mer_ADP;
            this.dnaARelease_releasedDnaAIndexs(this.enzymeIndexs_DnaA_Nmer_ATP(2:end)) = this.enzymeIndexs_DnaA_1mer_ATP;
            this.dnaARelease_releasedDnaAIndexs(this.enzymeIndexs_DnaA_Nmer_ADP(2:end)) = this.enzymeIndexs_DnaA_1mer_ADP;
            
            this.dnaARelease_releasedDnaAMatrix = zeros(numel(this.enzymeWholeCellModelIDs)); %released X original bound
            this.dnaARelease_releasedDnaAMatrix(this.enzymeIndexs_DnaA_1mer_ATP, this.enzymeIndexs_DnaA_1mer_ATP) = 1;
            this.dnaARelease_releasedDnaAMatrix(this.enzymeIndexs_DnaA_1mer_ADP, this.enzymeIndexs_DnaA_1mer_ADP) = 1;
            this.dnaARelease_releasedDnaAMatrix(this.enzymeIndexs_DnaA_1mer_ATP, this.enzymeIndexs_DnaA_Nmer_ATP(2:end)) = 1;
            this.dnaARelease_releasedDnaAMatrix(this.enzymeIndexs_DnaA_1mer_ADP, this.enzymeIndexs_DnaA_Nmer_ADP(2:end)) = 1;
        end
    end

    %model
    methods
        %Calculate
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, ~, states)
            %substrate and byproducts
            bmProd = zeros(size(this.substrateWholeCellModelIDs));
            byProd = zeros(size(this.substrateWholeCellModelIDs));
            
            nATPInc = (1:7) * states.complexProductions(this.enzymeGlobalIndexs(this.enzymeIndexs_DnaA_Nmer_ATP)); %activated newly produced DnaA
            nADPInc = (1:7) * states.complexProductions(this.enzymeGlobalIndexs(this.enzymeIndexs_DnaA_Nmer_ADP)); %activated newly produced DnaA
            nATPHyd = 4 * numel(this.enzymeIndexs_DnaA_Nmer_ATP);                                                  %regenerate polymerized DnaA-ADP released by replisome from OriC DnaA boxes (Note: the model makes the simplification that only polymerized DnaA-ATP is deactivated)
            bmProd(this.substrateIndexs_atp)       = nATPHyd + nATPInc;
            bmProd(this.substrateIndexs_water)     = nATPHyd;
            byProd(this.substrateIndexs_adp)       = nATPHyd - nADPInc;
            byProd(this.substrateIndexs_phosphate) = nATPHyd;
            byProd(this.substrateIndexs_hydrogen)  = nATPHyd;
            
            %current level of enzymes, so that parameters don't need to be refit
            %to be consistent with replication initiation duration
            minEnzExp = this.enzymes;
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
        end
        
        %initialization
        function initializeState(this)
            import edu.stanford.covert.util.ComputationUtil;
            
            %unbind all DnaA-AxP
            this.enzymes(this.enzymeIndexs_DnaA_Nmer_ATP) = this.enzymes(this.enzymeIndexs_DnaA_Nmer_ATP) + this.boundEnzymes(this.enzymeIndexs_DnaA_Nmer_ATP);
            this.enzymes(this.enzymeIndexs_DnaA_Nmer_ADP) = this.enzymes(this.enzymeIndexs_DnaA_Nmer_ADP) + this.boundEnzymes(this.enzymeIndexs_DnaA_Nmer_ADP);
            this.boundEnzymes(this.enzymeIndexs_DnaA_Nmer_ATP) = 0;
            this.boundEnzymes(this.enzymeIndexs_DnaA_Nmer_ADP) = 0;
            
            %break down any DnaA-AxP
            this.substrates(this.substrateIndexs_atp) = ...
                + this.substrates(this.substrateIndexs_atp) ...
                + (1:7) * this.enzymes(this.enzymeIndexs_DnaA_Nmer_ATP);
            this.substrates(this.substrateIndexs_adp) = ...
                + this.substrates(this.substrateIndexs_adp) ...
                + (1:7) * this.enzymes(this.enzymeIndexs_DnaA_Nmer_ADP);
            this.enzymes(this.enzymeIndexs_DnaA) = ...
                + this.enzymes(this.enzymeIndexs_DnaA) ...
                + (1:7) * this.enzymes(this.enzymeIndexs_DnaA_Nmer_ATP) ...
                + (1:7) * this.enzymes(this.enzymeIndexs_DnaA_Nmer_ADP);
            this.enzymes(this.enzymeIndexs_DnaA_Nmer_ATP) = 0;
            this.enzymes(this.enzymeIndexs_DnaA_Nmer_ADP) = 0;
            
            %initialize to ATP/ADP and bound/unbound state
            this.initializeStateBasedOnFinalConditions();
            %this.initializeStateBasedOnTheory();
        end
        
        function initializeStateBasedOnFinalConditions(this)
            DnaA_total = this.enzymes(this.enzymeIndexs_DnaA);
            
            isOriComplexFormed = (DnaA_total >= 50) || this.randStream.stochasticRound(0.5 * (DnaA_total >= 47));
            nonOriCDnaAATP_total = DnaA_total - 29 * isOriComplexFormed;
            f = this.kb1ATP * numel(this.dnaABoxIndexs_9mer) / (this.kb1ATP * numel(this.dnaABoxIndexs_9mer) + this.kb2ATP * numel(this.dnaABoxIndexs_8mer));
            nonOriCDnaAATP_9mer = this.randStream.stochasticRound(f * nonOriCDnaAATP_total);
            nonOriCDnaAATP_8mer = nonOriCDnaAATP_total - nonOriCDnaAATP_9mer;
            
            %update DnaA
            this.enzymes(:) =  0;            
            this.enzymes(this.enzymeIndexs_DnaA_1mer_ATP) = nonOriCDnaAATP_total + isOriComplexFormed;
            this.enzymes(this.enzymeIndexs_DnaA_Nmer_ATP(end)) = 4 * isOriComplexFormed;
            this.boundEnzymes(:) = 0;
            
            %bind DnaA to chromosome
            positionsStrands_7mer = [this.dnaABoxStartPositions(this.dnaABoxIndexs_R1234) ones(size(this.dnaABoxIndexs_R1234))];
            positionsStrands_1mer = [this.dnaABoxStartPositions(this.dnaABoxIndexs_R5) ones(size(this.dnaABoxIndexs_R5))];
            tfs7 = this.bindProteinToChromosome(positionsStrands_7mer, this.enzymeIndexs_DnaA_Nmer_ATP(end), 4 * isOriComplexFormed, [], true, false, 1, false, [], true);
            tfs1 = this.bindProteinToChromosome(positionsStrands_1mer, this.enzymeIndexs_DnaA_1mer_ATP, isOriComplexFormed, [], true, false, 1, false, [], true);
            if sum(tfs7) ~= 4 * isOriComplexFormed || sum(tfs1) ~= isOriComplexFormed
                throw(MException('ReplicationInitiation:error', 'DnaA-ATP should be bound to oriC sites'));
            end
            
            positionsStrands_9mer = [this.dnaABoxStartPositions(this.dnaABoxIndexs_9mer) ones(size(this.dnaABoxIndexs_9mer))];
            positionsStrands_8mer = [this.dnaABoxStartPositions(this.dnaABoxIndexs_8mer) ones(size(this.dnaABoxIndexs_8mer))];
            tfs9 = this.bindProteinToChromosome(positionsStrands_9mer, this.enzymeIndexs_DnaA_1mer_ATP, nonOriCDnaAATP_9mer, [], true, false, 1, false, [], true);
            tfs8 = this.bindProteinToChromosome(positionsStrands_8mer, this.enzymeIndexs_DnaA_1mer_ATP, nonOriCDnaAATP_8mer, [], true, false, 1, false, [], true);
            if sum(tfs9) ~= nonOriCDnaAATP_9mer || sum(tfs8) ~= nonOriCDnaAATP_8mer
                throw(MException('ReplicationInitiation:error', 'DnaA-ATP should be bound to 8- and 9-mer sites'));
            end
        end
        
        function initializeStateBasedOnTheory(this)
            %compute steady state
            DnaA_total = this.enzymes(this.enzymeIndexs_DnaA);
            x0 = [
                DnaA_total; %free DnaA-ATP
                0;          %9mer DnaA-ATP
                0];         %8mer DnaA-ATP
            if any(x0)
                x = this.initializeState_cubicRoot(x0);
            else
                x = x0;
            end
            
            %round to nearest integer-value solution
            x(2) = ComputationUtil.roundHalfUp(x(2)); 
            x(3) = ComputationUtil.roundHalfDown(x(3));
            x(3) = floor(x(3)) + (rem(x(3),1) >  0.5); %round half down
            x(1) = DnaA_total - sum(x(2:end));
            if any(x < 0)
                throw(MException('ReplicationInitiation:error', 'all values must be positive'));
            end            
            DnaA_ATP_total = sum(x([1 2 3]));
            DnaA_ADP_total = 0;
            
            %update DnaA
            this.enzymes(this.enzymeIndexs_DnaA) =  DnaA_total - DnaA_ATP_total - DnaA_ADP_total;
            this.enzymes(this.enzymeIndexs_DnaA_1mer_ATP) = DnaA_ATP_total;
            this.enzymes(this.enzymeIndexs_DnaA_1mer_ADP) = DnaA_ADP_total;
            
            %bind DnaA to chromosome
            positionsStrands_9mer = [this.dnaABoxStartPositions(this.dnaABoxIndexs_9mer) ones(size(this.dnaABoxIndexs_9mer))];
            positionsStrands_8mer = [this.dnaABoxStartPositions(this.dnaABoxIndexs_8mer) ones(size(this.dnaABoxIndexs_8mer))];
            tfs9 = this.bindProteinToChromosome(positionsStrands_9mer, this.enzymeIndexs_DnaA_1mer_ATP, x(2), [], true, false, 1, false, [], true);
            tfs8 = this.bindProteinToChromosome(positionsStrands_8mer, this.enzymeIndexs_DnaA_1mer_ATP, x(3), [], true, false, 1, false, [], true);
            if sum(tfs9) ~= x(2) || sum(tfs8) ~= x(3)
                throw(MException('ReplicationInitiation:error', 'DnaA-ATP should be bound to 8- and 9-mer sites'));
            end
            
            this.bindProteinToChromosome(positionsStrands_9mer(~tfs9, :), this.enzymeIndexs_DnaA_1mer_ADP, 0, true, false, 1, false, [], true);
            this.bindProteinToChromosome(positionsStrands_8mer(~tfs8, :), this.enzymeIndexs_DnaA_1mer_ADP, 0, true, false, 1, false, [], true);
        end
        
        function x = initializeState_cubicRoot(this, x)            
            Navo = edu.stanford.covert.util.ConstantUtil.nAvogadro;
            V = this.geometry.volume;
            
            k8 = this.kb2ATP * 1e9 / 3600 / Navo / V;
            k9 = this.kb1ATP * 1e9 / 3600 / Navo / V;
            kd = this.kd1ATP / 3600;
            N8 = numel(this.dnaABoxIndexs_8mer);
            N9 = numel(this.dnaABoxIndexs_9mer);
            
            Xt = sum(x);
                        
            a = -k8*k9;
            b = Xt*k8*k9 - kd*(k8+k9) - k8*k9*(N8+N9);
            c = Xt*kd*(k8+k9) - kd^2 - kd*(k9*N9+k8*N8);
            d = Xt*kd^2; 
            
            X0 = edu.stanford.covert.util.ComputationUtil.cubicfcn(a, b, c, d);
            
            X0 = X0(1);
            X8 = k8*N8*X0 / (k8*X0 + kd);
            X9 = k9*N9*X0 / (k9*X0 + kd);
            
            x = [X0; X9; X8];
        end
        
        %resource requirements
        function result = calcResourceRequirements_Current(this)
            result = zeros(size(this.substrates));
            result(this.substrateIndexs_atp) = ...
                + this.enzymes(this.enzymeIndexs_DnaA) ...          %number of free DnaA
                + this.enzymes(this.enzymeIndexs_DnaA_1mer_ADP) ... %number of DnaA-ADP to be regenerated
                + this.boundEnzymes(this.enzymeIndexs_DnaA_1mer_ADP);
            result(this.substrateIndexs_water) = ...
                (2:7) * this.enzymes(this.enzymeIndexs_DnaA_polymer_ATP);
        end

        %simulation
        function evolveState(this)
            [polATPs, polADPs] = this.calculateDnaAR1234Polymerization();
            
            %Activate free DnaA to DnaA-ATP
            this.activateFreeDnaA();
            
            %Dissociate free DnaA-ATP polymers into monomers and hydrolyze ATPs
            this.inactivateFreeDnaAATP();
            
            %DnaA binding is reduced when DNA isn't supercoiled (Margulies et
            %al, 1996)
            if collapse(this.chromosome.supercoiled)
                %bind and polymerized DnaA-ATP
                [polATPs, polADPs] = this.bindAndPolymerizeDnaAATP(polATPs, polADPs);
                
                %bind polymerized DnaA-ADP
                [polATPs, polADPs] = this.bindAndPolymerizeDnaAADP(polATPs, polADPs);
            end
            
            %Stochastically release bound DnaA-AxP with uniform probability
            this.releaseDnaAAxP(polATPs, polADPs);
            
            %Reactivate free DnaA from free DnaA-ADPs
            this.reactivateFreeDnaAADP();
        end
    end

    %model helper functions
    methods       
        %Activate free DnaA to DnaA-ATP
        function activateFreeDnaA(this)
            nActivations = min(...
                this.substrates(this.substrateIndexs_atp), ...
                this.enzymes(this.enzymeIndexs_DnaA));
            
            %stop early if no activation
            if nActivations == 0
                return; 
            end
            
            %update DnaA, DnaA-ATP, ATP
            this.substrates(this.substrateIndexs_atp)     = this.substrates(this.substrateIndexs_atp)     - nActivations;
            this.enzymes(this.enzymeIndexs_DnaA)          = this.enzymes(this.enzymeIndexs_DnaA)          - nActivations;
            this.enzymes(this.enzymeIndexs_DnaA_1mer_ATP) = this.enzymes(this.enzymeIndexs_DnaA_1mer_ATP) + nActivations;
        end
        
        %Dissociate free DnaA-ATP polymers into monomers and hydrolyze ATPs
        %
        %These polymers are normally only found bound to functional DnaA
        %boxes near the ORI. However if they are released from the
        %chromosome by another protein, for example by DNA polymerase at the
        %start of replication, they could become free. These freed polymers
        %must dissociate into individual DnaA-ADP molecules before they can
        %again bind to the chromosome.
        function inactivateFreeDnaAATP(this)
            nDissociatingPolymers = this.enzymes(this.enzymeIndexs_DnaA_polymer_ATP);
            
            %limit by water availability
            while true
                nDissociatingMonomers = (2:7) * nDissociatingPolymers;
                if nDissociatingMonomers == 0
                    return;
                end
                if this.substrates(this.substrateIndexs_water) >= nDissociatingMonomers
                    break;
                end
                idx = this.randStream.randsample(6, 1, false, nDissociatingPolymers);
                nDissociatingPolymers(idx) = nDissociatingPolymers(idx) - 1;
            end

            %update DnaA-ATP polymers, DnaA-ADP, water, phosphate, hydrogen
            this.substrates(this.substrateIndexs_water)      = this.substrates(this.substrateIndexs_water)      - nDissociatingMonomers;
            this.substrates(this.substrateIndexs_phosphate)  = this.substrates(this.substrateIndexs_phosphate)  + nDissociatingMonomers;
            this.substrates(this.substrateIndexs_hydrogen)   = this.substrates(this.substrateIndexs_hydrogen)   + nDissociatingMonomers;
            this.enzymes(this.enzymeIndexs_DnaA_1mer_ADP)    = this.enzymes(this.enzymeIndexs_DnaA_1mer_ADP)    + nDissociatingMonomers;
            this.enzymes(this.enzymeIndexs_DnaA_polymer_ATP) = this.enzymes(this.enzymeIndexs_DnaA_polymer_ATP) - nDissociatingPolymers;
        end
        
        function [polATPs, polADPs] = bindAndPolymerizeDnaAATP(this, polATPs, polADPs)
            %% number of free DnaA-ATP
            numFreeDnaAATP = this.enzymes(this.enzymeIndexs_DnaA_1mer_ATP);
            if numFreeDnaAATP == 0
                return;
            end
            
            %% upper bound of binding
            %rates of binding each DnaA box type
            [positionsStrands, bindingRates, avgBindingRate] = this.calculateDnaAATPBindingRates(polATPs, polADPs);
            
            %Estimate number of free DnaA boxes
            [~, complexs] = find(this.chromosome.complexBoundSites);
            numFreeBindingSites = numel(bindingRates) ...
                - sum(complexs == this.enzymeGlobalIndexs(this.enzymeIndexs_DnaA_1mer_ATP)) ...
                - sum(complexs == this.enzymeGlobalIndexs(this.enzymeIndexs_DnaA_1mer_ADP));
            
            %Maximum number of DnaA-ATP binding
            maxBinding = this.randStream.stochasticRound(min([
                avgBindingRate * numFreeBindingSites * numFreeDnaAATP + nnz(bindingRates(this.dnaABoxIndexs_7mer, :));
                numFreeBindingSites;
                numFreeDnaAATP]));
            
            %% upper bound of polymerization
            %polymerization state of each function DnaA box
            maxPol = 0;
            polATPRates = zeros(size(polATPs));
            if any(polATPs(:))
                
                %Rates of polymerizing each DnaA box
                polATPRates = this.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs);
                
                %number of site available for polymerization
                numFreePolSites = nnz(polATPRates);
                if numFreePolSites > 0
                    
                    %Calculate total number of DnaA-ATP molecules that bind an
                    %accessible DnaA box
                    maxPol = this.randStream.stochasticRound(min([
                        sum(polATPRates(:)) * numFreeDnaAATP;
                        numFreePolSites;
                        numFreeDnaAATP]));
                end
            end
            
            %% allocate DnaA-ATP
            totBindingRate = avgBindingRate * numFreeBindingSites;
            totPolRate = sum(polATPRates(:));
            
            maxBinding = min(maxBinding, this.randStream.stochasticRound(numFreeDnaAATP * totBindingRate / (totBindingRate + totPolRate)));
            
            %% bind, polymerize
            [polATPs, polADPs, nBound] = this.bindDnaAATP(polATPs, polADPs, maxBinding, positionsStrands, bindingRates);
            [polATPs, polADPs] = this.polymerizeDnaAATP(polATPs, polADPs, min(maxPol, numFreeDnaAATP - nBound), polATPRates);
        end
        
        function [polATPs, polADPs] = bindAndPolymerizeDnaAADP(this, polATPs, polADPs)
            %% number of free DnaA-ADP
            numFreeDnaAADP = this.enzymes(this.enzymeIndexs_DnaA_1mer_ADP);
            if numFreeDnaAADP == 0
                return;
            end
            
            %% binding
            %rates of binding each DnaA box type
            [positionsStrands, bindingRates, avgBindingRate] = this.calculateDnaAADPBindingRates(polATPs, polADPs);
            
            %Estimate number of free DnaA boxes
            [~, complexs] = find(this.chromosome.complexBoundSites);
            numFreeBindingSites = numel(bindingRates) ...
                - sum(complexs == this.enzymeGlobalIndexs(this.enzymeIndexs_DnaA_1mer_ATP)) ...
                - sum(complexs == this.enzymeGlobalIndexs(this.enzymeIndexs_DnaA_1mer_ADP));
            
            %Maximum number of DnaA-ADP binding
            maxBinding = this.randStream.stochasticRound(min([
                avgBindingRate * numFreeBindingSites * numFreeDnaAADP + nnz(bindingRates(this.dnaABoxIndexs_7mer, :));
                numFreeBindingSites;
                numFreeDnaAADP]));
            
            %% polymerization
            %polymerization state of each function DnaA box
            maxPol = 0;
            polADPRates = zeros(size(polADPs));
            if any(polATPs(:))
                
                %Rates of polymerizing each DnaA box
                polADPRates = this.calculateDnaAR1234ADPPolymerizationRates(polATPs, polADPs);
                
                %number of site available for polymerization
                numFreePolSites = nnz(polADPRates);
                if numFreePolSites > 0
                    
                    %Calculate total number of DnaA-ADP molecules that bind an
                    %accessible DnaA box
                    maxPol = this.randStream.stochasticRound(min([
                        sum(polADPRates(:)) * numFreeDnaAADP;
                        numFreePolSites;
                        numFreeDnaAADP]));
                end
            end
            
            %% allocation
            totBindingRate = avgBindingRate * numFreeBindingSites;
            totPolRate = sum(polADPRates(:));
            
            maxBinding = min(maxBinding, this.randStream.stochasticRound(numFreeDnaAADP * totBindingRate / (totBindingRate + totPolRate)));
            
            %% bind, polymerize
            [polATPs, polADPs, nBound] = this.bindDnaAADP(polATPs, polADPs, maxBinding, positionsStrands, bindingRates);
            [polATPs, polADPs] = this.polymerizeDnaAADP(polATPs, polADPs, min(maxPol, numFreeDnaAADP - nBound), polADPRates);
        end
                
        %polymerized DnaA-ATP
        %
        %Choose binding sites by sampling from CDF.
        %
        %Calculate maximum amount of DnaA-ATP that could form
        %Browning et al 2004, say that if there is enough ATP, all new
        %DnaA is assumed to be in the ATP form. The consequence of this is
        %that for most of the simulation there is pretty much only DnaA-ATP
        %and no DnaA-ADP. Our framework however doesn't distinguish newly
        %constructed DnaA-ATPs from those recently released from the
        %chromosome. Consequently, because we would like to model the
        %hydrolysis of ATP bound to recently released DnaA molecules, we
        %must model the binding of ATP to DnaA as occuring just prior to
        %binding to the chromosome.
        %
        %Note: to implement the Browning et al model, the hydrolysis of ATP
        %bound to recently released DnaA would have to be handled by the
        %setRegionProteinUnbound method of the Chromosome class.
        function [polATPs, polADPs] = polymerizeDnaAATP(this, polATPs, polADPs, maxPol, polATPRates)
            if maxPol == 0
                return;
            end
            
            %stochastically select DnaA boxes to polymerize according to rate distribution
            idxs = this.randStream.randsample(numel(polATPRates), maxPol, false, polATPRates(:));
            
            %update ATP polymerization status
            polATPs(idxs) = polATPs(idxs) + 1;
            
            %update chromosome state
            positionsStrands = [
                this.dnaABoxStartPositions(this.dnaABoxIndexs_R1234(mod(idxs - 1, 4) + 1), 1) ...
                2*ceil(idxs/4)-1];
            
            this.modifyProteinOnChromosome(positionsStrands, this.enzymeIndexs_DnaA_Nmer_ATP(polATPs(idxs)));
            
            counts = histc(polATPs(idxs), 2:7);
            this.enzymes(this.enzymeIndexs_DnaA_1mer_ATP)      = this.enzymes(this.enzymeIndexs_DnaA_1mer_ATP)      - maxPol;
            this.enzymes(this.enzymeIndexs_DnaA_Nmer_ATP(1:6)) = this.enzymes(this.enzymeIndexs_DnaA_Nmer_ATP(1:6)) - counts(:);
            this.enzymes(this.enzymeIndexs_DnaA_Nmer_ATP(2:7)) = this.enzymes(this.enzymeIndexs_DnaA_Nmer_ATP(2:7)) + counts(:);
        end
        
        %polymerized DnaA-ADP
        function [polATPs, polADPs] = polymerizeDnaAADP(this, polATPs, polADPs, maxPol, polADPRates)
            if maxPol == 0
                return;
            end
            
            %stochastically select DnaA boxes to polymerize according to rate distribution
            idxs = this.randStream.randsample(numel(polADPRates), maxPol, false, polADPRates(:));
            
            %update AxP polymerization status
            polADPs(idxs) = polATPs(idxs) + 1;
            polATPs(idxs) = 0;
            
            %update chromosome state
            positionsStrands = [
                this.dnaABoxStartPositions(this.dnaABoxIndexs_R1234(mod(idxs - 1, 4) + 1), 1) ...
                2*ceil(idxs/4)-1];
            
            this.modifyProteinOnChromosome(positionsStrands, this.enzymeIndexs_DnaA_Nmer_ADP(polADPs(idxs)));
            
            counts = histc(polADPs(idxs), 2:7);
            this.enzymes(this.enzymeIndexs_DnaA_1mer_ADP)      = this.enzymes(this.enzymeIndexs_DnaA_1mer_ADP)      - maxPol;
            this.enzymes(this.enzymeIndexs_DnaA_Nmer_ATP(1:6)) = this.enzymes(this.enzymeIndexs_DnaA_Nmer_ATP(1:6)) - counts(:);
            this.enzymes(this.enzymeIndexs_DnaA_Nmer_ADP(2:7)) = this.enzymes(this.enzymeIndexs_DnaA_Nmer_ADP(2:7)) + counts(:);
        end
        
        %Bind DnaA-ATP
        function [polATPs, polADPs, nBound] = bindDnaAATP(this, polATPs, polADPs, maxBinding, positionsStrands, bindingRates)
            if maxBinding == 0;
                nBound = 0;
                return;
            end
            
            %Bind DnaA-ATP stochastically to accessible DnaA boxes according to rate distribution
            [tfs, ~, ~, nBound] = this.bindProteinToChromosome(positionsStrands, ...
                this.enzymeIndexs_DnaA_1mer_ATP, maxBinding, bindingRates, true, false, 1, false, [], true);
            polATPs(tfs(this.dnaABoxIndexs_R1234, 1), 1) = 1;
        end
        
        %Bind DnaA-ADP
        function [polATPs, polADPs, nBound] = bindDnaAADP(this, polATPs, polADPs, maxBinding, positionsStrands, bindingRates)
            if maxBinding == 0;
                nBound = 0;
                return;
            end
            
            %Bind DnaA-ADP stochastically to accessible DnaA boxes according to rate distribution
            [tfs, ~, ~, nBound] = this.bindProteinToChromosome(positionsStrands, ...
                this.enzymeIndexs_DnaA_1mer_ADP, maxBinding, bindingRates, true, false, 1, false, [], true);
            polADPs(tfs(this.dnaABoxIndexs_R1234, 1), 1) = 1;
        end
                
        %Stochastically release/depolymerize bound DnaA-AxP with uniform probability
        function [polATPs, polADPs] = releaseDnaAAxP(this, polATPs, polADPs)
            %get bound DnaA complexes
            c = this.chromosome;
            [posStrnds, complexGblIdxs] = find(c.complexBoundSites);
            complexLclIdxs = ismembc2(complexGblIdxs, this.enzymeComplexGlobalIndexs);
            tfs = complexLclIdxs > 0;
            
            %protect completed polymerized layers at R1-4 positions
            minPol = this.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            
            if minPol(1) > 0
                protectedPositions1 = this.dnaABoxStartPositions(...
                    this.dnaABoxIndexs_R1234(polATPs(:, 1) == minPol(1)));
                if minPol(1) >= 7
                    protectedPositions1 = [protectedPositions1;
                        this.dnaABoxStartPositions(this.dnaABoxIndexs_R5, :)];
                end
                if ~isempty(protectedPositions1)
                    tfs(tfs) = ~edu.stanford.covert.util.SparseMat.ismember_subs(...
                        posStrnds(tfs, :), [protectedPositions1 ones(size(protectedPositions1))], [c.sequenceLen c.nCompartments]);
                end
            end
            
            if minPol(2) > 0
                protectedPositions2 = this.dnaABoxStartPositions(...
                    this.dnaABoxIndexs_R1234(polATPs(:, 2) == minPol(2)));
                if minPol(2) >= 7
                    protectedPositions2 = [protectedPositions2;
                        this.dnaABoxStartPositions(this.dnaABoxIndexs_R5, :)];
                end
                if ~isempty(protectedPositions2)
                    tfs(tfs) = ~edu.stanford.covert.util.SparseMat.ismember_subs(...
                        posStrnds(tfs, :), [protectedPositions2 3*ones(size(protectedPositions2))], [c.sequenceLen c.nCompartments]);                   
                end
            end
            
            %choose complexes to release
            tfs(tfs) = this.randStream.rand(sum(tfs), 1) < this.kd1ATP / 3600;
            if ~any(tfs)
                return;
            end
            complexLclIdxs = complexLclIdxs(tfs, 1);
            posStrnds = posStrnds(tfs, :);
            
            %update enzymes, bound enzymes
            remainingComplexLclIdxs = this.dnaARelease_remainingDnaAIndexs(complexLclIdxs);
            
            this.releaseProteinFromSites(posStrnds(remainingComplexLclIdxs == 0, :), false);
            
            posStrnds = posStrnds(remainingComplexLclIdxs ~= 0, :);
            complexLclIdxs = complexLclIdxs(remainingComplexLclIdxs ~= 0, 1);
            remainingComplexLclIdxs = remainingComplexLclIdxs(remainingComplexLclIdxs ~= 0, 1);
            releasedComplexLclIdxs = this.dnaARelease_releasedDnaAIndexs(complexLclIdxs);
            this.modifyProteinOnChromosome(posStrnds, remainingComplexLclIdxs);
            
            count = histc(complexLclIdxs, 1:numel(this.enzymes));
            remainingCount = histc(remainingComplexLclIdxs, 1:numel(this.enzymes));
            releasedCount = histc(releasedComplexLclIdxs, 1:numel(this.enzymes));
            this.enzymes = this.enzymes -count(:) + remainingCount(:) + releasedCount(:);
        end
    
        %Reactivate free DnaA-ATPs from free DnaA-ADPs, upto available DnaA-ADP and ATP
        %
        %Incorporating DnaA-ADP rejuvenation to DnaA-ATP as in Atlas 2008.
        %This reaction is promoted by the aciding phospholipids cadiolipin
        %and phosphatidylglycerol (yat-ming 1998, crooke 1992, sekimizu
        %1987).
        function reactivateFreeDnaAADP(this)
            %DnaA-ADP available for regeneration
            numFreeDnaAADP = this.enzymes(this.enzymeIndexs_DnaA_1mer_ADP);
            if numFreeDnaAADP == 0
                return; 
            end
            
            %upper bound on DnaA-ADP regeneration to DnaA-ATP
            membraneConc = this.mass.metaboliteWt(1, this.mass.compartment.membraneIndexs) / this.geometry.volume; %g/L
            numRegenerations = min([
                this.substrates(this.substrateIndexs_atp);
                numFreeDnaAADP;
                this.randStream.stochasticRound(numFreeDnaAADP * ...
                (this.k_Regen / 3600 * membraneConc) / (this.K_Regen_P4 + membraneConc) * ...
                this.stepSizeSec)]);
            
            %stop early if no regeneration
            if numRegenerations == 0
                return; 
            end
            
            %update DnaA-ATP, DnaA-ADP, ATP, ADP
            this.enzymes(this.enzymeIndexs_DnaA_1mer_ADP) = this.enzymes(this.enzymeIndexs_DnaA_1mer_ADP) - numRegenerations;
            this.enzymes(this.enzymeIndexs_DnaA_1mer_ATP) = this.enzymes(this.enzymeIndexs_DnaA_1mer_ATP) + numRegenerations;
            this.substrates(this.substrateIndexs_atp)     = this.substrates(this.substrateIndexs_atp)     - numRegenerations;
            this.substrates(this.substrateIndexs_adp)     = this.substrates(this.substrateIndexs_adp)     + numRegenerations;
        end
    end
    
    %model helper helper methods
    methods
        function polRates = calculateDnaAR1234ATPPolymerizationRates(this, polATPs, polADPs)
            %initialize rates
            polRates = zeros(4, 2);
            
            %complex size
            pol = this.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            
            %stop if no polymerization possible
            if ~any(pol)
                return;
            end
            
            %complex size range for cooperativity to apply
            polRange = min(6, max(1, pol));
            
            %Rates of polymerizing each DnaA box
            nAvo = edu.stanford.covert.util.ConstantUtil.nAvogadro;
            polRates(4,   :) = this.kb1ATP * 1e9 / 3600 / nAvo / this.geometry.volume * this.stepSizeSec;
            polRates(1:3, :) = this.kb2ATP * 1e9 / 3600 / nAvo / this.geometry.volume * this.stepSizeSec;
            polRates = polRates .* (polATPs == [polRange; polRange; polRange; polRange]);
            if ~any(polRates(:))
                return;
            end
            
            %Include effect of cooperativity
            polRates = polRates .* this.calculateDnaAR1234ATPPolymerizationCooperativity(polATPs, polADPs, pol);
        end
        
        function polRates = calculateDnaAR1234ADPPolymerizationRates(this, polATPs, polADPs)
            %initialize rates
            polRates = zeros(4, 2);
            
            %complex size
            pol = this.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            
            %stop if no polymerization possible
            if ~any(pol)
                return;
            end
            
            %complex size range for cooperativity to apply
            polRange = min(6, max(1, pol));
            
            %Rates of polymerizing each DnaA box
            nAvo = edu.stanford.covert.util.ConstantUtil.nAvogadro;
            polRates(4,   :) = this.kb1ADP * 1e9 / 3600 / nAvo / this.geometry.volume * this.stepSizeSec;
            polRates(1:3, :) = this.kb2ADP * 1e9 / 3600 / nAvo / this.geometry.volume * this.stepSizeSec;
            polRates = polRates .* (polATPs == [polRange; polRange; polRange; polRange]);
        end
        
        function cooperativity = calculateDnaAR1234ATPPolymerizationCooperativity(this, polATPs, polADPs, pol)
            %site cooperativity
            cooperativity = this.siteCooperativity * [
                polATPs(4, :) >  pol;
                all(polATPs([1 4], :) > [pol; pol], 1);
                all(polATPs([1 4], :) > [pol; pol], 1);
                any(polATPs > [pol; pol; pol; pol], 1);
                ];
            
            %additive state cooperativity
            %- we considered multiplicative and exponential models, but these
            %  give rise to higher variance in the replication initiation time
            cooperativity(4, :) = cooperativity(4, :) + this.stateCooperativity * pol; 
            
            %set cooperativity to 1 for sites which can't polymerize
            cooperativity = max(1, (polATPs == [pol; pol; pol; pol] & polADPs == 0) .* cooperativity);
        end
        
        %Rate which DnaA-ATP binds to each DnaA box
        function [positionsStrands, bindingRates, avgBindingRate] = calculateDnaAATPBindingRates(this, polATPs, polADPs)
            %8-mer, 9-mer, average binding rates
            nAvo = edu.stanford.covert.util.ConstantUtil.nAvogadro;
            rate9mer = this.kb1ATP * 1e9 / 3600 / nAvo / this.geometry.volume * this.stepSizeSec;
            rate8mer = this.kb2ATP * 1e9 / 3600 / nAvo / this.geometry.volume * this.stepSizeSec;
            avgBindingRate = (rate9mer*numel(this.dnaABoxIndexs_9mer) + rate8mer*numel(this.dnaABoxIndexs_8mer)) / ...
                numel(this.dnaABoxStartPositions);
            
            %functional DnaA complex size
            complexSize = this.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            
            %cooperativity
            cooperativity = this.calculateDnaAR1234ATPPolymerizationCooperativity(polATPs, polADPs, complexSize);
            
            %First chromosome positions, strands, binding rates
            %set rate of R5 to Inf if R1-4 fully bound
            %baseline rate of binding 7mer site low, and binding is HIGHLY
            %cooperative with complete occupancy of R1-4 sites
            positionsStrands = [this.dnaABoxStartPositions ones(size(this.dnaABoxStartPositions))];
            bindingRates = zeros(numel(this.dnaABoxStartPositions), 1);
            bindingRates(this.dnaABoxIndexs_9mer, 1) = rate9mer;
            bindingRates(this.dnaABoxIndexs_8mer, 1) = rate8mer;
            bindingRates(this.dnaABoxIndexs_R5,   1) = (complexSize(1) == 7) * realmax;
            
            bindingRates(this.dnaABoxIndexs_R1234) = ...
                bindingRates(this.dnaABoxIndexs_R1234) .* ...
                cooperativity(:, 1);
            
            %Second chromosome positions, strands, binding rates
            %- We ingore setting the R5 box rate for computational efficiency
            %  because a complex should never be form on the second chromosome
            %- We also ignore cooperativity on the second chromosome
            if collapse(this.chromosome.polymerizedRegions) > 2 * this.chromosome.sequenceLen
                positionsStrands2 = [this.dnaABoxStartPositions 3*ones(size(this.dnaABoxStartPositions))];
                bindingRates2 = zeros(numel(this.dnaABoxStartPositions), 1);
                bindingRates2(this.dnaABoxIndexs_9mer, 1) = rate9mer;
                bindingRates2(this.dnaABoxIndexs_8mer, 1) = rate8mer;
                bindingRates2(this.dnaABoxIndexs_R5,   1) = (complexSize(2) == 7) * realmax;
                
                bindingRates2(this.dnaABoxIndexs_R1234) = ...
                    bindingRates2(this.dnaABoxIndexs_R1234) .* ...
                    cooperativity(:, 2);
                
                if collapse(this.chromosome.polymerizedRegions) < 4 * this.chromosome.sequenceLen
                    [pos, len] = find(this.chromosome.doubleStrandedRegions);
                    len = len(pos(:, 2) == 4, 1);
                    pos = pos(pos(:, 2) == 4, 1);
                    tfs = this.dnaABoxStartPositions < pos(1) + len(1) | this.dnaABoxStartPositions > pos(end);
                    positionsStrands2 = positionsStrands2(tfs, :);
                    bindingRates2 = bindingRates2(tfs);
                end
                
                bindingRates = [bindingRates; bindingRates2];
                positionsStrands = [positionsStrands; positionsStrands2];
            end
        end
        
        %Rate which DnaA-ADP binds to each DnaA box
        function [positionsStrands, bindingRates, avgBindingRate] = calculateDnaAADPBindingRates(this, polATPs, polADPs)
            %8-mer, 9-mer, average binding rates
            nAvo = edu.stanford.covert.util.ConstantUtil.nAvogadro;
            rate9mer = this.kb1ADP * 1e9 / 3600 / nAvo / this.geometry.volume * this.stepSizeSec;
            rate8mer = this.kb2ADP * 1e9 / 3600 / nAvo / this.geometry.volume * this.stepSizeSec;
            avgBindingRate = (rate9mer*numel(this.dnaABoxIndexs_9mer) + rate8mer*numel(this.dnaABoxIndexs_8mer)) / numel(this.dnaABoxStartPositions);
            
            %functional DnaA complex size
            complexSize = this.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            
            %First chromosome positions, strands, binding rates
            %set rate of R5 to Inf if R1-4 fully bound
            %baseline rate of binding 7mer site low, and binding is HIGHLY
            %cooperative with complete occupancy of R1-4 sites
            positionsStrands = [this.dnaABoxStartPositions ones(size(this.dnaABoxStartPositions))];
            bindingRates = zeros(numel(this.dnaABoxStartPositions), 1);
            bindingRates(this.dnaABoxIndexs_9mer, 1) = rate9mer;
            bindingRates(this.dnaABoxIndexs_8mer, 1) = rate8mer;
            bindingRates(this.dnaABoxIndexs_R5,   1) = (complexSize(1) == 7) * realmax;
            
            %Second chromosome positions, strands, binding rates
            %We ingore setting the R5 box rate for computational efficiency
            %because a complex should never be form on the second chromosome
            if collapse(this.chromosome.polymerizedRegions) > 2 * this.chromosome.sequenceLen
                positionsStrands2 = [this.dnaABoxStartPositions 3*ones(size(this.dnaABoxStartPositions))];
                bindingRates2 = zeros(numel(this.dnaABoxStartPositions), 1);
                bindingRates2(this.dnaABoxIndexs_9mer, 1) = rate9mer;
                bindingRates2(this.dnaABoxIndexs_8mer, 1) = rate8mer;
                bindingRates2(this.dnaABoxIndexs_R5,   1) = (complexSize(2) == 7) * realmax;
                
                if collapse(this.chromosome.polymerizedRegions) < 4 * this.chromosome.sequenceLen
                    [pos, len] = find(this.chromosome.doubleStrandedRegions);
                    len = len(pos(:, 2) == 4, 1);
                    pos = pos(pos(:, 2) == 4, 1);
                    tfs = this.dnaABoxStartPositions < pos(1) + len(1) | this.dnaABoxStartPositions > pos(end);
                    positionsStrands2 = positionsStrands2(tfs, :);
                    bindingRates2 = bindingRates2(tfs);
                end
                
                bindingRates = [bindingRates; bindingRates2];
                positionsStrands = [positionsStrands; positionsStrands2];
            end
        end
        
        function [polATP, polADP] = calculateDnaAR1234Polymerization(this)
            polATP = zeros(4, 2);
            polADP = zeros(4, 2);
            c = this.chromosome;
            [posStrnds, complexs] = find(c.complexBoundSites);
            
            %ATP
            idxsATP = ismembc2(complexs, this.enzymeGlobalIndexs(this.enzymeIndexs_DnaA_Nmer_ATP));
            tfsATP = idxsATP > 0;
            idxsATP = idxsATP(tfsATP, 1);
            [tfs2, idxs2] = edu.stanford.covert.util.SparseMat.ismember_subs(...
                [this.dnaABoxStartPositions(this.dnaABoxIndexs_R1234) ones(4, 1)], ...
                posStrnds(tfsATP, :), [c.sequenceLen c.nCompartments]);
            polATP(tfs2, 1) = idxsATP(idxs2(tfs2, 1), 1);
            
            %ADP
            if nargout >= 2
                idxsADP = ismembc2(complexs, this.enzymeGlobalIndexs(this.enzymeIndexs_DnaA_Nmer_ADP));
                tfsADP = idxsADP > 0;
                idxsADP = idxsADP(tfsADP, 1);
                [tfs2, idxs2] = edu.stanford.covert.util.SparseMat.ismember_subs(...
                    [this.dnaABoxStartPositions(this.dnaABoxIndexs_R1234) ones(4, 1)], ...
                    posStrnds(tfsADP, :), [c.sequenceLen c.nCompartments]);
                polADP(tfs2, 1) = idxsADP(idxs2(tfs2, 1), 1);
            end
            
            %second chromosome
            if nnz(this.chromosome.polymerizedRegions) ~= 2
                %ATP
                [tfs2, idxs2] = edu.stanford.covert.util.SparseMat.ismember_subs(...
                    [this.dnaABoxStartPositions(this.dnaABoxIndexs_R1234) 3*ones(4, 1)], ...
                    posStrnds(tfsATP, :), [c.sequenceLen c.nCompartments]);
                polATP(tfs2, 2) = idxsATP(idxs2(tfs2, 1), 1);
                
                %ADP
                if nargout >= 2
                    [tfs2, idxs2] = edu.stanford.covert.util.SparseMat.ismember_subs(...
                        [this.dnaABoxStartPositions(this.dnaABoxIndexs_R1234) 3*ones(4, 1)], ...
                        posStrnds(tfsADP, :), [c.sequenceLen c.nCompartments]);
                    polADP(tfs2, 2) = idxsADP(idxs2(tfs2, 1), 1);
                end
            end
        end
        
        function siz = calculateDnaAR1234ComplexSize(~, polATPs, polADPs)
            siz = min(max(polATPs, polADPs-1), [], 1);
        end
        
        function tf = calcuateIsDnaAR5Occupied(this)
            tf = this.chromosome.complexBoundSites([
                this.dnaABoxStartPositions(this.dnaABoxIndexs_R5) 1
                this.dnaABoxStartPositions(this.dnaABoxIndexs_R5) 3])' == ...
                this.enzymeGlobalIndexs(this.enzymeIndexs_DnaA_1mer_ATP);
        end
        
        %Returns true if all functional DnaA boxes (R1-5) are maximally
        %occupied:
        %- R4   (9mer, high affinity)   is  bound by DnaA-ATP 7mer
        %- R1-3 (8mer, medium affinity) are bound by DnaA-ATP 7mers
        %- R5   (7mer, low affinity)    is  bound by DnaA-ATP 1mer
        function tf = calculateIsDnaAORIComplexAssembled(this)
            tf = all(this.calculateDnaAR1234Polymerization() >= [7 7 7 7; 7 7 7 7]', 1) & ...
                this.calcuateIsDnaAR5Occupied();
        end
        
        %approximate status of DnaA boxes (approximate to speed calculation;
        %approximate because doesn't use Chromosome isRegionPolymerized
        %method)
        %
        %A more exact calculation would replace the second and third sections
        %below with:
        %  status = repmat([this.dnaABoxStatus_NotBound this.dnaABoxStatus_NotExist], size(this.dnaABoxStartPositions, 1), 1);
        %  status(reshape(this.chromosome.isRegionPolymerized(positionsStrands, this.enzymeDNAFootprints(this.enzymeIndexs_DnaA_1mer_ATP)), [], 2)) = this.dnaABoxStatus_NotBound;
        function status = calculateDnaABoxStatus(this)
            %starting positions and strands of DnaA boxes
            positionsStrands = [repmat(this.dnaABoxStartPositions, 2, 1) reshape(repmat([1 3], size(this.dnaABoxStartPositions, 1), 1), [], 1)];
                        
            %initialize status
            status = repmat(this.dnaABoxStatus_NotBound, size(this.dnaABoxStartPositions, 1), 2);
            
            %find polymerized sites
            [starts1, lengths1] = find(this.chromosome.polymerizedRegions(:, 3));
            starts2 = find(this.chromosome.polymerizedRegions(:, 4));
            if isempty(starts1)
                if isempty(starts2)
                    status(:, 2) = this.dnaABoxStatus_NotExist;
                else
                    status(this.dnaABoxStartPositions < starts2(end), 2) = this.dnaABoxStatus_NotExist;
                end
            else
                if isempty(starts2)
                    status(this.dnaABoxStartPositions > starts1(1), 2) = this.dnaABoxStatus_NotExist;
                else
                    status(this.dnaABoxStartPositions > starts1(1) + lengths1(1) & this.dnaABoxStartPositions < starts2(end), 2) = this.dnaABoxStatus_NotExist;
                end
            end
            
            %find bound sites
            boundComplexs = reshape(this.chromosome.complexBoundSites(positionsStrands), [], 2);
            status(ismembc(boundComplexs, this.enzymeGlobalIndexs(this.enzymeIndexs_DnaA_Nmer_ATP))) = this.dnaABoxStatus_DnaAATPBound;
            status(ismembc(boundComplexs, this.enzymeGlobalIndexs(this.enzymeIndexs_DnaA_Nmer_ADP))) = this.dnaABoxStatus_DnaAADPBound;
        end
    end
end