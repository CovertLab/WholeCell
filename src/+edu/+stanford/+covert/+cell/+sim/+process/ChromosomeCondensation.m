%ChromosomeCondensation
%
% @wholeCellModelID Process_ChromosomeCondensation
% @name             Chromosome Condensation
% @description
%   Biology
%   =====================================================================
%   Chromosome segregation requires that the DNA be highly compacted.
%   Structural maintenance of chromosome (SMC) complexes are "V" shaped
%   proteins (with an head and two legs) that induce positive supercoils in
%   double stranded DNA (Porter et al., 2004). The complexes are believed to
%   work with a lock and key mechanism in which first DNA is looped around
%   the legs of the SMC complex, and then an ATP is bound between the two
%   tails to lock the SMC complex in place. The complexes bind and clamp the
%   DNA causing many loops in the DNA and compacting it. The loops around
%   each leg occupy 90bp. A loop of about 450bp forms between the two SMC
%   complex legs (Jensen and Shapiro, 2003, Strick and Kawaguchi, 2004).
%   Further, it has been inferred that there is about 1 SMC complex per every
%   7000bp (Jensen and Shapiro, 2003).
%
%   The initial chromosome is bound by SMC complexes averaging 7000bp spacing.
%   As the replication loop proceeds, the SMC complexes that it encounters
%   are displaced.  Once the DNA has been replicated, SMC complexes are
%   randomly bound to the DNA such that their spacing averages 7000bp. Each
%   SMC complex occupies 630bp, and no two SMC complexes can occupy the same
%   space. Both the decision of whether an SMC complex will bind at a given
%   time point and the binding location are random. SMC complexes do not fall
%   off the chromosomes unless due to the force of the replication loop. The
%   two chromosomes are tracked separately.
%
%   Knowledge Base
%   =====================================================================
%   The knowledge base contains the values of the parameters: smcSepNt and
%   smcSepProbCenter. The knowledge base also contains the measured DNA
%   footprint of each SMC complex. These footprints are loaded by the
%   Chromosome class from the knowledge base. This class retrieves the
%   footprint from the Chromosome class.
%
%   Representation
%   =====================================================================
%   substrates, enzymes, and boundEnzymes represent the counts of free
%   metabolites, free SMC complexes, and chromosome-bound SMC complexes. The
%   chromosomes property represents the specific base positions where the
%   chromosome-bound SMC complexes are located.
%
%   enzymeDNAFootprints represents the experimentally measured DNA footprint of each
%   SMC complex. smcSepNt represents the experimentally observed average SMC
%   complex spacing  [PUB_0517]. smcSepProbCenter is a parameter which controls
%   the SMC binding probability transfer function. smcSepProbCenter is not an
%   experimentally measured quantity. Rather, its value is pinned by several
%   constraints: implemented in the ChromosomeCondensation_Test class
%   testCalculateBindingPositionWithinRegion and testInitializeStateConverged
%   methods. These constraints are:
%   - Consistent with an SMC density of approximately 1/smSepNt
%   - Consistent with fast binding of multiple SMCs to large unbound regions
%   - Consistent with slow binding of SMCs to small unbound regions
%
%   Algorithm
%   =====================================================================
%   1. Calculate expected number of SMC complexes that should bind chromosomes
%      a. Calculate regions where SMCs can bind (regions either between SMC
%         complexes or between SMC complexes and replication bubble)
%      b. Compute the expected number of binding complexes as the ratio of the
%         sum of the lengths of the regions to the average SMC spacing
%   2. For 1 to minimum of free SMC complexs, ATP, and expected number of SMCs
%      that that should bind chromosome
%      A. Calculate regions where SMCs can bind (regions either between SMC
%         complexes or between SMC complexes and replication bubble)
%         a. Starting coordinate
%         b. Chromosomes
%         c. Lengths
%         d. Probability of SMC binding each region
%      B. Pick a region for SMC complex to bind
%      C. Pick a position within region for SMC complex to bind
%      D. Form SMC-ADP complex:
%         a. Decrement SMC. Increment SMC-ADP
%         b. Decrement ATP, H2O. Increment PI, H.
%      E. Bind SMC-ADP complex to chromosome
%         a. Decrement free SMC-ADP complex
%         b. Increment bound SMC-ADP complex
%
%   The probability that an SMC complex binds a region of length L is given by
%   the step function
%     p(L) = 1/smcSepNt * max(0, L/2-smcSepProbCenter)
%
%   The conditional probability that an SMC complex binds a position within
%   region assuming that the SMC complex is binding the region is given by
%    p(x) = 1/(L-2*smcSepProbCenter)  if x>smcSepProbCenter and x<L-smcSepProbCenter,
%           0                         otherwise
%
%   References
%   =====================================================================
%   1. Ullsperger, C., Cozzarelli, N.R. (1996). Contrasting enzymatic activities
%      of topoisomerase IV and DnA gyrase from Escherichia coli. Journal of Bio
%      Chem 271: 31549-31555.
%   2. Dekker, N.H., Viard, T., Bouthier de la Tour, C., Duguet, M., Bensimon,
%      D., Croquette, V. (2003). Thermophilic Topoisomerase I on a single DNA
%      molecule. Journal of molecular biology 329: 271-282.
%   3. Gore, J., Bryant, Z., Stone, M.D., Nollmann, M., Cozzarelli, N.R.,
%      Bustamante, C. (2006). Mechanochemical analysis of DNA gyrase using rotor
%      bead tracking. Nature 439: 100-104.
%   4. Bates, A. (2006). DNA Topoisomerases: Single Gyrase Caught in the Act.
%      Current Biology 16: 204-206.
%   5. Jensen, R.B, Shapiro, L. (2003). Cell-Cycle-Regulated Expression and
%      Subcellular Localization of the Caulobacter crescentus SMC Chromosome
%      Structural Protein. Journal of Bacteriology 185: 3068-3075. [PUB_0517]
%   6. Strick, T.R., Kawaguchi, T. (2004). Real-time detection of single-molecule
%      DNA compaction by condensing I. Current biology 14: 874-880.
%   7. Tadesse, S., Mascarenhas, J., Kosters, B., Hasilik, A., Graumann, P.L.
%      (2005). Genetic interaction of the SMC complex with topoisomerase IV in
%      Bacillus subtilis. Microbiology 151: 3729-3737.
%   8. Bloom, K., Joglekar, A. (2010). Towards building a chromosome segregation
%      machine. Nature 463: 446-456.
%   9. Porter, I.M., Khoudoli, G.A., Swedlow, J.R. (2004). Chromosome
%      condensation: DNA compaction in real time. Current Biology 14: 554-556.
%
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 6/9/2010
classdef ChromosomeCondensation < edu.stanford.covert.cell.sim.Process & edu.stanford.covert.cell.sim.ChromosomeProcessAspect
    
    %indices
    properties (Constant)
        optionNames__              = {}; %names of properties that are options
        fixedConstantNames__       = {   %names of fixed constant properties
            'smcSepNt';
            'smcSepProbCenter';
            };
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {}; %names of properties that are simulation state owned by the simulation or other processes
    end
    
    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs = {};
        
        substrateWholeCellModelIDs = {
            'ATP'; 'ADP'; 'PI'; 'H2O'; 'H'};
        substrateIndexs_atp       = 1;
        substrateIndexs_adp       = 2;
        substrateIndexs_phosphate = 3;
        substrateIndexs_water     = 4;
        substrateIndexs_hydrogen  = 5;
        
        enzymeWholeCellModelIDs = {
            'MG_213_214_298_6MER';      %Chromosome Segregation Protein SMC with SCP Proteins
            'MG_213_214_298_6MER_ADP';  %Chromosome Segregation Protein SMC with SCP Proteins-ADP
            };
        enzymeIndexs_SMC     = 1;
        enzymeIndexs_SMC_ADP = 2;
    end
   
    %fixed biological constants
    properties
        smcSepNt           %1 SMC complex per 7130 nt [PUB_0517]
        smcSepProbCenter   %where SMC binding probability transitions from 0 to 1, on the order of smcSepNt
    end
    
    methods
        %constructor
        function this = ChromosomeCondensation(wholeCellModelID, name)
            this@edu.stanford.covert.cell.sim.Process(wholeCellModelID, name);
        end
    end
    
    %communication between process/simulation
    methods
        %set references to state objects
        function storeObjectReferences(this, simulation)
            this.storeObjectReferences@edu.stanford.covert.cell.sim.Process(simulation);
            this.storeObjectReferences@edu.stanford.covert.cell.sim.ChromosomeProcessAspect(simulation);
        end
        
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)
            this.initializeConstants@edu.stanford.covert.cell.sim.Process(...
                knowledgeBase, simulation, varargin{:});
            this.initializeConstants@edu.stanford.covert.cell.sim.ChromosomeProcessAspect(...
                knowledgeBase, simulation, varargin{:});
        end       
        
        %Calculate 
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, constants, states)
            %% initialize
            bmProd = zeros(size(this.substrateWholeCellModelIDs));
            byProd = zeros(size(this.substrateWholeCellModelIDs));
            minEnzExp = zeros(size(this.enzymeWholeCellModelIDs));
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
            
            %% substrate and byproducts
            %ATP required for SMCs to bind chromosome
            %- 1 ATP for each new SMCs for new chromosome
            %- 1 ATP to rebind SMCs following displacement by DNA polymerase
            %  during replication
            %- 1 ATP to rebind SMCs following displacement by RNA polymerase
            %  during transcription throughout cell cycle
            %  Estimated by multiplying fraction of chromosome covered by
            %  SMCs by number of RNA polymerases and integrating over time.
            %  Factor of 0.2 added based on observing simulations
            nSMCs = sum(states.complexProductions(this.enzymeGlobalIndexs(this.enzymeIndexs_SMC_ADP)));
            nATP = 0.2 * nSMCs * this.enzymeDNAFootprints(this.enzymeIndexs_SMC_ADP) / this.chromosome.sequenceLen * ...
                sum(states.complexs0(this.complex.rnaPolymeraseIndexs)) * constants.states.Time.cellCycleLength / log(2);
            
            %ATP and hydrolysis byproducts (ATP + H2O ==> ADP + PI + H)
            bmProd(this.substrateIndexs_atp)       = nATP;
            bmProd(this.substrateIndexs_water)     = nATP;
            byProd(this.substrateIndexs_adp)       = nATP - nSMCs;
            byProd(this.substrateIndexs_phosphate) = nATP;
            byProd(this.substrateIndexs_hydrogen)  = nATP;
            
            %% enzymes: sufficient SMC proteins to predict spacing consistent with experimental observations
            minEnzExp(this.enzymeIndexs_SMC) = ...
                ceil(this.chromosome.sequenceLen / this.smcSepNt);
        end       
        
        %Initialize to steady SMC-bound state. A test class verifies convergence
        %to a steady-state.
        function initializeState(this)
            substrates = this.substrates;
            this.substrates(this.substrateIndexs_atp) = Inf;
            this.substrates(this.substrateIndexs_water) = Inf;
            for i = 1:20
                this.evolveState();
            end
            this.substrates = substrates;
        end
        
        %resource requirements
        function result = calcResourceRequirements_Current(this)
            result = zeros(size(this.substrates));
            result(this.substrateIndexs_atp)   = this.enzymes(this.enzymeIndexs_SMC);
            result(this.substrateIndexs_water) = this.enzymes(this.enzymeIndexs_SMC);
        end
        
        %simulation
        function evolveState(this)
            %% dissociate free SMC-ADP complexes
            this.enzymes(this.enzymeIndexs_SMC)       = this.enzymes(this.enzymeIndexs_SMC)       + this.enzymes(this.enzymeIndexs_SMC_ADP);
            this.substrates(this.substrateIndexs_adp) = this.substrates(this.substrateIndexs_adp) + this.enzymes(this.enzymeIndexs_SMC_ADP);
            this.enzymes(this.enzymeIndexs_SMC_ADP)   = 0;
            
            %% bind SMCs to DNA up to limit of available SMCs and energy
            nBindingMax = min([...
                this.substrates(this.substrateIndexs_atp) ...
                this.substrates(this.substrateIndexs_water) ...
                this.enzymes(this.enzymeIndexs_SMC)]);
            if nBindingMax == 0
                return;
            end
            
            c = this.chromosome;
            [posStrnds, lens] = find(c.polymerizedRegions);
            smcPosStrands = find(c.complexBoundSites == this.enzymeGlobalIndexs(this.enzymeIndexs_SMC_ADP));
            smcPosStrands = [
                mod(smcPosStrands(:,1)-this.smcSepNt/2-this.smcSepProbCenter/2+this.enzymeDNAFootprints(this.enzymeIndexs_SMC_ADP)/2 -1, c.sequenceLen)+1  2*ceil(smcPosStrands(:,2)/2)-1;
                mod(smcPosStrands(:,1)-this.smcSepNt/2-this.smcSepProbCenter/2+this.enzymeDNAFootprints(this.enzymeIndexs_SMC_ADP)/2 -1, c.sequenceLen)+1  2*ceil(smcPosStrands(:,2)/2)];
            [posStrnds, lens] = c.excludeRegions(posStrnds, lens, smcPosStrands, this.smcSepNt(ones(size(smcPosStrands, 1), 1), 1) + this.smcSepProbCenter);
            if isempty(posStrnds)
                return;
            end
            
            nBound = this.bindProteinToChromosomeStochastically(...
                this.enzymeIndexs_SMC_ADP,...
                nBindingMax, posStrnds, lens, [], [], @this.calcNewRegions);
            if nBound == 0
                return;
            end
            
            %% update molecule counts
            %- covalently bind ATP to SMC
            %- hydrolyze ATP, leaving ADP bound to SMC and releasing PI
            %- bind SMC-ADP to chromosome (already handled by bindProtein method)
            this.enzymes(this.enzymeIndexs_SMC)             = this.enzymes(this.enzymeIndexs_SMC)             - nBound;
            this.enzymes(this.enzymeIndexs_SMC_ADP)         = this.enzymes(this.enzymeIndexs_SMC_ADP)         + nBound;
            
            this.substrates(this.substrateIndexs_atp)       = this.substrates(this.substrateIndexs_atp)       - nBound;
            this.substrates(this.substrateIndexs_water)     = this.substrates(this.substrateIndexs_water)     - nBound;
            this.substrates(this.substrateIndexs_hydrogen)  = this.substrates(this.substrateIndexs_hydrogen)  + nBound;
            this.substrates(this.substrateIndexs_phosphate) = this.substrates(this.substrateIndexs_phosphate) + nBound;
        end
        
        function [rgnPosStrnds, rgnLens, rgnProbs] = calcNewRegions(this, rgnPosStrnds, rgnLens, ~, rgnIdx, offset)
            c = this.chromosome;
            [rgnPosStrnds, rgnLens] = c.excludeRegions(rgnPosStrnds, rgnLens, ...
                [mod(rgnPosStrnds(rgnIdx,1)+offset -this.smcSepNt/2-this.smcSepProbCenter/2+this.enzymeDNAFootprints(this.enzymeIndexs_SMC_ADP)/2 -1, c.sequenceLen)+1 ...
                rgnPosStrnds(rgnIdx,2)], this.smcSepNt + this.smcSepProbCenter);
            rgnProbs = max(0, rgnLens - this.enzymeDNAFootprints(this.enzymeIndexs_SMC_ADP) + 1);
        end
    end
end
