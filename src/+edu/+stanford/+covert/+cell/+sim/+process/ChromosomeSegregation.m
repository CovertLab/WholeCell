%Chromosome Segregation
%
% @wholeCellModelID Process_ChromosomeSegregation
% @name             Chromosome Segregation
% @description
%   Biology
%   ======================
%   Chromosome segregation in Mycoplasma genitalium has not been well described.
%   Chromosome segregation in M. genitalium is believed to result from two
%   factors
%   - entropically favorable segregation throughout DNA replication
%   - four chromosome segregation proteins which act to decatenate the
%     chromosomes following the completion of replication
%   Little is known about the times at which these two factors contribute to
%   chromosome segregation.
%
%   Bacterial chromosome segregation is believed to be entropically favorable
%   and occur during replication (Bloom et al. 2010, Jun et al. 2006). That is,
%   its believed to be entropically unfavorable for two chromosomes to be
%   located near each other, and so as the replication fork moves, the
%   replicated DNA molecules migrate away from each other toward the poles. We
%   assume that by the end of replication the chromosomes have already migrated
%   toward the poles.
%
%   Four proteins are believed to be involved in chromosome decatenation in M.
%   genitalium:
%   - MG_470_MONOMER  nucleotide binding domain protein CobQ/CobB/MinD/ParA
%   - MG_221_OCTAMER  MraZ
%   - MG_387_MONOMER  GTP binding protein Era
%   - MG_384_MONOMER  GTPase Obg
%   The specific function, kinetics, and metabolic costs of these proteins are
%   not known. Furthemore, because M. genitalium contains a reduced complement
%   of segregation proteins, it is difficult to infer their functions from
%   studies of other bacterial species.
%
%   Knowledge Base
%   ======================
%   The knowledge base contains the value of the gtpCost parameter.
%
%   Representation
%   ======================
%   The substrates and enzymes properties represent the counts of available
%   metabolites and segregation enzymes. gtpCost represents the energetic
%   (GTP) cost of chromosome segregation. The segregation status of the
%   chromosomes is represented as a boolean properties (segregated) of the
%   Chromosome class.
%
%   Initialization
%   ======================
%   Chromosome initializes to a state with 1 chromosome (which obviously
%   hasn't yet segregated from the second chromosome which will be later
%   produced).
%
%   Simulation
%   ======================
%   Because little is known about the molecular biology of chromosome
%   segregation, we have chosen to implement this process as a simple boolean
%   rule: chromosome segregation can occur if:
%   - the chromosome is replicated
%   - the chromosome is properly supercoiled
%   - there is at least one free molecule of each segregation proteins, and
%   - there is at least gtpCost free GTP molecules
%
%   Note, although the Glass et al. gene essentiality study suggests that the
%   cobQ/cobB/minD/parA gene is non-essential, we model this gene as essential
%   because we don't know its specific function and how the other segregation
%   proteins compensate in its absence.
%
%   References
%   ======================
%   1. Bloom, K., Joglekar, A. (2010). Towards building a chromosome segregation
%      machine. Nature 463: 446-456.
%   2. Jun, S., Mulder, B. (2006). Entropy-driven special organization of highly
%      confined polymers: Lessons for the bacterial chromosome. PNAS 103:
%      12388-12393.
%   3. Glass, J.I., Assad-Garcia, N., Alperovich, N., Yooseph, S., Lewis, M.R.,
%      Maruf, M., Hutchison III, C.A., Smith, H.O., Venter, J.C. (2006).
%      Essential genes of a minimal bacterium. PNAS 103: 425-430.
%
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/7/2010
classdef  ChromosomeSegregation < edu.stanford.covert.cell.sim.Process & edu.stanford.covert.cell.sim.ChromosomeProcessAspect
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
            'gtpCost'};
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {}; %names of simulation state properties redundant with timecourses in this or other processes or the simulation
    end
    
    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs   = {}; %whole cell model IDs of stimuli
        
        substrateWholeCellModelIDs = {   %whole cell model IDs of substrates
            'GTP';'GDP';'H';'H2O';'PI'};
        substrateIndexs_gtp       = 1; %index of ATP within substrateWholeCellModelIDs
        substrateIndexs_gdp       = 2; %index of ADP within substrateWholeCellModelIDs
        substrateIndexs_hydrogen  = 3; %index of hydrogen within substrateWholeCellModelIDs
        substrateIndexs_water     = 4; %index of water within substrateWholeCellModelIDs
        substrateIndexs_phosphate = 5; %index of inorganic phosphate within substrateWholeCellModelIDs
        
        enzymeWholeCellModelIDs = {      %whole cell model IDs of enzymes
            'MG_470_MONOMER'       %CobQ/CobB/MinD/ParA nucleotide binding domain
            'MG_221_OCTAMER'       %mraZ protein
            'MG_387_MONOMER'       %GTP-binding protein Era
            'MG_384_MONOMER'       %GTPase1 Obg
            'MG_203_204_TETRAMER'  %DNA topoisomerase IV
            };
        enzymeIndexs_cobQ   = 1;
        enzymeIndexs_mraZ   = 2;
        enzymeIndexs_obg    = 3;
        enzymeIndexs_era    = 4;
        enzymeIndexs_topoIV = 5;
    end
    
    %fixed biological constants
    properties
        gtpCost %number of GTP required for chromosome segregation (1); true value unknown
    end
    
    %constructor
    methods
        function this = ChromosomeSegregation(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.Process(wholeCellModelID, name);
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
    end
    
    %model
    methods
        %Calculate
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, ~, ~)
            %% initialize
            bmProd = zeros(size(this.substrateWholeCellModelIDs));
            byProd = zeros(size(this.substrateWholeCellModelIDs));
            minEnzExp = zeros(size(this.enzymeWholeCellModelIDs));
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
            
            %% substrate and byproducts: GTP required to decatenate chromosomes, and its hydrolysis (GTP + H2O ==> GDP + PI + H)
            bmProd(this.substrateIndexs_gtp)       = this.gtpCost;
            bmProd(this.substrateIndexs_water)     = this.gtpCost;
            byProd(this.substrateIndexs_gdp)       = this.gtpCost;
            byProd(this.substrateIndexs_phosphate) = this.gtpCost;
            byProd(this.substrateIndexs_hydrogen)  = this.gtpCost;
            
            %% Segregation requires at least 1 copy of every enzyme
            minEnzExp(:) = 2;
        end
        
        %initialization: chromosome state initializes chromosome to state with 1
        %chromosome (which trivially is not segregated)
        function initializeState(~)
        end
        
        %resource requirements
        function result = calcResourceRequirements_Current(this)
            c = this.chromosome;
            result = zeros(size(this.substrates));
            if ...
                    ~c.segregated && ...
                    collapse(c.polymerizedRegions) == c.nCompartments * c.sequenceLen && ...
                    all(this.enzymes)
                
                result(this.substrateIndexs_gtp) = this.gtpCost;
                result(this.substrateIndexs_water) = result(this.substrateIndexs_gtp);
            end
        end
        
        %simulation
        function evolveState(this)
            c = this.chromosome;
            
            if ...
                    ~c.segregated && ...
                    collapse(c.polymerizedRegions) == c.nCompartments * c.sequenceLen && ...
                    collapse(c.supercoiled) == c.nCompartments && ...
                    all(this.enzymes) && ...
                    this.substrates(this.substrateIndexs_gtp) >= this.gtpCost && ...
                    this.substrates(this.substrateIndexs_water) >= this.gtpCost
                
                c.segregated = true;
                
                this.substrates(this.substrateIndexs_gtp)       = this.substrates(this.substrateIndexs_gtp)       - this.gtpCost;
                this.substrates(this.substrateIndexs_water)     = this.substrates(this.substrateIndexs_water)     - this.gtpCost;
                this.substrates(this.substrateIndexs_gdp)       = this.substrates(this.substrateIndexs_gdp)       + this.gtpCost;
                this.substrates(this.substrateIndexs_phosphate) = this.substrates(this.substrateIndexs_phosphate) + this.gtpCost;
                this.substrates(this.substrateIndexs_hydrogen)  = this.substrates(this.substrateIndexs_hydrogen)  + this.gtpCost;
            end
        end
    end
end
