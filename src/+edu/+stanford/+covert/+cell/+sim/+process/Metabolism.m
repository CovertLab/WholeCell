%Metabolism
%
% @wholeCellModelID Process_Metabolism
% @name             Metabolism
% @description
%   Biology
%   ===============
%   To grow and replicate cells must uptake or produce their building blocks and
%   intermediate energy stores, particularly nucleic and amino acids and lipids,
%   as well as secrete metabolites such as modified nucleic acids which it
%   cannot catabolize and extract usable material or energy from. Furthermore,
%   the import/synthesis and export/breakdown of metabolites must be carefully
%   matched to the metabolic demands of the cell to ensure that cellular
%   processes aren't limited by too few metabolites or poisoned by too many
%   metabolites.
%
%   This process simulates
%   - the uptake of nutrients from the external environment
%   - the processing of imported nutrients to intermediate energy forms (ATP,
%     GTP) and macromolecule building blocks (dNTPs, NTPs, and amino acids)
%   - the assembly of lipids, their insertion into the membrane, and their
%     maturation within the membrane
%   - the catabolism of byproducts of macromolecule assembly and degradation
%   - the export to the external environment of uncatabolizable chemicals such
%     as modified nucleobases
%
%   Knowledge Base
%   ===============
%   M. genitalium metabolism was reconstructed from a variety of sources,
%   including flux-balance analysis metabolic models of other bacterial species
%   and the reaction kinetics database SABIO-RK, and was organized into 641
%   reactions in the knowledge base. These reactions are loaded into this process
%   by the initializeConstants method.
%
%      Object       No.
%      ==========   ===
%      substrates   580
%      enzymes      100
%      reactions    641
%        chemical   433
%        transport  208
%
%   Representation
%   ===============
%   The properties substrates and enzymes represent the counts of metabolites
%   and metabolic enzymes.
%
%   fbaReactionStoichiometryMatrix represents the stoichiometry and compartments
%   of metabolites and biomass in each of the 641 chemical/transport reactions,
%   exchange pseudoreactions, and biomass production pseudoreaction.
%   fbaReactionCatalysisMatrix represents the enzyme which catalyzes each
%   reaction. fbaEnzymeBounds represents the foward and backward kcat of the
%   catalyzing enzyme of each reaction. fbaReactionBounds represents the maximal
%   import and export rates of each  metabolite. fbaObjective indicates which
%   reaction represents the biomass production pseudoreaction. fbaRightHandSide
%   is a vector of zeros representing the change in concentration over time of
%   each metabolite and biomass. metabolismProduction is redundant with the
%   biomass production reaction in fbaReactionStoichiometryMatrix.
%   metabolismProduction is calculated by summing the metabolic demands of all
%   the other processes over the entire cell cycle. The table below lists the
%   units of several properties of this process.
%
%      Property                       Units
%      ===========================    ==============================
%      fbaEnzymeBounds                molecules/enzyme/s
%      fbaReactionBounds              molecules/(gram dry biomass)/s
%      metabolites                    molecules
%      enzymes                        molecules
%      stepSizeSec                    s
%      lowerBounds                    reactions/s
%      upperBounds                    reactions/s
%      growth                         cell/s
%      biomassComposition             molecules/cell
%      metabolismProduction           molecules/cell
%      chamberVolume                  L
%      setValues                      molecules/chamber
%      growthAssociatedMaintanence
%
%   Initialization
%   ===============
%   The simulation is initialized with 1 cell weight of macromolecules and
%   water, and few free metabolites by the simulation class' initializeState
%   method. In addition this process is initialized to a positive growth rate of
%   approximately log(2)/cellCycleLength computing using the flux-balance
%   analysis model implemented in evolveState.
%
%   Simulation
%   ===============
%   Transport and metabolism are modeled using the constraint-based method
%   flux-balance analysis (FBA). Briefly, FBA assumes that bacteria have evolved
%   to maximize growth, and poses transport and metabolism as the optimization
%   of cellular building block and energy production subject to available
%   nutrients, and allowed chemical reactions. FBA models are typically
%   constrained by experimentally measured transport and diffusion rates. In
%   addition to these constraints, we constrain constrain our FBA process by our
%   model's predicted enzyme abundances, and by experimentally measured kinetic
%   parameters. These additional constraints yield a more accurate model of
%   transport and metabolism. Furthermore, we use our other processes to obtain a
%   more extensive and accurate metabolic objective than used in earlier FBA
%   models. The optimization problem specific by FBA is posed as a linear
%   optimization problem, and solving using one of several publically available
%   linear programming packages.
%
%   Algorithm
%   ++++++++++++++++
%   1. Compute reaction bounds based on
%     - enzyme kinetics,
%     - enzyme availability
%     - maximal metabolite exchange rates
%     - external metabolite availability
%     - protein availability
%   2. Computes optimal reaction fluxes which maximize biomass production
%   3. Computes integer-valued production of biomass components, and export of
%      byproducts
%   4. Updates amounts of biomass components and byproducts
%
%   References
%   ===============
%   1. Orth JD, Thiele I, Palsson BO (2010). What is flux balance analysis?
%      Nat Biotechnol. 28(3):245-8 [PUB_0687].
%   2. Thiele I, Palsson BO (2010). A protocol for generating a
%      high-quality genome-scale metabolic reconstruction. Nat Protoc.
%      5(1): 93-121. [PUB_0686].
%   3. Covert MW, Xiao N, Chen TJ, Karr JR (2008). Integrating metabolic,
%      transcriptional regulatory and signal transduction models in
%      Escherichia coli. Bioinformatics. 24(18):2044-50. [PUB_0684].
%   4. Covert MW, Knight EM, Reed JL, Herrgard MJ, Palsson BO (2004).
%      Integrating high-throughput and computational data elucidates
%      bacterial networks. Nature. 429 (6987): 92-6. [PUB_0618].
%   5. Covert MW, Palsson BO (2003). Constraints-based models: regulation
%      of gene expression reduces the steady-state solution space. J Theor
%      Biol. 221(3): 309-25. [PUB_0685].
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Markus Covert, mcovert@stanford.edu
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 3/22/2011
classdef  Metabolism < edu.stanford.covert.cell.sim.ReactionProcess
    %property annotations
    properties (Constant)
        optionNames__              = {   %names of option properties
            'linearProgrammingOptions';
            'macromoleculeStateInitialization';
            'tolerance'
            'realmax'
            };
        fixedConstantNames__       = {   %names of fixed constant properties
            'growthAssociatedMaintenance';
            'nonGrowthAssociatedMaintenance';
            'cellCycleLength';
            'exchangeRateUpperBound_carbon';
            'exchangeRateUpperBound_noncarbon';
            'fbaObjective';
            'fbaRightHandSide';
            'fbaReactionStoichiometryMatrix';
            'fbaReactionCatalysisMatrix';
            'fbaReactionBounds';
            'fbaEnzymeBounds';
            'metabolismProduction';
            'monomerCompartments';
            };
        fittedConstantNames__      = {   %names of fitted constant properties
            'macromoleculeStateInitializationGrowthFactor'
            'unaccountedEnergyConsumption';
            };
        localStateNames__          = {}; %names of simulation state properties redundant with timecourses in this or other processes or the simulation
    end
    
    %options
    properties
        %linear programming
        linearProgrammingOptions = struct(...
            'solver', 'glpk',...
            'solverOptions', struct(...
                'glpk', struct('lpsolver', 1, 'presol', 1, 'scale', 1, 'msglev', 0, 'tolbnd', 10e-7),...
                'linprog', struct('Display','off'),...
                'lp_solve', struct('verbose', 0, 'scaling', 3 + 64 + 128, 'presolve', 0), ...
                'qsopt', struct()));
        
        %initial macromolecule state distribution
        macromoleculeStateInitialization = 'multinomial';
        
        %fitting tolerance
        tolerance = 1e-3;
        
        %realmax
        % Set to several orders of magnitude higher than non-cycle highest
        % fluxs. In tests gplk will tolerates as high as 10^9. Higher than
        % 10^9 glpk and other solvers have trouble finding feasible
        % solutions.
        realmax = 1e6;
    end
    
    %IDs, names, and local indices
    properties
        substrateIndexs_amp                   %index of AMP within reactionStoichiometryMatrix
        substrateIndexs_adp                   %index of ADP within reactionStoichiometryMatrix
        substrateIndexs_atp                   %index of ATP within reactionStoichiometryMatrix
        substrateIndexs_diphosphate           %index of PPi within reactionStoichiometryMatrix
        substrateIndexs_phosphate             %index of Pi  within reactionStoichiometryMatrix
        substrateIndexs_water                 %index of H2O within reactionStoichiometryMatrix
        substrateIndexs_hydrogen              %index of H   within reactionStoichiometryMatrix
        substrateIndexs_atpHydrolysis         %index of ATP, H2O, ADP, PI, H within reactionStoichiometryMatrix
        substrateIndexs_energy                %index of NxPs, PPI, PI, H within reactionStoichiometryMatrix
        
        substrateIndexs_limitableProteins     %indices of limitable proteins within reactionStoichiometryMatrix
        substrateIndexs_fba                   %indices of FBA substrates within reactionStoichiometryMatrix
        substrateIndexs_externalExchangedMetabolites  %indices of exchanged metabolites within reactionStoichiometryMatrix
        substrateIndexs_internalExchangedMetabolites  %indices of exchanged metabolites within reactionStoichiometryMatrix
        substrateIndexs_internalExchangedLimitedMetabolites %indices of limited, exchanged metabolites within reactionStoichiometryMatrix
        
        fbaSubstrateIndexs_substrates         %indices of substrates (metabolites, RNAs, proteins) within fbaReactionStoichiometryMatrix
        fbaSubstrateIndexs_metaboliteInternalExchangeConstraints %indices of substrates with are internally exchanged
        fbaSubstrateIndexs_biomass            %index of biomass within fbaReactionStoichiometryMatrix
        
        reactionIndexs_chemical               %indices of chemical reactions within reactionStoichiometryMatrix
        reactionIndexs_transport              %indices of transport reactions within reactionStoichiometryMatrix
        reactionIndexs_fba                    %indices of FBA reactions within reactionStoichiometryMatrix
        
        fbaReactionIndexs_metabolicConversion %indices of metabolic reactions (chemical, transport) within fbaReactionStoichiometryMatrix
        fbaReactionIndexs_metaboliteExternalExchange %indices of exchange reactions within fbaReactionStoichiometryMatrix
        fbaReactionIndexs_metaboliteInternalExchange %indices of exchange reactions within fbaReactionStoichiometryMatrix
        fbaReactionIndexs_metaboliteInternalLimitedExchange %indices of exchange reactions within fbaReactionStoichiometryMatrix
        fbaReactionIndexs_metaboliteInternalUnlimitedExchange %indices of exchange reactions within fbaReactionStoichiometryMatrix
        fbaReactionIndexs_biomassProduction   %index of biomass production reaction within fbaReactionStoichiometryMatrix
        fbaReactionIndexs_biomassExchange     %index of biomass efflux reaction within fbaReactionStoichiometryMatrix
        
        compartmentIndexs_cytosol             %index of cytosol compartment within reactionStoichiometryMatrix
        compartmentIndexs_extracellular       %index of extracellular compartment within reactionStoichiometryMatrix
        compartmentIndexs_membrane            %index of membrane compartment within reactionStoichiometryMatrix
    end
    
    %fixed biological constants
    properties
        macromoleculeStateInitializationGrowthFactor %set by running fitting once and seeing what mean growth rate is predicted; see fitEnzymes
        
        growthAssociatedMaintenance         %mmol ATP/gDCW [59.810; PUB_0558]
        nonGrowthAssociatedMaintenance      %mmol ATP/gDCW [8.390; PUB_0558]
        cellCycleLength                     %cell cycle length (s)
        exchangeRateUpperBound_carbon       %exchange rate upper bound -- carbon-containing metabolites (mmol/gDCW/h) [12]
        exchangeRateUpperBound_noncarbon    %exchange rate upper bound -- non-carbon-containing metabolites (mmol/gDCW/h) [20]
        substrateExternalExchangeBounds     %substrate exchange bounds                                                          [substrates X 2]
        substrateExchangeBounds             %substrate exchange bounds                                                          [substrates X 2]
        fbaObjective                        %FBA objective function corresponding to maximize growth                            [reactions;growth;biomass efflux] X 1
        fbaRightHandSide                    %right-hand side of equation Sv=0                                                   [substrates;biomass] X 1
        fbaReactionStoichiometryMatrix      %FBA stoichiometry matrix                                                           [substrates;biomass] X [reactions;growth;biomass efflux] X compartments
        fbaReactionCatalysisMatrix          %Enzymes catalyzing each FBA reaction                                               [reactions;growth;biomass efflux] X enzymes
        fbaReactionBounds                   %Exchange bounds associated with each FBA reaction (molecules/gram dry biomass/s)   [reactions;growth;biomass efflux] X 2
        fbaEnzymeBounds                     %Kinetics bounds associated with each FBA reaction (molecules/enzyme/s)             [reactions;growth;biomass efflux] X 2
        metabolismProduction                %biomass composition + energy + byproducts (molecules/cell)                         substrates X compartments
        metabolismNewProduction             %metabolism output represented by biomass pseudoreaction
        metabolismRecyclingProduction       %metabolism output represented by internal exchange reactions
        unaccountedEnergyConsumption        %unaccounted energy consumption (ATP molecules / cell life cycle)                   ATP X 1
        monomerCompartments                 %Indexs of compartments of protein monomers
        proteinLimitableProteinComposition  %protein complex composition                                                        [monomers; complexs] X [monomers; limitable complexs]
    end
    
    %global state (referenced locally for convenience)
    properties
        mass                      %cell mass
        metabolicReaction         %metabolic flux
    end
    
    %constructor
    methods
        function this = Metabolism(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.ReactionProcess(wholeCellModelID, name);
        end
    end
    
    %communication between process/simulation
    methods
        %set references to state objects
        function storeObjectReferences(this, simulation)
            this.storeObjectReferences@edu.stanford.covert.cell.sim.ReactionProcess(simulation);
            
            this.mass = simulation.state('Mass');
            this.metabolicReaction = simulation.state('MetabolicReaction');
            this.states = [this.states; {this.mass; this.metabolicReaction}];
        end
        
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)
            import edu.stanford.covert.util.ConstantUtil;
            
            %% options
            options = struct('retainSubstrateCompartments', true);
            
            %% call super class method
            this.initializeConstants@edu.stanford.covert.cell.sim.ReactionProcess(...
                knowledgeBase, simulation, options);
            
            %% cell cycle length
            this.cellCycleLength = simulation.state('Time').cellCycleLength;
            
            %% compartments
            this.compartmentIndexs_cytosol       = 1;
            this.compartmentIndexs_extracellular = 2;
            this.compartmentIndexs_membrane      = 3;
            
            compartmentIndexs = [
                this.compartment.cytosolIndexs;
                this.compartment.extracellularIndexs;
                this.compartment.membraneIndexs];
            
            this.substrateStimulusGlobalIndexs        = this.substrateStimulusGlobalIndexs(  :, compartmentIndexs);
            this.substrateMetaboliteGlobalIndexs      = this.substrateMetaboliteGlobalIndexs(:, compartmentIndexs);
            this.substrateRNAGlobalIndexs             = this.substrateRNAGlobalIndexs(       :, compartmentIndexs);
            this.substrateMonomerGlobalIndexs         = this.substrateMonomerGlobalIndexs(   :, compartmentIndexs);
            this.substrateComplexGlobalIndexs         = this.substrateComplexGlobalIndexs(   :, compartmentIndexs);
            
            this.substrateStimulusCompartmentIndexs   = this.substrateStimulusCompartmentIndexs(  :, compartmentIndexs);
            this.substrateMetaboliteCompartmentIndexs = this.substrateMetaboliteCompartmentIndexs(:, compartmentIndexs);
            this.substrateRNACompartmentIndexs        = this.substrateRNACompartmentIndexs(       :, compartmentIndexs);
            this.substrateMonomerCompartmentIndexs    = this.substrateMonomerCompartmentIndexs(   :, compartmentIndexs);
            this.substrateComplexCompartmentIndexs    = this.substrateComplexCompartmentIndexs(   :, compartmentIndexs);
            
            this.substrateStimulusGlobalCompartmentIndexs = sub2ind(...
                [numel(this.stimulus.wholeCellModelIDs) this.compartment.count],...
                this.substrateStimulusGlobalIndexs,...
                this.substrateStimulusCompartmentIndexs);
            this.substrateMetaboliteGlobalCompartmentIndexs = sub2ind(...
                [numel(this.metabolite.wholeCellModelIDs) this.compartment.count],...
                this.substrateMetaboliteGlobalIndexs,...
                this.substrateMetaboliteCompartmentIndexs);
            this.substrateRNAGlobalCompartmentIndexs = sub2ind(...
                [numel(this.rna.wholeCellModelIDs) this.compartment.count],...
                this.rna.matureIndexs(this.substrateRNAGlobalIndexs),...
                this.substrateRNACompartmentIndexs);
            this.substrateMonomerGlobalCompartmentIndexs = sub2ind(...
                [numel(this.monomer.wholeCellModelIDs) this.compartment.count],...
                this.monomer.matureIndexs(this.substrateMonomerGlobalIndexs),...
                this.substrateMonomerCompartmentIndexs);
            this.substrateComplexGlobalCompartmentIndexs = sub2ind(...
                [numel(this.complex.wholeCellModelIDs) this.compartment.count],...
                this.complex.matureIndexs(this.substrateComplexGlobalIndexs),...
                this.substrateComplexCompartmentIndexs);
            
            this.reactionStoichiometryMatrix = this.reactionStoichiometryMatrix(:, :, compartmentIndexs);
            this.reactionCoenzymeMatrix = this.reactionCoenzymeMatrix(:, :, compartmentIndexs);
            
            %% metabolites
            this.substrateIndexs_amp           = this.substrateIndexs({'AMP'});
            this.substrateIndexs_adp           = this.substrateIndexs({'ADP'});
            this.substrateIndexs_atp           = this.substrateIndexs({'ATP'});
            this.substrateIndexs_diphosphate   = this.substrateIndexs({'PPI'});
            this.substrateIndexs_phosphate     = this.substrateIndexs({'PI'});
            this.substrateIndexs_water         = this.substrateIndexs({'H2O'});
            this.substrateIndexs_hydrogen      = this.substrateIndexs({'H'});
            this.substrateIndexs_atpHydrolysis = this.substrateIndexs({'ATP'; 'H2O'; 'ADP'; 'PI'; 'H'});
            this.substrateIndexs_energy        = sort(this.substrateIndexs({'AMP', 'CMP', 'GMP', 'UMP', 'ADP', 'CDP', 'GDP', 'UDP', 'ATP', 'CTP', 'GTP', 'UTP', 'PI', 'PPI', 'H'}));
            
            %% reactions
            this.reactionIndexs_chemical  = find(strcmp(this.reactionTypes, 'chemical'));
            this.reactionIndexs_transport = find(strcmp(this.reactionTypes, 'transport'));
            
            %% FBA
            %comlexes
            complexFormationProcesses = [knowledgeBase.proteinComplexs(this.substrateComplexGlobalIndexs(:, 1)).complexFormationProcess];
            this.substrateIndexs_limitableProteins = [...
                this.substrateMonomerLocalIndexs
                this.substrateComplexLocalIndexs(strcmp({complexFormationProcesses.wholeCellModelID}, 'Process_MacromolecularComplexation'), :)];
            
            substratesByLimitableSubstrates = zeros(numel(this.substrateWholeCellModelIDs), numel(this.substrateWholeCellModelIDs));
            substratesByLimitableSubstrates(this.substrateComplexLocalIndexs, this.substrateMonomerLocalIndexs) = ...
                any(knowledgeBase.proteinComplexMonomerComposition(this.gene.mRNAIndexs(this.substrateMonomerGlobalIndexs(:, 1)), this.substrateComplexGlobalIndexs(:, 1), :), 3)';
            substratesByLimitableSubstrates(this.substrateComplexLocalIndexs, this.substrateComplexLocalIndexs) =  ...
                any(knowledgeBase.proteinComplexComplexComposition(this.substrateComplexGlobalIndexs(:, 1), this.substrateComplexGlobalIndexs(:, 1), :), 3)';
            substratesByLimitableSubstrates(:, setdiff(1:end, this.substrateIndexs_limitableProteins)) = 0;
            substratesByLimitableSubstrates(this.substrateIndexs_limitableProteins, :) = 0;
            substratesByLimitableSubstrates(sub2ind(...
                size(substratesByLimitableSubstrates), ....
                this.substrateIndexs_limitableProteins, ...
                this.substrateIndexs_limitableProteins)) = 1;
            
            this.proteinLimitableProteinComposition = substratesByLimitableSubstrates([this.substrateMonomerLocalIndexs; this.substrateComplexLocalIndexs], this.substrateIndexs_limitableProteins);
            
            %exchange bounds
            mets = knowledgeBase.metabolites(this.substrateMetaboliteGlobalIndexs(:, 1));
            this.substrateExternalExchangeBounds = zeros(numel(this.substrateWholeCellModelIDs), 2);
            for i = 1:numel(mets)
                met = mets(i);
                if isfield(met.empiricalFormula, 'C') && met.empiricalFormula.C > 0
                    this.substrateExternalExchangeBounds(this.substrateMetaboliteLocalIndexs(i), 1) = max(-this.exchangeRateUpperBound_carbon, met.exchangeLowerBound);
                    this.substrateExternalExchangeBounds(this.substrateMetaboliteLocalIndexs(i), 2) = min( this.exchangeRateUpperBound_carbon, met.exchangeUpperBound);
                else
                    this.substrateExternalExchangeBounds(this.substrateMetaboliteLocalIndexs(i), 1) = max(-this.exchangeRateUpperBound_noncarbon, met.exchangeLowerBound);
                    this.substrateExternalExchangeBounds(this.substrateMetaboliteLocalIndexs(i), 2) = min( this.exchangeRateUpperBound_noncarbon, met.exchangeUpperBound);
                end
            end
            this.substrateExternalExchangeBounds = this.substrateExternalExchangeBounds / ...
                ConstantUtil.secondsPerHour * (ConstantUtil.nAvogadro / 1000);
            
            %% Monomers
            this.monomerCompartments = knowledgeBase.proteinMonomerCompartments;
        end
        
        %Formulate FBA model for current biomass composition stored in
        %metabolite state and this module's reactions
        %- objective
        %- stoichiometry matrix
        %- bounds
        %
        %To decrease linear programming run-time and increase linear
        %programming accuracy this method automatically reduces the size of
        %the linear programming problem by removing
        %- reaction and substrates disconnected from biomass
        %- leaf metabolites (and reactions involving them) which cannot be
        %  balanced with non-zero fluxes
        %- metabolites (and reactions involving them) which can't move in
        %  both the forward and backward directions and thus cannot be
        %  balanced
        %- other reactions (and substrates unique to them) whose flux is
        %  zero in all solutions of the matrix problem Sv=0. That is
        %  reactions which have only zero components in the null space of
        %  S.
        function [subCmpRemoved, rxnRemoved] = formulateFBA(this, metabolismProduction, unaccountedEnergyConsumption, simplifyNetwork)
            import edu.stanford.covert.util.ComputationUtil;
            import edu.stanford.covert.util.ConstantUtil;
            import edu.stanford.covert.util.findNonInteractingRowsAndColumns;
            
            %% metabolism output
            if nargin >= 2 && ~isempty(metabolismProduction)
                metabolismProduction(abs(metabolismProduction) < 1e-4) = 0; %zero out metabolites that almost zero up to rounding error
                this.metabolismProduction = zeros(numel(this.substrateWholeCellModelIDs), size(this.substrateMetaboliteGlobalIndexs, 2));
                this.metabolismProduction(this.substrateMetaboliteLocalIndexs, :) = metabolismProduction(sub2ind(...
                    size(metabolismProduction), ...
                    this.substrateMetaboliteGlobalIndexs, ...
                    this.substrateMetaboliteCompartmentIndexs));
            end
            if nargin >= 3 && ~isempty(unaccountedEnergyConsumption)
                this.unaccountedEnergyConsumption = unaccountedEnergyConsumption;
            end
            assert(~any(sum(this.metabolismProduction ~= 0, 2) > 1));
            
            this.metabolismNewProduction = max(0, this.metabolismProduction);
            this.metabolismRecyclingProduction = min(0, this.metabolismProduction);
            
            this.metabolismRecyclingProduction(this.substrateIndexs({'ATP'; 'CTP'; 'GTP'; 'UTP'}), :) = ...
                + this.metabolismRecyclingProduction(this.substrateIndexs({'ATP'; 'CTP'; 'GTP'; 'UTP'}), :) ...
                - this.metabolismRecyclingProduction(this.substrateIndexs({'ADP'; 'CDP'; 'GDP'; 'UDP'}), :) ...
                - this.metabolismRecyclingProduction(this.substrateIndexs({'AMP'; 'CMP'; 'GMP'; 'UMP'}), :);
            this.metabolismNewProduction(this.substrateIndexs({'ATP'; 'CTP'; 'GTP'; 'UTP'}), :) = ...
                + this.metabolismNewProduction(this.substrateIndexs({'ATP'; 'CTP'; 'GTP'; 'UTP'}), :) ...
                + this.metabolismRecyclingProduction(this.substrateIndexs({'ADP'; 'CDP'; 'GDP'; 'UDP'}), :) ...
                + this.metabolismRecyclingProduction(this.substrateIndexs({'AMP'; 'CMP'; 'GMP'; 'UMP'}), :);
            
            this.metabolismRecyclingProduction(this.substrateIndexs({'H2O'}), :) = ...
                + this.metabolismRecyclingProduction(this.substrateIndexs({'H2O'}), :) ...
                - sum(this.metabolismRecyclingProduction(this.substrateIndexs({'ADP'; 'CDP'; 'GDP'; 'UDP'}), :)) ...
                - sum(this.metabolismRecyclingProduction(this.substrateIndexs({'AMP'; 'CMP'; 'GMP'; 'UMP'}), :));
            this.metabolismNewProduction(this.substrateIndexs({'H2O'}), :) = ...
                + this.metabolismNewProduction(this.substrateIndexs({'H2O'}), :) ...
                + sum(this.metabolismRecyclingProduction(this.substrateIndexs({'ADP'; 'CDP'; 'GDP'; 'UDP'}), :)) ...
                + sum(this.metabolismRecyclingProduction(this.substrateIndexs({'AMP'; 'CMP'; 'GMP'; 'UMP'}), :));
            
            this.metabolismRecyclingProduction(this.substrateIndexs({'ATP'; 'ADP'; 'PI'; 'H2O'; 'H'}), this.compartmentIndexs_cytosol) = ...
                + this.metabolismRecyclingProduction(this.substrateIndexs({'ATP'; 'ADP'; 'PI'; 'H2O'; 'H'}), this.compartmentIndexs_cytosol) ...
                + this.unaccountedEnergyConsumption * [-1; 1; 1; -1; 1];
            this.metabolismNewProduction(this.substrateIndexs({'ATP'; 'ADP'; 'PI'; 'H2O'; 'H'}), this.compartmentIndexs_cytosol) = ...
                + this.metabolismNewProduction(this.substrateIndexs({'ATP'; 'ADP'; 'PI'; 'H2O'; 'H'}), this.compartmentIndexs_cytosol) ...
                - this.unaccountedEnergyConsumption * [-1; 1; 1; -1; 1];
            
            this.metabolismNewProduction(this.substrateIndexs({'MET'; 'H'; 'H2O'}), this.compartmentIndexs_cytosol) = ...
                + this.metabolismNewProduction(this.substrateIndexs({'MET'; 'H'; 'H2O'}), this.compartmentIndexs_cytosol) ...
                + [1; 1; -1] .* this.metabolismRecyclingProduction(this.substrateIndexs({'FMET'}), this.compartmentIndexs_cytosol);
            this.metabolismRecyclingProduction(this.substrateIndexs({'MET'; 'H'; 'H2O'}), this.compartmentIndexs_cytosol) = ...
                + this.metabolismRecyclingProduction(this.substrateIndexs({'MET'; 'H'; 'H2O'}), this.compartmentIndexs_cytosol) ...
                - [1; 1; -1] .* this.metabolismRecyclingProduction(this.substrateIndexs({'FMET'}), this.compartmentIndexs_cytosol);
            
            tmp = ...
                + this.metabolismNewProduction(this.substrateIndexs({'FTHF10'; 'FOR'; 'THF'}), this.compartmentIndexs_cytosol) ...
                + this.metabolismRecyclingProduction(this.substrateIndexs({'FTHF10'; 'FOR'; 'THF'}), this.compartmentIndexs_cytosol);
            this.metabolismNewProduction(this.substrateIndexs({'FTHF10'; 'FOR'; 'THF'}), this.compartmentIndexs_cytosol) = max(0, ...
                + tmp ...
                - [1; -1; -1] .* max(0, min(tmp ./ [1; -1; -1])));
            this.metabolismRecyclingProduction(this.substrateIndexs({'FTHF10'; 'FOR'; 'THF'}), this.compartmentIndexs_cytosol) = ...
                + tmp ...
                - this.metabolismNewProduction(this.substrateIndexs({'FTHF10'; 'FOR'; 'THF'}), this.compartmentIndexs_cytosol);
            
            nInternalExchangeReactions = nnz(this.metabolismRecyclingProduction);
            nInternalExchangeReactionConstraints = 7;
            
            %% numbers and indices of components
            nCompartments = size(this.reactionStoichiometryMatrix, 3);
            nSubstrates   = size(this.reactionStoichiometryMatrix, 1) + nInternalExchangeReactionConstraints + 1;
            nEnzymes      = length(this.enzymeWholeCellModelIDs);
            
            subIdxs_substrates = (1:size(this.reactionStoichiometryMatrix, 1))';
            subIdxs_metaboliteInternalExchangeConstraints = subIdxs_substrates(end) + (1:nInternalExchangeReactionConstraints)';
            subIdxs_biomass = nSubstrates;
            
            subCmpIdxs_metaboliteInternalExchangeConstraints = sub2ind([nSubstrates nCompartments], ...
                subIdxs_metaboliteInternalExchangeConstraints, repmat(this.compartmentIndexs_cytosol, size(subIdxs_metaboliteInternalExchangeConstraints)));
            subCmpIdxs_biomass = sub2ind([nSubstrates nCompartments], subIdxs_biomass, this.compartmentIndexs_cytosol);
            
            nReactions = ...
                + size(this.reactionStoichiometryMatrix, 2) ...
                + length(subIdxs_substrates) ...
                + nInternalExchangeReactions ...
                + 2;
            
            rxnIdxs_metabolicConversion = (1:size(this.reactionStoichiometryMatrix, 2))';
            rxnIdxs_metaboliteExternalExchange = rxnIdxs_metabolicConversion(end) + (1:length(subIdxs_substrates))';
            rxnIdxs_metaboliteInternalExchange = rxnIdxs_metaboliteExternalExchange(end) + (1:nInternalExchangeReactions)';
            rxnIdxs_biomassProduction = nReactions - 1;
            rxnIdxs_biomassExchange = nReactions;
            
            %% reactions
            %initialize
            fbaSMat = zeros(nSubstrates, nReactions, nCompartments);
            
            %chemical and transport reactions
            fbaSMat(subIdxs_substrates, rxnIdxs_metabolicConversion, :) = ...
                this.reactionStoichiometryMatrix;
            
            %metabolite external exchange
            fbaSMat(sub2ind(...
                [nSubstrates, nReactions, nCompartments],...
                subIdxs_substrates,...
                rxnIdxs_metaboliteExternalExchange,...
                repmat(this.compartmentIndexs_extracellular, length(subIdxs_substrates), 1))) = ...
                1;
            
            %metabolite internal exchange
            [i, j] = find(this.metabolismRecyclingProduction);
            fbaSMat(sub2ind(...
                [nSubstrates, nReactions, nCompartments],...
                i, rxnIdxs_metaboliteInternalExchange, j)) = -1;
            rxnIdxs_metaboliteInternalLimitedExchange = ...
                rxnIdxs_metaboliteInternalExchange(...
                this.metabolismRecyclingProduction(this.metabolismRecyclingProduction ~= 0) < 0);
            rxnIdxs_metaboliteInternalUnlimitedExchange = ...
                rxnIdxs_metaboliteInternalExchange(...
                this.metabolismRecyclingProduction(this.metabolismRecyclingProduction ~= 0) > 0);
            
            [~, ntpIdxs] = ismember(this.substrateIndexs({'ATP'; 'CTP'; 'GTP'; 'UTP'}), i);
            [~, ndpIdxs] = ismember(this.substrateIndexs({'ADP'; 'CDP'; 'GDP'; 'UDP'}), i);
            [~, nmpIdxs] = ismember(this.substrateIndexs({'AMP'; 'CMP'; 'GMP'; 'UMP'}), i);
            [~, fthf10Idxs] = ismember(this.substrateIndexs({'FTHF10'; 'THF'; 'FOR'}), i);
            [~, metIdxs] = ismember(this.substrateIndexs({'MET'; 'FMET'}), i);
            tmp = zeros(nInternalExchangeReactionConstraints, nInternalExchangeReactions);
            idxs = [ntpIdxs(1) ndpIdxs(1) nmpIdxs(1)]; tmp(1, idxs(idxs ~= 0)) = 1;
            idxs = [ntpIdxs(2) ndpIdxs(2) nmpIdxs(2)]; tmp(2, idxs(idxs ~= 0)) = 1;
            idxs = [ntpIdxs(3) ndpIdxs(3) nmpIdxs(3)]; tmp(3, idxs(idxs ~= 0)) = 1;
            idxs = [ntpIdxs(4) ndpIdxs(4) nmpIdxs(4)]; tmp(4, idxs(idxs ~= 0)) = 1;
            idxs = fthf10Idxs(1:2); tmp(5, idxs(idxs ~= 0)) = 1;
            idxs = [fthf10Idxs([1 3]); metIdxs(2)]; tmp(6, idxs(idxs ~= 0)) = 1;
            idxs = metIdxs; tmp(7, idxs(idxs ~= 0)) = 1;
            fbaSMat(subIdxs_metaboliteInternalExchangeConstraints, rxnIdxs_metaboliteInternalExchange, this.compartmentIndexs_cytosol) = ...
                + fbaSMat(subIdxs_metaboliteInternalExchangeConstraints, rxnIdxs_metaboliteInternalExchange, this.compartmentIndexs_cytosol) ...
                + tmp;
            
            %biomass production
            fbaSMat(subIdxs_substrates, rxnIdxs_biomassProduction, :) = ...
                -permute(this.metabolismNewProduction, [1 3 2]);
            fbaSMat(...
                subIdxs_biomass, ...
                rxnIdxs_biomassProduction, ...
                this.compartmentIndexs_cytosol) = 1;
            
            %biomass exchange
            fbaSMat(...
                subIdxs_biomass,...
                rxnIdxs_biomassExchange,...
                this.compartmentIndexs_cytosol) = 1;
            
            %reshape
            fbaSMatFull = fbaSMat;
            fbaSMatFull(subIdxs_substrates, rxnIdxs_biomassProduction, :) = ...
                -permute(this.metabolismProduction, [1 3 2]);
            fbaSMat = reshape(...
                permute(fbaSMat, [2 1 3]), ...
                size(fbaSMat, 2), [])';
            fbaSMatFull = reshape(...
                permute(fbaSMatFull, [2 1 3]), ...
                size(fbaSMatFull, 2), [])';
            
            %% reaction catalysis
            fbaCatMat = zeros(nReactions, nEnzymes);
            fbaCatMat(rxnIdxs_metabolicConversion, :) = this.reactionCatalysisMatrix;
            
            %% reaction bounds
            fbaRxnBnds = repmat([-Inf Inf], nReactions, 1);
            
            this.reactionBounds(this.reactionBounds < 0) = -Inf;
            this.reactionBounds(this.reactionBounds > 0) =  Inf;
            fbaRxnBnds(rxnIdxs_metabolicConversion, :) = this.reactionBounds;
            
            exchangeLBs = this.substrateExternalExchangeBounds(:, 1);
            exchangeUBs = this.substrateExternalExchangeBounds(:, 2);
            fbaRxnBnds(rxnIdxs_metaboliteExternalExchange(~isnan(exchangeLBs)), 1) = exchangeLBs(~isnan(exchangeLBs));
            fbaRxnBnds(rxnIdxs_metaboliteExternalExchange(~isnan(exchangeUBs)), 2) = exchangeUBs(~isnan(exchangeUBs));
            
            fbaRxnBnds(rxnIdxs_metaboliteInternalLimitedExchange, 1) = -Inf;
            fbaRxnBnds(rxnIdxs_metaboliteInternalLimitedExchange, 2) = 0;
            fbaRxnBnds(rxnIdxs_metaboliteInternalUnlimitedExchange, 1) = 0;
            fbaRxnBnds(rxnIdxs_metaboliteInternalUnlimitedExchange, 2) = Inf;
            
            fbaRxnBnds(rxnIdxs_biomassProduction, :) = [0 Inf];
            fbaRxnBnds(rxnIdxs_biomassExchange, :) = [-Inf 0];
            
            %% enzyme kinetics
            fbaEnzBnds = NaN(nReactions, 2);
            fbaEnzBnds(rxnIdxs_metabolicConversion, :) = this.enzymeBounds;
            
            %% right hand side
            fbaRHS = zeros(nSubstrates * nCompartments, 1);
            
            %% objective
            fbaObj = zeros(nReactions, 1);
            fbaObj(rxnIdxs_biomassProduction) = 1e3;
            fbaObj(rxnIdxs_metaboliteInternalLimitedExchange) = 1 / sum(min(0, this.metabolismNewProduction(:)));
            
            %% simplify stoichiometry matrix
            subCmpIdxs_fba = (1:size(fbaSMat, 1))';
            rxnIdxs_fba = (1:size(fbaSMat, 2))';
            
            subCmpRemoved = zeros(size(fbaSMat, 1), 3);
            rxnRemoved = zeros(size(fbaSMat, 2), 3);
            
            if ~exist('simplifyNetwork', 'var') || simplifyNetwork
                iter = 0;
                while true
                    iter = iter + 1;
                    
                    %factor
                    [subCmpBlockIdxs, rxnBlockIdxs] = findNonInteractingRowsAndColumns(fbaSMatFull(subCmpIdxs_fba, rxnIdxs_fba));
                    if all(subCmpBlockIdxs == 1)
                        break;
                    end
                    
                    blockIdx = subCmpBlockIdxs(subCmpIdxs_fba == subCmpIdxs_biomass);
                    subCmpRemoved(subCmpIdxs_fba(subCmpBlockIdxs ~= blockIdx), 1) = iter;
                    rxnRemoved(rxnIdxs_fba(rxnBlockIdxs ~= blockIdx), 1) = iter;
                    subCmpRemoved(subCmpIdxs_fba(subCmpBlockIdxs ~= blockIdx), 3) = 1;
                    rxnRemoved(rxnIdxs_fba(rxnBlockIdxs ~= blockIdx), 3) = 1;
                    subCmpIdxs_fba = subCmpIdxs_fba(subCmpBlockIdxs == blockIdx);
                    rxnIdxs_fba = rxnIdxs_fba(rxnBlockIdxs == blockIdx);
                    
                    %prune
                    flags = true(3, 1);
                    iter2 = 0;
                    while any(flags)
                        iter2 = iter2 + 1;
                        flags = true(3, 1);
                        
                        %remove leaf metabolites (and reactions involving them)
                        %which cannot be balanced with non-zero fluxes
                        while flags(1)
                            tmpSubCmpTfs = sum(fbaSMatFull(subCmpIdxs_fba, rxnIdxs_fba) ~= 0, 2) == 1;
                            tmpRxnTFs = any(fbaSMatFull(subCmpIdxs_fba(tmpSubCmpTfs), rxnIdxs_fba), 1);
                            
                            subCmpRemoved(subCmpIdxs_fba(tmpSubCmpTfs), 1) = iter;
                            rxnRemoved(rxnIdxs_fba(tmpRxnTFs), 1) = iter;
                            subCmpRemoved(subCmpIdxs_fba(tmpSubCmpTfs), 2) = iter2;
                            rxnRemoved(rxnIdxs_fba(tmpRxnTFs), 2) = iter2;
                            subCmpRemoved(subCmpIdxs_fba(tmpSubCmpTfs), 3) = 2;
                            rxnRemoved(rxnIdxs_fba(tmpRxnTFs), 3) = 2;
                            
                            subCmpIdxs_fba = subCmpIdxs_fba(~tmpSubCmpTfs);
                            rxnIdxs_fba = rxnIdxs_fba(~tmpRxnTFs);
                            
                            flags(1) = any(tmpSubCmpTfs);
                        end
                        
                        %remove metabolites (and reactions involving them)
                        %which can't move in both the forward and backward
                        %directions and thus cannot be balanced
                        tmpSubCmpTfs = ...
                            (any(fbaSMatFull(subCmpIdxs_fba, rxnIdxs_fba) .* repmat(fbaRxnBnds(rxnIdxs_fba, 1)', numel(subCmpIdxs_fba), 1) < 0, 2) | ...
                            any(fbaSMatFull(subCmpIdxs_fba, rxnIdxs_fba) .* repmat(fbaRxnBnds(rxnIdxs_fba, 2)', numel(subCmpIdxs_fba), 1) < 0, 2)) & ...
                            (any(fbaSMatFull(subCmpIdxs_fba, rxnIdxs_fba) .* repmat(fbaRxnBnds(rxnIdxs_fba, 1)', numel(subCmpIdxs_fba), 1) > 0, 2) | ...
                            any(fbaSMatFull(subCmpIdxs_fba, rxnIdxs_fba) .* repmat(fbaRxnBnds(rxnIdxs_fba, 2)', numel(subCmpIdxs_fba), 1) > 0, 2));
                        tmpRxnTFs = ~any(fbaSMatFull(subCmpIdxs_fba(~tmpSubCmpTfs), rxnIdxs_fba), 1);
                        
                        subCmpRemoved(subCmpIdxs_fba(~tmpSubCmpTfs), 1) = iter;
                        rxnRemoved(rxnIdxs_fba(~tmpRxnTFs), 1) = iter;
                        subCmpRemoved(subCmpIdxs_fba(~tmpSubCmpTfs), 2) = iter2;
                        rxnRemoved(rxnIdxs_fba(~tmpRxnTFs), 2) = iter2;
                        subCmpRemoved(subCmpIdxs_fba(~tmpSubCmpTfs), 3) = 3;
                        rxnRemoved(rxnIdxs_fba(~tmpRxnTFs), 3) = 3;
                        
                        subCmpIdxs_fba = subCmpIdxs_fba(tmpSubCmpTfs);
                        rxnIdxs_fba = rxnIdxs_fba(tmpRxnTFs);
                        
                        flags(2) = ~all(tmpSubCmpTfs);
                        
                        %remove other reactions (and substrates unique to them)
                        %whose flux is zero in all solutions of the matrix
                        %problem Sv=0. That is reactions which have only zero
                        %components in the null space of S.
                        tmpRxnTFs = any(null(fbaSMatFull(subCmpIdxs_fba, rxnIdxs_fba), 'r'), 2);
                        tmpSubCmpTfs = any(fbaSMatFull(subCmpIdxs_fba, rxnIdxs_fba(tmpRxnTFs)), 2);
                        
                        subCmpRemoved(subCmpIdxs_fba(~tmpSubCmpTfs), 1) = iter;
                        rxnRemoved(rxnIdxs_fba(~tmpRxnTFs), 1) = iter;
                        subCmpRemoved(subCmpIdxs_fba(~tmpSubCmpTfs), 2) = iter2;
                        rxnRemoved(rxnIdxs_fba(~tmpRxnTFs), 2) = iter2;
                        subCmpRemoved(subCmpIdxs_fba(~tmpSubCmpTfs), 3) = 4;
                        rxnRemoved(rxnIdxs_fba(~tmpRxnTFs), 3) = 4;
                        
                        subCmpIdxs_fba = subCmpIdxs_fba(tmpSubCmpTfs);
                        rxnIdxs_fba = rxnIdxs_fba(tmpRxnTFs);
                        
                        flags(3) = ~all(tmpSubCmpTfs);
                    end
                end
            end
            
            subCmpIdxs_fba = unique([subCmpIdxs_fba; subCmpIdxs_metaboliteInternalExchangeConstraints]);
            subCmpRemoved(subCmpIdxs_metaboliteInternalExchangeConstraints, :) = 0;
            
            subCmpRemoved = subCmpRemoved(reshape(sub2ind([nSubstrates nCompartments], ...
                repmat(subIdxs_substrates, 1, nCompartments), ...
                repmat(1:nCompartments, numel(subIdxs_substrates), 1)), [], 1), :);
            rxnRemoved = rxnRemoved(rxnIdxs_metabolicConversion, :);
            
            %substrate indices
            [subIdxs_fba, cmpIdxs_fba] = ind2sub([nSubstrates, nCompartments], subCmpIdxs_fba);
            this.fbaSubstrateIndexs_substrates = find(ismember(subIdxs_fba, subIdxs_substrates));
            this.fbaSubstrateIndexs_metaboliteInternalExchangeConstraints = find(ismember(subCmpIdxs_fba, subCmpIdxs_metaboliteInternalExchangeConstraints));
            this.substrateIndexs_fba = sub2ind(...
                [numel(this.substrateWholeCellModelIDs) nCompartments], ...
                subIdxs_fba(this.fbaSubstrateIndexs_substrates), ...
                cmpIdxs_fba(this.fbaSubstrateIndexs_substrates));
            [tfs, idxs] = ismember(rxnIdxs_fba, rxnIdxs_metaboliteExternalExchange);
            this.substrateIndexs_externalExchangedMetabolites = idxs(tfs);
            this.substrateIndexs_internalExchangedMetabolites = ...
                find(this.metabolismRecyclingProduction);
            this.substrateIndexs_internalExchangedLimitedMetabolites = ...
                find(this.metabolismRecyclingProduction < 0);
            
            this.fbaSubstrateIndexs_biomass = find(ismember(subCmpIdxs_fba, subCmpIdxs_biomass));
            
            %reaction indices
            assert(all(ismember(rxnIdxs_metaboliteInternalLimitedExchange, rxnIdxs_fba)))
            
            this.reactionIndexs_fba = find(ismember(rxnIdxs_metabolicConversion, rxnIdxs_fba));
            this.fbaReactionIndexs_metabolicConversion = find(ismember(rxnIdxs_fba, rxnIdxs_metabolicConversion));
            this.fbaReactionIndexs_metaboliteExternalExchange = find(ismember(rxnIdxs_fba, rxnIdxs_metaboliteExternalExchange));
            this.fbaReactionIndexs_metaboliteInternalExchange = find(ismember(rxnIdxs_fba, rxnIdxs_metaboliteInternalExchange));
            this.fbaReactionIndexs_metaboliteInternalLimitedExchange = find(ismember(rxnIdxs_fba, rxnIdxs_metaboliteInternalLimitedExchange));
            this.fbaReactionIndexs_metaboliteInternalUnlimitedExchange = find(ismember(rxnIdxs_fba, rxnIdxs_metaboliteInternalUnlimitedExchange));
            this.fbaReactionIndexs_biomassProduction = find(ismember(rxnIdxs_fba, rxnIdxs_biomassProduction));
            this.fbaReactionIndexs_biomassExchange = find(ismember(rxnIdxs_fba, rxnIdxs_biomassExchange));
            
            fbaSMat = fbaSMat(subCmpIdxs_fba, rxnIdxs_fba);
            fbaCatMat = fbaCatMat(rxnIdxs_fba, :);
            fbaRHS = fbaRHS(subCmpIdxs_fba, :);
            fbaRxnBnds = fbaRxnBnds(rxnIdxs_fba, :);
            fbaEnzBnds = fbaEnzBnds(rxnIdxs_fba, :);
            fbaObj = fbaObj(rxnIdxs_fba, :);
            
            %% set properties
            this.fbaReactionStoichiometryMatrix = fbaSMat;
            this.fbaReactionCatalysisMatrix = fbaCatMat;
            this.fbaRightHandSide = fbaRHS;
            this.fbaReactionBounds = fbaRxnBnds;
            this.fbaEnzymeBounds = fbaEnzBnds;
            this.fbaObjective = fbaObj;
        end
    end
    
    %model
    methods
        %Calculate
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, ~, states)
            %% import classes
            import edu.stanford.covert.util.ComputationUtil;
            
            %% substrate and byproducts: none used by metabolism
            bmProd = zeros(size(this.substrateWholeCellModelIDs));
            byProd = zeros(size(this.substrateWholeCellModelIDs));
            
            %% fit enzyme expression consistent with growth rate
            if strcmp(this.macromoleculeStateInitialization, 'expected')
                growth0 = this.metabolicReaction.growth0;
            else
                growth0 = this.metabolicReaction.growth0 * this.macromoleculeStateInitializationGrowthFactor;
            end
            
            substrates = zeros(size(this.substrates));
            substrates(this.substrateMetaboliteLocalIndexs, :) = this.substrates(this.substrateMetaboliteLocalIndexs, :);
            substrates(this.substrateIndexs_internalExchangedLimitedMetabolites) = ...
                -growth0 *  this.metabolismRecyclingProduction(this.substrateIndexs_internalExchangedLimitedMetabolites);
            substrates(this.substrateMonomerLocalIndexs, :) = repmat(states.monomers0(this.substrateMonomerGlobalIndexs(:, 1)), 1, size(substrates, 2));
            substrates(this.substrateComplexLocalIndexs, :) = repmat(states.complexs0(this.substrateComplexGlobalIndexs(:, 1)), 1, size(substrates, 2));
            
            enzymes = zeros(size(this.enzymeWholeCellModelIDs));
            enzymes(this.enzymeMonomerLocalIndexs) = states.monomers0(this.enzymeMonomerGlobalIndexs);
            enzymes(this.enzymeComplexLocalIndexs) = states.complexs0(this.enzymeComplexGlobalIndexs);
            
            enzymes = this.fitEnzymes(growth0, substrates, enzymes);
            
            %% metabolic enzyme levels constistent with observed growth rate
            [minEnzExp, maxEnzExp] = this.calcMinMaxEnzymes(substrates, enzymes);
        end
        
        %calculate enzyme expression consistent with growth rate
        %- by minimizing the fractional change in fluxes by linear
        %  optimization (similar to MOMA)
        %- if linear optimization fails, increase/decrease all flux
        %  bounds by a factor of growth/growth0 and check if this yields
        %  a growth rate of growth0
        function enzymes = fitEnzymes(this, growth0, substrates, enzymes)
            % import classes
            import edu.stanford.covert.util.ComputationUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            %constants
            nEnzymes = numel(enzymes);
            enzMWs = this.enzymeMolecularWeights;
            
            minAvgExp = zeros(size(enzymes));
            minAvgExp(this.enzymeMonomerLocalIndexs) = this.monomer.minimumAverageExpression;
            minAvgExp(this.enzymeComplexLocalIndexs) = this.complex.minimumAverageExpression;
            
            loIdxs = find(enzymes < minAvgExp & any(this.fbaReactionCatalysisMatrix, 1)');
            hiIdxs = find(enzymes > minAvgExp);
            if (minAvgExp(loIdxs)-enzymes(loIdxs))' * enzMWs(loIdxs) < (enzymes(hiIdxs) - minAvgExp(hiIdxs))' * enzMWs(hiIdxs)
                enzymes(hiIdxs) = minAvgExp(hiIdxs) + (enzymes(hiIdxs) - minAvgExp(hiIdxs)) * ...
                    ((enzymes(hiIdxs)-minAvgExp(hiIdxs))' * enzMWs(hiIdxs) - ...
                    (minAvgExp(loIdxs)-enzymes(loIdxs))' * enzMWs(loIdxs)) / ...
                    ((enzymes(hiIdxs)-minAvgExp(hiIdxs))' * enzMWs(hiIdxs));
                enzymes(loIdxs) = minAvgExp(loIdxs);
            else
                warning('WholeCell:warning', 'Cannot ensure minimum enzyme expression');
            end
            
            initEnzymes = enzymes;
            
            % calculate current growth and fluxes
            [fbaObj, fbaSMat, fbaRxnBounds] = this.calcEffectiveFBANetwork();
            fluxBounds = this.calcFluxBounds(substrates, enzymes, fbaRxnBounds, this.fbaEnzymeBounds);
            [growth, ~, oldFbaReactionFluxs] = this.calcGrowthRate(fluxBounds, fbaObj, fbaSMat);
            
            if abs(growth - growth0) / growth0 < this.tolerance || abs(growth - growth0) < 1e-12
                return;
            end
            
            if growth < growth0 && growth0 > this.calcGrowthRate(this.calcFluxBounds(...
                    substrates, enzymes, fbaRxnBounds, this.fbaEnzymeBounds, false), fbaObj, fbaSMat)
                throw(MException('Metabolism:error', 'Cannot fit enzyme expression to match growth rate'));
            end
            
            % reactions constrainted by enzyme expression
            [forIdxs, revIdxs] = this.calcEnzymeKineticallyLimitedReactions(substrates, enzymes, fbaObj, fbaSMat, fbaRxnBounds);
            unconstrainedEnzIdxs = find(~any(this.fbaReactionCatalysisMatrix([forIdxs; revIdxs], :), 1)');
            
            % minimize difference from current enzyme expression
            objectiveFunc = zeros(size(fbaObj));
            objectiveFunc([forIdxs; revIdxs]) = 1 ./ oldFbaReactionFluxs([forIdxs; revIdxs]);
            
            % minimum enzyme expression
            loFluxBounds = fluxBounds(:, 1);
            upFluxBounds = fluxBounds(:, 2);
            
            if growth < growth0
                tmpFluxBounds = this.calcFluxBounds(...
                    substrates, enzymes, fbaRxnBounds, this.fbaEnzymeBounds, false);
                loFluxBounds(forIdxs) = oldFbaReactionFluxs(forIdxs);
                upFluxBounds(forIdxs) = min(tmpFluxBounds(forIdxs, 2), oldFbaReactionFluxs(forIdxs) * (1 + this.tolerance) * growth0 / growth);
                loFluxBounds(revIdxs) = max(tmpFluxBounds(revIdxs, 1), oldFbaReactionFluxs(revIdxs) * (1 + this.tolerance) * growth0 / growth);
                upFluxBounds(revIdxs) = oldFbaReactionFluxs(revIdxs);
                
                objective = 'minimize';
            else
                loFluxBounds(forIdxs) = 0;
                upFluxBounds(forIdxs) = oldFbaReactionFluxs(forIdxs);
                loFluxBounds(revIdxs) = oldFbaReactionFluxs(revIdxs);
                upFluxBounds(revIdxs) = 0;
                
                objective = 'maximize';
            end
            
            loFluxBounds(this.fbaReactionIndexs_biomassProduction) = growth0;
            upFluxBounds(this.fbaReactionIndexs_biomassProduction) = growth0;
            
            loFluxBounds = max(loFluxBounds, -this.realmax);
            upFluxBounds = min(upFluxBounds,  this.realmax);
            [fbaReactionFluxs, ~, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                objective, objectiveFunc, fbaSMat, ...
                this.fbaRightHandSide, loFluxBounds, upFluxBounds, ...
                'S', 'C', this.linearProgrammingOptions);
            if errFlag
                throw(MException('Metabolism:error', 'Unable to optimize fluxes: %s', errMsg));
            end
            
            %adjust enzyme expression
            fbaReactionFluxs = max(min(fbaReactionFluxs, upFluxBounds), loFluxBounds);
            
            reactionEnzymes = zeros(size(fbaReactionFluxs));
            reactionEnzymes(forIdxs) = fbaReactionFluxs(forIdxs) ./ this.fbaEnzymeBounds(forIdxs, 2);
            reactionEnzymes(revIdxs) = fbaReactionFluxs(revIdxs) ./ this.fbaEnzymeBounds(revIdxs, 1);
            
            enzymes = max(this.fbaReactionCatalysisMatrix([forIdxs; revIdxs], :) .* ...
                repmat(reactionEnzymes([forIdxs; revIdxs]), 1, nEnzymes), [], 1)';
            enzymes(unconstrainedEnzIdxs) = initEnzymes(unconstrainedEnzIdxs);
            enzymes = this.adjustEnzymeExpressionSparsely(growth0, substrates, enzymes, initEnzymes);
            enzymes = max(minAvgExp, enzymes);
            
            if growth < growth0
                incEnzymes = find(enzymes > initEnzymes);
                decEnzymes = setdiff(unique([find(enzymes < initEnzymes); unconstrainedEnzIdxs]), incEnzymes);
                enzymes = max(initEnzymes, enzymes);
                enzymes(decEnzymes) = minAvgExp(decEnzymes) + (initEnzymes(decEnzymes) - minAvgExp(decEnzymes)) * ...
                    ((initEnzymes(decEnzymes) - minAvgExp(decEnzymes))' * enzMWs(decEnzymes) - ...
                    (enzymes(incEnzymes) - initEnzymes(incEnzymes))' * enzMWs(incEnzymes)) / ...
                    ((initEnzymes(decEnzymes) - minAvgExp(decEnzymes))' * enzMWs(decEnzymes));
            else
                incEnzymes = unique([find(enzymes > initEnzymes); unconstrainedEnzIdxs]);
                decEnzymes = setdiff(find(enzymes < initEnzymes), incEnzymes);
                enzymes = min(initEnzymes, enzymes);
                enzymes(incEnzymes) = minAvgExp(incEnzymes) + (initEnzymes(incEnzymes) - minAvgExp(incEnzymes)) * ...
                    ((initEnzymes(incEnzymes) - minAvgExp(incEnzymes))' * enzMWs(incEnzymes) - ...
                    (enzymes(decEnzymes) - initEnzymes(decEnzymes))' * enzMWs(decEnzymes)) / ...
                    ((initEnzymes(incEnzymes) - minAvgExp(incEnzymes))' * enzMWs(incEnzymes));
            end
            
            %check growth rate fit
            if abs(this.calcGrowthRate(this.calcFluxBounds(...
                    substrates, enzymes, fbaRxnBounds, this.fbaEnzymeBounds), fbaObj, fbaSMat) ...
                    - growth0) / growth0 > this.tolerance
                throw(MException('Metabolism:error', 'Growth rate should not change'));
            end
            if abs((enzymes' * enzMWs - initEnzymes' * enzMWs) / (initEnzymes' * enzMWs)) > this.tolerance
                throw(MException('Metabolism:error', 'Protein mass should not change'));
            end
        end
        
        function enzymes = adjustEnzymeExpressionSparsely(this, growth0, substrates, enzymes, initEnzymes)
            [fbaObj, fbaSMat, fbaRxnBounds] = this.calcEffectiveFBANetwork();
            
            tfs = this.adjustEnzymeExpressionSparselyHelper(growth0, substrates, enzymes, initEnzymes, true(size(enzymes)));
            enzymes(tfs) = initEnzymes(tfs);
            idxs = find(~tfs);
            
            %one at a time
            tf = false;
            for i = 1:numel(idxs)
                tmpEnzymes = enzymes;
                tmpEnzymes(idxs(i)) = initEnzymes(idxs(i));
                
                if abs(this.calcGrowthRate(this.calcFluxBounds(...
                        substrates, tmpEnzymes, fbaRxnBounds, this.fbaEnzymeBounds), fbaObj, fbaSMat) ...
                        - growth0) / growth0 <= this.tolerance
                    tf = true;
                    break;
                end
            end
            if ~tf
                return;
            end
            
            %combinations
            try
                tmpTfs = dec2bin(0:2^numel(idxs)-1) == '1';
            catch %#ok<CTCH>
                return;
            end
            [~, order] = sort(sum(tmpTfs, 2), 1, 'descend');
            tmpTfs = tmpTfs(order, :);
            
            for i = 1:2^numel(idxs)
                tmpEnzymes = enzymes;
                tmpEnzymes(idxs(tmpTfs(i, :))) = initEnzymes(idxs(tmpTfs(i, :)));
                
                if abs(this.calcGrowthRate(this.calcFluxBounds(...
                        substrates, tmpEnzymes, fbaRxnBounds, this.fbaEnzymeBounds), fbaObj, fbaSMat) ...
                        - growth0) / growth0 <= this.tolerance
                    enzymes = tmpEnzymes;
                    return;
                end
            end
        end
        
        function tfs = adjustEnzymeExpressionSparselyHelper(this, growth0, substrates, enzymes, initEnzymes, tfs)
            [fbaObj, fbaSMat, fbaRxnBounds] = this.calcEffectiveFBANetwork();
            
            tmpEnzymes = enzymes;
            tmpEnzymes(tfs) = initEnzymes(tfs);
            
            if sum(tfs) == 0 || abs(this.calcGrowthRate(this.calcFluxBounds(...
                    substrates, tmpEnzymes, fbaRxnBounds, this.fbaEnzymeBounds), fbaObj, fbaSMat) ...
                    - growth0) / growth0 <= this.tolerance
                return;
            elseif sum(tfs) == 1
                tfs(:) = false;
                return;
            end
            
            idxs = find(tfs);
            tfs1 = tfs;
            tfs2 = tfs;
            tfs1(idxs(1:ceil(end/2))) = false;
            tfs2(idxs(ceil(end/2))+1:end) = false;
            tfs = ...
                this.adjustEnzymeExpressionSparselyHelper(growth0, substrates, enzymes, initEnzymes, tfs1) | ...
                this.adjustEnzymeExpressionSparselyHelper(growth0, substrates, enzymes, initEnzymes, tfs2);
        end
        
        %calculate minimum and maximum enzyme weight consistent with growth
        %rate
        %- by minimizing/maximizing the total enzyme weight consistent with
        %  the current growth rate by linear optimization
        %- if linear optimization fails, minEnzymes/maxEnzymes is set to
        %  the current enyzme level
        function [minEnzymes, maxEnzymes] = calcMinMaxEnzymes(this, substrates, enzymes)
            % import classes
            import edu.stanford.covert.util.ComputationUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            %constants
            nEnzymes = numel(enzymes);
            
            minAvgExp = zeros(size(enzymes));
            minAvgExp(this.enzymeMonomerLocalIndexs) = this.monomer.minimumAverageExpression;
            minAvgExp(this.enzymeComplexLocalIndexs) = this.complex.minimumAverageExpression;
            
            % calculate current growth and fluxes
            [fbaObj, fbaSMat, fbaRxnBounds] = this.calcEffectiveFBANetwork();
            fluxBounds = this.calcFluxBounds(...
                substrates, enzymes, fbaRxnBounds, this.fbaEnzymeBounds);
            [growth, ~, oldFbaReactionFluxs] = this.calcGrowthRate(fluxBounds, fbaObj, fbaSMat);
            
            % reactions constrainted by enzyme expression
            [forIdxs, revIdxs] = this.calcEnzymeKineticallyLimitedReactions(substrates, enzymes, fbaObj, fbaSMat, fbaRxnBounds);
            partlyUnconstrainedEnzIdxs = find(any(this.fbaReactionCatalysisMatrix(setdiff(1:end, [forIdxs; revIdxs]), :), 1));
            unconstrainedEnzIdxs = find(~any(this.fbaReactionCatalysisMatrix([forIdxs; revIdxs], :), 1));
            
            if isempty(forIdxs) && isempty(revIdxs)
                warning('WholeCell:warning', 'Unable to fit growth rate');
                minEnzymes = enzymes;
                maxEnzymes = enzymes;
                return;
            end
            
            % enzyme expression range
            objectiveFunc = zeros(size(fbaObj));
            objectiveFunc(forIdxs) = (this.fbaReactionCatalysisMatrix(forIdxs, :) * this.enzymeMolecularWeights) ./ this.fbaEnzymeBounds(forIdxs, 2);
            objectiveFunc(revIdxs) = (this.fbaReactionCatalysisMatrix(revIdxs, :) * this.enzymeMolecularWeights) ./ this.fbaEnzymeBounds(revIdxs, 1);
            
            % minimum enzyme expression
            loFluxBounds = fluxBounds(:, 1);
            upFluxBounds = fluxBounds(:, 2);
            
            loFluxBounds(this.fbaReactionIndexs_biomassProduction) = growth;
            upFluxBounds(this.fbaReactionIndexs_biomassProduction) = growth;
            
            loFluxBounds(forIdxs) = 0;
            upFluxBounds(forIdxs) = oldFbaReactionFluxs(forIdxs);
            loFluxBounds(revIdxs) = oldFbaReactionFluxs(revIdxs);
            upFluxBounds(revIdxs) = 0;
            
            loFluxBounds = max(loFluxBounds, -this.realmax);
            upFluxBounds = min(upFluxBounds,  this.realmax);
            [fbaReactionFluxs, ~, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                'minimize', objectiveFunc, fbaSMat, ...
                this.fbaRightHandSide, loFluxBounds, upFluxBounds, ...
                'S', 'C', this.linearProgrammingOptions);
            if errFlag
                warning('WholeCell:warning', errMsg);
                minEnzymes = enzymes;
            else
                fbaReactionFluxs = max(min(fbaReactionFluxs, upFluxBounds), loFluxBounds);
                
                reactionEnzymes = zeros(size(oldFbaReactionFluxs));
                reactionEnzymes(forIdxs) = fbaReactionFluxs(forIdxs) ./ this.fbaEnzymeBounds(forIdxs, 2);
                reactionEnzymes(revIdxs) = fbaReactionFluxs(revIdxs) ./ this.fbaEnzymeBounds(revIdxs, 1);
                
                minEnzymes = max(this.fbaReactionCatalysisMatrix([forIdxs; revIdxs], :) .* ...
                    repmat(reactionEnzymes([forIdxs; revIdxs]), 1, nEnzymes), [], 1)';
                minEnzymes(unconstrainedEnzIdxs) = 1;
                minEnzymes(partlyUnconstrainedEnzIdxs) = max(minEnzymes(partlyUnconstrainedEnzIdxs), 1);
                minEnzymes = min(enzymes, minEnzymes);
                minEnzymes = max(minAvgExp, minEnzymes);
            end
            
            if abs(this.calcGrowthRate(this.calcFluxBounds(...
                    substrates, minEnzymes, fbaRxnBounds, this.fbaEnzymeBounds), fbaObj, fbaSMat) ...
                    - growth) / growth > this.tolerance || ...
                    any(minEnzymes > enzymes)
                warning('WholeCell:warning', 'Unable to calculate minimum enzymes');
                minEnzymes = enzymes;
            end
            
            % maximum enzyme expression
            loFluxBounds = fluxBounds(:, 1);
            upFluxBounds = fluxBounds(:, 2);
            
            loFluxBounds(this.fbaReactionIndexs_biomassProduction) = growth;
            upFluxBounds(this.fbaReactionIndexs_biomassProduction) = growth;
            
            tmpFluxBounds = this.calcFluxBounds(...
                substrates, enzymes, fbaRxnBounds, this.fbaEnzymeBounds, false);
            loFluxBounds(forIdxs) = oldFbaReactionFluxs(forIdxs);
            upFluxBounds(forIdxs) = tmpFluxBounds(forIdxs, 2);
            loFluxBounds(revIdxs) = tmpFluxBounds(revIdxs, 1);
            upFluxBounds(revIdxs) = oldFbaReactionFluxs(revIdxs);
            
            loFluxBounds = max(loFluxBounds, -this.realmax);
            upFluxBounds = min(upFluxBounds,  this.realmax);
            [fbaReactionFluxs, ~, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                'maximize', objectiveFunc, fbaSMat, ...
                this.fbaRightHandSide, loFluxBounds, upFluxBounds, ...
                'S', 'C', this.linearProgrammingOptions);
            if errFlag
                warning('WholeCell:warning', errMsg);
                maxEnzymes = enzymes;
            else
                fbaReactionFluxs = max(min(fbaReactionFluxs, upFluxBounds), loFluxBounds);
                
                reactionEnzymes = zeros(size(oldFbaReactionFluxs));
                reactionEnzymes(forIdxs) = fbaReactionFluxs(forIdxs) ./ this.fbaEnzymeBounds(forIdxs, 2);
                reactionEnzymes(revIdxs) = fbaReactionFluxs(revIdxs) ./ this.fbaEnzymeBounds(revIdxs, 1);
                
                tmpReactionCatalysisMatrix = this.fbaReactionCatalysisMatrix;
                tmpReactionCatalysisMatrix(tmpReactionCatalysisMatrix == 0) = NaN;
                maxEnzymes = 0*min(tmpReactionCatalysisMatrix([forIdxs; revIdxs], :) .* ...
                    repmat(reactionEnzymes([forIdxs; revIdxs]), 1, nEnzymes), [], 1)';
                maxEnzymes(unconstrainedEnzIdxs) = realmax; %#ok<CPROP,PROP>
                maxEnzymes(partlyUnconstrainedEnzIdxs) = max(maxEnzymes(partlyUnconstrainedEnzIdxs), 1);
                maxEnzymes = max(enzymes, maxEnzymes);
            end
            
            if abs(this.calcGrowthRate(this.calcFluxBounds(...
                    substrates, maxEnzymes, fbaRxnBounds, this.fbaEnzymeBounds), fbaObj, fbaSMat) ...
                    - growth) / growth > this.tolerance || ...
                    any(maxEnzymes < enzymes)
                warning('WholeCell:warning', 'Unable to calculate maximum enzymes');
                maxEnzymes = enzymes;
            end
        end
        
        %calculate indices of reactions whose flux is actively limited by
        %the product of the expression of the catalyzing enzyme and its
        %reaction kinetics
        function [forIdxs, revIdxs] = calcEnzymeKineticallyLimitedReactions(this, substrates, enzymes, fbaObj, fbaSMat, fbaRxnBounds)
            fluxBounds = this.calcFluxBounds(...
                substrates, enzymes, fbaRxnBounds, this.fbaEnzymeBounds);
            fluxBounds_nokinetics = this.calcFluxBounds(...
                substrates, enzymes, fbaRxnBounds, this.fbaEnzymeBounds, false);
            [~, ~, fbaReactionFluxs] = this.calcGrowthRate(fluxBounds, fbaObj, fbaSMat);
            
            rxnIdxs = this.fbaReactionIndexs_metabolicConversion;
            forTfs = ...
                fbaReactionFluxs(rxnIdxs) > 0 & ...
                fluxBounds(rxnIdxs, 2) < fluxBounds_nokinetics(rxnIdxs, 2);
            revTfs = ...
                fbaReactionFluxs(rxnIdxs) < 0 & ...
                fluxBounds(rxnIdxs, 1) > fluxBounds_nokinetics(rxnIdxs, 1);
            
            forIdxs = rxnIdxs(forTfs);
            revIdxs = rxnIdxs(revTfs);
        end
        
        function [fbaObj, fbaSMat, fbaRxnBounds] = calcEffectiveFBANetwork(this, n)
            if nargin < 2
                n = 4;
            end
                
            fbaObj = this.fbaObjective;
            fbaObj(:) = 0;
            fbaObj(this.fbaReactionIndexs_biomassProduction) = 1;
            
            fbaSMat = this.fbaReactionStoichiometryMatrix;
            fbaSMat(this.fbaSubstrateIndexs_substrates, this.fbaReactionIndexs_biomassProduction) = ...
                - this.metabolismNewProduction(this.substrateIndexs_fba) ...
                - n * this.metabolismRecyclingProduction(this.substrateIndexs_fba);
            
            fbaRxnBounds = this.fbaReactionBounds;
            fbaRxnBounds(this.fbaReactionIndexs_metaboliteInternalExchange, :) = 0;
        end
        
        %initialization
        %- Compute growth rate, reaction fluxs, reduced costs, duals
        %- Don't update metabolite counts
        function initializeState(this)
            [this.metabolicReaction.growth, this.metabolicReaction.fluxs] = ...
                this.calcGrowthRate(this.calcFluxBounds(this.substrates, this.enzymes, this.fbaReactionBounds, this.fbaEnzymeBounds));
        end
        
        %resource requirements
        function result = calcResourceRequirements_Current(this)
            result = zeros(size(this.substrates));
            result(this.substrateIndexs_externalExchangedMetabolites, this.compartmentIndexs_extracellular) = ...
                this.fbaReactionBounds(this.fbaReactionIndexs_metaboliteExternalExchange, 2) * sum(this.mass.cellDry) * this.stepSizeSec;
            result = ...
                + result ...
                - this.metabolismRecyclingProduction * 2 * max(...
                    sum(this.mass.cellDry) / this.mass.cellInitialDryWeight * log(2) / this.cellCycleLength, ...
                    this.metabolicReaction.growth);
        end
        
        %simulation
        function evolveState(this)
            %% Calculate growth rate
            [this.metabolicReaction.growth, this.metabolicReaction.fluxs, fbaReactionFluxs] = ...
                this.calcGrowthRate(this.calcFluxBounds(this.substrates, this.enzymes, this.fbaReactionBounds, this.fbaEnzymeBounds));
            
            %% Compute real-valued amount of nutrient imported, biomass (DNA,
            %RNA, protein, membrane, etc) precursors produced and consumed,
            %and monomers modified. Stochastically round to:
            %1. Instantaneously, maintain integer-valued amounts of metabolites
            %2. Over time, maintain experimentally determined ratios of biomass
            %   components (stored in this.metabolismProduction)
            
            %nutrient uptake
            this.substrates(this.substrateIndexs_externalExchangedMetabolites, this.compartmentIndexs_extracellular) = ...
                + this.substrates(this.substrateIndexs_externalExchangedMetabolites, this.compartmentIndexs_extracellular) ...
                - this.randStream.stochasticRound(fbaReactionFluxs(this.fbaReactionIndexs_metaboliteExternalExchange) * this.stepSizeSec);
            
            %recycled metabolites
            this.substrates(this.substrateIndexs_internalExchangedMetabolites) = ...
                + this.substrates(this.substrateIndexs_internalExchangedMetabolites) ...
                + this.randStream.stochasticRound(fbaReactionFluxs(this.fbaReactionIndexs_metaboliteInternalExchange));
            
            %new metabolites
            this.substrates = ...
                + this.substrates ...
                + this.randStream.stochasticRound(this.metabolismNewProduction * this.metabolicReaction.growth * this.stepSizeSec);
            
            %unaccounted energy consumption
            this.substrates(this.substrateIndexs_atpHydrolysis) = ...
                + this.substrates(this.substrateIndexs_atpHydrolysis) ...
                + [-1; -1; 1; 1; 1] * this.randStream.stochasticRound(...
                this.unaccountedEnergyConsumption * this.metabolicReaction.growth * this.stepSizeSec);
            
            %make metabolites counts positive (in case stochastic rounding
            %made them slightly negative eg -1) except H2O, H+
            if any(any(this.substrates(this.substrateMetaboliteLocalIndexs(:, 1), :) < -1))
                [i, j] = find(this.substrates(this.substrateMetaboliteLocalIndexs(:, 1), :) < -1);
                compIDs = cell(3, 1);
                compIDs(this.compartmentIndexs_cytosol) = {'c'};
                compIDs(this.compartmentIndexs_extracellular) = {'e'};
                compIDs(this.compartmentIndexs_membrane) = {'m'};
                n = max(cellfun(@length, this.substrateWholeCellModelIDs(this.substrateMetaboliteLocalIndexs(i, 1))));
                msg = cellfun(@(id, comp, val) sprintf(['%-' num2str(n) 's  %-3s  %10d'], id, comp, val), ...
                    this.substrateWholeCellModelIDs(this.substrateMetaboliteLocalIndexs(i, 1)), compIDs(j), ...
                    num2cell(this.substrates(sub2ind(size(this.substrates), this.substrateMetaboliteLocalIndexs(i, 1), j))), ...
                    'UniformOutput', false);
                warning('WholeCell:warning:negativeMetabolites', ...
                    ['%d metabolites have negative counts\n%-' num2str(n) 's  %-3s  %10s\n%-' num2str(n) 's  %-3s  %10s\n%s'], ...
                    numel(i), ...
                    'Met', 'Cmp', 'Val', ...
                    repmat('=', 1, n), '===', repmat('=', 1, 10), ...
                    strjoin(sprintf('\n'), msg{:}));
            end
            this.substrates(this.substrateMetaboliteLocalIndexs(:, 1), :) = max(0, this.substrates(this.substrateMetaboliteLocalIndexs(:, 1), :));
            
            %% update cell volume
            this.mass.calcMass();
            this.geometry.calculateVolume();
        end
    end
    
    %evolve state helper methods
    methods
        %calculates the growth rate given the supplied flux bounds
        function [growth, reactionFluxs, fbaReactionFluxs, ...
                reducedCosts, fbaReducedCosts, ...
                shadowPrices, fbaShadowPrices] = ...
                calcGrowthRate(this, fluxBounds, fbaObj, fbaSMat)
            % import classes
            import edu.stanford.covert.util.ComputationUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            %inputs
            if nargin < 3
                fbaObj = this.fbaObjective;
            end
            if nargin < 4
                fbaSMat = this.fbaReactionStoichiometryMatrix;
            end
            
            %flux bounds
            loFluxBounds = fluxBounds(:, 1);
            upFluxBounds = fluxBounds(:, 2);
            
            %real-valued linear programming
            loFluxBounds = max(loFluxBounds, -this.realmax);
            upFluxBounds = min(upFluxBounds,  this.realmax);
            [fbaReactionFluxs, lambda, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                'maximize', fbaObj, fbaSMat, ...
                this.fbaRightHandSide, loFluxBounds, upFluxBounds, ...
                'S', 'C', this.linearProgrammingOptions);
            if errFlag
                warning('WholeCell:warning', 'Linear programming error: %s. Returning feasible, but possibly non-optimal solution x=0.', errMsg);
                fbaReactionFluxs = zeros(size(loFluxBounds));
            end
            
            fbaReactionFluxs = max(min(fbaReactionFluxs, upFluxBounds), loFluxBounds);
            growth = fbaReactionFluxs(this.fbaReactionIndexs_biomassProduction);
            reactionFluxs = zeros(size(this.reactionStoichiometryMatrix, 2), 1);
            reactionFluxs(this.reactionIndexs_fba) = fbaReactionFluxs(this.fbaReactionIndexs_metabolicConversion);
            
            if nargout > 3
                fbaReducedCosts = lambda.reducedCosts;
                reducedCosts = zeros(size(this.reactionStoichiometryMatrix, 2), 1);
                reducedCosts(this.reactionIndexs_fba) = fbaReducedCosts(this.fbaReactionIndexs_metabolicConversion);
                
                fbaShadowPrices = lambda.shadowPrices;
                shadowPrices = zeros(size(this.substrates));
                shadowPrices(this.substrateIndexs_fba) = fbaShadowPrices(this.fbaSubstrateIndexs_substrates);
            end
        end
        
        %Compute reaction flux upper and lower bounds based on
        % 1. Enzyme kinetics         (fbaEnzymeBounds)
        % 2. Enzyme availability     (enzymes)
        % 3. Transport rates         (fbaReactionBounds)
        % 4. Metabolite availability (substrates)
        % 5. Protein availability    (substrates)
        function bounds = calcFluxBounds(this, substrates, enzymes, fbaReactionBounds, fbaEnzymeBounds, ...
                applyEnzymeKineticBounds, applyEnzymeBounds, applyDirectionalityBounds, ...
                applyExternalMetaboliteBounds, applyInternalMetaboliteBounds, applyProteinBounds)
            %import classes
            import edu.stanford.covert.util.ConstantUtil;
            
            %options
            if nargin < 6,  applyEnzymeKineticBounds      = true; end
            if nargin < 7,  applyEnzymeBounds             = true; end
            if nargin < 8,  applyDirectionalityBounds     = true; end
            if nargin < 9,  applyExternalMetaboliteBounds = true; end
            if nargin < 10, applyInternalMetaboliteBounds = true; end
            if nargin < 11, applyProteinBounds            = true; end
            
            %numbers
            nReactions = size(this.fbaReactionStoichiometryMatrix, 2);
            
            %initialize
            lowerBounds = -inf(nReactions, 1);
            upperBounds =  inf(nReactions, 1);
            
            %numbers of enzymes catalyzing each reaction, enzyme kinetics
            rxnEnzymes = this.fbaReactionCatalysisMatrix * enzymes;
            if applyEnzymeKineticBounds
                lowerBounds = max(lowerBounds, fbaEnzymeBounds(:, 1) .* rxnEnzymes);
                upperBounds = min(upperBounds, fbaEnzymeBounds(:, 2) .* rxnEnzymes);
            end
            
            %numbers of enzymes catalyzing each reaction, unkown enzyme kinetics
            if applyEnzymeBounds
                lowerBounds(any(this.fbaReactionCatalysisMatrix, 2) & rxnEnzymes <= 0) = 0;
                upperBounds(any(this.fbaReactionCatalysisMatrix, 2) & rxnEnzymes <= 0) = 0;
            end
            
            %reaction directionality / thermodynamics
            if applyDirectionalityBounds
                lowerBounds(this.fbaReactionIndexs_metabolicConversion) = max(lowerBounds(this.fbaReactionIndexs_metabolicConversion), fbaReactionBounds(this.fbaReactionIndexs_metabolicConversion, 1));
                upperBounds(this.fbaReactionIndexs_metabolicConversion) = min(upperBounds(this.fbaReactionIndexs_metabolicConversion), fbaReactionBounds(this.fbaReactionIndexs_metabolicConversion, 2));
                
                lowerBounds(this.fbaReactionIndexs_metaboliteInternalExchange) = max(lowerBounds(this.fbaReactionIndexs_metaboliteInternalExchange), fbaReactionBounds(this.fbaReactionIndexs_metaboliteInternalExchange, 1));
                upperBounds(this.fbaReactionIndexs_metaboliteInternalExchange) = min(upperBounds(this.fbaReactionIndexs_metaboliteInternalExchange), fbaReactionBounds(this.fbaReactionIndexs_metaboliteInternalExchange, 2));
                
                lowerBounds(this.fbaReactionIndexs_biomassExchange) = max(lowerBounds(this.fbaReactionIndexs_biomassExchange), fbaReactionBounds(this.fbaReactionIndexs_biomassExchange, 1));
                upperBounds(this.fbaReactionIndexs_biomassExchange) = min(upperBounds(this.fbaReactionIndexs_biomassExchange), fbaReactionBounds(this.fbaReactionIndexs_biomassExchange, 2));
                
                lowerBounds(this.fbaReactionIndexs_biomassProduction) = max(lowerBounds(this.fbaReactionIndexs_biomassProduction), fbaReactionBounds(this.fbaReactionIndexs_biomassProduction, 1));
                upperBounds(this.fbaReactionIndexs_biomassProduction) = min(upperBounds(this.fbaReactionIndexs_biomassProduction), fbaReactionBounds(this.fbaReactionIndexs_biomassProduction, 2));
            end
            
            %external metabolite availability
            if applyExternalMetaboliteBounds
                upperBounds(this.fbaReactionIndexs_metaboliteExternalExchange) = min(...
                    upperBounds(this.fbaReactionIndexs_metaboliteExternalExchange), ...
                    substrates(this.substrateIndexs_externalExchangedMetabolites, this.compartmentIndexs_extracellular) / this.stepSizeSec);
                
                cellDryMass = sum(this.mass.cellDry);
                lowerBounds(this.fbaReactionIndexs_metaboliteExternalExchange) = ...
                    max(lowerBounds(this.fbaReactionIndexs_metaboliteExternalExchange), ...
                    fbaReactionBounds(this.fbaReactionIndexs_metaboliteExternalExchange, 1) * cellDryMass);
                upperBounds(this.fbaReactionIndexs_metaboliteExternalExchange) = ...
                    min(upperBounds(this.fbaReactionIndexs_metaboliteExternalExchange), ...
                    fbaReactionBounds(this.fbaReactionIndexs_metaboliteExternalExchange, 2) * cellDryMass);
            end
            
            %internal metabolite availability
            if applyInternalMetaboliteBounds
                lowerBounds(this.fbaReactionIndexs_metaboliteInternalLimitedExchange) = max(...
                    lowerBounds(this.fbaReactionIndexs_metaboliteInternalLimitedExchange), ...
                    -substrates(this.substrateIndexs_internalExchangedLimitedMetabolites) / ...
                    this.stepSizeSec);
            end
            
            %protein monomers and complexes
            if applyProteinBounds
                limitedReactions = any(any(...
                    this.reactionStoichiometryMatrix([this.substrateMonomerLocalIndexs; this.substrateComplexLocalIndexs], this.reactionIndexs_fba, :) & ...
                    ~permute(repmat(this.proteinLimitableProteinComposition * substrates(this.substrateIndexs_limitableProteins, :), [1 1 numel(this.reactionIndexs_fba)]), [1 3 2]), ...
                    3), 1);
                lowerBounds(this.fbaReactionIndexs_metabolicConversion(limitedReactions)) = 0;
                upperBounds(this.fbaReactionIndexs_metabolicConversion(limitedReactions)) = 0;
            end
            
            %sanitize linear programming input: replace Infs with large numbers
            bounds = [lowerBounds upperBounds];
        end
    end
end