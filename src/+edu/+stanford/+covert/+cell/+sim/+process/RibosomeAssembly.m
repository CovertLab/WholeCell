%Assemble Ribosomes
%
% @wholeCellModelID Process_RibosomeAssembly
% @name             Ribosomal Assembly
% @description
%   Biology
%   ===============
%   Assembly of the 30S and 50S ribosomal particles is a special case of
%   macromolecular complexation in which subunits are incorporated in a
%   stereotyped pattern [PUB_0660, PUB_0661], and which requires six
%   energy-dependent GTPases – EngA, EngB, Era, Obg, RbfA, and RbgA. The
%   exact energy requirement of each GTPase is unknown. Here we assume that
%   each GTPase that has been reported to be required to form each
%   ribosomal particle requires 1 GTP per particle. 
%
%   Knowledge Base
%   ===============
%   The 30S and 50S ribosomal particles are represented as complexes in the
%   knowledge base. Their RNA and protein subunit composition was curated from
%   the literature, and is stored in association with the complexes in the
%   knowledge base. This composition is loaded by initializeConstants into this
%   process's proteinComplexRNAComposition and proteinComplexMonomerComposition
%   properties.
%
%   Representation
%   ===============
%   substrates, enzymes, RNAs, monomers, and complexs represent the counts of
%   free metabolites (eg. ATP, ADP, Pi, H2O, H+), the ribosomal assembly
%   GTPases, the amounts of RNA and protein monomer ribosomal subunits, and the
%   30S and 50S ribosomal particles. The process doesn't represent any
%   intermediate state of ribosomal particle assembly; ribosomal particle
%   assembly is assumed to be an all-or-nothing process on the time scale of
%   this process. That is, we make the simplifying assumption that either a
%   ribosomal particle completely forms within a single time step, or no
%   progress in assembly of that particle is made during that time step.
%
%   proteinComplexRNAComposition and proteinComplexMonomerComposition represent
%   the RNA and protein monomer composition of the 30S and 50S ribosomal
%   subunits. complexationCatalysisMatrix represents the GTPases required to
%   form each ribosomal particle.
%
%   Initialization
%   ===============
%   The process is initialized to a state with the maximal number of formed
%   ribosomal particles given the amounts of initialized RNA and protein
%   monomers. That is, the process is initialized to a state where insufficient
%   RNA and protein monomers are available to form additional ribosomal subunits.
%
%   Simulation
%   ===============
%   In a randomized order over particles, for each ribosomal particle:
%   1. Calculate the maximum number of particles that can form based on
%      available RNA and protein monomer subunits, GTPases, and GTP.
%   2. Increment the number of ribosomal particles. Decrement the numbers of RNA
%      and protein monomer subunits, and GTP and water. Increment the counts of
%      the byproducts of GTP hydrolysis (GDP, Pi, H).
%
%   This makes the simplying assumption that ribosomal assembly is fast compared
%   to the 1s time scale of this process, and is energetically favorable such
%   that in the presence of saturating enzymes and energy, ribosomal assembly is
%   limited by the subunit availability.
%
%   References
%   ===============
%   1. Nierhaus KH (1991). The assembly of prokaryotic ribosomes.
%      Biochimie. 76(3):739-55.[PUB_0660]
%   2. Culver GM (2003). Assembly of the 30S ribosomal subunit.
%      Biopolymers. 68(2):234-49. [PUB_0661]
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/4/2010
classdef RibosomeAssembly < edu.stanford.covert.cell.sim.Process
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
            'complexationCatalysisMatrix'
            'proteinComplexRNAComposition'
            'proteinComplexMonomerComposition'
            };
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {   %names of simulation state properties redundant with timecourses in this or other processes or the simulation
            'RNAs'
            'monomers'
            'complexs'
            };
    end
    
    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs = {}; %whole cell model IDs of stimuli

        substrateWholeCellModelIDs = { %whole cell model IDs of substrates
            'GTP';'GDP';'PI';'H2O';'H'};
        substrateIndexs_gtp       = 1; %index within substrates of GTP
        substrateIndexs_gdp       = 2; %index within substrates of GDP
        substrateIndexs_phosphate = 3; %index within substrates of inorganic phosphate
        substrateIndexs_water     = 4; %index within substrates of water
        substrateIndexs_hydrogen  = 5; %index within substrates of hydrogen

        enzymeWholeCellModelIDs = {    %whole cell model IDs of enzymes
            'MG_329_MONOMER';          %GTP-binding protein engA
            'MG_335_MONOMER';          %GTP-binding protein engB, putative
            'MG_387_MONOMER';          %GTP-binding protein Era
            'MG_384_MONOMER';          %GTPase1 Obg
            'MG_143_MONOMER';          %ribosome-binding factor A
            'MG_442_MONOMER'};         %ribosomal biogenesis GTPase
        enzymeIndexs_engA = 1;         %index within enzymes of EngA
        enzymeIndexs_engB = 2;         %index within enzymes of EngB
        enzymeIndexs_era  = 3;         %index within enzymes of Era
        enzymeIndexs_obg  = 4;         %index within enzymes of Obg
        enzymeIndexs_rbfA = 5;         %index within enzymes of RbfA
        enzymeIndexs_rbgA = 6;         %index within enzymes of RbgA
        enzymeIndexs_30S_assembly_gtpase = [3;5];     %index within enzymes of 30 S ribosome assembly GTPases
        enzymeIndexs_50S_assembly_gtpase = [1;2;4;6]; %index within enzymes of 50 S ribosome assembly GTPases

        rnaWholeCellModelIDs           %whole cell model IDs of ribosome RNA subunits
        monomerWholeCellModelIDs       %whole cell model IDs of ribosome protein monomer subunits
        complexWholeCellModelIDs = {   %whole cell model IDs of ribosomal particles
            'RIBOSOME_30S';            %30S ribosomal particle
            'RIBOSOME_50S'};           %50S ribosomal particle
        complexIndexs_30S_ribosome = 1;%index within complexs of 30S ribosomal particle
        complexIndexs_50S_ribosome = 2;%index within complexs of 50S ribosomal particle

        rnaGlobalIndexs                %indices within simulation.matureRNAIndexs of RNAs
        monomerGlobalIndexs            %indices within simulation.matureMonomerIndexs of protein monomers
        complexGlobalIndexs            %indices within simulation.matureComplexIndexs of ribosomal particles
    end

    %fixed biological constants
    properties
        complexationCatalysisMatrix      %assembly GTPases required to assemble each ribosomal particle [GTPases X ribosomal particles]
        proteinComplexRNAComposition     %RNA subunit composition of ribosomal particles [RNAs X ribosomal particles]
        proteinComplexMonomerComposition %protein monomer subunit composition of ribosomal particles [protein monomers X ribosomal particles]
    end

    %global state (stored locally for convenience)
    properties
        RNAs     %counts of ribosome RNA subunits
        monomers %counts of ribosome protein monomer subunits
        complexs %counts of ribosomal particles
    end

    %constructor
    methods
        function this = RibosomeAssembly(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.Process(wholeCellModelID, name);
        end
    end

    %communication between process/simulation
    methods
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)
            this.initializeConstants@edu.stanford.covert.cell.sim.Process(knowledgeBase, simulation, varargin{:});

            %complexes
            [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, this.complexGlobalIndexs] = ...
                this.initializeConstantsHelper(this.complexWholeCellModelIDs);

            %RNAs
            rnaComposition = knowledgeBase.proteinComplexRNAComposition(:, :, this.compartment.cytosolIndexs);

            this.rnaGlobalIndexs = ...
                find(any(rnaComposition(this.gene.rRNAIndexs, this.complexGlobalIndexs),2));

            this.rnaWholeCellModelIDs = this.rna.wholeCellModelIDs(...
                this.rna.matureIndexs(this.rna.matureRRNAIndexs(this.rnaGlobalIndexs)));

            this.proteinComplexRNAComposition = rnaComposition(...
                this.gene.rRNAIndexs(this.rnaGlobalIndexs), this.complexGlobalIndexs);
            assertEqual([0; 1], unique(this.proteinComplexRNAComposition), 'Stoichiometry of subunits in complexes must be 1');
            this.proteinComplexRNAComposition = this.proteinComplexRNAComposition == 1;

            %monomers
            monomerComposition = knowledgeBase.proteinComplexMonomerComposition(:, :, this.compartment.cytosolIndexs);

            this.monomerGlobalIndexs = ...
                find(any(monomerComposition(this.gene.mRNAIndexs, this.complexGlobalIndexs),2));

            this.monomerWholeCellModelIDs = this.monomer.wholeCellModelIDs(...
                this.monomer.matureIndexs(this.monomerGlobalIndexs));

            this.proteinComplexMonomerComposition = monomerComposition(...
                this.gene.mRNAIndexs(this.monomerGlobalIndexs), this.complexGlobalIndexs);
            assertEqual([0; 1], unique(this.proteinComplexMonomerComposition), 'Stoichiometry of subunits in complexes must be 1');
            this.proteinComplexMonomerComposition = this.proteinComplexMonomerComposition == 1;                      

            %GTPases
            this.complexationCatalysisMatrix = false(length(this.enzymeWholeCellModelIDs), length(this.complexWholeCellModelIDs));
            this.complexationCatalysisMatrix(this.enzymeIndexs_30S_assembly_gtpase, this.complexIndexs_30S_ribosome) = true;
            this.complexationCatalysisMatrix(this.enzymeIndexs_50S_assembly_gtpase, this.complexIndexs_50S_ribosome) = true;
        end

        %retrieve state from simulation
        function copyFromState(this)
            this.copyFromState@edu.stanford.covert.cell.sim.Process();

            this.RNAs = this.rna.counts(this.rna.matureIndexs(this.rna.matureRRNAIndexs(this.rnaGlobalIndexs)), this.compartment.cytosolIndexs, :);
            this.monomers = this.monomer.counts(this.monomer.matureIndexs(this.monomerGlobalIndexs),  this.compartment.cytosolIndexs, :);
            this.complexs = this.complex.counts(this.complex.nascentIndexs(this.complexGlobalIndexs), this.compartment.cytosolIndexs, :);
        end

        %send state to simulation
        function copyToState(this)
            this.copyToState@edu.stanford.covert.cell.sim.Process();

            this.rna.counts(this.rna.matureIndexs(this.rna.matureRRNAIndexs(this.rnaGlobalIndexs)), this.compartment.cytosolIndexs, :) = this.RNAs;
            this.monomer.counts(this.monomer.matureIndexs(this.monomerGlobalIndexs),  this.compartment.cytosolIndexs, :) = this.monomers;
            this.complex.counts(this.complex.nascentIndexs(this.complexGlobalIndexs), this.compartment.cytosolIndexs, :) = this.complexs;
        end
    end

    %allocate memory for state
    methods
        function allocateMemoryForState(this, numTimePoints)
            this.allocateMemoryForState@edu.stanford.covert.cell.sim.Process(numTimePoints);

            this.RNAs     = zeros(length(this.rnaGlobalIndexs),     1, numTimePoints);
            this.monomers = zeros(length(this.monomerGlobalIndexs), 1, numTimePoints);
            this.complexs = zeros(length(this.complexGlobalIndexs), 1, numTimePoints);
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
            
            %% substrate and byproducts: GTP hydrolysis by ribosome assembly GTPases
            gtpCost = 0;
            for i = 1:numel(this.complexWholeCellModelIDs)
                gtpCost = gtpCost + this.getGtpPerComplex(i) * min([
                    states.rnaProductions(this.rna.matureRRNAIndexs(this.rnaGlobalIndexs)) ./ this.proteinComplexRNAComposition(:, i);
                    states.monomerProductions(this.monomerGlobalIndexs)                    ./ this.proteinComplexMonomerComposition(:, i)]);
            end
            
            bmProd(this.substrateIndexs_gtp)       = gtpCost;
            bmProd(this.substrateIndexs_water)     = gtpCost;
            byProd(this.substrateIndexs_gdp)       = gtpCost;
            byProd(this.substrateIndexs_phosphate) = gtpCost;
            byProd(this.substrateIndexs_hydrogen)  = gtpCost;      
            
            %% enzymes
            minEnzExp([
                this.enzymeIndexs_30S_assembly_gtpase;
                this.enzymeIndexs_50S_assembly_gtpase]) = 2;
        end
        
        %initialization: complexs initialized by
        %1. Macromolecular complexation
        %2. Ribosome assembly
        %3. Protein folding
        %4. Protein activation
        %Here we push monomers to the ribosome-complexed state
        function initializeState(this)
            nCplxs = zeros(size(this.complexWholeCellModelIDs));
            for i = 1:numel(this.complexWholeCellModelIDs)
                nCplxs(i) = floor(min([
                    this.RNAs     ./ this.proteinComplexRNAComposition(:, i);
                    this.monomers ./ this.proteinComplexMonomerComposition(:, i)]));
            end
            this.RNAs     = this.RNAs     - this.proteinComplexRNAComposition     * nCplxs;
            this.monomers = this.monomers - this.proteinComplexMonomerComposition * nCplxs;
            this.complexs = this.complexs + nCplxs;
        end

        %resource requirements
        function result = calcResourceRequirements_Current(this)
            result = zeros(size(this.substrates));
            
            %assumes:
            %- protein produced slowly
            %- each RNA/protein used in only 1 complex
            %- stoichiometry of each RNA/protein in each complex is 1
            %- 1 GTP used per assembly GTPase per complex
            %
            %these assumptions are verified in the unit test and
            %initializeConstants
            result(this.substrateIndexs_gtp) = ...
                + this.getGtpPerComplex(2) * all(this.enzymes(this.complexationCatalysisMatrix(:, 1))) * (1 + min([
                    this.RNAs(this.proteinComplexRNAComposition(:, 1))
                    this.monomers(this.proteinComplexMonomerComposition(:, 1))
                    ])) ...
                + this.getGtpPerComplex(2) * all(this.enzymes(this.complexationCatalysisMatrix(:, 2))) * (1 + min([
                    this.RNAs(this.proteinComplexRNAComposition(:, 2))
                    this.monomers(this.proteinComplexMonomerComposition(:, 2))
                    ]));
            result(this.substrateIndexs_water) = result(this.substrateIndexs_gtp);
        end
        
        %simulation
        function evolveState(this)
            %stop if no GTP
            if ~this.substrates(this.substrateIndexs_gtp)
                return;
            end
            
            %randomize order in which to try forming complexes
            randOrder = this.randStream.randperm(numel(this.complexWholeCellModelIDs));
            for j = 1:numel(this.complexWholeCellModelIDs)
                i = randOrder(j);
                
                %number of GTPs needed for GTPases to form complex
                gtpPerComplex = this.getGtpPerComplex(i);
                
                %number of complexes that can form, limited by:
                %- enzymes
                %- GTP
                %- RNA subunits
                %- protein monomer subunits
                %
                %assumes stoichiometry of each subunit is 1; this is verified in
                %initializeConstants
                newComplexs = floor(min([
                    this.substrates(this.substrateIndexs_gtp) / gtpPerComplex;
                    this.substrates(this.substrateIndexs_water) / gtpPerComplex;
                    this.RNAs(this.proteinComplexRNAComposition(:,i));
                    this.monomers(this.proteinComplexMonomerComposition(:,i))]));
                if newComplexs == 0 || ~all(this.enzymes(this.complexationCatalysisMatrix(:,i)))
                    continue;
                end
                
                %update counts of complexes, RNAs, protein monomers, metabolites
                this.complexs(i) = this.complexs(i) + newComplexs;
                
                this.RNAs     = this.RNAs     - this.proteinComplexRNAComposition(:,i)     * newComplexs;
                this.monomers = this.monomers - this.proteinComplexMonomerComposition(:,i) * newComplexs;
                
                this.substrates = this.substrates + [-1; 1; 1; -1; 1] * newComplexs * gtpPerComplex; %substrate order is [GTP, GDP, Pi, H2O, H]
            end
        end
    end

    %get methods of dependent local state
    methods
        function result = getDryWeight(this)
            rnaMolecularWeights = this.rna.molecularWeights(...
                this.rna.matureIndexs(this.rna.matureRRNAIndexs(this.rnaGlobalIndexs))) / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            monomerMolecularWeights = this.monomer.molecularWeights(...
                this.monomer.matureIndexs(this.monomerGlobalIndexs)) / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            complexMolecularWeights = this.complex.molecularWeights(...
                this.complex.nascentIndexs(this.complexGlobalIndexs)) / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            
            if size(this.RNAs, 3) == 1
                result = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + ...
                    rnaMolecularWeights'     * this.RNAs     + ...
                    monomerMolecularWeights' * this.monomers + ...
                    complexMolecularWeights' * this.complexs;
            else
                result = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + ...
                    permute(rnaMolecularWeights'     * permute(this.RNAs,    [1 3 2]),[1 3 2]) + ...
                    permute(monomerMolecularWeights' * permute(this.monomers,[1 3 2]),[1 3 2]) + ...
                    permute(complexMolecularWeights' * permute(this.complexs,[1 3 2]),[1 3 2]);
            end
        end

        function result = getGtpPerComplex(this, i)
            result = sum(this.complexationCatalysisMatrix(:,i));
        end
    end
end
