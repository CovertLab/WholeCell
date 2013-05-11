%Macromolecular complexation
%
% @wholeCellModelID Process_MacromolecularComplexation
% @name             Macromolecular complexation
% @description
%   Biology
%   ===============
%   An important step in the synthesis of functional enzymes is the
%   stochiometric formation of macromolecular complexes. Macromolecular
%   complexation is kinetically fast, and energetically favorable.
%   Consequently, the model assumes that complexation is limited only by
%   subunit availability, and proceeds to completion rapidly. In addition,
%   we assume that each complex forms with the same specific rate.
%
%   Knowledge Base
%   ===============
%   As of 8/17/2010 the M. genitalium knowledge base contained 155
%   macromolecular complexes involving 262 protein monomer species and 5 RNA
%   species, and 277 links between protein/RNA/complex subunit species and
%   complexes. Of these 155 complexes, detailed information is available on
%   the formation of 6, and these six complexes are formed in other processes.
%   The remaining 149 are formed by this process. The following table lists
%   several statistics about the macromolecular complexation network.
%
%      Statistic                                      Value
%      ===========================================    ===========
%      Mean no. subunit species per complex           1.8 +/- 3.2
%      Min, Max no. subunits species per complex      1 - 34
%      Mean no. subunits per complex                  5 +/- 15.6
%      Min, Max no. subunits per complex              2 - 192
%      Mean no. complexes per monomer species         1 +/ 0.1
%      Min, Max no. complexes per monomer species     1 - 2
%      No. subunits participating in >1 complexes     5
%      No. subunits which are themselves complexes    5
%      No. complexes formed in this process            149
%      No. complexes formed in other processes          6
%
%   Representation
%   ===============
%   Two properties are used to represent the counts of the 149 macromolecular
%   complexes formed by this process (complexs), and of free macromolecular
%   complex subunits (substrates).
%
%   Four properties are used to represent the structure of the macromolecular
%   complex - subunit network: complexComposition, complexNetworks,
%   substrates2complexNetworks, and complexs2complexNetworks. complexComposition
%   is an adjancency matrix between subunits and complexes populated from the
%   knowledge base; entries contain the number of subunits of each type in each
%   complex. complexNetworks is cell array containing the clustering of
%   complexComposition built by findNonInteractingRowsAndColumns called by
%   initializeConstants; each entry represents a disconnected part of the
%   macromolecular complex - subunit network, which can be simulated separate from
%   all other disconnected parts. substrates2complexNetworks and
%   complexs2complexNetworks represent mappings between substrates and complexs
%   and the disconnected parts of the macromolecular complex - subunit network
%   stored in complexNetworks.
%
%   Initialization
%   ===============
%   Macromolecular complexes are initialized up to the amounts of RNA and
%   protein subunits initialized by other processes.
%
%   Simulation
%   ===============
%   Macromolecular complexes are formed assuming:
%   1) Complexation is highly energetically favorable and
%   2) Complexation is fast, and thus
%   3) Macromolecular complexes are formed to completion; that is until
%      there are insufficient free monomers to form additional complexes.
%
%   Complexes are formed according to Monte Carlo simulation for
%   each independent protein complex network (previously established by
%   initializeConstants and stored in complexNetworks,
%   substrates2complexNetworks, and complexs2complexNetworks). First,
%   we use mass-action kinetics to compute the relative formation rate of
%   each complex. Specifically we compute the relative formation rate as
%   the product of the concentration of all monomers raised to the power of
%   their stoichiometries within the complex. Second we stochastically form
%   complexes according the computed formation rates. This is repeated until
%   no further complexes can form.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Markus Covert, mcovert@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/4/2010
classdef  MacromolecularComplexation < edu.stanford.covert.cell.sim.Process
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
            'complexComposition';
            'complexNetworks';
            'substrates2complexNetworks';
            'complexs2complexNetworks';
			};
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {   %names of simulation state properties redundant with timecourses in this or other processes or the simulation
            'complexs'};
    end

    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs   = {}; %whole cell model IDs of stimuli
        substrateWholeCellModelIDs = {}; %whole cell model IDs of subunits
        enzymeWholeCellModelIDs    = {}; %enzyme whole cell model ids
        complexWholeCellModelIDs   = {}; %whole cell model IDs of macromolecular complexes

        complexGlobalIndexs              %complex row indices within simulation.complexs
        complexCompartmentIndexs         %complex compartment indices within simulation.compartment.wholeCellModelIDs
        complexGlobalCompartmentIndexs   %complex row & column indices within simulation.complexs
    end

    %fixed biological constants
    properties
        complexComposition           %protein complex composition (monomers X complexes X compartments)
        complexNetworks              %decomposition of complexComposition into independent networks
        substrates2complexNetworks   %network index of each subunit
        complexs2complexNetworks     %network index of each protein complex, index of zero indicates not to model complexation using this process
    end

    %global state (stored locally for convenience)
    properties
        complexs                     %numbers of macromolecular complexes
    end

    %constructor
    methods
        function this = MacromolecularComplexation(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.Process(wholeCellModelID, name);
        end
    end

    %communication between process/simulation
    methods
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)
            %subunits
            this.substrateWholeCellModelIDs = cell(knowledgeBase.numGenes,1);
            this.substrateWholeCellModelIDs(this.gene.mRNAIndexs) = this.monomer.wholeCellModelIDs(this.monomer.matureIndexs);
            this.substrateWholeCellModelIDs(this.gene.rRNAIndexs) = this.rna.wholeCellModelIDs(this.rna.matureIndexs(this.rna.matureRRNAIndexs));
            this.substrateWholeCellModelIDs(this.gene.sRNAIndexs) = this.rna.wholeCellModelIDs(this.rna.matureIndexs(this.rna.matureSRNAIndexs));
            this.substrateWholeCellModelIDs(this.gene.tRNAIndexs) = this.rna.wholeCellModelIDs(this.rna.matureIndexs(this.rna.matureTRNAIndexs));

            %complexs
            process = findobj(knowledgeBase.processes, 'wholeCellModelID', this.wholeCellModelID);
            complexProcesses = [knowledgeBase.proteinComplexs.complexFormationProcess]';
            complexIdxs = find(complexProcesses == process);
            this.complexGlobalIndexs = this.monomer.nascentIndexs(complexIdxs);
            this.complexCompartmentIndexs = this.complex.compartments(this.complexGlobalIndexs);
            this.complexGlobalCompartmentIndexs = sub2ind(...
                [numel(this.complex.wholeCellModelIDs) numel(this.compartment.wholeCellModelIDs)], ...
                this.complexGlobalIndexs, ...
                this.complexCompartmentIndexs);

            this.complexWholeCellModelIDs = this.complex.wholeCellModelIDs(this.complexGlobalIndexs);

            this.complexComposition = sum(...
                knowledgeBase.proteinComplexRNAComposition(:, complexIdxs, :) + ...
                knowledgeBase.proteinComplexMonomerComposition(:, complexIdxs, :), 3);

            %1. Reshape monomer-complex network: [monomers-Compartment-1 x Complexs; monomers-Compartment-2 x Complexs; ... ; monomers-Compartment-n x Complexs]
            %   a. permute dimensions
            %   b. reshape
            %   c. transpose
            %2. Ignore complexes formed by other processes
            %3. Break monomer-complex network up into disjoint sets and build complexes for each set
            this.initializeComplexNetworks();
            
            %call super class method to compute substrate mapping onto
            %simulation and get substrate molecular weights
            this.initializeConstants@edu.stanford.covert.cell.sim.Process(knowledgeBase, simulation);
        end
        
        function initializeComplexNetworks(this)
            [substrates2complexNetworks, complexs2complexNetworks, complexNetworks] = ...
                edu.stanford.covert.util.findNonInteractingRowsAndColumns(this.complexComposition); %#ok<*PROP>
            
            this.substrates2complexNetworks = zeros(size(substrates2complexNetworks));
            this.complexs2complexNetworks = zeros(size(complexs2complexNetworks));
            this.complexNetworks = cell(1, 1);
            j = 1;
            for i = 1:length(complexNetworks)
                if size(complexNetworks{i}, 2) == 0
                elseif size(complexNetworks{i}, 2) == 1
                    this.substrates2complexNetworks(substrates2complexNetworks == i) = 1;
                    this.complexs2complexNetworks(complexs2complexNetworks == i) = 1;
                else
                    j = j + 1;
                    
                    this.substrates2complexNetworks(substrates2complexNetworks == i) = j;
                    this.complexs2complexNetworks(complexs2complexNetworks == i) = j;
                    this.complexNetworks{j, 1} = complexNetworks{i};
                end
            end
            this.complexNetworks{1} = this.complexComposition(this.substrates2complexNetworks == 1, this.complexs2complexNetworks == 1);
            
            this.substrateWholeCellModelIDs = this.substrateWholeCellModelIDs(this.substrates2complexNetworks > 0);
            
            this.complexWholeCellModelIDs = this.complexWholeCellModelIDs(this.complexs2complexNetworks > 0);
            this.complexGlobalIndexs = this.complexGlobalIndexs(this.complexs2complexNetworks > 0);
            this.complexCompartmentIndexs = this.complexCompartmentIndexs(this.complexs2complexNetworks > 0);
            this.complexGlobalCompartmentIndexs = this.complexGlobalCompartmentIndexs(this.complexs2complexNetworks > 0);
            
            this.complexComposition = this.complexComposition(this.substrates2complexNetworks>0, this.complexs2complexNetworks>0);
            this.substrates2complexNetworks = this.substrates2complexNetworks(this.substrates2complexNetworks > 0);
            this.complexs2complexNetworks = this.complexs2complexNetworks(this.complexs2complexNetworks > 0);
        end

        %retrieve state from simulation
        function copyFromState(this)
            this.copyFromState@edu.stanford.covert.cell.sim.Process();

            numTime = size(this.complexs, 3);

            if numTime == 1
                this.complexs = this.complex.counts(this.complexGlobalCompartmentIndexs);
            else
                this.complexs = permute(this.complex.counts(sub2ind(...
                    size(this.complex.counts), ...
                    repmat(this.complexGlobalIndexs, 1, numTime), ...
                    repmat(this.complexCompartmentIndexs, 1, numTime),...
                    repmat(1:numTime, length(this.complexGlobalIndexs), 1))),...
                    [1 3 2]);
            end
        end

        %send state to simulation
        function copyToState(this)
            this.copyToState@edu.stanford.covert.cell.sim.Process();

            %complexs
            numTime = size(this.complexs, 3);
            
            if numTime == 1
                this.complex.counts(this.complexGlobalCompartmentIndexs) = this.complexs;
            else
                this.complex.counts(sub2ind(...
                    size(this.complex.counts), ...
                    repmat(this.complexGlobalIndexs, 1, numTime), ...
                    repmat(this.complexCompartmentIndexs, 1, numTime),...
                    repmat(1:numTime, length(this.complexGlobalIndexs), 1))) = ...
                    permute(this.complexs, [1 3 2]);
            end
        end
    end

    %allocate memory for state
    methods
        function allocateMemoryForState(this, numTimePoints)
            this.allocateMemoryForState@edu.stanford.covert.cell.sim.Process(numTimePoints);

            this.complexs = zeros(length(this.complexWholeCellModelIDs), 1, numTimePoints);
        end
    end

    %model
    methods
        %Calculate
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, ~, ~)
            %no substrates or enzymes required for passive macromolecular
            %complexation
            bmProd = zeros(size(this.substrateWholeCellModelIDs));
            byProd = zeros(size(this.substrateWholeCellModelIDs));
            minEnzExp = zeros(size(this.enzymeWholeCellModelIDs));
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
        end
        
        %initialization
        %- monomers initialized to mature/bound/inactivated state
        %  by simulation initializeState method
        %- complexes initialized by
        %  1. Macromolecular complexation
        %  2. Ribosome assembly
        %  3. Protein folding
        %  4. Protein activation
        %
        %Here we use evolveState to form protein complexes
        function initializeState(this)
            this.evolveState();
        end

        %resource requirements
        function result = calcResourceRequirements_Current(this)
            result = zeros(size(this.substrates));
        end

        %simulation
        function evolveState(this)
            newComplexs = zeros(size(this.complexs));
            
            %subunits only involved in one complex (i.e. no competition)
            newComplexs(this.complexs2complexNetworks == 1) = buildProteinComplexs_bounds(...
                this.substrates(this.substrates2complexNetworks == 1, 1),...
                this.complexNetworks{1});
            
            %subunits involved in multiple complexes (i.e. competition): Run
            %Monte Carlo simulation for each independent protein complex network
            for i = 2:length(this.complexNetworks)
                newComplexs(this.complexs2complexNetworks == i) = ...
                    buildProteinComplexs_montecarlokinetic(...
                    this.substrates(this.substrates2complexNetworks == i), ...
                    this.complexNetworks{i}, this.randStream);
            end
            
            %stop if no new complexes
            if ~any(newComplexs)
                return;
            end
            
            this.complexs = this.complexs + newComplexs;
            this.substrates = this.substrates - this.complexComposition * newComplexs;
        end
    end
    
    %get methods of dependent local state
    methods
        function value = getDryWeight(this)
            if size(this.complexs, 3) == 1
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + ...
                    this.complex.molecularWeights(this.complexGlobalIndexs)' * this.complexs / ...
                    edu.stanford.covert.util.ConstantUtil.nAvogadro;
            else
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + ...
                    permute(...
                    this.complex.molecularWeights(this.complexGlobalIndexs)' * permute(this.complexs, [1 3 2]),...
                    [1 3 2]) / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            end
        end
    end
end

function proteinComplexs = buildProteinComplexs_montecarlokinetic(...
        totalProteinMonomers, proteinComplexMatrix, randStream)
    nComplexs = size(proteinComplexMatrix, 2);
    proteinComplexs = zeros(nComplexs, 1);

    %build complexes stochastically
    while true
        %compute rate of building each complex
        cumprob = buildProteinComplexs_rates_collisionTheory(...
            totalProteinMonomers, proteinComplexMatrix, 'cumulative probability');

        %stop when we can't make any more complexes
        if isnan(cumprob(1)); break; end;

        %select complex
        selectedComplex = find(randStream.rand() < cumprob, 1, 'first');
        if isempty(selectedComplex)
            selectedComplex = find(cumprob == 1, 1, 'first');
        end

        %mass balance
        proteinComplexs(selectedComplex) = proteinComplexs(selectedComplex) + 1;
        totalProteinMonomers = totalProteinMonomers - proteinComplexMatrix(:, selectedComplex);
    end
end

function rates = buildProteinComplexs_rates_collisionTheory(...
        totalProteinMonomers, proteinComplexMatrix, normalization)
    %compute rate of building each complex
    %Assumes that complexes form by all subunits simultaneously colliding.
    %However, protein complexes do not assemble this way. Rather subunits are
    %added sequentially. Thus this method will will underestimate the rate of
    %formation of larger protein complexes relative to smaller complexes.
    %Note: to get correct units subtract one from the sum(coefficients)
    %term. Its not necessary to get the units correct because this factor
    %is later normalized out.
    %TODO: replace mean(totalProteinMonomers) with better estimate of volume
    rates = prod((totalProteinMonomers(:, ones(size(proteinComplexMatrix,2), 1)) / mean(totalProteinMonomers)) .^ proteinComplexMatrix, 1)';
    
    %set rate to zero if upper zero
    ub = buildProteinComplexs_bounds(totalProteinMonomers, proteinComplexMatrix);
    rates(ub == 0) = 0;
    
    if nargin >= 3
        switch normalization(1)
            case 'c' %cumulative probability
                rates = cumsum(rates);
                rates = rates / rates(end);
            case 'p' %probability
                rates = rates / sum(rates);
            otherwise
                throw(MException('MacromolecularComplexation:error', 'Invalid normalization'));
        end
    end
end

function ub = buildProteinComplexs_bounds(totalProteinMonomers, proteinComplexMatrix)
    ub = floor(min(totalProteinMonomers(:, ones(1, size(proteinComplexMatrix, 2))) ./ proteinComplexMatrix, [], 1))';
end
