 %Protein Activation
%
% @wholeCellModelID Process_ProteinActivation
% @name             Protein Activation
% @description
%   Biology
%   ================================
%   The activity of proteins and other macromolecules can be modulated by other
%   molecules such as free metabolites (and pseudo metabolites we've represented
%   as stimuli) both at the enzymatically active site, and at more distant
%   allosteric sites. This regulation helps proteins, and thereby the cell,
%   respond to changes in the internal and external environments, and maintain
%   homeostasis. From a network perspective, this regulation can give rise
%   to positive and negative feedback loops.
%
%   This process models protein regulation by transitioning proteins between
%   enzymatically active (simulation properties matureMonomers and
%   matureComplexs) and inactive states (simulation properties inactiveMonomers
%   and inactiveComplexs) according to boolean regulatory rules built from the
%   primary literature. Regulatory rules are evaluated independently for each
%   compartment. Proteins for which we have not implemented a regulatory
%   rule are assumed to always remain in the enyzmatically active state.
%
%   Knowledge Base
%   ================================
%   As of 8/9/2010 the M. genitalium knowledge base includes regulatory rules
%   for six proteins:
%   - MG_085_HEXAMER  HPr(Ser) kinase/phosphatase
%   - MG_101_MONOMER  Uncharacterized HTH-type transcriptional regulator
%   - MG_127_MONOMER  Spx subfamily protein
%   - MG_205_DIMER    heat-inducible transcription repressor HrcA, putative
%   - MG_236_MONOMER  ferric uptake repressor
%   - MG_409_DIMER    phosphate transport system regulatory protein PhoU,
%                     putative
%
%   Boolean activation rule syntax
%   ++++++++++++++++++++++++++++++++
%   - Boolean activation rules must evaluate to true or false where true
%     indicates that the protein species is active, and false indicates
%     the protein species is inactive. That is protein species are
%     activated/inactivated in an all-or-nothing fashion for each
%     compartment.
%   - Using Whole Cell Model IDs boolean rules can contain references to
%     counts/concentrations of several kinds of objects:
%     - stimuli                    count                (simulation.stimuli)
%     - metabolites                concentration (mM)   (simulation.metabolites)
%     - mature protein monomers    concentration (mM)   (simulation.matureIndexs)
%     - mature protein complexes   concentration (mM)   (simulation.matureComplexIndexs)
%     That is, Whole Cell Model IDs of knowledge base objects contained
%     within boolean rules will replaced by the current value or
%     concentration in mM of the corresponding knowledge base object before
%     evaluation of the boolean rule
%   - Boolean activation rules permit the operators: |, &, +, -, !, <, >, <=, >=, ==
%   - Boolean activation rules permit parenthesis for grouping
%   - Boolean activation rules permit spaces
%
%   Examples
%   - objectID1 > val1
%   - objectID1 <= val1
%   - objectID1 == val1 | objectID2>conc2
%   - objectID1 == val1 | !(objectID2>conc2 & objectID3<=conc3 & (objectID4+object5)<conc5)
%
%   Representation
%   ================================
%   The substrates and inactivatedSubstrates properties represent the counts of
%   enzymatically active proteins in each compartment (all of the compartments
%   in the simulation are separately mapped to this process). stimuli represents
%   the  values of external perturbations and pseudo free metabolites
%   (simulation property stimuli), free metabolites, and macromolecules which
%   affect the activity of the proteins.
%
%   activationRules represents the boolean regulatory rules governing each of
%   the proteins in substrates. The activation rule for each protein is
%   evaluated independently for each compartment. Boolean rules are transcoded
%   to MATLAB syntax from the syntax described above during initializeConstants.
%   Boolean rules are evaluated using the MATLAB eval command in a workspace
%   where variables having the names of whole cell model ids of the stimuli are
%   defined and set to the value/concentration of the corresponding stimulus.
%
%   Initialization
%   ================================
%   The same boolean regulatory rules that are evaluated during the simulation
%   are evaluated during initialization. Proteins for which their regulatory
%   rule evaluates to false are initialized to the inactive state
%   (inactivatedSubstrates); proteins for which their regulatory rule evaluates to
%   true are initialized to the active state (substrates). This is achieved by
%   calling the evolveState method.
%
%   Simulation
%   ================================
%   For each compartment
%     1. Use scaleComponents to compute the concentrations of objects (except
%        stimuli objects) mapped to the process's stimuli property.
%     2. Assign to local variables named with the stimuli Whole Cell Model IDs
%        values equal that of the corresponding object (for stimuli objects) or
%        the concentration (for all other objects) of the corresponding object
%        computed in (1).
%     3. Evaluate boolean regulatory rule of each protein using the MATLAB eval
%        command.
%     4. Update substrates and inactivatedSubstrates properties based on (3).
%        Proteins whose boolean regulation rule evaluates to false are
%        transitioned from substrates to inactivatedSubstrates. Proteins whose
%        boolean regulation rule evaluates to true undergo the opposite
%        transition.
%   End
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/9/2010

%TODO: extend boolean regulatory rules to permite RNAs
%      - edit knowledgebase/schema.xml
%      - edit knowledgebase/classes/KnowledgeBaseObject.php
%      - edit simulation/src/.../KnowledgeBase.m
%      - edit simulation/src/.../ProteinMonomer.m
%      - edit simulation/src/.../ProteinComplex.m

classdef  ProteinActivation < edu.stanford.covert.cell.sim.Process
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of properties that are options
        fixedConstantNames__       = {   %names of fixed constant properties
			'activationRules'};
        fittedConstantNames__      = {}; %names of properties that are fitted constants        
        localStateNames__          = {   %names of properties that are simulation state owned by the simulation or other processes
            'inactivatedSubstrates'};
    end

    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs     = {}; %whole cell models of stimuli, set by initializeConstants from activation rules in knowledge base
        substrateWholeCellModelIDs   = {}; %whole cell models of substrates, set by initializeConstants from activation rules in knowledge base
        enzymeWholeCellModelIDs      = {}; %whole cell models of enzymes, not used
                
        inactivatedSubstrateMonomerGlobalCompartmentIndexs %indices of inactivated monomer substrates within this.monomer.counts
        inactivatedSubstrateComplexGlobalCompartmentIndexs %indices of inactivated complex substrates within simulation.complexs
    end

    %fixed biological constants
    properties
        activationRules                       %boolean rules governing the (in)activation of each substrate
    end

    %local state
    properties
        inactivatedSubstrates                 %counts of inactivated substrates; these have same whole cell models ids and molecular weights as substrates
    end

    %constructor
    methods
        function this = ProteinActivation(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.Process(wholeCellModelID, name);
        end
    end

    %communication between process/simulation
    methods
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)           
            %find regulated monomers, complexs
            monomers = findobj(knowledgeBase.proteinMonomers, '-not', 'activationRule', '');
            complexs = findobj(knowledgeBase.proteinComplexs, '-not', 'activationRule', '');

            %whole cell model ids of regulated proteins
            this.substrateWholeCellModelIDs = {
                monomers.wholeCellModelID ...
                complexs.wholeCellModelID}';

            %boolean regulatory rules of regulated proteins
            this.activationRules = this.transcodeActivationRules({
                monomers.activationRule ...
                complexs.activationRule}');
            
            %whole cell model ids of stimuli, metabolites, and proteins which
            %regulate the regulated proteins -- these links were created by parsing the
            %activation rules when they were saved to the knowledge base,
            %and loading into the MATLAB representation of the knowledge
            %base from the MySQL database representation
            regulatorIDs = cell(0, 1);
            for i = 1:numel(monomers)
                if ~isempty(monomers(i).stimuliRegulators)
                    regulatorIDs = [
                        regulatorIDs
                        {monomers(i).stimuliRegulators.wholeCellModelID}'
                        ]; %#ok<*AGROW>
                end
                if ~isempty(monomers(i).metaboliteRegulators)
                    regulatorIDs = [
                        regulatorIDs
                        {monomers(i).metaboliteRegulators.wholeCellModelID}'
                        ];
                end
                if ~isempty(monomers(i).proteinMonomerRegulators)
                    regulatorIDs = [
                        regulatorIDs
                        {monomers(i).proteinMonomerRegulators.wholeCellModelID}'
                        ];
                end
                if ~isempty(monomers(i).proteinComplexRegulators)
                    regulatorIDs = [
                        regulatorIDs
                        {monomers(i).proteinComplexRegulators.wholeCellModelID}'
                        ];
                end
            end
            for i = 1:numel(complexs)
                if ~isempty(complexs(i).stimuliRegulators)
                    regulatorIDs = [
                        regulatorIDs
                        {complexs(i).stimuliRegulators.wholeCellModelID}'
                        ];
                end
                if ~isempty(complexs(i).metaboliteRegulators)
                    regulatorIDs = [
                        regulatorIDs
                        {complexs(i).metaboliteRegulators.wholeCellModelID}'
                        ];
                end
                if ~isempty(complexs(i).proteinMonomerRegulators)
                    regulatorIDs = [
                        regulatorIDs
                        {complexs(i).proteinMonomerRegulators.wholeCellModelID}'
                        ];
                end
                if ~isempty(complexs(i).proteinComplexRegulators)
                    regulatorIDs = [
                        regulatorIDs
                        {complexs(i).proteinComplexRegulators.wholeCellModelID}'
                        ];
                end
            end
            this.stimuliWholeCellModelIDs = unique(regulatorIDs);

            %call super class method to construct mapping between process stimuli
            %and substrates and the simulation
            this.initializeConstants@edu.stanford.covert.cell.sim.Process(...
                knowledgeBase, simulation, struct(...
                    'retainStimuliCompartments', true,...
                    'retainSubstrateCompartments', true,...
                    'retainEnzymeCompartments', true));

            %inactivated substrate indices
            this.inactivatedSubstrateMonomerGlobalCompartmentIndexs = sub2ind(...
                [numel(this.monomer.wholeCellModelIDs) this.compartment.count], ...
                this.monomer.inactivatedIndexs(this.substrateMonomerGlobalIndexs), ...
                this.substrateMonomerCompartmentIndexs);
            this.inactivatedSubstrateComplexGlobalCompartmentIndexs = sub2ind(...
                [numel(this.complex.wholeCellModelIDs) this.compartment.count], ...
                this.complex.inactivatedIndexs(this.substrateComplexGlobalIndexs), ...
                this.substrateComplexCompartmentIndexs);
        end

        %translate activation rules into MATLAB expressions
        %- replace negation operator '!' with MATLAB negation operator '~'
        function value = transcodeActivationRules(~, value)
            for i = 1:length(value)
                value{i} = strrep(value{i},'!','~');
            end
        end
        
        %retrieve state from simulation
        function copyFromState(this)
            this.copyFromState@edu.stanford.covert.cell.sim.Process();
            
            numTimePoints = size(this.inactivatedSubstrates, 3);
            if numTimePoints == 1
                this.inactivatedSubstrates(this.substrateMonomerLocalIndexs, :) = ...
                    this.monomer.counts(this.inactivatedSubstrateMonomerGlobalCompartmentIndexs);
                
                this.inactivatedSubstrates(this.substrateComplexLocalIndexs, :) = ...
                    this.complex.counts(this.inactivatedSubstrateComplexGlobalCompartmentIndexs);
            else
                this.inactivatedSubstrates(this.substrateMonomerLocalIndexs, :, :) = this.monomer.counts(...
                    sub2ind(size(this.monomer.counts), ...
                    repmat(this.monomer.inactivatedIndexs(this.substrateMonomerGlobalIndexs), [1 1 numTimePoints]), ...
                    repmat(this.substrateMonomerCompartmentIndexs, [1 1 numTimePoints]), ...
                    permute(repmat((1:numTimePoints)',[1 size(this.substrateMonomerGlobalIndexs)]),[2 3 1])));
                
                this.inactivatedSubstrates(this.substrateComplexLocalIndexs, :, :) = this.complex.counts(...
                    sub2ind(size(this.complex.counts), ...
                    repmat(this.complex.inactivatedIndexs(this.substrateComplexGlobalIndexs), [1 1 numTimePoints]), ...
                    repmat(this.substrateComplexCompartmentIndexs, [1 1 numTimePoints]), ...
                    permute(repmat((1:numTimePoints)',[1 size(this.substrateComplexGlobalIndexs)]),[2 3 1])));
            end
        end
        
        %send state to simulation
        function copyToState(this)
            this.copyToState@edu.stanford.covert.cell.sim.Process();                     
            
            numTimePoints = size(this.inactivatedSubstrates, 3);
            if numTimePoints == 1
                this.monomer.counts(this.inactivatedSubstrateMonomerGlobalCompartmentIndexs) = ...
                    this.inactivatedSubstrates(this.substrateMonomerLocalIndexs, :);
                
                this.complex.counts(this.inactivatedSubstrateComplexGlobalCompartmentIndexs) = ...
                    this.inactivatedSubstrates(this.substrateComplexLocalIndexs, :);
            else
                this.monomer.counts(...
                    sub2ind(size(this.monomer.counts), ...
                    repmat(this.monomer.inactivatedIndexs(this.substrateMonomerGlobalIndexs), [1 1 numTimePoints]), ...
                    repmat(this.substrateMonomerCompartmentIndexs, [1 1 numTimePoints]), ...
                    permute(repmat((1:numTimePoints)',[1 size(this.substrateMonomerGlobalIndexs)]),[2 3 1]))) = ...
                    this.inactivatedSubstrates(this.substrateMonomerLocalIndexs, :, :);
                
                this.complex.counts(...
                    sub2ind(size(this.complex.counts), ...
                    repmat(this.complex.inactivatedIndexs(this.substrateComplexGlobalIndexs), [1 1 numTimePoints]), ...
                    repmat(this.substrateComplexCompartmentIndexs, [1 1 numTimePoints]), ...
                    permute(repmat((1:numTimePoints)',[1 size(this.substrateComplexGlobalIndexs)]),[2 3 1]))) = ...
                    this.inactivatedSubstrates(this.substrateComplexLocalIndexs, :, :);
            end
        end
    end

    %allocate memory for state
    methods
        function allocateMemoryForState(this, numTimePoints)
            this.allocateMemoryForState@edu.stanford.covert.cell.sim.Process(numTimePoints);

            this.inactivatedSubstrates = zeros(size(this.substrates));
        end
    end

    %model
    methods
        %Calculate
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, ~, ~)
            %no metabolites or enzymes needed
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
        %Here we use evolveState to evaluate and apply all protein activation
        %rules to both monomers and complexes
        function initializeState(this)
            this.evolveState();
        end

        %resource requirements
        function result = calcResourceRequirements_Current(this)
            result = zeros(size(this.substrates));
        end

        %simulation
        function evolveState(this)
            i = this.evaluateActivationRules();
            this.substrates(i) = this.substrates(i) + this.inactivatedSubstrates(i);
            this.inactivatedSubstrates(i) = 0;

            j = ~i;
            this.inactivatedSubstrates(j) = this.inactivatedSubstrates(j) + this.substrates(j);
            this.substrates(j) = 0;
        end

        function activatedSubstrates = evaluateActivationRules(this)            
            %convert units of stimuli to mM
            stimuli = this.stimuli(:, this.compartment.extracellularIndexs) / ...
                (edu.stanford.covert.util.ConstantUtil.nAvogadro/1000) / ...
                this.geometry.volume;
            stimuli(this.stimulusStimulusLocalIndexs) = ...
                this.stimuli(this.stimulusStimulusLocalIndexs, this.compartment.extracellularIndexs);

            %evaluate regulatory rules
            activatedSubstrates = false(size(this.substrates));
            for j = 1:size(this.substrates, 2)
                if ~any(this.substrates(:, j)) && ~any(this.inactivatedSubstrates(:, j))
                    continue;
                end
                
                for i = 1:size(this.stimuli, 1)
                    eval([this.stimuliWholeCellModelIDs{i} '=' num2str(stimuli(i)) ';']);
                end
                
                for i = 1:size(this.substrates, 1)
                    if ~(this.substrates(i, j) || this.inactivatedSubstrates(i, j))
                        continue;
                    end
                    activatedSubstrates(i, j) = logical(eval(this.activationRules{i}));
                end
            end
        end
    end

    %get methods of dependent local state
    methods
        function value = getDryWeight(this)
            if size(this.inactivatedSubstrates, 3) == 1
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + ...
                    this.substrateMolecularWeights' * sum(this.inactivatedSubstrates, 2) ...
                     / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            else
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + ...
                    permute(this.substrateMolecularWeights' * permute(sum(this.inactivatedSubstrates, 2), [1 3 2]), [1 3 2]) ...
                     / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            end
        end
    end
end
