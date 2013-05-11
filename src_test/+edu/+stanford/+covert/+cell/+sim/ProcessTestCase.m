%ProcessTestCase
% Base class for whole cell process test classes. This purposes of this
% class are to provide functions for
% - loading the process test fixture prior to the execution of each test.
%   This is automatically called by the setUp function.
% - testing for run-time errors in these process methods:
%   - allocateMemoryForState
%   - computeBiomassProductionByproducts
%   - computeMinimumEnzymeExpression
%   - resourceRequirement
%   - initializeState
%   - evolveState
% - testing that the process obeys mass balance
% - testing gene essentiality with regard to the process's primary function
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 9/10/2010
classdef ProcessTestCase < TestCase
    properties
        process               %process object
        processClass          %full class name of process
        processClassName      %short class name of process
        processBaseLocation   %base directory of process
    end
    
    %constructor
    methods
        function this = ProcessTestCase(name)
            %parent class constructor
            this = this@TestCase(name);
            
            %get class name of process
            this.processClass     = edu.stanford.covert.util.substr(class(this), 1, -5);
            this.processClassName = this.processClass(find(this.processClass == '.', 1, 'last') + 1:end);
            
            %find parent folder of process -- used to store test fixture and
            %expected outputs
            locationParts = strfind(this.Location, filesep);
            if this.Location(locationParts(end - 1) + 1) == '@'
                this.processBaseLocation = this.Location(1:locationParts(end-1));
            else
                this.processBaseLocation = this.Location(1:locationParts(end));
            end
        end
    end
    
    %setUp and tearDown
    methods
        %Called prior to the execution of each test case.
        %Performs three functions:
        % - instantiate process
        % - initialize process using test fixture
        % - seed random number generator
        function setUp(this)
            constructor = str2func(this.processClass);
            this.process = constructor([], this.processClassName);
            this.loadTestFixture();
            this.process.verbosity = 0;
            this.process.seed = 1;
            this.process.seedRandStream();
        end
        
        %Called following the execution of each test case.
        function tearDown(this)
            this.process = []; %free memory
        end
        
        %Loads knowledge base test fixture into process. Called by setUp prior
        %to the execution of each test. You may wish to override this method to
        %set additional properties.
        function loadTestFixture(this)
            this.process = edu.stanford.covert.cell.sim.ProcessFixture.load(this.process);
        end
    end
    
    %Suggested tests for subclasses to implement:
    % - Demonstrate main purpose of process using minimal, hard-coded fake data.
    % - Withhold one or all stimuli.
    % - Withhold one or all substrates.
    % - Withhold one or all enzymes.
    % - Predict which genes are essential to the process (testGeneEssentiality).
    
    %tests for absence of run-time errors
    methods (Sealed)
        function testAllocateMemoryForState(this)
            this.process.allocateMemoryForState(1);  %no run-time errors
        end
        
        function testCalcResourceRequirements_LifeCycle(this)
            %references
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], true);
            m = sim.process(this.process.wholeCellModelID(9:end));                        
            g = sim.gene;
            t = sim.state('Time');
            r = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            %process local state
            m.stimuli = this.process.stimuli;
            m.substrates = this.process.substrates;
            m.enzymes = this.process.enzymes;
            m.boundEnzymes = this.process.boundEnzymes;
            this.process = m;
            
            %constants
            constants = edu.stanford.covert.util.StructUtil.catstruct(...
                sim.getFixedConstants(), sim.getFittedConstants());
            constants.processes.DNADamage.substrateGlobalIndexs = sim.process('DNADamage').substrateGlobalIndexs;
            constants.processes.DNADamage.substrateStimulusGlobalIndexs = sim.process('DNADamage').substrateStimulusGlobalIndexs;
            constants.processes.DNADamage.substrateStimulusCompartmentIndexs = sim.process('DNADamage').substrateStimulusCompartmentIndexs;
            constants.processes.DNADamage.substrateStimulusGlobalCompartmentIndexs = sim.process('DNADamage').substrateStimulusGlobalCompartmentIndexs;
            constants.processes.DNADamage.substrateStimulusLocalIndexs = sim.process('DNADamage').substrateStimulusLocalIndexs;
            
            %RNA and protein
            rnas = sum( r.counts( r.matureIndexs, :), 2) + ...
                sum(r.counts(r.boundIndexs, :), 2) + ...
                sum(r.counts(r.aminoacylatedIndexs, :), 2) + ...
                sum(r.counts(r.damagedIndexs, :), 2);
            mons = ...
                sum(pm.counts(pm.matureIndexs, :), 2) + ...
                sum(pm.counts(pm.boundIndexs, :), 2) + ...
                sum(pm.counts(pm.damagedIndexs, :), 2) + ...
                sum(pm.counts(pm.inactivatedIndexs, :), 2);
            cpxs = sum(pc.counts(pc.matureIndexs, :), 2) + ...
                sum(pc.counts(pc.boundIndexs, :), 2) + ...
                sum(pc.counts(pc.damagedIndexs, :), 2) + ...
                sum(pc.counts(pc.inactivatedIndexs, :), 2);
            
            rnaInCpxs = zeros(size(rnas));
            rnaInCpxs(setdiff(1:end, r.matureMRNAIndexs)) = ...
                sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * cpxs;
            monInCpxs = sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * cpxs;
            
            cpxDecays = cpxs .* pc.decayRates(pc.matureIndexs) * t.cellCycleLength / log(2);
            rnaDecays = rnas .*  r.decayRates( r.matureIndexs) * t.cellCycleLength / log(2);
            rnaDecays(setdiff(1:end, r.matureMRNAIndexs)) = ...
                rnaDecays(setdiff(1:end, r.matureMRNAIndexs)) + ...
                sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * cpxDecays;
            monDecays = mons .* pm.decayRates(pm.matureIndexs) * t.cellCycleLength / log(2) + ...
                sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * cpxDecays;
            
            states = struct;
            
            states.rnas0               = rnas;
            states.rnas                = rnas * t.cellCycleLength / log(2);
            states.rnaDecays           = rnaDecays;
            states.rnaDecays0          = rnaDecays * log(2) / t.cellCycleLength;
            states.rnaProductions      = rnas + rnaInCpxs + rnaDecays;
            states.rnaProductions0     = (rnas + rnaInCpxs) * t.cellCycleLength / log(2) + (rnas .* r.decayRates(r.matureIndexs));
            states.rnaProductions0(setdiff(1:end, r.matureMRNAIndexs)) = ...
                states.rnaProductions0(setdiff(1:end, r.matureMRNAIndexs)) + ...
                sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * (cpxs .* pc.decayRates(pc.matureIndexs));
            
            states.monomers0           = mons;
            states.monomers            = mons * t.cellCycleLength / log(2);
            states.monomerDecays       = monDecays;
            states.monomerDecays0      = monDecays * log(2) / t.cellCycleLength;
            states.monomerProductions  = mons + monInCpxs + monDecays;
            states.monomerProductions0 = (mons + monInCpxs) * t.cellCycleLength / log(2) + (mons .* pm.decayRates(pm.matureIndexs));
            states.monomerProductions0 = ...
                states.monomerProductions0 + ...
                sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * (cpxs .* pc.decayRates(pc.matureIndexs));
            
            states.complexs0           = cpxs;
            states.complexs            = cpxs * t.cellCycleLength / log(2);
            states.complexDecays       = cpxDecays;
            states.complexDecays0      = cpxDecays * log(2) / t.cellCycleLength;
            states.complexProductions  = cpxs + cpxDecays;
            states.complexProductions0 = cpxs * t.cellCycleLength / log(2) + (cpxs .* pc.decayRates(pc.matureIndexs));            
            
            %calculate life cycle resource requirements
            [bmProd, byProd, minEnzExp, maxEnzExp] = m.calcResourceRequirements_LifeCycle(constants, states);
            
            %assertions
            assertEqual(size(m.substrateWholeCellModelIDs), size(bmProd));
            assertEqual(size(m.substrateWholeCellModelIDs), size(byProd));
            assertEqual(size(m.enzymeWholeCellModelIDs), size(minEnzExp));
            assertEqual(size(m.enzymeWholeCellModelIDs), size(maxEnzExp));
        end
        
        function testCalcResourceRequirements_Current(this)
            this.process.calcResourceRequirements_Current();  %no run-time errors
        end
        
        function testInitializeState(this)
            this.process.initializeState();  %no run-time errors
        end
        
        function testEvolveState(this)
            m = this.process;
            m.evolveState();  %no run-time errors
            validateattributes(m.substrates(setdiff(1:end, m.substrateStimulusLocalIndexs), :), {'numeric'}, {'real', 'nonnegative', 'integer'});
            validateattributes(m.enzymes, {'numeric'}, {'real', 'nonnegative', 'integer'});
            validateattributes(m.boundEnzymes, {'numeric'}, {'real', 'nonnegative', 'integer'});
        end
    end
    
    %other tests provided by this class
    methods (Sealed)
        %Test that process obeys mass balance
        function testMassBalance(this)
            import edu.stanford.covert.util.ConstantUtil;
            m = this.process;
            initialWeight = sum(m.dryWeight);
            
            m.evolveState();
            assertElementsAlmostEqual(...
                initialWeight * ConstantUtil.nAvogadro, sum(m.dryWeight) * ConstantUtil.nAvogadro, ...
                'relative', 1e-8);
        end
        
        %Test that process uses rand streams correctly
        %- only changes state of process (and possibly state) rand
        %  stream(s)
        %- doesn't change state of global rand stream
        function testRandStreamUsage(this)
            %% hold onto state of default stream
            defaultStream = RandStream.getDefaultStream;
            defaultStreamType = defaultStream.Type;
            defaultStreamSeed = defaultStream.Seed;
            defaultStreamNumStreams = defaultStream.NumStreams;
            defaultStreamStreamIndex = defaultStream.StreamIndex;
            defaultStreamState = defaultStream.State;
            defaultStreamSubstram = defaultStream.Substream;
            defaultStreamRandnAlg = defaultStream.RandnAlg;
            defaultStreamAntithetic = defaultStream.Antithetic;
            defaultStreamFullPrecision = defaultStream.FullPrecision;
            
            %% simulate
            %run-1
            this.setUp();
            m = this.process;
            m.evolveState();
            streams1 = cell(1+numel(m.states), 1);
            streams1{1} = this.process.randStream;
            for i = 1:numel(m.states)
                streams1{i+1} = m.states{i}.randStream;
            end
            this.tearDown();
            
            %run-2
            this.setUp();
            m = this.process;
            m.evolveState();
            streams2 = cell(1+numel(m.states), 1);
            streams2{1} = this.process.randStream;
            for i = 1:numel(m.states)
                streams2{i+1} = m.states{i}.randStream;
            end
            this.tearDown();
            
            %% assert state of process and state random streams updated the same in both runs
            for i = 1:numel(streams1)
                assertTrue(streams1{i} == streams2{i});
            end
            
            %% assert state of default stream unchanged
            defaultStream = RandStream.getDefaultStream;
            assertEqual(defaultStreamType, defaultStream.Type);
            assertEqual(defaultStreamSeed, defaultStream.Seed);
            assertEqual(defaultStreamNumStreams, defaultStream.NumStreams);
            assertEqual(defaultStreamStreamIndex, defaultStream.StreamIndex);
            assertEqual(defaultStreamState, defaultStream.State);
            assertEqual(defaultStreamSubstram, defaultStream.Substream);
            assertEqual(defaultStreamRandnAlg, defaultStream.RandnAlg);
            assertEqual(defaultStreamAntithetic, defaultStream.Antithetic);
            assertEqual(defaultStreamFullPrecision, defaultStream.FullPrecision);
        end
    end
    
    %test helpers
    methods (Sealed)
        function helpTestGeneEssentiality(...
                this, expectedEssentialGenes, isFunctioningOkay, options)
            %process options
            if ~exist('options','var')
                options = struct;
            end
            lengthSec = this.optionVal(options, 'lengthSec', 1);
            knockoutStimuli = this.optionVal(options, 'knockoutStimuli', 1);
            knockoutSubstrates = this.optionVal(options, 'knockoutSubstrates', 1);
            knockoutEnzymes = this.optionVal(options, 'knockoutEnzymes', 1);
            
            m = this.process;
            
            base_globalStates = struct;
            for i = 1:numel(m.states)
                base_globalStates.(m.states{i}.wholeCellModelID) = struct;
                for j = 1:numel(m.states{i}.stateNames)
                    base_globalStates.(m.states{i}.wholeCellModelID).(m.states{i}.stateNames{j}) = ...
                        m.states{i}.(m.states{i}.stateNames{j});
                end
            end
            
            base_localStates = struct;
            for i = 1:numel(m.localStateNames)
                base_localStates.(m.localStateNames{i}) = m.(m.localStateNames{i});
            end
            
            baseStates = edu.stanford.covert.util.StructUtil.catstruct(base_globalStates, base_localStates);
            
            %compute gene composition of stimuli, substrates, enzymes
            geneWholeCellModelIDs    = m.gene.wholeCellModelIDs;
            stimuliGeneComposition   = m.stimuliGeneComposition();
            substrateGeneComposition = m.substrateGeneComposition();
            enzymeGeneComposition    = m.enzymeGeneComposition();
            
            includedGenes = false(size(geneWholeCellModelIDs));
            if knockoutStimuli
                includedGenes = includedGenes | any(stimuliGeneComposition,2);
            end
            if knockoutSubstrates
                includedGenes = includedGenes | any(substrateGeneComposition,2);
            end
            if knockoutEnzymes
                includedGenes = includedGenes | any(enzymeGeneComposition,2);
            end
            includedGenes = find(includedGenes);
            
            geneWholeCellModelIDs    = geneWholeCellModelIDs(includedGenes);
            stimuliGeneComposition   = stimuliGeneComposition(includedGenes, :);
            substrateGeneComposition = substrateGeneComposition(includedGenes, :);
            enzymeGeneComposition    = enzymeGeneComposition(includedGenes, :);
            
            %verify all expected essential genes included in test
            assertTrue(...
                all(ismember(expectedEssentialGenes, geneWholeCellModelIDs)), ...
                'not all expected essential genes are included in this test');
            
            %verify first that the process functions properly with all genes intact
            for j = 1:lengthSec
                m.evolveState();
                if isFunctioningOkay(m, baseStates)
                    break;
                end
            end
            assertTrue(...
                isFunctioningOkay(m, baseStates),...
                'process is not functioning properly with all genes intact');
            
            mispredictedAsEssential = false(size(geneWholeCellModelIDs));
            mispredictedAsNotEssential = false(size(geneWholeCellModelIDs));
            
            %for each gene
            for i = 1:length(geneWholeCellModelIDs)
                %reset process state
                m.seedRandStream();
                
                for j = 1:numel(m.states)
                    for k = 1:numel(m.states{j}.stateNames)
                        m.states{j}.(m.states{j}.stateNames{k}) = ...
                            base_globalStates.(m.states{j}.wholeCellModelID).(m.states{j}.stateNames{k});
                    end
                end
                
                for j = 1:numel(m.localStateNames)
                    m.(m.localStateNames{j}) = base_localStates.(m.localStateNames{j});
                end
                
                %knock out gene
                m.stimuli(     stimuliGeneComposition(i,:) >0, :) = 0;
                m.substrates(substrateGeneComposition(i,:) >0, :) = 0;
                m.enzymes(      enzymeGeneComposition(i,:) >0, :) = 0;
                m.boundEnzymes( enzymeGeneComposition(i,:) >0, :) = 0;
                
                this.helpKnockoutGene(geneWholeCellModelIDs{i});
                
                initial_localStates = base_localStates;
                initial_localStates.stimuli      = m.stimuli;
                initial_localStates.substrates   = m.substrates;
                initial_localStates.enzymes      = m.enzymes;
                initial_localStates.boundEnzymes = m.boundEnzymes;
                
                initialStates = edu.stanford.covert.util.StructUtil.catstruct(base_globalStates, initial_localStates);
                
                expectedEssential = ...
                    any(strcmp(geneWholeCellModelIDs{i}, expectedEssentialGenes));
                
                %execute process
                for j = 1:lengthSec
                    m.evolveState();
                    if isFunctioningOkay(m, initialStates)
                        if expectedEssential
                            mispredictedAsEssential(i) = true;
                        end
                        break;
                    end
                end
                
                %determine whether process correctly predicted gene essentiality
                if isFunctioningOkay(m, initialStates)
                    if expectedEssential
                        mispredictedAsEssential(i) = true;
                    end
                else
                    if ~expectedEssential
                        mispredictedAsNotEssential(i) = true;
                    end
                end
            end
            
            assertTrue(all(~mispredictedAsEssential), ...
                ['genes are not essential: '...
                strjoin(', ', geneWholeCellModelIDs{mispredictedAsEssential})]);
            assertTrue(all(~mispredictedAsNotEssential), ...
                ['genes are essential: '...
                strjoin(', ', geneWholeCellModelIDs{mispredictedAsNotEssential})]);
        end
        
        function result = optionVal(~, options, name, defaultValue)
            if isfield(options, name)
                result = options.(name);
            else
                result = defaultValue;
            end
        end
    end
    
    %helper methods optionally overridden by subclasses
    methods
        function helpKnockoutGene(this, geneWholeCellModelID) %#ok<MANU,INUSD>
        end
    end
end
