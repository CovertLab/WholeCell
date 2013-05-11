%Simulation unit test cases
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 9/14/2010
classdef Simulation_Test < TestCase
    properties
        simulation
        outputDir
    end

    methods
        function this = Simulation_Test(name)
            this = this@TestCase(name);
        end

        function setUp(this)
            this.simulation = this.newSimulation();
            this.outputDir = '';
        end

        function tearDown(this)
            if ~isempty(this.outputDir) && exist(this.outputDir, 'dir')
                rmdir(this.outputDir, 's');
            end
        end
    end

    %tests
    methods
        function testConstants(this)
            sim = this.simulation;
            
            assertAllEqual(false, isnan(sim.state('Metabolite').molecularWeights));
            assertEqual(1, sim.compartment.cytosolIndexs);
        end
        
        function testAllocateMemoryForState(this)
            this.simulation.allocateMemoryForState(1);  %no run-time errors
        end
        
        function testApplyPerturbations(this)
            sim = this.simulation;
            g =  sim.gene;
            
            % 1 mRNA genes
            this.helper_testApplyPerturbations(g.mRNAIndexs(1));
            
            % 2 mRNA genes
            this.helper_testApplyPerturbations(g.mRNAIndexs([1 101]));
            
            % mRNA genes
            this.helper_testApplyPerturbations(g.mRNAIndexs);
            
            % rRNA genes
            this.helper_testApplyPerturbations(g.rRNAIndexs);
            
            % sRNA genes
            this.helper_testApplyPerturbations(g.sRNAIndexs);
            
            % tRNA genes
            this.helper_testApplyPerturbations(g.tRNAIndexs);
            
            % all genes
            this.helper_testApplyPerturbations((1:numel(g.wholeCellModelIDs))');
        end
        
        function helper_testApplyPerturbations(this, geneIdxs)
            %% references
            sim = this.simulation;
            g =  sim.gene;
            r = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            ftsZRing = sim.state('FtsZRing');
            c = sim.state('Chromosome');
            rnaPol = sim.state('RNAPolymerase');
            rib = sim.state('Ribosome');
            pol = sim.state('Polypeptide');
            transcript = sim.state('Transcript');
            
            %% initialize state
            r.decayRates(:) = 0;
            pm.decayRates(:) = 0;
            pc.decayRates(:) = 0;
            r.counts(:) = 1;
            pm.counts(:) = 1;
            pc.counts(:) = 1;
            
            r.counts(r.boundIndexs(r.matureSRNAIndexs), :) = 0;
            
            c.initialize();
            vals = setdiff((1:numel(pm.boundIndexs))', pm.translationFactorIndexs);
            monPosStrnds = [100 * (1:numel(vals))'  ones(size(vals))];
            c.monomerBoundSites(monPosStrnds) = vals;
            vals = setdiff((1:numel(pc.boundIndexs))', [
                pc.ribosome70SIndexs
                pc.translationFactorIndexs
                pc.ftsZGTPIndexs
                pc.ftsZGDPIndexs
                pc.dnaAPolymerIndexs
                pc.replisomeIndexs
                ]);
            cpxPosStrnds = [1000 * (1:numel(vals))' + monPosStrnds(end, 1)  ones(size(vals))];
            c.complexBoundSites(cpxPosStrnds) = vals;
            assertIn(cpxPosStrnds(end, 1), [1 c.sequenceLen-1000]);
            pc.counts(pc.boundIndexs([pc.dnaAPolymerIndexs; pc.replisomeIndexs]), :) = 0;
            
            rnaPol.states = [
                rnaPol.nonSpecificallyBoundValue
                rnaPol.specificallyBoundValue
                ];
            rnaPol.positionStrands = [
                find(c.complexBoundSites == pc.rnaPolymeraseIndexs(1))
                find(c.complexBoundSites == pc.rnaPolymeraseIndexs(2))
                ];
            transcript.boundTranscriptionUnits = [0; 1];
            transcript.boundTranscriptProgress = [0; 0];
            transcript.boundTranscriptChromosome = [0; 0];
            
            rib.states = rib.activeValue;
            rib.boundMRNAs = 1;
            rib.mRNAPositions = 1;
            rib.tmRNAPositions = 0;
            pol.boundMRNAs = rib.boundMRNAs;
            pol.nascentMonomerLengths = rib.mRNAPositions;
            pol.proteolysisTagLengths = rib.tmRNAPositions;
            
            pc.counts(pc.boundIndexs([pc.ftsZGTPIndexs; pc.ftsZGDPIndexs]), :) = 0;
            ftsZRing.numEdgesOneStraight = 0;
            ftsZRing.numEdgesTwoStraight = 0;
            ftsZRing.numEdgesTwoBent = 0;
            ftsZRing.numResidualBent = 0;
            
            %% apply perturbations
            monIdxs = find(ismember(g.mRNAIndexs, geneIdxs));
            matRnaIdxs = find(any(r.matureRNAGeneComposition(setdiff(geneIdxs, g.mRNAIndexs), :), 1));
            cpxIdxs = find(any(any(pc.proteinComplexComposition(geneIdxs, :, :), 3), 1));
            
            otherMatRnaIdxs = setdiff(1:numel(r.matureIndexs), matRnaIdxs);
            otherMonIdxs = setdiff(1:numel(pm.matureIndexs), monIdxs);
            otherCpxIdxs = setdiff(1:numel(pc.matureIndexs), cpxIdxs);
            
            rnaDecayRates = r.decayRates;
            monDecayRates = pm.decayRates;
            cpxDecayRates = pc.decayRates;
            
            rnaCounts = r.counts;
            monCounts = pm.counts;
            cpxCounts = pc.counts;
            
            sim.applyOptions('geneticKnockouts', g.wholeCellModelIDs(geneIdxs));
            sim.applyPerturbations();
            
            %% assertions
            %decay rates
            decayRate = log(2)/(log(2)/1e6);
            
            assertAllEqual(0, r.decayRates(r.nascentIndexs));
            assertAllEqual(0, r.decayRates(r.intergenicIndexs));
            if ~isempty(matRnaIdxs)
                assertAllEqual(decayRate, r.decayRates(r.processedIndexs(matRnaIdxs)));
                assertAllEqual(decayRate, r.decayRates(r.matureIndexs(matRnaIdxs)));
                assertAllEqual(decayRate, r.decayRates(r.boundIndexs(matRnaIdxs)));
                assertAllEqual(decayRate, r.decayRates(r.misfoldedIndexs(matRnaIdxs)));
                assertAllEqual(decayRate, r.decayRates(r.damagedIndexs(matRnaIdxs)));
                assertAllEqual(decayRate, r.decayRates(r.aminoacylatedIndexs(matRnaIdxs)));
            end
            assertEqual(rnaDecayRates(r.processedIndexs(otherMatRnaIdxs)), r.decayRates(r.processedIndexs(otherMatRnaIdxs)));
            assertEqual(rnaDecayRates(r.matureIndexs(otherMatRnaIdxs)), r.decayRates(r.matureIndexs(otherMatRnaIdxs)));
            assertEqual(rnaDecayRates(r.boundIndexs(otherMatRnaIdxs)), r.decayRates(r.boundIndexs(otherMatRnaIdxs)));
            assertEqual(rnaDecayRates(r.misfoldedIndexs(otherMatRnaIdxs)), r.decayRates(r.misfoldedIndexs(otherMatRnaIdxs)));
            assertEqual(rnaDecayRates(r.damagedIndexs(otherMatRnaIdxs)), r.decayRates(r.damagedIndexs(otherMatRnaIdxs)));
            assertEqual(rnaDecayRates(r.aminoacylatedIndexs(otherMatRnaIdxs)), r.decayRates(r.aminoacylatedIndexs(otherMatRnaIdxs)));
            
            if ~isempty(monIdxs)
                assertAllEqual(decayRate, pm.decayRates(pm.nascentIndexs(monIdxs)));
                assertAllEqual(decayRate, pm.decayRates(pm.signalSequenceIndexs(monIdxs)));
                assertAllEqual(decayRate, pm.decayRates(pm.processedIIndexs(monIdxs)));
                assertAllEqual(decayRate, pm.decayRates(pm.processedIIIndexs(monIdxs)));
                assertAllEqual(decayRate, pm.decayRates(pm.foldedIndexs(monIdxs)));
                assertAllEqual(decayRate, pm.decayRates(pm.matureIndexs(monIdxs)));
                assertAllEqual(decayRate, pm.decayRates(pm.inactivatedIndexs(monIdxs)));
                assertAllEqual(decayRate, pm.decayRates(pm.boundIndexs(monIdxs)));
                assertAllEqual(decayRate, pm.decayRates(pm.misfoldedIndexs(monIdxs)));
                assertAllEqual(decayRate, pm.decayRates(pm.damagedIndexs(monIdxs)));
            end
            assertEqual(monDecayRates(pm.nascentIndexs(otherMonIdxs)), pm.decayRates(pm.nascentIndexs(otherMonIdxs)));
            assertEqual(monDecayRates(pm.signalSequenceIndexs(otherMonIdxs)), pm.decayRates(pm.signalSequenceIndexs(otherMonIdxs)));
            assertEqual(monDecayRates(pm.processedIIndexs(otherMonIdxs)), pm.decayRates(pm.processedIIndexs(otherMonIdxs)));
            assertEqual(monDecayRates(pm.processedIIIndexs(otherMonIdxs)), pm.decayRates(pm.processedIIIndexs(otherMonIdxs)));
            assertEqual(monDecayRates(pm.foldedIndexs(otherMonIdxs)), pm.decayRates(pm.foldedIndexs(otherMonIdxs)));
            assertEqual(monDecayRates(pm.matureIndexs(otherMonIdxs)), pm.decayRates(pm.matureIndexs(otherMonIdxs)));
            assertEqual(monDecayRates(pm.inactivatedIndexs(otherMonIdxs)), pm.decayRates(pm.inactivatedIndexs(otherMonIdxs)));
            assertEqual(monDecayRates(pm.boundIndexs(otherMonIdxs)), pm.decayRates(pm.boundIndexs(otherMonIdxs)));
            assertEqual(monDecayRates(pm.misfoldedIndexs(otherMonIdxs)), pm.decayRates(pm.misfoldedIndexs(otherMonIdxs)));
            assertEqual(monDecayRates(pm.damagedIndexs(otherMonIdxs)), pm.decayRates(pm.damagedIndexs(otherMonIdxs)));
            
            if ~isempty(cpxIdxs)
                assertAllEqual(decayRate, pc.decayRates(pc.nascentIndexs(cpxIdxs)));
                assertAllEqual(decayRate, pc.decayRates(pc.matureIndexs(cpxIdxs)));
                assertAllEqual(decayRate, pc.decayRates(pc.inactivatedIndexs(cpxIdxs)));
                assertAllEqual(decayRate, pc.decayRates(pc.boundIndexs(cpxIdxs)));
                assertAllEqual(decayRate, pc.decayRates(pc.misfoldedIndexs(cpxIdxs)));
                assertAllEqual(decayRate, pc.decayRates(pc.damagedIndexs(cpxIdxs)));
            end
            assertEqual(cpxDecayRates(pc.nascentIndexs(otherCpxIdxs)), pc.decayRates(pc.nascentIndexs(otherCpxIdxs)));
            assertEqual(cpxDecayRates(pc.matureIndexs(otherCpxIdxs)), pc.decayRates(pc.matureIndexs(otherCpxIdxs)));
            assertEqual(cpxDecayRates(pc.inactivatedIndexs(otherCpxIdxs)), pc.decayRates(pc.inactivatedIndexs(otherCpxIdxs)));
            assertEqual(cpxDecayRates(pc.boundIndexs(otherCpxIdxs)), pc.decayRates(pc.boundIndexs(otherCpxIdxs)));
            assertEqual(cpxDecayRates(pc.misfoldedIndexs(otherCpxIdxs)), pc.decayRates(pc.misfoldedIndexs(otherCpxIdxs)));
            assertEqual(cpxDecayRates(pc.damagedIndexs(otherCpxIdxs)), pc.decayRates(pc.damagedIndexs(otherCpxIdxs)));
            
            %counts
            assertAllEqual(1, r.counts(r.nascentIndexs, :));
            assertAllEqual(1, r.counts(r.intergenicIndexs, :));
            if ~isempty(matRnaIdxs)
                assertAllEqual(0, r.counts(r.processedIndexs(matRnaIdxs), :));
                assertAllEqual(0, r.counts(r.matureIndexs(matRnaIdxs), :));
                assertAllEqual(0, r.counts(r.boundIndexs(matRnaIdxs), :));
                assertAllEqual(0, r.counts(r.misfoldedIndexs(matRnaIdxs), :));
                assertAllEqual(0, r.counts(r.damagedIndexs(matRnaIdxs), :));
                assertAllEqual(0, r.counts(r.aminoacylatedIndexs(matRnaIdxs), :));
            end
            assertEqual(rnaCounts(r.processedIndexs(otherMatRnaIdxs, :)), r.counts(r.processedIndexs(otherMatRnaIdxs, :)));
            assertEqual(rnaCounts(r.matureIndexs(otherMatRnaIdxs, :)), r.counts(r.matureIndexs(otherMatRnaIdxs, :)));
            assertEqual(rnaCounts(r.boundIndexs(otherMatRnaIdxs, :)), r.counts(r.boundIndexs(otherMatRnaIdxs, :)));
            assertEqual(rnaCounts(r.misfoldedIndexs(otherMatRnaIdxs, :)), r.counts(r.misfoldedIndexs(otherMatRnaIdxs, :)));
            assertEqual(rnaCounts(r.damagedIndexs(otherMatRnaIdxs, :)), r.counts(r.damagedIndexs(otherMatRnaIdxs, :)));
            assertEqual(rnaCounts(r.aminoacylatedIndexs(otherMatRnaIdxs, :)), r.counts(r.aminoacylatedIndexs(otherMatRnaIdxs, :)));
            
            if ~isempty(monIdxs)
                assertAllEqual(0, pm.counts(pm.nascentIndexs(monIdxs), :));
                assertAllEqual(0, pm.counts(pm.signalSequenceIndexs(monIdxs), :));
                assertAllEqual(0, pm.counts(pm.processedIIndexs(monIdxs), :));
                assertAllEqual(0, pm.counts(pm.processedIIIndexs(monIdxs), :));
                assertAllEqual(0, pm.counts(pm.foldedIndexs(monIdxs), :));
                assertAllEqual(0, pm.counts(pm.matureIndexs(monIdxs), :));
                assertAllEqual(0, pm.counts(pm.inactivatedIndexs(monIdxs), :));
                assertAllEqual(0, pm.counts(pm.boundIndexs(monIdxs), :));
                assertAllEqual(0, pm.counts(pm.misfoldedIndexs(monIdxs), :));
                assertAllEqual(0, pm.counts(pm.damagedIndexs(monIdxs), :));
            end
            assertEqual(monCounts(pm.nascentIndexs(otherMonIdxs, :)), pm.counts(pm.nascentIndexs(otherMonIdxs, :)));
            assertEqual(monCounts(pm.signalSequenceIndexs(otherMonIdxs, :)), pm.counts(pm.signalSequenceIndexs(otherMonIdxs, :)));
            assertEqual(monCounts(pm.processedIIndexs(otherMonIdxs, :)), pm.counts(pm.processedIIndexs(otherMonIdxs, :)));
            assertEqual(monCounts(pm.processedIIIndexs(otherMonIdxs, :)), pm.counts(pm.processedIIIndexs(otherMonIdxs, :)));
            assertEqual(monCounts(pm.foldedIndexs(otherMonIdxs, :)), pm.counts(pm.foldedIndexs(otherMonIdxs, :)));
            assertEqual(monCounts(pm.matureIndexs(otherMonIdxs, :)), pm.counts(pm.matureIndexs(otherMonIdxs, :)));
            assertEqual(monCounts(pm.inactivatedIndexs(otherMonIdxs, :)), pm.counts(pm.inactivatedIndexs(otherMonIdxs, :)));
            assertEqual(monCounts(pm.boundIndexs(otherMonIdxs, :)), pm.counts(pm.boundIndexs(otherMonIdxs, :)));
            assertEqual(monCounts(pm.misfoldedIndexs(otherMonIdxs, :)), pm.counts(pm.misfoldedIndexs(otherMonIdxs, :)));
            assertEqual(monCounts(pm.damagedIndexs(otherMonIdxs, :)), pm.counts(pm.damagedIndexs(otherMonIdxs, :)));
            
            if ~isempty(cpxIdxs)
                assertAllEqual(0, pc.counts(pc.nascentIndexs(cpxIdxs), :));
                assertAllEqual(0, pc.counts(pc.matureIndexs(cpxIdxs), :));
                assertAllEqual(0, pc.counts(pc.inactivatedIndexs(cpxIdxs), :));
                assertAllEqual(0, pc.counts(pc.boundIndexs(cpxIdxs), :));
                assertAllEqual(0, pc.counts(pc.misfoldedIndexs(cpxIdxs), :));
                assertAllEqual(0, pc.counts(pc.damagedIndexs(cpxIdxs), :));
            end
            assertEqual(cpxCounts(pc.nascentIndexs(otherCpxIdxs, :)), pc.counts(pc.nascentIndexs(otherCpxIdxs, :)));
            assertEqual(cpxCounts(pc.matureIndexs(otherCpxIdxs, :)), pc.counts(pc.matureIndexs(otherCpxIdxs, :)));
            assertEqual(cpxCounts(pc.inactivatedIndexs(otherCpxIdxs, :)), pc.counts(pc.inactivatedIndexs(otherCpxIdxs, :)));
            assertEqual(cpxCounts(pc.boundIndexs(otherCpxIdxs, :)), pc.counts(pc.boundIndexs(otherCpxIdxs, :)));
            assertEqual(cpxCounts(pc.misfoldedIndexs(otherCpxIdxs, :)), pc.counts(pc.misfoldedIndexs(otherCpxIdxs, :)));
            assertEqual(cpxCounts(pc.damagedIndexs(otherCpxIdxs, :)), pc.counts(pc.damagedIndexs(otherCpxIdxs, :)));
            
            %% test membrane protein decay
            %monomer
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], true);
            comp = sim.compartment;
            g = sim.gene;
            pm = sim.state('ProteinMonomer');
            protDcy = sim.process('ProteinDecay');
            
            idx = find(pm.compartments(pm.matureIndexs) == comp.membraneIndexs, 1, 'first');
            
            sim.applyOptions('geneticKnockouts', g.wholeCellModelIDs(g.mRNAIndexs(idx)));
            sim.applyPerturbations();
            
            monIdx = pm.matureIndexs(idx);
            pm.counts(monIdx, :) = 0;
            pm.counts(monIdx, comp.membraneIndexs) = 1;
            
            protDcy.copyFromState();
            protDcy.enzymes(:) = 1e6;
            protDcy.substrates(:) = 1e6;
            protDcy.evolveState();
            protDcy.copyToState();
            
            assertEqual(0, pm.counts(monIdx, comp.membraneIndexs));
            
            %complex
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], true);
            comp = sim.compartment;
            g = sim.gene;
            pc = sim.state('ProteinComplex');
            protDcy = sim.process('ProteinDecay');
            
            idx = find(pc.compartments(pc.matureIndexs) == comp.membraneIndexs, 1, 'first');
            geneIdx = find(any(pc.proteinComplexComposition(:, idx, :), 3), 1, 'first');
            
            sim.applyOptions('geneticKnockouts', g.wholeCellModelIDs(geneIdx));
            sim.applyPerturbations();
            
            cpxIdx = pc.matureIndexs(idx);
            pc.counts(cpxIdx, :) = 0;
            pc.counts(cpxIdx, comp.membraneIndexs) = 1;
            
            protDcy.copyFromState();
            protDcy.enzymes(:) = 1e6;
            protDcy.substrates(:) = 1e6;
            protDcy.evolveState();
            protDcy.copyToState();
            
            assertEqual(0, pc.counts(cpxIdx, comp.membraneIndexs));
        end
        
        function testEvolveState(this)
            sim = this.simulation;
            met = sim.state('Metabolite');
            
            t0 = sim.state('Time').values;
            
            %mechanism for recording a call sequence
            calls = {};
            function r = recorder(s, value)
                if ~exist('value', 'var')
                    value = [];
                end
                function result = onCall(varargin)
                    calls{end+1} = s;
                    if isa(value, 'function_handle')
                        value(varargin{:}); %#ok<NOEFF>
                        result = [];
                    else
                        result = value;
                    end
                end
                r = @onCall;
            end

            %mechanism for using substrates
            function u = user(substrates)
                function use(process)
                    process.substrates = process.substrates - substrates;
                end
                u = @use;
            end
            
            %metabolite indices of interest (no special significance)
            M1 = 7; M2 = 123; M3 = 432;
            
            %three fake processes
            sim.setForTest('processes', {
                %simple case: all metabolites in same compartment, all
                %substrates are metabolites
                edu.stanford.covert.test.mock(...
                    'substrates', [0 0]',...
                    'substrateMetaboliteLocalIndexs', [1 2]',...
                    'substrateMetaboliteGlobalIndexs', [M1 M3]',...
                    'substrateMetaboliteCompartmentIndexs', [1 1]',...
                    'substrateMetaboliteGlobalCompartmentIndexs', sub2ind(size(met.counts), [M1 M3]', [1 1]'), ...
                    'calcResourceRequirements_Current', recorder('r1', [10 10]'),...
                    'copyFromState', recorder('f1'),...
                    'copyToState', recorder('t1'),...
                    'evolveState', recorder('e1', user([3 4]')))
                %more complex: metabolites in various compartments, some
                %substrates are not metabolites
                edu.stanford.covert.test.mock(...
                    'substrates', [0 0 0 0 0]',...
                    'substrateMetaboliteLocalIndexs', [1 2 5]',...
                    'substrateMetaboliteGlobalIndexs', [M1 M2 M3]',...
                    'substrateMetaboliteCompartmentIndexs', [3 1 1]',...
                    'substrateMetaboliteGlobalCompartmentIndexs', sub2ind(size(met.counts), [M1 M2 M3]', [3 1 1]'), ...
                    'calcResourceRequirements_Current', recorder('r2', [5 5 0 0 5]'),...
                    'copyFromState', recorder('f2'),...
                    'copyToState', recorder('t2'),...
                    'evolveState', recorder('e2', user([1 0 0 0 -1]')))
                %each metabolite in multiple compartments
                edu.stanford.covert.test.mock(...
                    'substrates', [0 0 0 0; 0 0 0 0]',...
                    'substrateMetaboliteLocalIndexs', [1 3 4]',...
                    'substrateMetaboliteGlobalIndexs', [M2 M3 M1; M2 M3 M1]',...
                    'substrateMetaboliteCompartmentIndexs', [5 4 5; 3 1 3]',...
                    'substrateMetaboliteGlobalCompartmentIndexs', sub2ind(size(met.counts), [M2 M3 M1; M2 M3 M1]', [5 4 5; 3 1 3]'), ...
                    'calcResourceRequirements_Current', recorder('r3', [5 0 5 0; 0 0 5 5]'),...
                    'copyFromState', recorder('f3'),...
                    'copyToState', recorder('t3'),...
                    'evolveState', recorder('e3'))
                });
            sim.setForTest('processEvalOrderIndexs', [3 2 1]');
            sim.setForTest('processesInEvalOrder', sim.processes(sim.processEvalOrderIndexs));
            met.counts(M1,:) = 20;
            met.counts(M2,:) = 0;
            met.counts(M3,:) = 10;

            sim.evolveState();
            assertEqual(t0 + sim.stepSizeSec, sim.state('Time').values);
            assertEqual(...
                'f1 r1 f2 r2 f3 r3 f3 e3 t3 f2 e2 t2 f1 e1 t1',...
                strjoin(' ', calls{:}));

            %processes receive proporational to their allocation
            assertEqual(...
                [17 1],...
                sim.processes{1}.substrates');
            assertEqual(...
                [9 0 0 0 3],...
                sim.processes{2}.substrates');
            assertEqual(...
                [0 0 10 0;
                 0 0 2 10],...
                sim.processes{3}.substrates');
            
            %verify that metabolite/substrate usage/production propagate back
            assertEqual([17 20 19 20 20 20], met.counts(M1,:));
            assertEqual([ 0  0  0  0  0  0], met.counts(M2,:));
            assertEqual([ 7 10 10 10 10 10], met.counts(M3,:));
        end
        
        function testRunLogging(this)
            sim = this.simulation;
            sim.applyOptions('lengthSec', 10);
            sim.state('MetabolicReaction').initialGrowthFilterWidth = Inf;
            
            %1 logger
            logger = edu.stanford.covert.cell.sim.util.SummaryLogger(2, 0);
            sim.run(logger);
            
            %2 loggers
            logger1 = edu.stanford.covert.cell.sim.util.SummaryLogger(2, 0);
            logger2 = edu.stanford.covert.cell.sim.util.SummaryLogger(2, 0);
            sim.run({logger1 logger2});
        end
        
        function testMassBalance(this)
            sim = this.simulation;
            sim.applyOptions('lengthSec', 100);
            sim.state('MetabolicReaction').initialGrowthFilterWidth = Inf;
            
            sim.allocateMemoryForState(1);
            sim.initializeState();
            sim.applyPerturbations();
            
            mass = sim.state('Mass');
            initialMass = sum(mass.cell);
            
            sim.run();
            
            assertElementsAlmostEqual(initialMass, sum(mass.cell));
        end

        function testPlots(this)
            sim = this.simulation;

            %create figure, add axis
            handle = figure();
            axesHandle = subplot(1,1,1,'Parent',handle);
            legend(axesHandle,'off');
            set(axesHandle,'Box','on');

            %inputs
            time = 1;
            compartments = 1;

            %call all plotting functions
            simulationMetaData = metaclass(sim);
            failedPlotMethodNames = {};
            for i = 1:length(simulationMetaData.Methods)
                if length(simulationMetaData.Methods{i}.Name) > 4 && ...
                   strcmp(simulationMetaData.Methods{i}.Access, 'public') && ...
                   strcmp(simulationMetaData.Methods{i}.Name(1:4), 'plot') && ...
                   numel(simulationMetaData.Methods{i}.InputNames) == 4
                    cla(axesHandle);
                    try
                        sim.(simulationMetaData.Methods{i}.Name)(axesHandle, time, compartments)
                    catch %#ok<CTCH>
                        failedPlotMethodNames{end+1}=simulationMetaData.Methods{i}.Name; %#ok<AGROW>
                    end
                end
            end

            %close figure
            close(handle);

            assert(isempty(failedPlotMethodNames), ...
                ['The following plot methods failed: ' strjoin(', ',failedPlotMethodNames{:})]);
        end
        
        function testOptions(this)
            sim = this.simulation;
            assertEqual(0, sim.getOptions.verbosity);
            assertEqual(1, sim.getOptions.stepSizeSec);

            sim.applyOptions();
            assertEqual(0, sim.getOptions.verbosity);
            assertEqual(1, sim.getOptions.stepSizeSec);

            sim.applyOptions('verbosity', 3, 'stepSizeSec', 0.5);
            assertEqual(3, sim.getOptions.verbosity);
            assertEqual(0.5, sim.getOptions.stepSizeSec);

            sim.applyOptions('verbosity', 4);
            assertEqual(4, sim.getOptions.verbosity);
            assertEqual(0.5, sim.getOptions.stepSizeSec);
        end
    end

    %helpers
    methods (Access = private)
        function sim = newSimulation(~)
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load();
        end

        function prepareOutputDir(this)
            this.outputDir = ['tmp/Simulation_Integrated_Test_' this.MethodName];
            if ~exist(this.outputDir, 'dir')
                mkdir(this.outputDir);
            end
            delete([this.outputDir '/*']);
        end

        function result = createMetadata(this)
            result = struct;
            result.shortDescription = 'Wild-type simulation';
            result.longDescription  = 'Wild-type simulation';
            result.email            = 'jkarr@stanford.edu';
            result.firstName        = 'Jonathan';
            result.lastName         = 'Karr';
            result.affiliation      = 'Covert Lab, Department of Bioengineering, Stanford University';
            result.knowledgeBaseWID = this.simulation.knowledgeBaseWID;
            result.downsampleStep   = 5000;
            result.outputDirectory  = this.outputDir;
            result.startTime        = datestr(clock,'yyyy-mm-dd HH:MM:SS');
            result.endTime          = datestr(clock,'yyyy-mm-dd HH:MM:SS');
            [result.userName, result.hostName, result.ipAddress] = ...
                edu.stanford.covert.util.computerInfo();
            result.revision = 234;
            result.differencesFromRevision = '';
        end
    end

    %static helpers
    methods (Static, Access = private)
        function assertPropertyValuesEqual(...
                metadata, options, parameters, fittedConstants, timeCourses, sim)
            %meta data
            metadataNames = fieldnames(sim.metadata);
            for i = 1:length(metadataNames)
                name = metadataNames{i};
                assertEqual(metadata.(name), sim.metadata.(name), name);
            end
            
            %options
            optionNames = sim.optionNames;
            for i = 1:length(optionNames)
                name = optionNames{i};
                switch class(options.Core__.(name))
                    case {'cell','struct'}
                        assertEqual(options.Core__.(name), sim.(name));
                    otherwise
                        assertElementsAlmostEqual(...
                            options.Core__.(name), sim.(name),...
                            'relative', 1e-4, name);
                end
                    
            end
            for i = 1:length(sim.processes)
                process = sim.processes{i};
                optionNames = process.optionNames;
                for j = 1:length(optionNames)
                    name = optionNames{j};
                    switch class(options.Core__.(name))
                        case {'cell','struct'}
                            assertEqual(options.(process.wholeCellModelID).(name), process.(name));
                        otherwise
                            assertElementsAlmostEqual(...
                                options.(process.wholeCellModelID).(name),...
                                process.(name),...
                                'relative', 1e-4, name);
                    end
                end
            end

            %parameters
            parameterNames = sim.parameterNames;
            simulationfieldnames = fieldnames(sim);
            for i = 1:length(parameterNames)
                name = parameterNames{i};
                if ismember(name, simulationfieldnames)
                    if isnumeric(parameters.Core__.(name))
                        assertElementsAlmostEqual(...
                            parameters.Core__.(name), sim.(name),...
                            'relative', 1e-4, name);
                    else
                        assertEqual(parameters.Core__.(name), sim.(name), name);
                    end
                end
            end
            for i = 1:length(sim.processes)
                process = sim.processes{i};
                parameterNames = process.parameterNames;
                processFieldnames = fieldnames(process);
                for j = 1:length(parameterNames)
                    name = parameterNames{j};
                    if ismember(name, processFieldnames)
                        if isnumeric(parameters.(process.wholeCellModelID).(name))
                            assertElementsAlmostEqual(...
                                parameters.(process.wholeCellModelID).(name), ...
                                process.(name),...
                                'relative', 1e-4, name);
                        else
                            assertEqual(...
                                parameters.(process.wholeCellModelID).(name), ...
                                process.(name), name);
                        end
                    end
                end
            end

            %fit constants
            fittedConstantNames = sim.fittedConstantNames;
            for i = 1:length(fittedConstantNames)
                name = fittedConstantNames{i};
                assertElementsAlmostEqual(...
                    fittedConstants.Core__.(name), sim.(name),...
                    'relative', 1e-4, name);
            end
            for i = 1:length(sim.processes)
                process = sim.processes{i};
                fittedConstantNames = process.fittedConstantNames;
                for j = 1:length(fittedConstantNames)
                    name = fittedConstantNames{j};
                    assertElementsAlmostEqual(...
                        fittedConstants.(process.wholeCellModelID).(name),...
                        process.(name),...
                        'relative', 1e-4, name);
                end
            end

            %time courses
            timeCourseNames = sim.timeCourseNames;

            for i = 1:length(timeCourseNames)
                name = timeCourseNames{i};                
                assertElementsAlmostEqual(...
                    timeCourses.Core__.(name), sim.(name),...
                    'relative', 1e-4, name);
            end

            for i = 1:length(sim.processes)
                process = sim.processes{i};
                timeCourseNames = process.timeCourseNames;
                for j = 1:length(timeCourseNames)
                    name = timeCourseNames{j};
                    if isa(timeCourses.(sim.processWholeCellModelIDs{i}).(name), 'edu.stanford.covert.util.CircularSparseMat')
                        assertTrue(...
                            isalmostequal(timeCourses.Core__.(name), sim.(name)),...
                            name);
                    else
                        assertElementsAlmostEqual(...
                            timeCourses.(process.wholeCellModelID).(name), ...
                            process.(name),...
                            'relative', 1e-4,...
                            name);
                    end
                end
            end
        end
    end
    
    %test division
    methods
        function testDivideState(~)
            import edu.stanford.covert.cell.sim.constant.Condition;
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;      
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil
            
            %load simulation
            %sim = CachedSimulationObjectUtil.load();
            
            %baseOutDir = [SimulationDiskUtil.getBaseDir() filesep '2012_10_24_00_49_53'];
            %simDir = [baseOutDir filesep '2'];
            %md = DiskLogger.loadMetadata(simDir);
            %sim = loadSimulation(simDir, md.lengthSec, md.lengthSec, 1, 'extract');
            
            tmp = load('src_test/+edu/+stanford/+covert/+cell/+sim/fixtures/Simulation_EndOfCellCycle.mat');
            sim = tmp.simulation;
            
            %references to constants, states, processes
            c = sim.compartment;
            chr = sim.state('Chromosome');
            met = sim.state('Metabolite');
            r = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            trl = sim.process('Translation');
            
            %divide cell into 2 daughters
            daughters = sim.divideState();
            daughter1 = daughters{1};
            daughter2 = daughters{2};
            
            %test options, parameters copied correctly
            assertEqual(sim.getOptions(), daughter1.getOptions());
            assertEqual(sim.getParameters(), daughter1.getParameters());
            
            %% test
            % mass
            assertElementsAlmostEqual(sim.state('Mass').cell, ...
                + daughter1.state('Mass').cell ...
                + daughter2.state('Mass').cell, 'relative', 1e-100);
            
            %bound proteins
            d1_chr = daughter1.state('Chromosome');
            d2_chr = daughter2.state('Chromosome');
            d1_pm = daughter1.state('ProteinMonomer');
            d2_pm = daughter2.state('ProteinMonomer');
            d1_pc = daughter1.state('ProteinComplex');
            d2_pc = daughter2.state('ProteinComplex');

            [~, vals] = find(chr.monomerBoundSites);
            [~, vals1] = find(d1_chr.monomerBoundSites);
            [~, vals2] = find(d2_chr.monomerBoundSites);
            cnts = histc(vals, (1:numel(pm.boundIndexs)));
            cnts1 = histc(vals1, (1:numel(pm.boundIndexs)));
            cnts2 = histc(vals2, (1:numel(pm.boundIndexs)));
            cnts = cnts(:);
            cnts1 = cnts1(:);
            cnts2 = cnts2(:);
            assertEqual(cnts, cnts1 + cnts2);
            trlFactorIdxs = setdiff(pm.getIndexs(trl.enzymeWholeCellModelIDs(trl.enzymeIndexs_elongationFactors)), 0);
            idxs = setdiff((1:numel(pm.boundIndexs))', trlFactorIdxs);
            assertEqual(cnts(idxs), pm.counts(pm.boundIndexs(idxs), c.cytosolIndexs));
            
            [~, vals] = find(chr.complexBoundSites);
            [~, vals1] = find(d1_chr.complexBoundSites);
            [~, vals2] = find(d2_chr.complexBoundSites);
            cnts = histc(vals, (1:numel(pc.boundIndexs)));
            cnts1 = histc(vals1, (1:numel(pc.boundIndexs)));
            cnts2 = histc(vals2, (1:numel(pc.boundIndexs)));
            cnts = cnts(:);
            cnts1 = cnts1(:);
            cnts2 = cnts2(:);
            assertEqual(cnts, cnts1 + cnts2);
            trlFactorIdxs = setdiff(pc.getIndexs(trl.enzymeWholeCellModelIDs(trl.enzymeIndexs_elongationFactors)), 0);
            idxs = setdiff((1:numel(pc.boundIndexs))', [
                trlFactorIdxs
                pc.getIndexs({'MG_224_9MER_GTP'; 'MG_224_9MER_GDP'})
                pc.ribosome70SIndexs
                ]);
            assertEqual(cnts(idxs), pc.counts(pc.boundIndexs(idxs), c.cytosolIndexs));
            
            %total metabolites, RNA, protein except, FtsZ
            d1_met = daughter1.state('Metabolite');
            d2_met = daughter2.state('Metabolite');
            d1_r = daughter1.state('Rna');
            d2_r = daughter2.state('Rna');
            
            assertEqual(met.counts(:, setdiff(1:end, c.extracellularIndexs)), ...
                + d1_met.counts(:, setdiff(1:end, c.extracellularIndexs)) ...
                + d2_met.counts(:, setdiff(1:end, c.extracellularIndexs)));
            
            mediaIdxs = met.setCounts(:, Condition.objectIndexs);
            nonMediaIdxs = setdiff((1:numel(met.wholeCellModelIDs))', mediaIdxs);
            assertEqual(met.counts(nonMediaIdxs, c.extracellularIndexs), ...
                + d1_met.counts(nonMediaIdxs, c.extracellularIndexs) ...
                + d2_met.counts(nonMediaIdxs, c.extracellularIndexs));
            assertEqual(...
                d1_met.counts(mediaIdxs, c.extracellularIndexs), ...
                d2_met.counts(mediaIdxs, c.extracellularIndexs));
            
            assertEqual(r.counts, ...
                + d1_r.counts ...
                + d2_r.counts);
            assertEqual(pm.counts, ...
                + d1_pm.counts ...
                + d2_pm.counts);
            
            tmp = pc.getIndexs({'MG_224_9MER_GTP'; 'MG_224_9MER_GDP'});
            ftsZIdxs = [pc.matureIndexs(tmp) pc.boundIndexs(tmp)];
            idxs = setdiff(1:numel(pc.wholeCellModelIDs), ftsZIdxs(:));
            assertEqual(pc.counts(idxs, :), ...
                + d1_pc.counts(idxs, :) ...
                + d2_pc.counts(idxs, :));
            assertEqual(sum(pc.counts(ftsZIdxs(1, :), :), 1), ...
                + sum(d1_pc.counts(ftsZIdxs(1, :), :), 1) ...
                + sum(d2_pc.counts(ftsZIdxs(1, :), :), 1));
            assertEqual(sum(pc.counts(ftsZIdxs(2, :), :), 1), ...
                + sum(d1_pc.counts(ftsZIdxs(2, :), :), 1) ...
                + sum(d2_pc.counts(ftsZIdxs(2, :), :), 1));
            assertAllEqual(0, d1_pc.counts(ftsZIdxs(:, 2), :));
            assertAllEqual(0, d2_pc.counts(ftsZIdxs(:, 2), :));
            
            %metabolite, RNA, protein counts non-negative
            validateattributes(daughter1.state('Metabolite').counts, {'numeric'}, {'nonnegative'});
            validateattributes(daughter2.state('Metabolite').counts, {'numeric'}, {'nonnegative'});
            validateattributes(d1_r.counts, {'numeric'}, {'nonnegative'});
            validateattributes(d2_r.counts, {'numeric'}, {'nonnegative'});
            validateattributes(d1_pm.counts, {'numeric'}, {'nonnegative'});
            validateattributes(d2_pm.counts, {'numeric'}, {'nonnegative'});
            validateattributes(d1_pc.counts, {'numeric'}, {'nonnegative'});
            validateattributes(d2_pc.counts, {'numeric'}, {'nonnegative'});
            
            %total transcripts, RNA polymerase
            rnaPol = sim.state('RNAPolymerase');
            d1_rnaPol = daughter1.state('RNAPolymerase');
            d2_rnaPol = daughter2.state('RNAPolymerase');
            transcript = sim.state('Transcript');
            d1_transcript = daughter1.state('Transcript');
            d2_transcript = daughter2.state('Transcript');
            
            tmp1 = [rnaPol.states rnaPol.positionStrands transcript.boundTranscriptionUnits transcript.boundTranscriptProgress transcript.boundTranscriptChromosome];
            tmp2 = [
                d1_rnaPol.states d1_rnaPol.positionStrands(:, 1) d1_rnaPol.positionStrands(:, 2)                                                 d1_transcript.boundTranscriptionUnits d1_transcript.boundTranscriptProgress d1_transcript.boundTranscriptChromosome
                d2_rnaPol.states d2_rnaPol.positionStrands(:, 1) d2_rnaPol.positionStrands(:, 2)+2*(d2_rnaPol.states ~= rnaPol.freeValue) d2_transcript.boundTranscriptionUnits d2_transcript.boundTranscriptProgress 2*d2_transcript.boundTranscriptChromosome
                ];
            tmp1 = tmp1(tmp1(:, 1) ~= rnaPol.notExistValue, :);
            tmp2 = tmp2(tmp2(:, 1) ~= rnaPol.notExistValue, :);
            tmp1 = sortrows(tmp1, 1:6);
            tmp2 = sortrows(tmp2, 1:6);
            assertEqual(tmp1, tmp2);
            
            for i = 1:numel(daughters)
                d_pc = daughters{i}.state('ProteinComplex');
                d_rnaPol = daughters{i}.state('RNAPolymerase');
                d_transcript = daughters{i}.state('Transcript');
                
                assertEqual(...
                    sum(d_pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs), c.cytosolIndexs)), ...
                    d_rnaPol.nActive + d_rnaPol.nNonSpecificallyBound + d_rnaPol.nSpecificallyBound);
                assertEqual(...
                    d_pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs(2), c.cytosolIndexs)), ...
                    sum(d_rnaPol.states == rnaPol.activelyTranscribingValue) + d_rnaPol.nSpecificallyBound);
                assertEqual(...
                    d_pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs(1), c.cytosolIndexs)), ...
                    sum(d_rnaPol.states > rnaPol.activelyTranscribingValue) + d_rnaPol.nNonSpecificallyBound);
                
                assertEqual(...
                    sum(d_pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), c.cytosolIndexs)), ...
                    d_rnaPol.nFree);
                
                assertEqual(double(d_rnaPol.states >= rnaPol.activelyTranscribingValue | d_rnaPol.states == rnaPol.specificallyBoundValue)', ...
                    d_transcript.boundTranscriptChromosome');
            end
            
            %total polypeptides, ribosomes
            rib = sim.state('Ribosome');
            d1_rib = daughter1.state('Ribosome');
            d2_rib = daughter2.state('Ribosome');
            
            tmp1 = [rib.states rib.boundMRNAs rib.mRNAPositions rib.tmRNAPositions];
            tmp2 = [
                d1_rib.states d1_rib.boundMRNAs d1_rib.mRNAPositions d1_rib.tmRNAPositions
                d2_rib.states d2_rib.boundMRNAs d2_rib.mRNAPositions d2_rib.tmRNAPositions
                ];
            tmp1 = tmp1(tmp1(:, 1) ~= rib.notExistValue, :);
            tmp2 = tmp2(tmp2(:, 1) ~= rib.notExistValue, :);
            tmp1 = sortrows(tmp1, 1:4);
            tmp2 = sortrows(tmp2, 1:4);
            assertEqual(tmp1, tmp2)
            
            assertEqual(d1_pc.counts(pc.boundIndexs(pc.ribosome70SIndexs), c.cytosolIndexs), d1_rib.nActive + d1_rib.nStalled);
            assertEqual(d2_pc.counts(pc.boundIndexs(pc.ribosome70SIndexs), c.cytosolIndexs), d2_rib.nActive + d2_rib.nStalled);
            
            assertEqual(rib.nStalled, r.counts(r.boundIndexs(trl.tmRNAGlobalIndex), c.cytosolIndexs));
            assertEqual(d1_rib.nStalled, d1_r.counts(r.boundIndexs(trl.tmRNAGlobalIndex), c.cytosolIndexs));
            assertEqual(d2_rib.nStalled, d2_r.counts(r.boundIndexs(trl.tmRNAGlobalIndex), c.cytosolIndexs));
            
            %chromosome
            for i = 1:numel(chr.stateNames)
                sName = chr.stateNames{i};
                if isequal(sName, 'segregated')
                    continue;
                end
                assertEqual(chr.(sName), ...
                    [daughter1.state('Chromosome').(sName)(:, [1 2]) ...
                    daughter2.state('Chromosome').(sName)(:, [1 2])]);
            end
            assertEqual(false, daughter1.state('Chromosome').segregated)
            assertEqual(false, daughter2.state('Chromosome').segregated)
            
            %stimulus
            assertEqual(sim.state('Stimulus').values, daughter1.state('Stimulus').values)
            assertEqual(sim.state('Stimulus').values, daughter2.state('Stimulus').values)
            
            %geometry
            assertElementsAlmostEqual(sim.state('Geometry').volume, ...
                + daughter1.state('Geometry').volume ...
                + daughter2.state('Geometry').volume, ...
                'relative', 1e-24);
            
            %FtsZ ring
            assertEqual(0, daughter1.state('FtsZRing').numEdgesOneStraight);
            assertEqual(0, daughter1.state('FtsZRing').numEdgesTwoStraight);
            assertEqual(0, daughter1.state('FtsZRing').numEdgesTwoBent);
            assertEqual(0, daughter1.state('FtsZRing').numResidualBent);
            
            assertEqual(0, daughter2.state('FtsZRing').numEdgesOneStraight);
            assertEqual(0, daughter2.state('FtsZRing').numEdgesTwoStraight);
            assertEqual(0, daughter2.state('FtsZRing').numEdgesTwoBent);
            assertEqual(0, daughter2.state('FtsZRing').numResidualBent);
               
            %host
            sNames = sim.state('Host').stateNames;
            for i = 1:numel(sNames)
                assertEqual(sim.state('Host').(sNames{i}), daughter1.state('Host').(sNames{i}));
            end
            
            %time
            assertEqual(0, daughter1.state('Time').values)
            assertEqual(0, daughter2.state('Time').values)
        end
        
        function testInitializeFromDaughter(~)
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            
            %sim = CachedSimulationObjectUtil.load();
            tmp = load('src_test/+edu/+stanford/+covert/+cell/+sim/fixtures/Simulation_EndOfCellCycle.mat');
            sim = tmp.simulation;
            
            %divide cell
            [~, state_daughters] = sim.divideState();
            
            %setup daughter cell simulation
            tmpDirName = 'output/runMediumTests/initializeFromDaughter';
            if exist(tmpDirName, 'dir')
                rmdir(tmpDirName, 's')
            end
            
            mkdir(tmpDirName)
            
            state_daughter1 = state_daughters(1);
            save([tmpDirName filesep 'initialConditions.mat'], '-struct', 'state_daughter1');
            copyfile('src_test/+edu/+stanford/+covert/+cell/+sim/fixtures/initializeFromDaughter.xml', ...
                [tmpDirName filesep 'conditions.xml']);
            
            %run simulation
            daughter = runSimulation('outDir', tmpDirName, 'logToDisk', true);
            
            %test options, parameters copied correctly
            sim.setForTest('lengthSec', 100);
            sim.setForTest('seed', daughter.seed);
            sim.setForTest('media', reshape(sim.media, [], 6));
            sim.setForTest('stimulus', reshape(sim.stimulus, [], 6));
            sim.setForTest('geneticKnockouts', reshape(sim.geneticKnockouts, [], 1));
            
            for i = 1:numel(sim.states)
                s = sim.states{i};
                s.seed = daughter.states{i}.seed;
            end
            for i = 1:numel(sim.processes)
                p = sim.processes{i};
                p.seed = daughter.processes{i}.seed;
            end
            
            assertEqual(sim.getOptions(), daughter.getOptions());
            assertEqual(sim.getParameters(), daughter.getParameters());
            
            %test state copied correctly
            tmp = load([tmpDirName filesep 'state-0.mat']);
            assertEqual(state_daughter1, tmp)

            %clean up
            rmdir(tmpDirName, 's')
        end
    end
end
