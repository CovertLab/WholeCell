%Analysis_Test
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/12/2010
classdef Analysis_Test < TestCase
    properties
        simulation
    end
    
    methods
        function this = Analysis_Test(name)
            this = this@TestCase(name);
        end
    end
    
    methods
        function setUp(this)
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], true);
            sim.applyOptions('lengthSec', 10);
            this.simulation = sim;
        end
    end
    
    methods
        function testBiomassCompositionProductionAnalysis(this)
            %check for no runtime errors
            edu.stanford.covert.cell.sim.analysis.BiomassCompositionProduction.run(this.simulation, ...
                'output/runAnalysisTests/BiomassCompositionProduction');
        end
        
        function disabled_testCellGeometryAnimation(~)
            %check for no runtime errors
            edu.stanford.covert.cell.sim.analysis.CellGeometryAnimation().run();
        end
        
        function testCellOverviewAnalysis(~)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            %check for no runtime errors
            simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            simIdxs = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            simIdxs = simIdxs(1:min(2, end));
            edu.stanford.covert.cell.sim.analysis.CellOverview.run(simBatchDir, 'output/runAnalysisTests/CellOverview', simIdxs);
        end
        
        function testCellStateAnalysis(this)
            %check for no runtime errors
            edu.stanford.covert.cell.sim.analysis.CellState.run(this.simulation, 'output/runAnalysisTests/CellState');
        end
        
        function testChromosomeAnimationAnalysis(~)
            %check for no runtime errors
            edu.stanford.covert.cell.sim.analysis.ChromosomeAnimation(...
                [], [], ['output' filesep 'runAnalysisTests' filesep 'ChromosomeAnimation.avi'], 1000, [], [], 1).run();
        end
        
        function testChromosomePositionHistogram(~)
            %check for no runtime errors
            edu.stanford.covert.cell.sim.analysis.ChromosomePositionHistogram.run(...
                [], 'output/runAnalysisTests/ChromosomePositionHistogram');
        end
        
        function testChromosomeSpaceTimePlot(~)
            %check for no runtime errors
            edu.stanford.covert.cell.sim.analysis.ChromosomeSpaceTimePlot.run(...
                [], [], 'output/runAnalysisTests/ChromosomeSpaceTimePlot');
        end
        
        function testConstantsAnalysis(this)
            %check for no runtime errors
            edu.stanford.covert.cell.sim.analysis.Constants.run(this.simulation, 'output/runAnalysisTests/Constants');
        end
        
        function testDNADamageAnalysis(this)
            %check for no runtime errors
            edu.stanford.covert.cell.sim.analysis.DNADamage.run(this.simulation, 'output/runAnalysisTests/DNADamage');
        end
        
        function testFBAAnalysis(this)
            %check for no runtime errors
            edu.stanford.covert.cell.sim.analysis.FBA.run(this.simulation, 'output/runAnalysisTests/FBA');
        end
        
        function testFlipbookAnimationAnalysis(~)
            import edu.stanford.covert.cell.sim.analysis.FlipbookAnimation;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            %find recent wild-type simulation
            simGroup = SimulationDiskUtil.getLatestWildTypeSimulationGroup();
            simIdxs = 1;
            
            %check for no runtime errors
            animator = FlipbookAnimation(...
                simGroup, simIdxs, ['output' filesep 'runAnalysisTests' filesep 'FlipbookAnimation.avi'], 1000, [], [], 1);
            animator.run([], true, true, false, true);
        end
        
        function testMinimalGenes(~)
            import edu.stanford.covert.cell.sim.analysis.MinimalGenes;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            %% analyze simple network
            %--> A
            %--> B
            %g1, e1: A --> C
            %g2, e2: A --> C
            %g3, e3: A --> C
            %g4+g5, e4: B --> D
            %g4+g6, e5: B --> D
            %C + D --> Biomass
            %Biomass -->
            
            nEnz = 5;
            nGene = 6;
            nMet = 5;
            nRxn = 9;
            S = [
                1  0  -1  -1  -1   0    0   0   0
                0  1   0   0   0  -1   -1   0   0
                0  0   1   1   1   0    0  -1   0
                0  0   0   0   0   1    1  -1   0
                0  0   0   0   0   0    0   1  -1
                ];
            growthRxnIdx = nRxn - 1;
            
            subGeneComp = zeros(nGene, nMet);
            
            enzGeneComp = zeros(nGene, nEnz);
            enzGeneComp(1, 1) = 1;
            enzGeneComp(2, 2) = 1;
            enzGeneComp(3, 3) = 1;
            enzGeneComp([4 5], 4) = 2;
            enzGeneComp([4 6], 5) = 2;
            
            rxnCatMat = zeros(nRxn, nEnz);
            rxnCatMat(3, 1) = 1;
            rxnCatMat(4, 2) = 1;
            rxnCatMat(5, 3) = 1;
            rxnCatMat(6, 4) = 1;
            rxnCatMat(7, 5) = 1;
            
            fbaFluxBounds = [
                0 Inf
                0 Inf
                0 1
                0 1.5
                0 3
                0 50
                0 50
                0 Inf
                0 Inf
                ];
            
            opts = struct('solver', 'cplex');
            
            [growths, incGenes, fluxs] = MinimalGenes.calcMaxGeneSetGrowthRates(...
                S, fbaFluxBounds, subGeneComp, enzGeneComp, rxnCatMat, growthRxnIdx, [], opts);
            
            geneFirstAppearances = zeros(nGene, 1);
            for i = 1:nGene
                geneFirstAppearances(i) = find(growths > 0 & incGenes(:, i), 1, 'first');
            end
            assertEqual(3, geneFirstAppearances(3));
            assertTrue(geneFirstAppearances(1) > geneFirstAppearances(2));
            assertTrue(geneFirstAppearances(2) > geneFirstAppearances(3));
            assertEqual(3, geneFirstAppearances(4));
            assertEqual(geneFirstAppearances(4), min(geneFirstAppearances(5:6)));
            
            [axesHandle, figHandle] = PlotUtil.newAxesHandle();
            plot(axesHandle, 0:nGene, growths)
            saveas(figHandle, 'output/runAnalysisTests/MinimalGenes/ToyNetwork.pdf');
            close(figHandle);
            
            %% full analysis
            MinimalGenes.run('output/runAnalysisTests/MinimalGenes');
        end
        
        function testMultiGenerationsAnalysis(~)
            %% import classes
            import edu.stanford.covert.cell.sim.analysis.MultiGenerations;
           
            %% test ancestry
            assertEqual([
                1  0  1 NaN  3 4
                2  0  2 NaN  5 6
                3  1  1  1  7 8
                4  1  1  1  9 10
                5  1  2  2  11 12
                6  1  2  2  13 14 
                7  2  1  3  NaN NaN
                8  2  1  3  NaN NaN
                9  2  1  4  NaN NaN
                10 2  1  4  NaN NaN
                11 2  2  5  NaN NaN
                12 2  2  5  NaN NaN
                13 2  2  6  NaN NaN
                14 2  2  6  NaN NaN
                ], MultiGenerations.calcAncestry(3, 2));
            
            relations = cat(3, [
                0   1  1  2  2  2  2 3 3 3 3 3 3 3 3
                -1  0  0  1  1  1  1 2 2 2 2 2 2 2 2
                -1  0  0  1  1  1  1 2 2 2 2 2 2 2 2
                -2 -1 -1  0  0  0  0 1 1 1 1 1 1 1 1
                -2 -1 -1  0  0  0  0 1 1 1 1 1 1 1 1
                -2 -1 -1  0  0  0  0 1 1 1 1 1 1 1 1
                -2 -1 -1  0  0  0  0 1 1 1 1 1 1 1 1
                -3 -2 -2 -1 -1 -1 -1 0 0 0 0 0 0 0 0
                -3 -2 -2 -1 -1 -1 -1 0 0 0 0 0 0 0 0
                -3 -2 -2 -1 -1 -1 -1 0 0 0 0 0 0 0 0
                -3 -2 -2 -1 -1 -1 -1 0 0 0 0 0 0 0 0
                -3 -2 -2 -1 -1 -1 -1 0 0 0 0 0 0 0 0
                -3 -2 -2 -1 -1 -1 -1 0 0 0 0 0 0 0 0
                -3 -2 -2 -1 -1 -1 -1 0 0 0 0 0 0 0 0
                -3 -2 -2 -1 -1 -1 -1 0 0 0 0 0 0 0 0
                ], [
                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                0 0 1 0 0 1 1 0 0 0 0 1 1 1 1
                0 1 0 1 1 0 0 1 1 1 1 0 0 0 0 
                0 0 1 0 1 2 2 0 0 1 1 2 2 2 2
                0 0 1 1 0 2 2 1 1 0 0 2 2 2 2
                0 1 0 2 2 0 1 2 2 2 2 0 0 1 1
                0 1 0 2 2 1 0 2 2 2 2 1 1 0 0
                0 0 1 0 1 2 2 0 1 2 2 3 3 3 3
                0 0 1 0 1 2 2 1 0 2 2 3 3 3 3
                0 0 1 1 0 2 2 2 2 0 1 3 3 3 3
                0 0 1 1 0 2 2 2 2 1 0 3 3 3 3
                0 1 0 2 2 0 1 3 3 3 3 0 1 2 2
                0 1 0 2 2 0 1 3 3 3 3 1 0 2 2
                0 1 0 2 2 1 0 3 3 3 3 2 2 0 1
                0 1 0 2 2 1 0 3 3 3 3 2 2 1 0
                ]);
            tmp = MultiGenerations.calcAncestryRelations(4, 1);
            assertEqual(relations(:, :, 1), tmp(:, :, 1));
            assertEqual(relations(:, :, 2), tmp(:, :, 2));
            
            %% analysis
            MultiGenerations.run('output/runAnalysisTests/MultiGenerations');
        end
        
        function testOverviewAnimationAnalysis(~)
            import edu.stanford.covert.cell.sim.analysis.OverviewAnimation;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            %find recent wild-type simulation
            simGroup = SimulationDiskUtil.getLatestWildTypeSimulationGroup();
            simIdxs = 1:min(3, SimulationDiskUtil.getNumSimulations(simGroup));
            
            %check for no runtime errors
            animator = OverviewAnimation(...
                simGroup, simIdxs, ['output' filesep 'runAnalysisTests' filesep 'OverviewAnimation.avi'], 1000, [], [], 1);
            animator.run();
        end
        
        function disabled_testReplicationAnalysis(~)
            %check for no runtime errors
            edu.stanford.covert.cell.sim.analysis.Replication.run('output/runAnalysisTests/Replication');
        end
        
        function testReplicationInitiationDurationAnalysis(~)
            %check for no runtime errors
            edu.stanford.covert.cell.sim.analysis.ReplicationInitiationDuration.run(...
                'output/runSmallTests/ReplicationInitiation', ...
                'output/runAnalysisTests/ReplicationInitiation');
        end
        
        function testRNAAnalysis(this)
            %check for no runtime errors
            edu.stanford.covert.cell.sim.analysis.RNA.run(this.simulation, 3, 1000, 'output/runAnalysisTests/RNA');
        end
        
        function testPopulationAnalysis(~)
			import edu.stanford.covert.cell.sim.analysis.Population;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            %check for no runtime errors
            simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            selectedSims = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            selectedSims = selectedSims(1:min(end, 2));
            Population.run(simBatchDir, 'output/runAnalysisTests/Population', selectedSims);
        end
        
        function testProcessMetaboliteUsage(~)
            import edu.stanford.covert.cell.sim.analysis.ProcessMetaboliteUsage;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            %check for no runtime errors
            simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            selectedSims = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            selectedSims = selectedSims(1:min(end, 2));
            ProcessMetaboliteUsage.run(simBatchDir, selectedSims, 'output/runAnalysisTests/ProcessMetaboliteUsage');
        end
        
        function testRNAExpressionAnalysis(this)
            import edu.stanford.covert.cell.sim.analysis.RNAExpression;
            
            RNAExpression.run(this.simulation, 'output/runAnalysisTests/RNAExpression');
        end
        
		%disabled because test takes too long
        function disabled_testSingleGeneDeletionsAnalysis(~)
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            %check for no runtime errors
            simBatchDir = SimulationDiskUtil.getLatestSingleGeneDeletionSimulationGroup();
            SingleGeneDeletions.run(simBatchDir, false);
            SingleGeneDeletions.run(simBatchDir, true);
        end
        
        function testSingleCellAnalysis(~)
            %check for no runtime errors
            edu.stanford.covert.cell.sim.analysis.SingleCell.run([], [], 'output/runAnalysisTests/SingleCell');
        end
        
        function testSingleCellVariationAnalysis(this)
            %check for no runtime errors
            edu.stanford.covert.cell.sim.analysis.SingleCellVariation.run(this.simulation, 'output/runAnalysisTests/SingleCellVariation');
        end
        
        function testSimulationStructureAnalysis(this)
            %check for no runtime errors
            edu.stanford.covert.cell.sim.analysis.SimulationStructure.run(this.simulation, 'output/runAnalysisTests/SimulationStructure');
        end
        
        function disabled_testTimeCourseAnalysis(this)
            %check for no runtime errors
            edu.stanford.covert.cell.sim.analysis.TimeCourse.run(this.simulation, 'output/runAnalysisTests/TimeCourse');
        end

	function testTranslationAnalysis(~)
            %check for no runtime errors
            edu.stanford.covert.cell.sim.analysis.TranslationAnalysis.run([], 'output/runAnalysisTests/TranslationAnalysis');
        end
    end
end
