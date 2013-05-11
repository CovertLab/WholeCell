% Poster Figures
% Make figures for posters
%
% Author: Derek Macklin, macklin@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Created: 9/13/2011
classdef PosterFigures
    methods (Static = true);
        function run(eventName, outputDir)
            import edu.stanford.covert.cell.sim.analysis.PosterFigures
            
            PosterFigures.(eventName)(outputDir);
        end
        
        function pioneer2011(outputDir)
            import edu.stanford.covert.cell.sim.analysis.SingleCell;
            import edu.stanford.covert.cell.sim.analysis.ChromosomeSpaceTimePlot;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            import edu.stanford.covert.cell.sim.analysis.ProcessMetaboliteUsage;
            
            simBatchDir = '2011_09_09_19_17_51';
            simNum = 1;
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep num2str(simNum)]);
            
            %% Gene -> RNA -> Protein -> Growth Correlation
            [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            SingleCell.growth(sim, simBatchDir, simNum, figHandle)
            saveas(figHandle, [outputDir filesep 'Growth.eps'], 'epsc');
            close(figHandle);
            
            %% Chromosomal stuff
            [simDir, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep num2str(simNum)]);
            c = sim.state('Chromosome');
            
            simIdx = num2str(SimulationDiskUtil.getSimulationIndex(simDir));
            simTimeStamp = SimulationDiskUtil.getSimulationTimeStamp(simDir);
            stateNames = {'ProteinComplex' 'counts'
                'Chromosome' 'complexBoundSites'
                'Chromosome' 'monomerBoundSites'
                'Chromosome' 'polymerizedRegions'
                'RNAPolymerase' 'positionStrands'
                'RNAPolymerase' 'states'
                'RNAPolymerase' 'transcriptionFactorBindingProbFoldChange'
                'Time' 'values'};
            states = DiskLogger.load(simDir, stateNames, [], [], [], 'extract');
            
            height = 1000;
            width = 1000;
            
            RNAPSparseMat = ChromosomeSpaceTimePlot.makeRNAPSparseMat(states, sim);
            
            [ax, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            [densityBEMatrix, maxTime, maxSpace] = ChromosomeSpaceTimePlot.makeSpaceTimeDensityPlots(...
                ax, sim, states, RNAPSparseMat, ...
                ChromosomeSpaceTimePlot.chromProteins, [], width, height, simTimeStamp, simIdx);
            close(figHandle);
            [ax, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            set(figHandle, 'Color', [1 1 1]);
            ChromosomeSpaceTimePlot.plotCircularDensity(...
                ax, densityBEMatrix, maxTime, sim, false, ...
                sprintf('Circular Chromosome Density Plot\nSimulation: %s #%s', simTimeStamp, simIdx));
            set(ax, 'Position', [0.2 0.2 0.6 0.6])
            saveas(figHandle, [outputDir filesep 'CircularDensity.eps'], 'epsc');
            close(figHandle);
            
            selectedTime = 4 * 3600;
            [ax, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            ChromosomeSpaceTimePlot.makeRingPlot(...
                ax, sim, states, RNAPSparseMat, ...
                ChromosomeSpaceTimePlot.chromProteins, [], selectedTime, ...
                maxSpace);
            saveas(figHandle, [outputDir filesep 'RingPlot.eps'], 'epsc');
            close(figHandle);

            
            [ax, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            whatToShow = {'Topoisomerase IV', 'DnaA Complex', 'LuxR', 'HTH Regulator', 'Ferric uptake repressor', 'Helicase'};
            ChromosomeSpaceTimePlot.plotSpaceTimeOverlay(...
                ax, sim, states, RNAPSparseMat, ...
                ChromosomeSpaceTimePlot.chromProteins, whatToShow, ...
                sprintf('Chromosome 1\nSimulation: %s #%s', simTimeStamp, simIdx));
            set(figHandle, 'Renderer', 'Painters')
            saveas(figHandle, [outputDir filesep 'SpaceTimeOverlay.eps'], 'epsc');
            close(figHandle);
            clear('states', 'stateNames', 'RNAPSparseMat', 'whatToShow', 'height', 'width');
%             
            %% Metabolite Usages
            wIDsToShow = {'ATP', 'GTP'};
            desiredInfo = 'Usages';
            [ax, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            ProcessMetaboliteUsage.processSingleCellMetaboliteHandler(wIDsToShow, 'Usages', figHandle, '2011_09_09_19_17_51', simNum)
            saveas(figHandle, [outputDir filesep 'ProcessUsages.eps'], 'epsc');
        end
    end
end