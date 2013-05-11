% Figures 3-4
%
% Author: Derek Macklin, macklin@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Created: 1/12/2012
classdef Figures34
    methods (Static = true)
        function run(outDirectory)
            import edu.stanford.covert.cell.sim.analysis.Figures34;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            %output directory
            if nargin < 1
                outDirectory = 'documentation/paper/figures';
            end
            
            %select simulation
            simBatchDir = [SimulationDiskUtil.getBaseDir() filesep '2011_10_25_01_45_44' filesep];
            selectedSim = 10;
            
            %setup linux path
            if isunix()
                a = getenv('LD_LIBRARY_PATH');
                a = ['/usr/lib:' a];
                setenv('LD_LIBRARY_PATH', a);
            end
            
            % figure 2: fitting, validation
            Figures34.figure2(simBatchDir, selectedSim, outDirectory)
            
            % figure 3: chromosome
            Figures34.figure3(simBatchDir, selectedSim, outDirectory)
            
            % figure 4: cell cycle
            Figures34.figure4(simBatchDir, selectedSim, outDirectory)
            
            % figure 5: energy
            Figures34.figure5(simBatchDir, selectedSim, outDirectory)
            
            %growth sensitivity
            Figures34.growthSensitivity(simBatchDir, selectedSim, outDirectory)
        end
        
        %fitting, validation
        function figure2(simBatchDir, selectedSim, outDirectory)
            import edu.stanford.covert.cell.sim.analysis.Figures34;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            outDir = [outDirectory filesep 'fittingAndValidation' filesep];
            
            %% figure
            w = 11.4; h = 17; % Cell 1.5 column
            [~, figHandle] = PlotUtil.newAxesHandle();
            clf(figHandle);
            set(figHandle, 'PaperUnits', 'centimeters');
            set(figHandle, 'PaperSize', [w h]);
            set(figHandle, 'PaperPositionMode', 'manual');
            set(figHandle, 'PaperPosition', [0 0 get(figHandle, 'PaperSize')]);
            set(figHandle, 'Units', 'normalized');
            
            %A: growth measurement
            axesPosA = [0.10 0.82 0.415 0.165];
            PlotUtil.labelSubFigure('A', [0.128 0.992 -0.0286 -0.0208], figHandle);
            Figures34.growthMeasurement(figHandle, axesPosA, outDir);
            
            %B: simulated growth
            PlotUtil.labelSubFigure('B', [0.607 0.992 -0.0286 -0.0208], figHandle);
            Figures34.populationMass(simBatchDir, selectedSim, figHandle, [0.585 0.82 0.415 0.165], outDir);
            
            %C: biomass composition
            axesPosC = [0.10 0.595 0.415 0.165];
            PlotUtil.labelSubFigure('C', [0.128 0.772 -0.0286 -0.0208], figHandle);
            Figures34.dryBiomassComposition(simBatchDir, figHandle, axesPosC, outDir);
            
            %D: mass fraction dynamics
            PlotUtil.labelSubFigure('D', [0.607 0.772 -0.0286 -0.0208], figHandle);
            Figures34.normCellMass(simBatchDir, selectedSim, outDir, figHandle, [0.577 0.595 0.423 0.165]);
            
            %E: metabolic map
            axesPosE = [0.07 0.29 0.50 0.215];
            PlotUtil.labelSubFigure('E', [0.135 0.527 -0.0286 -0.0208], figHandle);
            Figures34.metabolicMap(simBatchDir, selectedSim, outDir, figHandle, axesPosE);
            
            %F: ratio of predicted to observed metabolite concentrations
            axesPosF = [0.11 0.057 0.58 0.18];
            PlotUtil.labelSubFigure('F', [0.135 0.255 -0.0286 -0.0208], figHandle);
            Figures34.metaboliteConcentrations(simBatchDir, selectedSim, outDir, figHandle, axesPosF)
            
            %G: mRNA, protein expression dynamics
            PlotUtil.labelSubFigure('G', [0.684 0.527 -0.0286 -0.0208], figHandle);
            Figures34.proteinBursting(simBatchDir, selectedSim, outDir, figHandle, [0.645 0.29 0.35 0.215]);
            
            %H: scatter plot of RNA vs protein expression with marginals
            PlotUtil.labelSubFigure('H', [0.799 0.255 -0.0286 -0.0208], figHandle);
            Figures34.RNAvProtein(simBatchDir, outDir, figHandle, [0.765 0.057 0.23 0.18]);
            
            % area labels
            axesHandleA = findobj(figHandle, 'position', axesPosA);
            axesHandleE = findobj(figHandle, 'position', axesPosE);
            text(-0.212, 0.5 * (axesPosC(2)+axesPosC(4)-axesPosA(2)) / axesPosA(4), 'Training Data', ...
                'FontSize', 10, 'rotation', 90, 'Parent', axesHandleA, ...
                'units', 'normalized', 'HorizontalAlign', 'center', 'VerticalAlign', 'middle', ...
                'FontWeight', 'bold');
            text((-0.212 - (axesPosE(1) - axesPosA(1))/axesPosA(3)) * axesPosA(3) / axesPosE(3), ...
                0.5 * (axesPosF(2)+axesPosF(4)-axesPosE(2)) / axesPosE(4), 'Independent Validation', ...
                'FontSize', 10, 'rotation', 90, 'Parent', axesHandleE, ...
                'units', 'normalized', 'HorizontalAlign', 'center', 'VerticalAlign', 'middle', ...
                'FontWeight', 'bold');
            annotation(figHandle, 'line', 0.035 * [1 1], [0.570 0.997], 'Color', 'k', 'LineWidth', 1)
            annotation(figHandle, 'line', 0.035 * [1 1], [0.003 0.530], 'Color', 'k', 'LineWidth', 1)
            annotation(figHandle, 'line', [0.035 0.05], 0.570 * [1 1], 'Color', 'k', 'LineWidth', 1)
            annotation(figHandle, 'line', [0.035 0.05], 0.997 * [1 1], 'Color', 'k', 'LineWidth', 1)
            annotation(figHandle, 'line', [0.035 0.05], 0.530 * [1 1], 'Color', 'k', 'LineWidth', 1)
            annotation(figHandle, 'line', [0.035 0.05], 0.003 * [1 1], 'Color', 'k', 'LineWidth', 1)
            
            %standarize tick lengths
            axesHandles = findobj(figHandle, 'type', 'axes', '-and', 'tag', '');
            axesPos = cell2mat(get(axesHandles, 'position'));
            tickLen = 0.01 * max(max(axesPos(:, 3:4), [], 1) .* [w h]);
            for i = 1:numel(axesHandles)
                maxVal = max(axesPos(i, 3:4) .* [w h]);
                set(axesHandles(i), 'TickLen', [tickLen / maxVal  0.0250]);
            end
            
            %save
            print(figHandle, [outDirectory filesep 'figure2.pdf'], '-dpdf', '-rgb');
            close(figHandle);
            return
            
            %% supplemental figures
            
            %gene expression
            w = 11.4; h = w * (1-0.09)/(1-0.075);
            [~, figHandle] = PlotUtil.newAxesHandle();
            set(figHandle, 'PaperUnits', 'centimeters');
            set(figHandle, 'PaperSize', [w h]);
            set(figHandle, 'PaperPositionMode', 'manual');
            set(figHandle, 'PaperPosition', [0 0 get(figHandle, 'PaperSize')]);
            set(figHandle, 'Units', 'normalized');
            figPos = get(figHandle, 'Position');
            set(figHandle, 'Position', [figPos(1:2) figPos(3) figPos(3) * h / w]);
            
            clf(figHandle);
            [exptGeneExpr, modelGeneExpr] = Figures34.geneExpression(simBatchDir, selectedSim, figHandle, [0.09 0.075 1-0.09 1-0.075], true);
            print(figHandle, [outDirectory filesep '../supplement/figures/rnaExpressionComparison_log.pdf'], '-dpdf', '-rgb');
            
            clf(figHandle);
            [exptGeneExpr, modelGeneExpr] = Figures34.geneExpression(simBatchDir, selectedSim, figHandle, [0.082 0.065 1-0.082 1-0.065], false);
            print(figHandle, [outDirectory filesep 'figureS1a.pdf'], '-dpdf', '-rgb');
            
            close(figHandle);
            
            %gene expression ratios
            [axesHandle, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            Figures34.geneExpressionRatios(exptGeneExpr, modelGeneExpr, axesHandle);
            print(figHandle, [outDir 'RNAExprRatios.pdf'], '-dpdf', '-rgb');
            close(figHandle);
            
            %dynamics of gene, RNA, protein copy number; reaction flux, growth rate
            [~, figHandle] = PlotUtil.newAxesHandle();
            Figures34.geneToPhenotype(simBatchDir, selectedSim, outDir, figHandle, [0.1 0.1 0.8 0.8]);
            print(figHandle, [outDir 'geneToPhenotype.pdf'], '-dpdf', '-rgb');
            close(figHandle);
            
            %metabolic map
            [~, figHandle] = PlotUtil.newAxesHandle();
            Figures34.metabolicMap(simBatchDir, selectedSim, outDir, figHandle, [0.1 0.2 0.8 0.6], true);
            print(figHandle, [outDirectory filesep 'figureS1b.pdf'], '-dpdf', '-rgb');
            close(figHandle);
        end
        
        %chromosome
        function figure3(simBatchDir, selectedSim, outDirectory)
            import edu.stanford.covert.cell.sim.analysis.Figures34;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            outDir = [outDirectory filesep 'chromosome' filesep];
            
            w = 11.4; h = 16.1; % Cell 1.5 column
            h2 = 0.2025;
            [~, figHandle] = PlotUtil.newAxesHandle([w h]);
            clf(figHandle);
            
            %A: DNA-bound protein density
            PlotUtil.labelSubFigure('A', [0.080 0.997 -0.0286 -0.0208], figHandle);
            Figures34.circosPlot(simBatchDir, selectedSim, outDir, figHandle, [0.00 0.66 0.56 0.3]);
            
            %B: time to DNA occupancy
            PlotUtil.labelSubFigure('B', [0.605 0.997 -0.0286 -0.0208], figHandle);
            Figures34.chromosomeSampling(simBatchDir, selectedSim, outDir, figHandle, [0.595 0.86 0.40 0.13]);
            
            %C: time to RNA expression
            PlotUtil.labelSubFigure('C', [0.605 0.817 -0.0286 -0.0208], figHandle);
            Figures34.rnaSampling(simBatchDir, selectedSim, outDir, figHandle, [0.595 0.68 0.40 0.13]);
            
            %D: space-time plot of DNA-bound protein
            PlotUtil.labelSubFigure('D', [0.080 0.647 -0.0286 -0.0208], figHandle);
            Figures34.spaceTimePlot(simBatchDir, selectedSim, outDir, figHandle, [0.06 0.36 0.925 0.28]);
            
            %E-G: collisions
            PlotUtil.labelSubFigure('E', [0.080 0.325 -0.0286 -0.0208], figHandle);
            PlotUtil.labelSubFigure('F', [0.633 0.325 -0.0286 -0.0208], figHandle);
            [posFreq, protDensity] = Figures34.analyzeDNABoundProteinDisplacements(simBatchDir, outDir, figHandle, ...
                [0.095 0.07 0.94 0.242]);
            
            %standarize tick lengths
            axesHandles = findobj(figHandle, 'type', 'axes', '-and', 'tag', '');
            axesPos = cell2mat(get(axesHandles, 'position'));
            tickLen = 0.0075 * max(max(axesPos(:, 3:4), [], 1) .* [w h]);
            for i = 1:numel(axesHandles)
                maxVal = max(axesPos(i, 3:4) .* [w h]);
                set(axesHandles(i), 'TickLen', [tickLen / maxVal  0.0250]);
            end
            
            %save
            print(figHandle, [outDirectory filesep 'figure3.pdf'], '-dpdf', '-rgb');
            close(figHandle);
            
            %% supplement
            %fit
            [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            clf(figHandle);
            Figures34.replicationInitiationVsReplicationDurationFit(simBatchDir, selectedSim, figHandle, [0.1 0.1 0.8 0.8]);
            print(figHandle, [outDir filesep 'ReplicationInitiationVsReplicationDurationFit.pdf'], '-dpdf', '-rgb');
            close(figHandle);
            
            %analysis of RNA polymerase exploration
            [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            clf(figHandle);
            Figures34.plotRNAPolSampling(simBatchDir, selectedSim, outDir, figHandle, [0.1 0.1 0.8 0.8]);
            print(figHandle, [outDir filesep 'RNAPolSampling.pdf'], '-dpdf', '-rgb');
            close(figHandle);
            
            %F: scatter plot of density vs collisions
            PlotUtil.labelSubFigure('F', [0.715 0.80 -0.0286 -0.0690], figHandle);
            Figures34.dnaBoundProteinDensityVsDisplacement(simBatchDir, selectedSim, outDir, figHandle, [0.71 0.755-h2 0.29 h2]);
            
            %more supplemental material            
            Figures34.cellShape(outDir, simBatchDir, selectedSim);
            Figures34.mrnaSampling(simBatchDir, selectedSim, outDir, figHandle, [0.71 0.755-0.12 0.29 0.12]);
        end
        
        %cell cycle
        function figure4(simBatchDir, selectedSim, outDirectory)
            import edu.stanford.covert.cell.sim.analysis.Figures34;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            outDir = [outDirectory filesep 'cellCycle' filesep];
            
            w = 11.4; h = 18; % Cell 1.5 column
            h2 = 0.2025;
            [~, figHandle] = PlotUtil.newAxesHandle([w h]);
            clf(figHandle);
            
            %F: histograms of cell cycle phase lengths
            PlotUtil.labelSubFigure('F', [0.420 0.36 -0.0286 -0.0208], figHandle, 8, 'normalized');
            Figures34.cellCycleLength(simBatchDir, selectedSim, figHandle, [0.385 0.16 0.27 0.18]);
            
            %G: scatter plot of replication vs. replication initiation durations
            PlotUtil.labelSubFigure('G', [0.750 0.36 -0.0286 -0.0208], figHandle, 8, 'normalized');
            Figures34.replicationInitiationVsReplicationDuration(simBatchDir, selectedSim, figHandle, [0.72 0.16 0.27 0.18]);
            
            %standarize tick lengths
            axesHandles = findobj(figHandle, 'type', 'axes', '-and', 'tag', '');
            axesPos = cell2mat(get(axesHandles, 'position'));
            tickLen = 0.0075 * max(max(axesPos(:, 3:4), [], 1) .* [w h]);
            for i = 1:numel(axesHandles)
                maxVal = max(axesPos(i, 3:4) .* [w h]);
                set(axesHandles(i), 'TickLen', [tickLen / maxVal  0.0250]);
            end
            
            %save
            print(figHandle, [outDirectory filesep 'figure4.pdf'], '-dpdf', '-rgb');
            close(figHandle);
        end
        
        %energy
        function figure5(simBatchDir, selectedSim, outDirectory, plotUsageRate)
            import edu.stanford.covert.cell.sim.analysis.Figures34;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SummaryLogger;
            import edu.stanford.covert.util.ConstantUtil;
            
            if nargin < 4
                plotUsageRate = false;
            end
            
            outDir = [outDirectory filesep 'energy' filesep];
            
            %% get data
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            c = sim.compartment;
            met = sim.state('Metabolite');
            metab = sim.process('Metabolism');
            
            nSims = SimulationDiskUtil.getNumSimulations(simBatchDir);
            
            if ~exist([outDir 'data.mat'], 'file')
                stateNames = {
                    'atpUsage' 'ATP'
                    'gtpUsage' 'GTP'
                    };
                ensemble = SimulationEnsemble(simBatchDir, stateNames, [], 1:nSims);
                simEndTimes = ensemble.stateData.simulationEndTimes;
                
                totalEnergyUsagePerCell = permute(nansum([
                    ensemble.stateData.values(ensemble.getPropertyIndices('atpUsage'), :, :, :)
                    ensemble.stateData.values(ensemble.getPropertyIndices('gtpUsage'), :, :, :)
                    ], 3), [1 4 2 3]) / ConstantUtil.nAvogadro * 1e3;
                
                metIdxs = sub2ind([numel(met.wholeCellModelIDs) numel(c.wholeCellModelIDs)], ...
                    met.getIndexs({'ATP'; 'GTP'}), c.cytosolIndexs([1 1])');
                stateNames = {
                    'Time'       'values'  ':'  ':'
                    'Geometry'   'volume'  ':'  ':'
                    'Metabolite' 'counts'  met.getIndexs({'ATP'; 'GTP'; 'FAD'; 'NAD'; 'NADP'; 'FADH2'; 'NADH'; 'NADPH'})  c.cytosolIndexs
                    'Metabolite' 'processUsages' metIdxs ':'
                    };
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSim);
                
                simEndTime = states.Time.values(end);
                
                energyCarrierDynamics = full(permute(states.Metabolite.counts, [1 3 2]));
                energyCarrierDynamics(3:5, :) = energyCarrierDynamics(3:5, :) + energyCarrierDynamics(6:8, :);
                energyCarrierDynamics = energyCarrierDynamics(1:5, :) ./ repmat(permute(states.Geometry.volume, [1 3 2]), [5 1]) / ConstantUtil.nAvogadro * 1e3;
                
                processEnergyUsageDynamics = full(states.Metabolite.processUsages);
                
                stateNames = {
                    'Metabolite' 'processUsages' [metIdxs(1) metIdxs(2)] ':'
                    };
                processEnergyUsages = zeros(numel(sim.processes), 1);
                for i = 1:nSims
                    states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', i);
                    processEnergyUsages = processEnergyUsages + full(sum(states.Metabolite.processUsages, 3)');
                end
                
                clear ensemble states;
                                
                metIdxs = sub2ind([numel(met.wholeCellModelIDs) numel(c.wholeCellModelIDs)], ...
                    met.getIndexs({'ATP'; 'CTP'; 'GTP'; 'UTP'; 'AMP'; 'CMP'; 'GMP'; 'UMP'; 'FAD'; 'NAD'; 'NADP'; 'FADH2'; 'NADH'; 'NADPH'}), c.cytosolIndexs(ones(14, 1)));
                stateNames = {
                    'MetabolicReaction' 'growth' ':' ':'
                    'Metabolite' 'counts' met.ntpIndexs' c.cytosolIndexs
                    'Metabolite' 'processUsages' metIdxs [sim.processIndex('Metabolism'); sim.processIndex('Transcription'); sim.processIndex('RNADecay')]
                    };
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSim);
                growthDynamics = full(permute(states.MetabolicReaction.growth, [1 3 2]));
                energyCarrierSynthesis = full(permute(states.Metabolite.processUsages([1 3 end-5:end], 1, :), [1 3 2]));
                energyCarrierSynthesis(3:5, :) = energyCarrierSynthesis(3:5, :) + energyCarrierSynthesis(6:8, :);
                energyCarrierSynthesis = energyCarrierSynthesis(1:5, :);
                energyCarrierSynthesis(1, :) = energyCarrierSynthesis(1, :) - metab.unaccountedEnergyConsumption * growthDynamics;
                
                ntpDynamics = full(permute(states.Metabolite.counts, [1 3 2]));
                
                transcriptionDynamics = full(permute(sum(states.Metabolite.processUsages(1:4, 2, :), 1), [1 3 2]));
                rnaDecayDynamics = -full(permute(sum(states.Metabolite.processUsages(5:8, 3, :), 1), [1 3 2]));
                
                stateNames = {
                    'mass' 'Mass'
                    'ntps' 'NTP'
                    'superhelicity' 'superhelicity'
                    };
                ensemble = SimulationEnsemble(simBatchDir, stateNames, [], 1:nSims);
                ntpCounts = permute(nanmean(...
                    ensemble.stateData.values(ensemble.getPropertyIndices('ntps'), :, :, :), ...
                    3), [1 4 2 3]) / ConstantUtil.nAvogadro * 1e3;
                superhelicity = permute(nanmean(...
                    ensemble.stateData.values(ensemble.getPropertyIndices('superhelicity'), :, :, :), ...
                    3), [1 4 2 3]);
                
                stateNames = {
                    'Mass' 'rnaWt'
                    'Mass' 'proteinWt'
                    'RNAPolymerase' 'nActive'
                    'Ribosome' 'nStalled'
                    'Geometry' 'volume'
                    };
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', 1:nSims);
                rnaWt = permute(nanmean(full(states.Mass.rnaWt(:, c.cytosolIndexs, :, :)), 3), [1 4 3 2]);
                proteinWt = permute(nanmean(full(states.Mass.proteinWt(:, c.cytosolIndexs, :, :)), 3), [1 4 3 2]);
                activeRNAPols =  permute(nanmean(full(states.RNAPolymerase.nActive), 3), [1 4 3 2]);
                stalledRibs = permute(nanmean(full(states.Ribosome.nStalled), 3), [1 4 3 2]);
                
                ntpConcs = permute(nanmean(...
                    ensemble.stateData.values(ensemble.getPropertyIndices('ntps'), :, 2:end, :) ./ ...
                    states.Geometry.volume, ...
                    3), [1 4 2 3]) / ConstantUtil.nAvogadro * 1e3;
                
                simStats = SummaryLogger.getSimulationStatistics(simBatchDir);
                repInitDuration = simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME);
                repDuration = diff(simStats(:, [SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME SummaryLogger.SIM_STATUS_INDEX_REP_TIME]), 1, 2);
                
                save([outDir 'data.mat'], ...
                    'simEndTimes', 'totalEnergyUsagePerCell', 'energyCarrierDynamics', ...
                    'processEnergyUsageDynamics', 'processEnergyUsages', ...
                    'growthDynamics', 'energyCarrierSynthesis', 'ntpDynamics', 'transcriptionDynamics', 'rnaDecayDynamics', ...
                    'ntpCounts', 'superhelicity', ...
                    'rnaWt', 'proteinWt', 'activeRNAPols', 'stalledRibs', 'ntpConcs', ...
                    'repInitDuration', 'repDuration');
            else
                load([outDir 'data.mat']);
            end            
            
            %% figure
            figW = 8.5; figH = 10; % Cell 1 column
            [~, figHandle] = PlotUtil.newAxesHandle([figW figH]);
            clf(figHandle);
            
            leftColX = 0.085;
            leftColW = 0.42;
            leftColLabelX = 0.105;
            
            window = 300;
            
            %% A: dynamics of energy carrier synthesis
            for i = 1:size(energyCarrierSynthesis, 1)
                energyCarrierSynthesis(i, :) = conv(energyCarrierSynthesis(i, :), ones(window, 1) / window, 'same');
            end
            
            PlotUtil.labelSubFigure('A', [leftColLabelX 0.993 -0.0286 -0.0208], figHandle, 8, 'normalized');
            axesHandle = subplot('Position', [leftColX 0.745 leftColW 0.23], 'Parent', figHandle);
            h = plot(axesHandle, (window:simEndTimes(selectedSim)-window+1) / 3600, -energyCarrierSynthesis(:, window:end-window+1) / ConstantUtil.nAvogadro * 1e21);
            set(h(1), 'Color', 'b');
            set(h(2), 'Color', 'g');
            set(h(3), 'Color', 'r')
            set(h(4), 'Color', 'c')
            set(h(5), 'Color', [1 0.5 0])
            
            text(4.1, 20.0, 'ATP',      'FontSize', 6, 'Color', get(h(1), 'Color'), 'HorizontalAlign', 'center');
            text(4.1, 0.50, 'GTP',      'FontSize', 6, 'Color', get(h(2), 'Color'), 'HorizontalAlign', 'center');
            text(4.1, 2e-5, 'FAD(H_2)', 'FontSize', 6, 'Color', get(h(3), 'Color'), 'HorizontalAlign', 'center');
            text(2.5, 5e-3, 'NAD(H)',   'FontSize', 6, 'Color', get(h(4), 'Color'), 'HorizontalAlign', 'center');
            text(5.5, 5e-3, 'NADP(H)',  'FontSize', 6, 'Color', get(h(5), 'Color'), 'HorizontalAlign', 'center');
            
            xlim(axesHandle, [0 simEndTimes(selectedSim) / 3600]);
            xlabel(axesHandle, 'Time (h)', 'FontSize', 7);
            ylabel(axesHandle, 'Synthesis (10^{-21} mol s^{-1})', 'FontSize', 7);
            set(axesHandle, 'FontSize', 5, 'Box', 'off', 'tickdir', 'out', ...
                'XTick', 0:4:8, 'yscale', 'log', 'ytick', 10.^(-6:2:2))
            
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', 1e-7, 'ytickoffset', 0.018);
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            for i = 1:numel(xTicks)
                xTickPos = get(xTicks(i), 'Position');
                set(xTicks(i), 'Position', [xTickPos(1) 3e-7 xTickPos(3:end)]);
            end
            set(xTicks, 'FontSize', 5);
            set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'right');
            set(yTicks(1), 'String', '10^{-6}')
            set(yTicks(2), 'String', '10^{-4}')
            set(yTicks(3), 'String', '10^{-2}')
            set(yTicks(4), 'String', '10^{0}')
            set(yTicks(5), 'String', '10^{2}')
            
            xLabelPos = get(get(axesHandle, 'xlabel'), 'position');
            yLabelPos = get(get(axesHandle, 'ylabel'), 'position');
            set(get(axesHandle, 'xlabel'), 'position', [xLabelPos(1) 1e-7 xLabelPos(end)]);
            set(get(axesHandle, 'ylabel'), 'position', [-1.0 yLabelPos(2:end)]);
            
            %% C energy dynamics of processes
            [~, pIdxs] = sort(processEnergyUsages);
            pIdxs = pIdxs(end:-1:end-14);
            
            for i = 1:size(processEnergyUsageDynamics, 1)
                for j = 1:size(processEnergyUsageDynamics, 2)
                    processEnergyUsageDynamics(i, j, :) = ...
                        conv(squeeze(processEnergyUsageDynamics(i, j, :)), ones(window, 1) / window, 'same');
                end
            end
            processEnergyUsageDynamics = processEnergyUsageDynamics / ConstantUtil.nAvogadro * 1e24;
            
            totEnergy = sum(max(0, processEnergyUsages)) + sim.process('Metabolism').unaccountedEnergyConsumption * nSims;
            
            nCols = 1;
            nRows = ceil(numel(pIdxs) / nCols);
            
            W = 0.27;
            H = 0.98;
            
            x = 0.65;
            y = 0.065;
            xMargin = 0.05;
            yMargin = 0.03;
            xAxisOffset = 0.01;
            yAxisOffset = xAxisOffset * figW / figH;
            w = (W - (nCols - 1) * xMargin) / nCols;
            h = (H - y - (nRows - 1) * yMargin) / nRows;
            
            PlotUtil.labelSubFigure('C', [x-0.01 0.993 -0.0286 -0.0208], figHandle, 8, 'normalized');
            
            for i = 1:numel(pIdxs)
                row = ceil(i / nCols);
                col = mod(i-1, nCols) + 1;
                axesHandle = subplot('Position', [x+(col-1)*(w+xMargin)  y+(nRows-row)*(h+yMargin) w  h], 'Parent', figHandle);
                handle = plot(axesHandle, (window:simEndTimes(selectedSim)-window+1) / 3600, ...
                    permute(processEnergyUsageDynamics(:, pIdxs(i), window:end-window+1), [1 3 2]));
                set(handle(1), 'Color', 'b');
                set(handle(2), 'Color', 'g');
                set(axesHandle, 'FontSize', 5, 'visible', 'off', 'yscale', 'linear');
                xlim(axesHandle, [0 simEndTimes(selectedSim) / 3600]);
                ylim(axesHandle, [0 max(max(processEnergyUsageDynamics(:, pIdxs(i), :)))]);
                annotation(figHandle, 'textbox', [x+(col-1)*(w+xMargin)  y+(nRows-row)*(h+yMargin)+h-0.02 w 1e-6], ...
                    'String', sim.processMetadata.names{pIdxs(i)}, ...
                    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                    'FontSize', 6);
                
                %y axis
                yAxesHandle = subplot('Position', [x-xAxisOffset y+(nRows-row)*(h+yMargin) 1e-6 h], 'Parent', figHandle);
                ylim(yAxesHandle, ylim(axesHandle));
                set(yAxesHandle, 'FontSize', 5, 'yscale', 'linear', 'box', 'off', 'tickdir', 'out');
                if row == ceil(nRows / 2) && col == 1
                    ylabel(yAxesHandle, {'NTP Use (10^{-24} mol s^{-1})'}, 'FontSize', 7);
                    yLabelPos = get(get(yAxesHandle, 'ylabel'), 'position');
                    set(get(yAxesHandle, 'ylabel'), 'position', [-140 yLabelPos(2:end)]);
                end
                
                tick2text(yAxesHandle, 'axis', 'y', 'ytickoffset', 20);
                yTicks = getappdata(yAxesHandle, 'YTickText');
                set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'right');
                
                delete(yTicks(2:end-1));
                
                %pie
                val = processEnergyUsages(pIdxs(i)) / totEnergy;
                
                scale = 1.4;
                axesHandle = subplot('Position', [1-scale*h*figH/figW-0.005  y+(nRows-row)*(h+yMargin)+h/2-scale*h/2  scale*h*figH/figW  scale*h], 'Parent', figHandle);                                
                pieHandle = pie(axesHandle, 100 * [val 1-val]);
                set(pieHandle(1), 'EdgeColor', 'k', 'FaceColor', 'r');
                set(pieHandle(3), 'EdgeColor', 'k', 'FaceColor', 0.75 * [1 1 1]);
                set(pieHandle(2:2:end), 'visible', 'off');
                set(axesHandle, 'visible', 'off', 'FontSize', 5);
                axis(axesHandle, 'equal');
                xlim(axesHandle, [-1.05 1.05])
                ylim(axesHandle, [-1.05 1.05])
            end
            
            %x axis
            xAxesHandle = subplot('Position', [x y-yAxisOffset w 1e-6], 'Parent', figHandle);
            xlim(xAxesHandle, [0 simEndTimes(selectedSim) / 3600]);
            xlabel(xAxesHandle, 'Time (h)', 'FontSize', 7);
            set(xAxesHandle, 'FontSize', 5, 'XTick', 0:4:8, 'box', 'off', 'tickdir', 'out');
            
            tick2text(xAxesHandle, 'axis', 'x', 'xtickoffset', 50);
            xTicks = getappdata(xAxesHandle, 'XTickText');
            set(xTicks, 'FontSize', 5);
            
            xLabelPos = get(get(xAxesHandle, 'xlabel'), 'position');
            set(get(xAxesHandle, 'xlabel'), 'position', [xLabelPos(1) -100 xLabelPos(3:end)]);
            
            %legend
            annotation(figHandle, 'ellipse', [1-0.06  0.015-0.01*figW/figH+0.02  0.01  0.01*figW/figH], ...
                'EdgeColor', 'none', 'FaceColor', 'b');
            annotation(figHandle, 'ellipse', [1-0.06  0.015-0.01*figW/figH  0.01  0.01*figW/figH], ...
                'EdgeColor', 'none', 'FaceColor', 'g');
            annotation(figHandle, 'textbox', [1 - 0.065, 0.005+0.02 1e-6 1e-6], 'String', 'ATP', ...
                'FontSize', 5, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'EdgeColor', 'none');
            annotation(figHandle, 'textbox', [1 - 0.065, 0.005 1e-6 1e-6], 'String', 'GTP', ...
                'FontSize', 5, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'EdgeColor', 'none');
            
            %% B: total energy usage vs cell cycle length
            PlotUtil.labelSubFigure('B', [leftColLabelX 0.673 -0.0286 -0.0208], figHandle, 8, 'normalized');
            axesHandle = subplot('Position', [leftColX 0.42 leftColW 0.23], 'Parent', figHandle);
            
            hold(axesHandle, 'on');
            if plotUsageRate
                ydata = totalEnergyUsagePerCell * 1e15 ./ simEndTimes(:, [1 1])' * 3600;
                ylim(axesHandle, [4.8 10.8]);
                ylabel(axesHandle, 'NTP Use (10^{-18} mol h^{-1})', 'FontSize', 7);
                set(axesHandle, 'YTick', 6:2:10);
                
                plot(axesHandle, 8, 10.5, 'b.', 'MarkerSize', 10);
                plot(axesHandle, 8, 9.8, 'g.', 'MarkerSize', 10);
                text(8.2, 10.5, 'ATP', 'FontSize', 6, 'VerticalAlignment', 'middle');
                text(8.2, 9.8, 'GTP', 'FontSize', 6, 'VerticalAlignment', 'middle');
            else
                ydata = totalEnergyUsagePerCell * 1e15;
                ylim(axesHandle, [50 max(totalEnergyUsagePerCell(:)) * 1e15]);
                ylabel(axesHandle, 'Tot NTP Use (10^{-18} mol)', 'FontSize', 7);
                set(axesHandle, 'YTick', 50:25:125);
                
                plot(axesHandle, 8, 125, 'b.', 'MarkerSize', 10);
                plot(axesHandle, 8, 115, 'g.', 'MarkerSize', 10);
                text(8.2, 125, 'ATP', 'FontSize', 6, 'VerticalAlignment', 'middle');
                text(8.2, 115, 'GTP', 'FontSize', 6, 'VerticalAlignment', 'middle');
            end
            h = plot(axesHandle, simEndTimes / 3600, ydata, '.');
            set(h(1), 'Color', 'b');
            set(h(2), 'Color', 'g');
            xlim(axesHandle, [7.7 14]);
            xlabel(axesHandle, 'Cell Cycle Length (h)', 'FontSize', 7);
            set(axesHandle, 'FontSize', 5, 'Box', 'off', 'TickDir', 'out', ...
                'XTick', 8:2:14)
            
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', 0.06, 'ytickoffset', 0.015);
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            set(xTicks, 'FontSize', 5);
            set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'right');
            
            xLabelPos = get(get(axesHandle, 'xlabel'), 'position');
            yLabelPos = get(get(axesHandle, 'ylabel'), 'position');
            if plotUsageRate
                set(get(axesHandle, 'xlabel'), 'position', [xLabelPos(1) 4 xLabelPos(end)]);
                set(get(axesHandle, 'ylabel'), 'position', [7.1 yLabelPos(2:end)]);
            else
                set(get(axesHandle, 'xlabel'), 'position', [xLabelPos(1) 41 xLabelPos(end)]);
                set(get(axesHandle, 'ylabel'), 'position', [7.02 yLabelPos(2:end)]);
            end
            
            %% D: energy allocation among processes
            processNames = sim.processMetadata.names;
            processNames(sim.processIndex('DNASupercoiling')) = {'Supercoil'};
            processNames(sim.processIndex('ProteinDecay')) = {'Prot Dcy'};
            processNames(sim.processIndex('tRNAAminoacylation')) = {'tRNA Acyl'};
            processNames(sim.processIndex('Transcription')) = {'Transcription'};
            processNames(sim.processIndex('Translation')) = {'Translation'};
            processNames(sim.processIndex('Replication')) = {'DNA Rep'};
            
            PlotUtil.labelSubFigure('D', [leftColLabelX 0.32 -0.0286 -0.0208], figHandle, 8, 'normalized');
            axesHandle = subplot('Position', [leftColLabelX-0.02 -0.085 leftColW 0.50], 'Parent', figHandle);
            [~, pIdxs] = sort(processEnergyUsages);
            pIdxs = pIdxs(end:-1:end-2);
            vals = [
                sum(max(0, processEnergyUsages(setdiff(1:end, pIdxs))))
                sim.process('Metabolism').unaccountedEnergyConsumption * nSims
                processEnergyUsages(pIdxs)                
                ];
            totEnergy = sum(vals);
            vals = vals / totEnergy * 100;
            labels = ['Other'; 'Unaccounted'; processNames(pIdxs)];
            
            colors = [                
                0 0 1
                1 0 0
                0 1 0
                0 1 1
                1 0.5 0
                ];
            explode = [1 0.5 0.5 0.7 0.55];
            shiftX = [0.25 0 -0.13 0.05 0];
            shiftY = [0 0 0 0 0];
            for i = 1:numel(vals)
                ptch = patch(...
                    [0 cos((linspace(sum(vals(1:i-1)), sum(vals(1:i)), 50)-vals(1)/2)/100 * 2 * pi)], ...
                    [0 sin((linspace(sum(vals(1:i-1)), sum(vals(1:i)), 50)-vals(1)/2)/100 * 2 * pi)], ...
                    1, ...
                    'Parent', axesHandle, ...
                    'EdgeColor', 'k', ...
                    'FaceColor', colors(i, :), ...
                    'FaceAlpha', 1, 'EdgeAlpha', 1);
            end
            for i = 1:numel(vals)
                text(...
                    explode(i) * cos((sum(vals(1:i-1))+vals(i)/2-vals(1)/2)/100 * 2 * pi) + shiftX(i), ...
                    explode(i) * sin((sum(vals(1:i-1))+vals(i)/2-vals(1)/2)/100 * 2 * pi) + shiftY(i) + 0.07, ...
                    labels{i}, 'FontSize', 6, ...
                    'HorizontalAlign', 'center', 'VerticalAlign', 'middle');
                text(...
                    explode(i) * cos((sum(vals(1:i-1))+vals(i)/2-vals(1)/2)/100 * 2 * pi) + shiftX(i), ...
                    explode(i) * sin((sum(vals(1:i-1))+vals(i)/2-vals(1)/2)/100 * 2 * pi) + shiftY(i) - 0.07, ...
                    sprintf('(%0.1f%%)', vals(i)), 'FontSize', 6, ...
                    'HorizontalAlign', 'center', 'VerticalAlign', 'middle');
            end
            
            xlim(axesHandle, [-1.05 1.05]);
            ylim(axesHandle, [-1.05 1.05]);
            axis(axesHandle, 'equal');
            set(axesHandle, 'FontSize', 5, 'visible', 'off');
            
            %% standarize tick lengths
            axesHandles = findobj(figHandle, 'type', 'axes', '-and', 'tag', '');
            axesPos = cell2mat(get(axesHandles, 'position'));
            tickLen = 0.01 * max(max(axesPos(:, 3:4), [], 1) .* [figW figH]);
            for i = 1:numel(axesHandles)
                maxVal = max(axesPos(i, 3:4) .* [figW figH]);
                set(axesHandles(i), 'TickLen', [tickLen / maxVal  0.0250]);
            end
            
            %% save
            print(figHandle, [outDirectory filesep 'figure5.pdf'], '-dpdf', '-rgb');
            print(figHandle, [outDirectory filesep 'figure5.eps'], '-deps', '-rgb');
            close(figHandle);
            
            %% supplemental figures
            
            %energy carrier production
            for i = 1:size(energyCarrierSynthesis, 1)
                energyCarrierSynthesis(i, :) = conv(energyCarrierSynthesis(i, :), ones(window, 1) / window, 'same');
            end
            [axesHandle, figHandle] = PlotUtil.newAxesHandle([8.5 8.5]);            
            h = plot(axesHandle, (window:simEndTimes(selectedSim)-window+1) / 3600, -energyCarrierSynthesis(:, window:end-window+1));
            xlim(axesHandle, [0 simEndTimes(selectedSim) / 3600]);
            xlabel(axesHandle, 'Time (h)', 'FontSize', 7);
            ylabel(axesHandle, 'Synthesis (molecule s^{-1})', 'FontSize', 7);
            set(axesHandle, 'box', 'off', 'tickdir', 'out', 'FontSize', 5, 'yscale', 'log', 'XTick', 0:2:8, 'YTick', 10.^(-2:2:4));
            legend(h, {'ATP'; 'GTP'; 'FAD(H_2)'; 'NAD(H)'; 'NADP(H)'}, 'Location', 'NorthWest')
            saveas(figHandle, [outDir 'energyCarrierSynthesis.pdf']);
            close(figHandle);            
            
            %RNA decay vs transcription dynamics
            [axesHandle, figHandle] = PlotUtil.newAxesHandle([8.5 4]);
            clf(figHandle);
            
            axesHandle = subplot('Position', [0.11 0.14 0.36 0.85], 'Parent', figHandle);
            x = -15:15;
            y = xcorr(transcriptionDynamics, ntpDynamics, x(end), 'coeff');
            plot(axesHandle, x, y);
            [~, idx] = max(y);
            xlim(axesHandle, [-15 15]);
            ylim(axesHandle, [min(y) max(y)])
            line(x(idx) * [1 1], ylim(axesHandle), 'Color', 'r', 'LineStyle', ':', 'Parent', axesHandle)
            text(2.75, min(y) + 0.95 * range(y), sprintf('%d s', x(idx)), 'Color', 'r', 'FontSize', 6, 'Parent', axesHandle);
            xlabel(axesHandle, 'Lag (s)', 'FontSize', 7);
            ylabel(axesHandle, {'NTP-Transcription X-Corr'}, 'FontSize', 7);
            set(axesHandle, 'box', 'off', 'tickdir', 'out', 'FontSize', 5, 'XTick', -10:10:10, 'YTick', 0.600:0.002:0.606);
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', 0.04, 'ytickoffset', 0.03);
            set(getappdata(axesHandle, 'XTickText'), 'FontSize', 5);
            set(getappdata(axesHandle, 'YTickText'), 'FontSize', 5, 'HorizontalAlign', 'right');
            set(get(axesHandle, 'ylabel'), 'Position', [-21 mean(ylim(axesHandle)) 0]);
            
            axesHandle = subplot('Position', [0.64 0.14 0.36 0.85], 'Parent', figHandle);
            x = -15:15;
            y = xcorr(transcriptionDynamics, rnaDecayDynamics, x(end), 'coeff');
            plot(axesHandle, x, y);
            [~, idx] = max(y);
            xlim(axesHandle, [-15 15]);
            ylim(axesHandle, [min(y) max(y)])
            line(x(idx) * [1 1], ylim(axesHandle), 'Color', 'r', 'LineStyle', ':', 'Parent', axesHandle)
            text(3.75, min(y) + 0.95 * range(y), sprintf('%d s', x(idx)), 'Color', 'r', 'FontSize', 6, 'Parent', axesHandle);
            xlabel(axesHandle, 'Lag (s)', 'FontSize', 7);
            ylabel(axesHandle, 'Transcription-Decay X-Corr', 'FontSize', 7);
            set(axesHandle, 'box', 'off', 'tickdir', 'out', 'FontSize', 5, 'XTick', -10:10:10, 'YTick', 0.4:0.2:0.8);
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', 0.04, 'ytickoffset', 0.03);
            set(getappdata(axesHandle, 'XTickText'), 'FontSize', 5);
            set(getappdata(axesHandle, 'YTickText'), 'FontSize', 5, 'HorizontalAlign', 'right')
            set(get(axesHandle, 'ylabel'), 'Position', [-19 mean(ylim(axesHandle)) 0]);
            
            saveas(figHandle, [outDir 'transcriptionVsRNADecayDynamics.pdf']);
            close(figHandle);
            
            %correlations with NTP concentrations
            [axesHandle, figHandle] = PlotUtil.newAxesHandle([17.4 17.4]);
            clf(figHandle);
            
            axesHandle = subplot(3, 3, 1, 'Parent', figHandle);
            plot(axesHandle, ntpCounts * 1e15, simEndTimes / 3600, '.');
            xlabel(axesHandle, 'Mean NTP Count (amol)', 'FontSize', 7);
            ylabel(axesHandle, 'Cell Cycle Length (h)', 'FontSize', 7);
            set(axesHandle, 'FontSize', 5, 'box', 'off', 'tickDir', 'out')
            
            axesHandle = subplot(3, 3, 2, 'Parent', figHandle);
            plot(axesHandle, ntpConcs, simEndTimes / 3600, '.');
            xlabel(axesHandle, 'Mean NTP Conc (mM)', 'FontSize', 7);
            ylabel(axesHandle, 'Cell Cycle Length (h)', 'FontSize', 7);
            set(axesHandle, 'FontSize', 5, 'box', 'off', 'tickDir', 'out')
            
            axesHandle = subplot(3, 3, 3, 'Parent', figHandle);
            plot(axesHandle, ntpConcs, repInitDuration / 3600, '.');
            xlabel(axesHandle, 'Mean NTP Conc (mM)', 'FontSize', 7);
            ylabel(axesHandle, 'Rep Init Duration (h)', 'FontSize', 7);
            set(axesHandle, 'FontSize', 5, 'box', 'off', 'tickDir', 'out')
            
            axesHandle = subplot(3, 3, 4, 'Parent', figHandle);
            plot(axesHandle, ntpConcs, repDuration / 3600, '.');
            xlabel(axesHandle, 'Mean NTP Conc (mM)', 'FontSize', 7);
            ylabel(axesHandle, 'Replication Duration (h)', 'FontSize', 7);
            set(axesHandle, 'FontSize', 5, 'box', 'off', 'tickDir', 'out')
            
            axesHandle = subplot(3, 3, 5, 'Parent', figHandle);
            plot(axesHandle, ntpConcs, rnaWt * 1e15, '.');
            xlabel(axesHandle, 'Mean NTP Conc (mM)', 'FontSize', 7);
            ylabel(axesHandle, 'RNA (fg)', 'FontSize', 7);
            set(axesHandle, 'FontSize', 5, 'box', 'off', 'tickDir', 'out')
            
            axesHandle = subplot(3, 3, 6, 'Parent', figHandle);
            plot(axesHandle, ntpConcs, proteinWt * 1e15, '.');
            xlabel(axesHandle, 'Mean NTP Conc (mM)', 'FontSize', 7);
            ylabel(axesHandle, 'Protein (fg)', 'FontSize', 7);
            set(axesHandle, 'FontSize', 5, 'box', 'off', 'tickDir', 'out')
            
            axesHandle = subplot(3, 3, 7, 'Parent', figHandle);
            plot(axesHandle, ntpConcs, superhelicity, '.');
            xlabel(axesHandle, 'Mean NTP Conc (mM)', 'FontSize', 7);
            ylabel(axesHandle, 'Superhelicity', 'FontSize', 7);
            set(axesHandle, 'FontSize', 5, 'box', 'off', 'tickDir', 'out')
            
            axesHandle = subplot(3, 3, 8, 'Parent', figHandle);
            plot(axesHandle, ntpConcs, activeRNAPols, '.');
            xlabel(axesHandle, 'Mean NTP Conc (mM)', 'FontSize', 7);
            ylabel(axesHandle, 'Active RNA Pols', 'FontSize', 7);
            set(axesHandle, 'FontSize', 5, 'box', 'off', 'tickDir', 'out')
            
            axesHandle = subplot(3, 3, 9, 'Parent', figHandle);
            plot(axesHandle, ntpConcs, stalledRibs, '.');
            xlabel(axesHandle, 'Mean NTP Conc (mM)', 'FontSize', 7);
            ylabel(axesHandle, 'Stalled Ribosomes', 'FontSize', 7);
            set(axesHandle, 'FontSize', 5, 'box', 'off', 'tickDir', 'out')
            
            saveas(figHandle, [outDir 'ntpCorrelations.pdf']);
            close(figHandle);
        end
        
        function growthSensitivity(simBatchDir, selectedSim, outDirectory)
            import edu.stanford.covert.cell.sim.analysis.Figures34;
            import edu.stanford.covert.cell.sim.util.FitConstants;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;            
            import edu.stanford.covert.util.ComputationUtil;
            
            if ~exist([outDirectory filesep 'growthSensitivity'], 'dir')
                mkdir([outDirectory filesep 'growthSensitivity']);
            end            
            outDir = [outDirectory filesep 'growthSensitivity' filesep];
            
            %% get data
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            
%             if ~exist([outDirectory filesep 'growthSensitivity.mat'], 'file')
%                 save([outDirectory filesep 'growthSensitivity.mat']);
%             else
%                 load([outDirectory filesep 'growthSensitivity.mat']);
%             end

            g = sim.gene;
            r = sim.state('Rna');
            mr = sim.state('MetabolicReaction');
            met = sim.process('Metabolism');
            
            %% path to optimal growth
            %get metabolic network
            fitter = FitConstants(sim);
            [~, freeMons, freeCpxs] = fitter.calcMacromolecularCounts(fitter.constructParameterVectorFromSimulation());
            
            growth0 = mr.growth0 * met.macromoleculeStateInitializationGrowthFactor;
            
            substrates = zeros(size(met.substrates));
            substrates(met.substrateMetaboliteLocalIndexs, :) = met.substrates(met.substrateMetaboliteLocalIndexs, :);
            substrates(met.substrateIndexs_internalExchangedLimitedMetabolites) = ...
                -growth0 *  met.metabolismRecyclingProduction(met.substrateIndexs_internalExchangedLimitedMetabolites);
            substrates(met.substrateMonomerLocalIndexs, :) = repmat(freeMons(met.substrateMonomerGlobalIndexs(:, 1)), 1, size(substrates, 2));
            substrates(met.substrateComplexLocalIndexs, :) = repmat(freeCpxs(met.substrateComplexGlobalIndexs(:, 1)), 1, size(substrates, 2));
            
            enzymes = zeros(size(met.enzymeWholeCellModelIDs));
            enzymes(met.enzymeMonomerLocalIndexs) = freeMons(met.enzymeMonomerGlobalIndexs);
            enzymes(met.enzymeComplexLocalIndexs) = freeCpxs(met.enzymeComplexGlobalIndexs);
            
            [fbaObj, fbaSMat, fbaRxnBounds] = met.calcEffectiveFBANetwork();
            fluxBounds = met.calcFluxBounds(substrates, enzymes, fbaRxnBounds, met.fbaEnzymeBounds);
            [initFbaReactionFluxs, ~, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                'maximize', fbaObj, fbaSMat, ...
                met.fbaRightHandSide, fluxBounds(:, 1), fluxBounds(:, 2), ...
                'S', 'C', met.linearProgrammingOptions);
            initGrowth = initFbaReactionFluxs(met.fbaReactionIndexs_biomassProduction);
            
            fbaEnzBounds = met.fbaEnzymeBounds;
            fbaRHS = met.fbaRightHandSide;
            fbaRxnCatMat = met.fbaReactionCatalysisMatrix;
            bmRxnTfs = false(size(fbaSMat, 2), 1);
            bmRxnTfs(met.fbaReactionIndexs_biomassProduction) = true;
            nonKineticFluxBounds = met.calcFluxBounds(...
                substrates, enzymes, fbaRxnBounds, met.fbaEnzymeBounds, ...
                false, true, true, false, true, true);
            
            kineticLimitedEnzymes = any(met.fbaReactionCatalysisMatrix(any(isfinite(met.fbaEnzymeBounds), 2), :), 1)';
            kineticUnlimitedEnzymes = any(met.fbaReactionCatalysisMatrix(any(isinf(met.fbaEnzymeBounds), 2), :), 1)';
            
            %split reversible reactions into positive and negative
            %components
            idxs = find(fluxBounds(:, 1) < 0);
            fluxBounds = [max(fluxBounds, 0); max(-fliplr(fluxBounds(idxs, :)), 0)];
            nonKineticFluxBounds = [max(nonKineticFluxBounds, 0); max(-fliplr(nonKineticFluxBounds(idxs, :)), 0)];
            fbaEnzBounds = [max(fbaEnzBounds, 0); max(-fliplr(fbaEnzBounds(idxs, :)), 0)];
            fbaSMat = [fbaSMat -fbaSMat(:, idxs)];
            fbaObj = [fbaObj; -fbaObj(idxs)];
            fbaRxnCatMat = [fbaRxnCatMat; fbaRxnCatMat(idxs, :)];
            bmRxnTfs = [bmRxnTfs; bmRxnTfs(idxs)];
            
            idxs = find(fluxBounds(:, 2) ~= 0);
            fluxBounds = fluxBounds(idxs, :);
            nonKineticFluxBounds = nonKineticFluxBounds(idxs, :);
            fbaEnzBounds = fbaEnzBounds(idxs, :);
            fbaSMat = fbaSMat(:, idxs);
            fbaObj = fbaObj(idxs);
            fbaRxnCatMat = fbaRxnCatMat(idxs, :);
            bmRxnTfs = bmRxnTfs(idxs);            
            
            nMet = size(fbaSMat, 1);
            nRxn = size(fbaSMat, 2);
            nEnz = size(fbaRxnCatMat, 2);
            nKinEnz = sum(kineticLimitedEnzymes);
            bmRxnIdx = find(bmRxnTfs);
            uFbaEnzBounds = unique(fbaEnzBounds(:, 2));
            kineticConstraints = (1 ./ min(uFbaEnzBounds(end-1), fbaEnzBounds(:, 2 * ones(nKinEnz, 1))))';
            kineticConstraints(fbaRxnCatMat(:, kineticLimitedEnzymes)' == 0) = 0;
            enzymeConstraints = zeros(nKinEnz, nEnz);
            enzymeConstraints(sub2ind(size(enzymeConstraints), 1:nKinEnz, find(kineticLimitedEnzymes)')) = 1;
            
            enzGeneComp = met.enzymeGeneComposition;
            
            %check that growth not affected by splitting
            [initFbaReactionFluxs_flux, ~, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                'maximize', fbaObj, fbaSMat, ...
                fbaRHS, ...
                fluxBounds(:, 1), fluxBounds(:, 2), ...
                'S', 'C', met.linearProgrammingOptions);
            if errFlag
                throw(MException('Metabolism:error', 'Unable to optimize growth: %s', errMsg));
            end
            assertElementsAlmostEqual(initGrowth, initFbaReactionFluxs_flux(bmRxnIdx));
            
            %check that growth not affected additional constraints on
            %enzyme mass
            [initFbaReactionFluxs_enzyme, ~, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                'maximize', [fbaObj; zeros(nEnz, 1)], [fbaSMat zeros(nMet, nEnz); zeros(1, nRxn) met.enzymeMolecularWeights'; kineticConstraints -enzymeConstraints], ...
                [fbaRHS; enzymes' * met.enzymeMolecularWeights; zeros(nKinEnz, 1)], ...
                [nonKineticFluxBounds(:, 1); zeros(nEnz, 1) + enzymes], [nonKineticFluxBounds(:, 2); inf(nEnz, 1)], ...
                [repmat('S', size(fbaRHS)); 'S'; repmat('U', nKinEnz, 1)], 'C', met.linearProgrammingOptions);
            if errFlag
                throw(MException('Metabolism:error', 'Unable to optimize growth: %s', errMsg));
            end
            assertElementsAlmostEqual(initGrowth, initFbaReactionFluxs_enzyme(bmRxnIdx));
            
            %find maximal growth: enzyme method
            [result, ~, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                'maximize', [fbaObj; zeros(nEnz, 1)], [fbaSMat zeros(nMet, nEnz); zeros(1, nRxn) met.enzymeMolecularWeights'; kineticConstraints -enzymeConstraints], ...
                [fbaRHS; enzymes' * met.enzymeMolecularWeights; zeros(nKinEnz, 1)], ...
                [nonKineticFluxBounds(:, 1); (kineticLimitedEnzymes & kineticUnlimitedEnzymes) + (~kineticLimitedEnzymes) .* enzymes], [nonKineticFluxBounds(:, 2); inf(nEnz, 1)], ...
                [repmat('S', size(fbaRHS)); 'S'; repmat('U', nKinEnz, 1)], 'C', met.linearProgrammingOptions);
            if errFlag
                throw(MException('Metabolism:error', 'Unable to optimize growth: %s', errMsg));
            end
            maxGrowth_enzyme = result(bmRxnIdx);
            assertElementsAlmostEqual(maxGrowth_enzyme, ...
                met.calcGrowthRate(met.calcFluxBounds(substrates, result(end-nEnz+1:end), met.fbaReactionBounds, met.fbaEnzymeBounds)));
            
            %find maximal growth: flux method
            totalEnzymeConstriant = 1 ./ fbaEnzBounds(:, 2) .* (fbaRxnCatMat * met.enzymeMolecularWeights) * 1e-7;
            totalEnzymeConstriant(isinf(fbaEnzBounds(:, 2)) | fbaEnzBounds(:, 2) == 0) = 0;
            [result, ~, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                'maximize', fbaObj, [fbaSMat; totalEnzymeConstriant'], ...
                [fbaRHS; enzymes(kineticLimitedEnzymes)' * met.enzymeMolecularWeights(kineticLimitedEnzymes) * 1e-7], ...
                nonKineticFluxBounds(:, 1), nonKineticFluxBounds(:, 2), ...
                [repmat('S', size(fbaRHS)); 'U'], 'C', met.linearProgrammingOptions);
            if errFlag
                throw(MException('Metabolism:error', 'Unable to optimize growth: %s', errMsg));
            end
            maxGrowth_flux = result(bmRxnIdx); 
            
            %find intermediate growths
            [optGrowths, optEnzymes, optEnzymeOrder] = Figures34.calcOptimalEnzymeExpressionChanges(fbaObj, fbaSMat, fbaRHS, ...
                kineticConstraints, enzymeConstraints, nonKineticFluxBounds, ...
                enzymes, kineticLimitedEnzymes, kineticUnlimitedEnzymes, bmRxnIdx, sim);
            
            %find intermediate growths
            nSteps = 40;
            growths_enzyme = [
                linspace(0, initGrowth, 20) ...
                linspace(initGrowth, maxGrowth_enzyme, 20)
                ];
            growths_flux = [
                linspace(0, initGrowth, 20) ...
                linspace(initGrowth, maxGrowth_flux, 20)
                ];
            distsL1 = zeros(nSteps, 2);
            distsL2 = zeros(nSteps, 2);
            distsInf = zeros(nSteps, 2);
            distsL1Reg = zeros(nSteps, 2);
            cplexopts = cplexoptimset('Display', 'off');
            for i = 1:nSteps
                %enzyme methods
                tmpNonKineticFluxBounds = nonKineticFluxBounds;
                tmpNonKineticFluxBounds(bmRxnIdx, :) = growths_enzyme(i);
                [~, ~, distsL1(i, 1), errFlag, errMsg] = ComputationUtil.linearProgramming(...
                    'minimize', [
                    zeros(nRxn, 1)
                    zeros(nEnz, 1)
                    ones(nEnz, 1)
                    ], [...
                    fbaSMat zeros(nMet, nEnz) zeros(nMet, nEnz)
                    zeros(1, nRxn) met.enzymeMolecularWeights' zeros(1, nEnz)
                    kineticConstraints -enzymeConstraints zeros(nKinEnz, nEnz)
                    zeros(nEnz, nRxn)  eye(nEnz) -eye(nEnz)
                    zeros(nEnz, nRxn) -eye(nEnz) -eye(nEnz)
                    ], [
                    fbaRHS
                    enzymes' * met.enzymeMolecularWeights
                    zeros(nKinEnz, 1)
                    initFbaReactionFluxs_enzyme(end-nEnz+1:end)
                    -initFbaReactionFluxs_enzyme(end-nEnz+1:end)
                    ], [
                    tmpNonKineticFluxBounds(:, 1)
                    (kineticLimitedEnzymes & kineticUnlimitedEnzymes) + (~kineticLimitedEnzymes) .* enzymes
                    -inf(nEnz, 1)
                    ], [
                    tmpNonKineticFluxBounds(:, 2)
                    inf(nEnz, 1)
                    inf(nEnz, 1)
                    ], ...
                    [
                    repmat('S', size(fbaRHS))
                    'S'
                    repmat('U', nKinEnz, 1)
                    repmat('U', 2 * nEnz, 1)
                    ], 'C', met.linearProgrammingOptions);
                if errFlag
                    throw(MException('Metabolism:error', 'Unable to optimize growth: %s', errMsg));
                end
                
                [~, distsL2(i, 1), errFlag, output] = cplexqp(...
                    2 * diag([zeros(nRxn, 1); ones(nEnz, 1)]), -2 * ([zeros(nRxn, 1); ones(nEnz, 1)]) .* initFbaReactionFluxs_enzyme, ...
                    [kineticConstraints -enzymeConstraints], zeros(nKinEnz, 1), ...
                    [fbaSMat zeros(nMet, nEnz); zeros(1, nRxn) met.enzymeMolecularWeights' * 1e-7], [fbaRHS; enzymes' * met.enzymeMolecularWeights * 1e-7], ...
                    [tmpNonKineticFluxBounds(:, 1); (kineticLimitedEnzymes & kineticUnlimitedEnzymes) + (~kineticLimitedEnzymes) .* enzymes], [tmpNonKineticFluxBounds(:, 2); inf(nEnz, 1)], ...
                    [], cplexopts);
                if errFlag ~= 1
                    throw(MException('Metabolism:error', 'Unable to optimize growth: %s', output.message));
                end
                
                [~, ~,  distsInf(i, 1), errFlag, errMsg] = ComputationUtil.linearProgramming(...
                    'minimize', [zeros(nRxn, 1); zeros(nEnz, 1); 1], [
                    fbaSMat zeros(nMet, nEnz) zeros(nMet, 1)
                    zeros(1, nRxn) met.enzymeMolecularWeights' 0
                    kineticConstraints -enzymeConstraints zeros(nKinEnz, 1)
                    zeros(nEnz, nRxn)  eye(nEnz) -ones(nEnz, 1)
                    zeros(nEnz, nRxn) -eye(nEnz) -ones(nEnz, 1)
                    ], [
                    fbaRHS
                    enzymes' * met.enzymeMolecularWeights
                    zeros(nKinEnz, 1)
                    initFbaReactionFluxs_enzyme(end-nEnz+1:end)
                    -initFbaReactionFluxs_enzyme(end-nEnz+1:end)
                    ], ...
                    [tmpNonKineticFluxBounds(:, 1); (kineticLimitedEnzymes & kineticUnlimitedEnzymes) + (~kineticLimitedEnzymes) .* enzymes; -inf], ...
                    [tmpNonKineticFluxBounds(:, 2); inf(nEnz, 1); inf], ...
                    [repmat('S', size(fbaRHS)); 'S'; repmat('U', nKinEnz, 1); repmat('U', 2 * nEnz, 1)], 'C', met.linearProgrammingOptions);
                if errFlag
                    throw(MException('Metabolism:error', 'Unable to optimize growth: %s', errMsg));
                end
                
                %flux methods
                tmpNonKineticFluxBounds = nonKineticFluxBounds;
                tmpNonKineticFluxBounds(bmRxnIdx, :) = growths_flux(i);
                
                [~, ~, distsL1(i, 2), errFlag, errMsg] = ComputationUtil.linearProgramming(...
                    'minimize', [zeros(nRxn, 1); ones(nRxn, 1)], ...
                    [[fbaSMat; totalEnzymeConstriant'] zeros(nMet+1, nRxn); eye(nRxn) -eye(nRxn); -eye(nRxn) -eye(nRxn)], ...
                    [fbaRHS; enzymes(kineticLimitedEnzymes)' * met.enzymeMolecularWeights(kineticLimitedEnzymes) * 1e-7; initFbaReactionFluxs_flux; -initFbaReactionFluxs_flux], ...
                    [tmpNonKineticFluxBounds(:, 1); -inf(nRxn, 1)], [tmpNonKineticFluxBounds(:, 2); inf(nRxn, 1)], ...
                    [repmat('S', nMet, 1); 'U'; repmat('U', 2 * nRxn, 1)], 'C', met.linearProgrammingOptions);
                if errFlag
                    throw(MException('Metabolism:error', 'Unable to optimize growth: %s', errMsg));
                end
                
                [~, distsL2(i, 2), errFlag, output] = cplexqp(...
                    2 * eye(nRxn), -2 * initFbaReactionFluxs_flux, ...
                    totalEnzymeConstriant', enzymes(kineticLimitedEnzymes)' * met.enzymeMolecularWeights(kineticLimitedEnzymes) * 1e-7, ...
                    fbaSMat, fbaRHS, ...
                    tmpNonKineticFluxBounds(:, 1), tmpNonKineticFluxBounds(:, 2), ...
                    [], cplexopts);
                if errFlag ~= 1
                    throw(MException('Metabolism:error', 'Unable to optimize growth: %s', output.message));
                end
                
                [~, ~, distsInf(i, 2), errFlag, errMsg] = ComputationUtil.linearProgramming(...
                    'minimize', [zeros(nRxn, 1); 1], ...
                    [[fbaSMat; totalEnzymeConstriant'] zeros(nMet+1, 1); eye(nRxn) -ones(nRxn, 1); -eye(nRxn) -ones(nRxn, 1)], ...
                    [fbaRHS; enzymes(kineticLimitedEnzymes)' * met.enzymeMolecularWeights(kineticLimitedEnzymes) * 1e-7; initFbaReactionFluxs_flux; -initFbaReactionFluxs_flux], ...
                    [tmpNonKineticFluxBounds(:, 1); -inf], [tmpNonKineticFluxBounds(:, 2); inf], ...
                    [repmat('S', nMet, 1); 'U'; repmat('U', 2 * nRxn, 1)], 'C', met.linearProgrammingOptions);
                if errFlag
                    throw(MException('Metabolism:error', 'Unable to optimize growth: %s', errMsg));
                end
            end
            
            % plot
            figHandle = figure;
            ylims = [initGrowth max(maxGrowth_enzyme, maxGrowth_flux)];
            
            axesHandle = subplot(1, 1, 1);
            h = plot(0:size(optGrowths, 1) - 1, optGrowths);
            set(axesHandle, 'FontSize', 5, 'TickDir', 'out', 'box', 'off');
            xlabel('Regularized Enzyme distance', 'FontSize', 7);
            ylabel('Growth (cell s^{-1})', 'FontSize', 7);
            ylim(ylims);
            xlim([-0.5 numel(optGrowths)-0.5]);
            for i = 1:size(optGrowths, 1)
                text(i-1, optGrowths(i), sprintf('%d ', optEnzymeOrder{i}), ...
                    'FontSize', 5, 'HorizontalAlign', 'left', 'VerticalAlign', 'bottom', 'interpreter', 'none', ...
                    'Color', get(h(2), 'Color'));
            end
            
            saveas(figHandle, [outDir 'growthEvolution-1 gene at a time.pdf']);
            close(figHandle);
            
            % plot all metrics
            figHandle = figure;
            ylims = [0 max(maxGrowth_enzyme, maxGrowth_flux)];
            
            axesHandle = subplot(2, 4, 1);
            h = plot(0:size(optGrowths, 1) - 1, optGrowths);
            set(axesHandle, 'FontSize', 5, 'TickDir', 'out', 'box', 'off');
            xlabel('Regularized Enzyme distance', 'FontSize', 7);
            ylabel('Growth (cell s^{-1})', 'FontSize', 7);
            ylim(ylims);
            xlim([-0.5 numel(optGrowths)-0.5]);
            for i = 1:size(optGrowths, 1)
                text(i-1, optGrowths(i), sprintf('%d ', optEnzymeOrder{i}), ...
                    'FontSize', 5, 'HorizontalAlign', 'left', 'VerticalAlign', 'bottom', 'interpreter', 'none', ...
                    'Color', get(h(2), 'Color'));
            end
            
            axesHandle = subplot(2, 4, 2);
            plot(distsL1(:, 1), growths_enzyme)
            xlabel('L_1 Enzyme distance', 'FontSize', 7)
            ylim(ylims);
            set(axesHandle, 'FontSize', 5, 'TickDir', 'out', 'box', 'off');
            
            axesHandle = subplot(2, 4, 3);
            plot(distsL2(:, 1), growths_enzyme)
            xlabel('L_2 Enzyme distance', 'FontSize', 7)
            ylim(ylims);
            set(axesHandle, 'FontSize', 5, 'TickDir', 'out', 'box', 'off');
            
            axesHandle = subplot(2, 4, 4);
            plot(distsInf(:, 1), growths_enzyme)
            xlabel('L_{\infty} Enzyme distance', 'FontSize', 7)
            ylim(ylims);
            set(axesHandle, 'FontSize', 5, 'TickDir', 'out', 'box', 'off');
            
            axesHandle = subplot(2, 4, 6);
            plot(distsL1(:, 2), growths_flux)
            xlabel('L_1 Flux distance', 'FontSize', 7)
            ylim(ylims);
            set(axesHandle, 'FontSize', 5, 'TickDir', 'out', 'box', 'off');
            
            axesHandle = subplot(2, 4, 7);
            plot(distsL2(:, 2), growths_flux)
            xlabel('L_2 Flux distance', 'FontSize', 7)
            ylim(ylims);
            set(axesHandle, 'FontSize', 5, 'TickDir', 'out', 'box', 'off');
            
            axesHandle = subplot(2, 4, 8);
            plot(distsInf(:, 2), growths_flux)
            xlabel('L_{\infty} Flux distance', 'FontSize', 7)
            ylim(ylims);
            set(axesHandle, 'FontSize', 5, 'TickDir', 'out', 'box', 'off');
            
            saveas(figHandle, [outDir 'growthEvolution-several metrics.pdf']);
            close(figHandle);
            
            %% compare optimal to suboptimal growth            
            paramVec0 = fitter.constructParameterVectorFromSimulation();
            [rnaExp, nmpComp, aaComp, rnaWtFracs, rnaDecayRates, ...
                biomassComposition, biomassProduction, byproducts, ...
                unaccountedEnergyConsumption] = fitter.extractParameterVector(paramVec0);
            
            import edu.stanford.covert.util.ConstantUtil;
            pm = sim.state('ProteinMonomer');
            t = sim.state('Time');
            mass = sim.state('Mass');
            rnaMWs = r.molecularWeights(r.matureIndexs) / ConstantUtil.nAvogadro;
            monMWs = pm.molecularWeights(pm.matureIndexs) / ConstantUtil.nAvogadro;
            
            rnaCnt = zeros(size(r.matureIndexs));
            rnaCnt(r.matureMRNAIndexs) = ...
                rnaExp(r.matureMRNAIndexs) * ...
                mass.cellInitialDryWeight * mass.initialFractionNTPsInRNAs * mass.dryWeightFractionRNA * ...
                rnaWtFracs(r.mRNAWeightFractionIndexs) / (rnaMWs(r.matureMRNAIndexs)' * rnaExp(r.matureMRNAIndexs));
            rnaCnt(r.matureRRNAIndexs) = ...
                rnaExp(r.matureRRNAIndexs) * ...
                mass.cellInitialDryWeight * mass.initialFractionNTPsInRNAs * mass.dryWeightFractionRNA * ...
                sum(rnaWtFracs(r.rRNAWeightFractionIndexs)) / (rnaMWs(r.matureRRNAIndexs)' * rnaExp(r.matureRRNAIndexs));
            rnaCnt(r.matureSRNAIndexs) = ...
                rnaExp(r.matureSRNAIndexs) * ...
                mass.cellInitialDryWeight * mass.initialFractionNTPsInRNAs * mass.dryWeightFractionRNA * ...
                rnaWtFracs(r.sRNAWeightFractionIndexs) / (rnaMWs(r.matureSRNAIndexs)' * rnaExp(r.matureSRNAIndexs));
            rnaCnt(r.matureTRNAIndexs) = ...
                rnaExp(r.matureTRNAIndexs) * ...
                mass.cellInitialDryWeight * mass.initialFractionNTPsInRNAs * mass.dryWeightFractionRNA * ...
                rnaWtFracs(r.tRNAWeightFractionIndexs) / (rnaMWs(r.matureTRNAIndexs)' * rnaExp(r.matureTRNAIndexs));
            
            monCnt = (r.matureRNAGeneComposition(g.mRNAIndexs, r.matureMRNAIndexs) * rnaCnt(r.matureMRNAIndexs)) ./ ...
                (log(2) / t.cellCycleLength + pm.decayRates(pm.matureIndexs));
            monCnt = mass.cellInitialDryWeight * mass.initialFractionAAsInMonomers * mass.dryWeightFractionProtein * monCnt / (monMWs' * monCnt);
            
            idx = find(diff(optGrowths(:, 1)) /  maxGrowth_enzyme > 1e-6, 1, 'last') + 1;
            metMonCnt = r.matureRNAGeneComposition(g.mRNAIndexs, :) * max(r.matureRNAGeneComposition .* repmat(enzGeneComp * squeeze(optEnzymes(idx, 1, :)), 1, size(r.matureRNAGeneComposition, 2)), [], 1)';
            tfs = metMonCnt ~= 0;
            monCnt1 = max(monCnt, metMonCnt);
            monCnt1(~tfs) = monCnt1(~tfs) * (monCnt' * monMWs - monCnt1(tfs)' * monMWs(tfs)) / (monCnt1(~tfs)' * monMWs(~tfs));            
            monCnt = monCnt1;
            
            mrnaExp = ComputationUtil.invertCompositionMatrix(r.matureRNAGeneComposition(g.mRNAIndexs, r.matureMRNAIndexs)) * ...
                (monCnt .* (log(2) / t.cellCycleLength + pm.decayRates(pm.matureIndexs)));
            rnaCnt(r.matureMRNAIndexs) = mrnaExp * ...
                mass.cellInitialDryWeight * mass.initialFractionNTPsInRNAs * mass.dryWeightFractionRNA * ...
                rnaWtFracs(r.mRNAWeightFractionIndexs) / (rnaMWs(r.matureMRNAIndexs)' * mrnaExp);
            rnaExp = rnaCnt / sum(rnaCnt);
            
            paramVec = fitter.constructParameterVector(...
                rnaExp, nmpComp, aaComp, rnaWtFracs, rnaDecayRates, ...
                biomassComposition, biomassProduction, byproducts, unaccountedEnergyConsumption);
            [~, freeMons, freeCpxs] = fitter.calcMacromolecularCounts(paramVec);
            tmpEnzymes = zeros(size(met.enzymeWholeCellModelIDs));
            tmpEnzymes(met.enzymeMonomerLocalIndexs) = freeMons(met.enzymeMonomerGlobalIndexs);
            tmpEnzymes(met.enzymeComplexLocalIndexs) = freeCpxs(met.enzymeComplexGlobalIndexs);
            assertElementsAlmostEqual(maxGrowth_enzyme, ...
                met.calcGrowthRate(met.calcFluxBounds(substrates, tmpEnzymes, met.fbaReactionBounds, met.fbaEnzymeBounds)), ...
                'relative', 1e-2);
            
%             %simulate growth rate distribution
%             paramVecs = {paramVec0; paramVec};
%             meanGrowthRates = mr.meanInitialGrowthRate * [1 maxGrowth_enzyme / initGrowth];
%             mr.initialGrowthFilterWidth = inf;
%             growths = zeros(numel(paramVecs), 100);
%             for i = 1:numel(paramVecs)
%                 fitter.applyParameterVectorToSimulation(paramVecs{i});
%                 mr.meanInitialGrowthRate = meanGrowthRates(i);
%                 for j = 1:100
%                     sim.applyOptions('seed', 100 * i + j);
%                     sim.initializeState();
%                     growths(i, j) = mr.growth;
%                 end
%             end
%             
%             %compare distribution
%             figHandle = figure();
%             
%             axesHandle = subplot(1, 1, 1);
%             edges = linspace(quantile(growths(:), 0.02), quantile(growths(:), 0.98), 10);
%             cnts = [
%                 histc(growths(1, :), edges)
%                 histc(growths(2, :), edges)
%                 ];
%             h = plot(edges, cnts);
%             legend(h, {'Baseline', 'Optimum'}, 'Location', 'NorthWest');
%             ylabel('Freq', 'FontSize', 7)
%             xlabel('Growth (cell s^{-1})', 'FontSize', 7);
%             set(axesHandle, 'FontSize', 5, 'tickDir', 'out', 'box', 'off');
%             
%             saveas(figHandle, [outDir 'growthDistributions.pdf']);
%             close(figHandle);
            
            %% metabolic / non-metabolic tradeoff            
            %max growth supported by various total metabolic enzyme masses
            totMetEnzMass = linspace(...
                ((kineticLimitedEnzymes & kineticUnlimitedEnzymes) + (~kineticLimitedEnzymes) .* enzymes)' * met.enzymeMolecularWeights, ...
                mass.cellInitialDryWeight * mass.dryWeightFractionProtein * ConstantUtil.nAvogadro, ...
                20);
            metEnzGrowths = zeros(numel(totMetEnzMass), 1);            
            for i = 1:numel(totMetEnzMass)
                [result, ~, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                    'maximize', [fbaObj; zeros(nEnz, 1)], [fbaSMat zeros(nMet, nEnz); zeros(1, nRxn) met.enzymeMolecularWeights'; kineticConstraints -enzymeConstraints], ...
                    [fbaRHS; totMetEnzMass(i); zeros(nKinEnz, 1)], ...
                    [nonKineticFluxBounds(:, 1); (kineticLimitedEnzymes & kineticUnlimitedEnzymes) + (~kineticLimitedEnzymes) .* enzymes], ...
                    [nonKineticFluxBounds(:, 2); inf(nEnz, 1)], ...
                    [repmat('S', size(fbaRHS)); 'S'; repmat('U', nKinEnz, 1)], 'C', met.linearProgrammingOptions);
                if errFlag
                    metEnzGrowths(i, 1) = NaN;
                else
                    metEnzGrowths(i, 1) = result(bmRxnIdx);
                end
            end
            
            %max growth supported by various total non-metabolic enzyme masses
            nonMetProts = ~any(enzGeneComp(g.mRNAIndexs, :), 2);
            [~, ~, ~, ~, ~, minMonCnt, ~, ~] = fitter.calcResourceRequirements(fitter.constructParameterVectorFromSimulation());
            totNonMetProtMass = minMonCnt(nonMetProts)' * pm.molecularWeights(pm.matureIndexs(nonMetProts));
            
            %max growth supported by media
            tmpFluxBounds = met.calcFluxBounds(...
                substrates, enzymes, met.fbaReactionBounds, met.fbaEnzymeBounds, ...
                false, true, true, true, true, true);
            [result, ~, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                'maximize', met.fbaObjective, met.fbaReactionStoichiometryMatrix, met.fbaRightHandSide, ...
                tmpFluxBounds(:, 1), tmpFluxBounds(:, 2), ...
                'S', 'C', met.linearProgrammingOptions);
            if errFlag
                throw(MException('Metabolism:error', 'Unable to optimize growth: %s', errMsg));
            end
            mediaGrowth = result(met.fbaReactionIndexs_biomassProduction);
            
            %plot
            figHandle = figure();
            
            axesHandle = subplot(1, 1, 1);
            h  = [
                plot(totMetEnzMass, metEnzGrowths, 'Color', 'b')
                line([0 totNonMetProtMass totNonMetProtMass totMetEnzMass(end)], [0 0 metEnzGrowths(end) metEnzGrowths(end)], 'Color', 'g')
                line(xlim(axesHandle), mediaGrowth * [1 1], 'Color', 'r', 'LineStyle', ':', 'Parent', axesHandle)
                ];
            xlabel(axesHandle, 'Tot Metabolic Enzyme Mass (Da)', 'FontSize', 7);
            ylabel(axesHandle, 'Growth (cell s^{-1}', 'FontSize', 7);
            xlim(axesHandle, [0 totMetEnzMass(end)]);
            ylim(axesHandle, [0 metEnzGrowths(end)]);
            legend(h, {'Metabolism', 'Non-metabolism', 'Media'});
            set(axesHandle, 'FontSize', 5, 'tickdir', 'out', 'box', 'off');
            
            saveas(figHandle, [outDir 'metabolism-nonmetabolism tradeoff.pdf']);
            close(figHandle);
            
            %% two gene tradeoff
            %calculate
            essMetGenes = g.getIndexs({
                'MG_006'; 'MG_013'; 'MG_023'; 'MG_034'; 'MG_037'; 'MG_038'; 'MG_041'; 
                'MG_042'; 'MG_043'; 'MG_044'; 'MG_045'; 'MG_047';
                'MG_053'; 'MG_058'; 'MG_066'; 'MG_069'; 'MG_071'; 
                'MG_077'; 'MG_078'; 'MG_079'; 'MG_080'; 'MG_102'; 'MG_107'; 
                'MG_111'; 'MG_112'; 'MG_114'; 'MG_118'; 'MG_124'; 'MG_128'; 
                'MG_137'; 'MG_145'; 'MG_171'; 'MG_179'; 'MG_180'; 'MG_181'; 
                'MG_212'; 'MG_215'; 'MG_216'; 'MG_228'; 'MG_229'; 'MG_230'; 
                'MG_231'; 'MG_245'; 'MG_270'; 'MG_271'; 'MG_272'; 'MG_273'; 
                'MG_274'; 'MG_275'; 'MG_276'; 'MG_278'; 'MG_287'; 'MG_299'; 
                'MG_300'; 'MG_301'; 'MG_302'; 'MG_303'; 'MG_304'; 'MG_321'; 
                'MG_322'; 'MG_323'; 'MG_330'; 'MG_342'; 'MG_351'; 'MG_357'; 
                'MG_368'; 'MG_382'; 'MG_383'; 'MG_394'; 'MG_396'; 'MG_407'; 
                'MG_429'; 'MG_430'; 'MG_431'; 'MG_434'; 'MG_437'; 'MG_453'; 
                'MG_458'; 'MG_517'});
            possEnz = find(any(enzGeneComp(essMetGenes, :), 1)' & kineticLimitedEnzymes);
            
            selectedEnzIdxs = find(any(enzGeneComp(g.getIndexs({'MG_107'; 'MG_111'}), :), 1));
            assert(all(ismember(selectedEnzIdxs, possEnz)))
            selectedEnz = false(nEnz, 1);
            selectedEnz(selectedEnzIdxs) = true;
            selectedEnz1 = false(nEnz, 1);
            selectedEnz1(selectedEnzIdxs(1)) = true;
            
            nSteps = 20;
            totEnzMass = linspace(0, 0.2 * enzymes(selectedEnz)' * met.enzymeMolecularWeights(selectedEnz), nSteps);
            optEnzExp = zeros(nSteps, 2);
            optGrowth = zeros(nSteps, 2);
            enzExpConstantGrowth = zeros(nSteps, 2, nSteps);
            for i = 1:nSteps
                [result, ~, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                    'maximize', [fbaObj; zeros(nEnz, 1)], ...
                    [fbaSMat zeros(nMet, nEnz); zeros(1, nRxn) selectedEnz'.*met.enzymeMolecularWeights'; kineticConstraints -enzymeConstraints], ...
                    [fbaRHS; totEnzMass(i); zeros(nKinEnz, 1)], ...
                    [nonKineticFluxBounds(:, 1); ~selectedEnz .* enzymes], ...
                    [nonKineticFluxBounds(:, 2); max(~selectedEnz .* enzymes, selectedEnz .* inf(nEnz, 1))], ...
                    [repmat('S', size(fbaRHS)); 'S'; repmat('U', nKinEnz, 1)], 'C', met.linearProgrammingOptions);
                if errFlag
                    throw(MException('Metabolism:error', 'Unable to optimize growth: %s', errMsg));
                end
                optEnzExp(i, :) = result(end - nEnz + selectedEnzIdxs);
                optGrowth(i, 1) = result(bmRxnIdx);

                tmpNonKineticFluxBounds = nonKineticFluxBounds;
                tmpNonKineticFluxBounds(bmRxnIdx, :) = optGrowth(i, 1);
                totEnzMass2 = linspace(totEnzMass(i), 2 * totEnzMass(end), nSteps);
                for j = 1:nSteps/2
                    [result, ~, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                        'minimize', [zeros(nRxn, 1); selectedEnz1], ...
                        [fbaSMat zeros(nMet, nEnz); zeros(1, nRxn) selectedEnz'.*met.enzymeMolecularWeights'; kineticConstraints -enzymeConstraints], ...
                        [fbaRHS; totEnzMass2(j); zeros(nKinEnz, 1)], ...
                        [tmpNonKineticFluxBounds(:, 1); ~selectedEnz .* enzymes], ...
                        [tmpNonKineticFluxBounds(:, 2); max(~selectedEnz .* enzymes, selectedEnz .* inf(nEnz, 1))], ...
                        [repmat('S', size(fbaRHS)); 'S'; repmat('U', nKinEnz, 1)], 'C', met.linearProgrammingOptions);
                    if errFlag
                        throw(MException('Metabolism:error', 'Unable to optimize growth: %s', errMsg));
                    end
                    enzExpConstantGrowth(i, :, nSteps/2 - j + 1) = result(end - nEnz + selectedEnzIdxs);

                    [result, ~, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                        'maximize', [zeros(nRxn, 1); selectedEnz1], ...
                        [fbaSMat zeros(nMet, nEnz); zeros(1, nRxn) selectedEnz'.*met.enzymeMolecularWeights'; kineticConstraints -enzymeConstraints], ...
                        [fbaRHS; totEnzMass2(j); zeros(nKinEnz, 1)], ...
                        [tmpNonKineticFluxBounds(:, 1); ~selectedEnz .* enzymes], ...
                        [tmpNonKineticFluxBounds(:, 2); max(~selectedEnz .* enzymes, selectedEnz .* inf(nEnz, 1))], ...
                        [repmat('S', size(fbaRHS)); 'S'; repmat('U', nKinEnz, 1)], 'C', met.linearProgrammingOptions);
                    if errFlag
                        throw(MException('Metabolism:error', 'Unable to optimize growth: %s', errMsg));
                    end
                    enzExpConstantGrowth(i, :, nSteps/2 + j) = result(end - nEnz + selectedEnzIdxs);
                end
            end
            
            %plot
            figHandle = figure();
            
            axesHandle = subplot(1, 1, 1);
            hold(axesHandle, 'on');
            h  = zeros(3, 1);
            h(1) = plot(axesHandle, optEnzExp(:, 1), optEnzExp(:, 2), 'r');
            h(2) = plot(axesHandle, ...
                reshape([permute(enzExpConstantGrowth(:, 1, :), [1 3 2]) NaN(nSteps, 1)]', 1, []), ...
                reshape([permute(enzExpConstantGrowth(:, 2, :), [1 3 2]) NaN(nSteps, 1)]', 1, []), ...
                'Color', 'b');
            h(3) = plot(axesHandle, ...
                reshape([zeros(nSteps, 1)  2 * totEnzMass' / met.enzymeMolecularWeights(selectedEnzIdxs(1))  NaN(nSteps, 1)]', 1, []), ...
                reshape([2 * totEnzMass' / met.enzymeMolecularWeights(selectedEnzIdxs(2))  zeros(nSteps, 1)  NaN(nSteps, 1)]', 1, []), ...
                'Color', [0.75 0.75 0.75]);
            uistack(h(2), 'top');
            uistack(h(1), 'top');
            xlabel(axesHandle, sprintf('Expression %s', met.enzymeNames{selectedEnzIdxs(1)}), 'FontSize', 7);
            ylabel(axesHandle, sprintf('Expression %s', met.enzymeNames{selectedEnzIdxs(2)}), 'FontSize', 7);
            xlim(axesHandle, [0 totEnzMass(end) / met.enzymeMolecularWeights(selectedEnzIdxs(1))]);
            ylim(axesHandle, [0 totEnzMass(end) / met.enzymeMolecularWeights(selectedEnzIdxs(2))]);
            legend(h, {'Optima', 'Constant Growth', 'Constant Mass'});
            set(axesHandle, 'FontSize', 5, 'tickdir', 'out', 'box', 'off');
            
            saveas(figHandle, [outDir 'enzyme A-B tradeoff.pdf']);
            close(figHandle);                        
            
            %% figure
            figW = 8.5; figH = 10; % Cell 1 column
            [~, figHandle] = PlotUtil.newAxesHandle([figW figH]);
            clf(figHandle);
            
            %%A
            PlotUtil.labelSubFigure('A', [leftColLabelX 0.993 -0.0286 -0.0208], figHandle, 8, 'normalized');
            axesHandle = subplot('Position', [leftColX 0.745 leftColW 0.23], 'Parent', figHandle);
            
            %% standarize tick lengths
            axesHandles = findobj(figHandle, 'type', 'axes', '-and', 'tag', '');
            axesPos = cell2mat(get(axesHandles, 'position'));
            tickLen = 0.01 * max(max(axesPos(:, 3:4), [], 1) .* [figW figH]);
            for i = 1:numel(axesHandles)
                maxVal = max(axesPos(i, 3:4) .* [figW figH]);
                set(axesHandles(i), 'TickLen', [tickLen / maxVal  0.0250]);
            end
            
            %% save
            print(figHandle, [outDirectory filesep 'growthSensitivity.pdf'], '-dpdf', '-rgb');
            close(figHandle);
        end
        
        function [optGrowths, optEnzymes, incEnzymes] = calcOptimalEnzymeExpressionChanges(fbaObj, fbaSMat, fbaRHS, ...
                kineticConstraints, enzymeConstraints, nonKineticFluxBounds, ...
                enzymes, kineticLimitedEnzymes, kineticUnlimitedEnzymes, bmRxnIdx, sim)
            import edu.stanford.covert.util.ComputationUtil;
            
            met = sim.process('Metabolism');
            nMet = size(fbaSMat, 1);
            nRxn = size(fbaSMat, 2);
            nEnz = numel(kineticLimitedEnzymes);
            nKinEnz = sum(kineticLimitedEnzymes);
            Emax = (enzymes' * met.enzymeMolecularWeights) ./ met.enzymeMolecularWeights;
            
            [result, ~, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                'maximize', [fbaObj; zeros(nEnz, 1)], [fbaSMat zeros(nMet, nEnz); zeros(1, nRxn) met.enzymeMolecularWeights'; kineticConstraints -enzymeConstraints], ...
                [fbaRHS; enzymes' * met.enzymeMolecularWeights; zeros(nKinEnz, 1)], ...
                [nonKineticFluxBounds(:, 1); enzymes], [nonKineticFluxBounds(:, 2); inf(nEnz, 1)], ...
                [repmat('S', size(fbaRHS)); 'S'; repmat('U', nKinEnz, 1)], 'C', met.linearProgrammingOptions);
            if errFlag
                throw(MException('Metabolism:error', 'Unable to optimize growth: %s', errMsg));
            end
            initGrowth = result(bmRxnIdx);
            
            [result, ~, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                'maximize', [fbaObj; zeros(nEnz, 1)], [fbaSMat zeros(nMet, nEnz); zeros(1, nRxn) met.enzymeMolecularWeights'; kineticConstraints -enzymeConstraints], ...
                [fbaRHS; enzymes' * met.enzymeMolecularWeights; zeros(nKinEnz, 1)], ...
                [nonKineticFluxBounds(:, 1); (kineticLimitedEnzymes & kineticUnlimitedEnzymes) + (~kineticLimitedEnzymes) .* enzymes], [nonKineticFluxBounds(:, 2); inf(nEnz, 1)], ...
                [repmat('S', size(fbaRHS)); 'S'; repmat('U', nKinEnz, 1)], 'C', met.linearProgrammingOptions);
            if errFlag
                throw(MException('Metabolism:error', 'Unable to optimize growth: %s', errMsg));
            end
            maxGrowth = result(bmRxnIdx);
            
            tmpNonKineticFluxBounds = nonKineticFluxBounds;
            tmpNonKineticFluxBounds(bmRxnIdx) = maxGrowth;
            [result, ~, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                'minimize', [
                zeros(nRxn, 1)
                zeros(nEnz, 1)
                ones(nEnz, 1)
                ], [
                fbaSMat             zeros(nMet, nEnz)             zeros(nMet, nEnz)
                zeros(1, nRxn)      met.enzymeMolecularWeights'   zeros(1, nEnz)
                kineticConstraints  -enzymeConstraints            zeros(nKinEnz, nEnz)
                zeros(nEnz, nRxn)   eye(nEnz)                     -diag(Emax - enzymes)
                zeros(1, nRxn)      zeros(1, nEnz)                ones(1, nEnz)
                ], [
                fbaRHS
                enzymes' * met.enzymeMolecularWeights
                zeros(nKinEnz, 1)
                enzymes
                nEnz
                ], [
                tmpNonKineticFluxBounds(:, 1)
                (kineticLimitedEnzymes & kineticUnlimitedEnzymes) + (~kineticLimitedEnzymes) .* enzymes
                zeros(nEnz, 1)
                ], [
                tmpNonKineticFluxBounds(:, 2)
                inf(nEnz, 1)
                ones(nEnz, 1)
                ], [
                repmat('S', size(fbaRHS))
                'S'
                repmat('U', nKinEnz, 1)
                repmat('U', nEnz, 1)
                'U'
                ], [
                repmat('C', nRxn, 1)
                repmat('C', nEnz, 1)
                repmat('I', nEnz, 1)
                ], met.linearProgrammingOptions);
            if errFlag
                throw(MException('Metabolism:error', 'Unable to optimize growth: %s', errMsg));
            end
            assertElementsAlmostEqual(maxGrowth, result(bmRxnIdx));
            optEnzyme = result((end-nEnz+1:end) - nEnz);
            
            maxModEnz = sum(result(end-nEnz+1:end));
            optGrowths = zeros(nEnz + 1, 2, 1);
            optEnzymes = zeros(nEnz + 1, 2, nEnz);
            incEnzymes = cell(nEnz + 1, 2);
            optGrowths(1, :) = initGrowth;
            optEnzymes(1, 1, :) = enzymes;
            optEnzymes(1, 2, :) = enzymes;
            nIncreases = 1;
            for i = 1:nEnz
                [result, ~, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                    'maximize', [
                    fbaObj
                    zeros(nEnz, 1)
                    zeros(nEnz, 1)
                    ], [
                    fbaSMat             zeros(nMet, nEnz)             zeros(nMet, nEnz)
                    zeros(1, nRxn)      met.enzymeMolecularWeights'   zeros(1, nEnz)
                    kineticConstraints  -enzymeConstraints            zeros(nKinEnz, nEnz)
                    zeros(nEnz, nRxn)   eye(nEnz)                     -diag(Emax - enzymes)
                    zeros(1, nRxn)      zeros(1, nEnz)                ones(1, nEnz)
                    ], [
                    fbaRHS
                    enzymes' * met.enzymeMolecularWeights
                    zeros(nKinEnz, 1)
                    enzymes
                    i
                    ], [
                    nonKineticFluxBounds(:, 1)
                    (kineticLimitedEnzymes & kineticUnlimitedEnzymes) + (~kineticLimitedEnzymes) .* enzymes
                    zeros(nEnz, 1)
                    ], [
                    nonKineticFluxBounds(:, 2)
                    inf(nEnz, 1)
                    ones(nEnz, 1)
                    ], [
                    repmat('S', size(fbaRHS))
                    'S'
                    repmat('U', nKinEnz, 1)
                    repmat('U', nEnz, 1)
                    'U'
                    ], [
                    repmat('C', nRxn, 1)
                    repmat('C', nEnz, 1)
                    repmat('I', nEnz, 1)
                    ], met.linearProgrammingOptions);
                if errFlag
                    throw(MException('Metabolism:error', 'Unable to optimize growth: %s', errMsg));
                end
                optGrowths(i+1, 1) = result(bmRxnIdx);
                optEnzymes(i+1, 1, :) = result((end-nEnz+1:end) - nEnz);
                incEnzymes{i+1, 1} = find(squeeze(optEnzymes(i+1, 1, :)) > enzymes);
                
                [result, ~, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                    'maximize', [
                    fbaObj
                    zeros(nEnz, 1)
                    zeros(nEnz, 1)
                    ], [
                    fbaSMat             zeros(nMet, nEnz)             zeros(nMet, nEnz)
                    zeros(1, nRxn)      met.enzymeMolecularWeights'   zeros(1, nEnz)
                    kineticConstraints  -enzymeConstraints            zeros(nKinEnz, nEnz)
                    zeros(nEnz, nRxn)   eye(nEnz)                     -diag(Emax - squeeze(optEnzymes(i, 2, :)))
                    zeros(1, nRxn)      zeros(1, nEnz)                ones(1, nEnz)
                    ], [
                    fbaRHS
                    enzymes' * met.enzymeMolecularWeights
                    zeros(nKinEnz, 1)
                    squeeze(optEnzymes(i, 2, :))
                    nIncreases
                    ], [
                    nonKineticFluxBounds(:, 1)
                    (kineticLimitedEnzymes & kineticUnlimitedEnzymes) + (~kineticLimitedEnzymes) .* enzymes
                    zeros(nEnz, 1)
                    ], [
                    nonKineticFluxBounds(:, 2)
                    inf(nEnz, 1)
                    ones(nEnz, 1)
                    ], [
                    repmat('S', size(fbaRHS))
                    'S'
                    repmat('U', nKinEnz, 1)
                    repmat('U', nEnz, 1)
                    'U'
                    ], [
                    repmat('C', nRxn, 1)
                    repmat('C', nEnz, 1)
                    repmat('I', nEnz, 1)
                    ], met.linearProgrammingOptions);
                if errFlag
                    throw(MException('Metabolism:error', 'Unable to optimize growth: %s', errMsg));
                end
                if result(bmRxnIdx) - optGrowths(i, 2) > 1e-6 * maxGrowth
                    optGrowths(i+1, 2) = result(bmRxnIdx);
                    optEnzymes(i+1, 2, :) = result((end-nEnz+1:end) - nEnz);
                    incEnzymes{i+1, 2} = find(optEnzymes(i+1, 2, :) > optEnzymes(i, 2, :));
                    nIncreases = 1;
                else
                    optGrowths(i+1, 2) = optGrowths(i, 2);
                    optEnzymes(i+1, 2, :) = optEnzymes(i, 2, :);
                    incEnzymes{i+1, 2} = [];
                    nIncreases = nIncreases + 1;
                end
            end
            optGrowths(maxModEnz+1:end, 1) = maxGrowth;
            optEnzymes(maxModEnz+1:end, 1, :) = repmat(permute(optEnzyme, [3 2 1]), nEnz - maxModEnz + 1, 1);
            
            idx = find(any(diff(optGrowths, 1, 1) > 1e-6 * maxGrowth, 2), 1, 'last') + 1;
            optGrowths = optGrowths(1:idx, :);
            optEnzymes = optEnzymes(1:idx, :, :);
            incEnzymes = incEnzymes(1:idx, :);
        end
        
        function RNAvProtein(simBatchDir, outDirectory, figHandle, position)
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.MacromoleculeUtil;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            selectedSimulations = 1:SimulationDiskUtil.getNumSimulations(simBatchDir);
            
            %% get data
            if exist([outDirectory 'RNAvProtein.mat'], 'file')
                load([outDirectory 'RNAvProtein.mat']);
            else
                [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
                g = sim.gene;
                pm = sim.state('ProteinMonomer');
                pc = sim.state('ProteinComplex');
                rna = sim.state('Rna');
                
                geneId = 'MG_218';
                [~, ~, ~, ~, ~, ...
                    ~, ~, ~, ~, ...
                    ~, ~, ~, ~, ...
                    matureRnaIdxs, ~, ~, ~, ...
                    monomerIdxs, ~, ~, ~, ...
                    ~, ~, ~, ~] = ...
                    MacromoleculeUtil.getMacromoleculeIndexsIDsNames(geneId, sim);
                
                stateNames = {
                    'Rna'               'counts' rna.matureIndexs(matureRnaIdxs)         '-sum'
                    'ProteinMonomer'    'counts' pm.matureIndexs(monomerIdxs)            '-sum'
                    };
                cr = zeros(numel(selectedSimulations), 1);
                cp = zeros(numel(selectedSimulations), 1);
                for i = 1:numel(selectedSimulations)
                    states = SimulationEnsemble.load(simBatchDir, {'Time' 'values'}, [], [], 1, 'extract', selectedSimulations(i));
                    
                    tIdx = randi(numel(states.Time.values), 1);
                    states = SimulationEnsemble.load(simBatchDir, stateNames, tIdx, tIdx, 1, 'extract', selectedSimulations(i));
                    cr(i) = states.Rna.counts;
                    cp(i) = states.ProteinMonomer.counts;
                end
                
                save([outDirectory 'RNAvProtein.mat'], 'cp', 'cr');
            end
            
            %% plot
            if nargin < 4
                [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
                position = [0.1 0.1 0.8 0.8];
            end
            x = position(1);
            y = position(2);
            w = position(3);
            h = position(4);
            figPos = get(figHandle, 'papersize');
            figW = figPos(1);
            figH = figPos(2);
            axesHandles = [
                subplot('Position', [x, y + h-0.15*w*figW/figH, 0.8*w, 0.15*w*figW/figH], 'Parent', figHandle)
                subplot('Position', [x + 0.85*w, y, 0.15*w, h-0.2*w*figW/figH], 'Parent', figHandle)
                subplot('Position', [x, y, 0.8*w, h-0.2*w*figW/figH], 'Parent', figHandle)
                ];
            
            xlims = [min(cr) - 0.5, max(cr) + 0.5];
            ylims = [0 60];
            
            %RNA marginal
            parmhat = gamfit(cr);
            
            axesHandle = axesHandles(1);
            hold(axesHandle, 'on');
            [n, edges] = hist(axesHandle, cr, min(cr):max(cr));
            n = n / sum(n) * 100;
            bar(axesHandle, edges, n, 'b', 'EdgeColor', 'none');
            x = linspace(min(cr), max(cr), 20);
            plot(axesHandle, x, gampdf(x, parmhat(1), parmhat(2)) * 100, 'r');
            ylim(axesHandle, [0 max(n)]);
            xlim(axesHandle, xlims);
            set(axesHandle, 'Box', 'off');
            set(axesHandle, 'TickDir', 'out');
            set(axesHandle, 'XTickLabel', '');
            set(axesHandle, 'XTick', []);
            set(axesHandle, 'YTick', 0:25:50);
            
            tick2text(axesHandle, 'axis', 'y', 'ytickoffset', 0.04);
            yTicks = getappdata(axesHandle, 'YTickText');
            set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'right');
            delete(yTicks(2));
            
            ylabel(axesHandle, 'Freq', 'FontSize', 7);
            ylabelPos = get(get(axesHandle, 'ylabel'), 'position');
            set(get(axesHandle, 'ylabel'), 'position', [-0.83 25 ylabelPos(3:end)]);
            
            %protein marginal
            parmhat = gamfit(cp);
            
            axesHandle = axesHandles(2);
            hold(axesHandle, 'on');
            [n, xout] = hist(cp, 50);
            n = n / sum(n) * 100;
            barh(axesHandle, xout, n, 'b', 'EdgeColor', 'none');
            x = linspace(0, 60, 100);
            plot(axesHandle, gampdf(x, parmhat(1), parmhat(2)) * 100, x, 'r');
            xlim(axesHandle, [0 max(n)]);
            ylim(axesHandle, ylims);
            set(axesHandle, 'Box', 'off');
            set(axesHandle, 'TickDir', 'out');
            set(axesHandle, 'YTickLabel', []);
            set(axesHandle, 'YTick', []);
            set(axesHandle, 'XTick', 0:2:6);
            
            tick2text(axesHandle, 'axis', 'x', 'xtickoffset', 0.062);
            xTicks = getappdata(axesHandle, 'XTickText');
            set(xTicks, 'FontSize', 5);
            delete(xTicks(2:3));
            
            xlabel(axesHandle, 'Freq', 'FontSize', 7);
            xlabelPos = get(get(axesHandle, 'xlabel'), 'position');
            set(get(axesHandle, 'xlabel'), 'position', [xlabelPos(1)-0.3 -6 xlabelPos(3:end)]);
            
            %RNA vs protein
            axesHandle = axesHandles(3);
            plot(axesHandle, cr, cp, 'b.');
            xlim(axesHandle, xlims);
            ylim(axesHandle, ylims);
            set(axesHandle, 'XTick', [0 1 2], 'YTick', 0:25:50);
            set(axesHandle, 'TickDir', 'out', 'box', 'off');
            
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', 0.062, 'ytickoffset', 0.04);
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            set(xTicks, 'FontSize', 5);
            set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'right');
            delete(yTicks(2));
            
            xlabel(axesHandle, 'mRNA count', 'FontSize', 7);
            ylabel(axesHandle, 'Protein count', 'FontSize', 7);
            xlabelPos = get(get(axesHandle, 'xlabel'), 'position');
            ylabelPos = get(get(axesHandle, 'ylabel'), 'position');
            set(get(axesHandle, 'xlabel'), 'position', [xlabelPos(1) -6 xlabelPos(3:end)]);
            set(get(axesHandle, 'ylabel'), 'position', [-0.83 ylabelPos(2:end)]);
        end
        
        function cellShape(outDirectory, simBatchDir, selectedSim)
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.analysis.Figures34;
            import edu.stanford.covert.cell.sim.state.CellGeometry;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if ~exist('selectedSim', 'var') || isempty(selectedSim)
                selectedSim = 1:SimulationDiskUtil.getNumSimulations(simBatchDir);
            end
            
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            
            stateNames={'Geometry' 'width'
                'Geometry' 'pinchedDiameter'
                'Chromosome' 'polymerizedRegions'
                'Chromosome' 'monomerBoundSites'
                'Chromosome' 'complexBoundSites'
                'Time' 'values'};
            
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSim);
            
            log = load([simBatchDir filesep num2str(selectedSim) filesep 'summary.mat']);
            
            cm = sim.state('Mass');
            geom = sim.state('Geometry');
            density = geom.density;
            mass = (log.mass(2:end) * cm.cellInitialDryWeight)/(1-cm.fractionWetWeight);
            width = squeeze(states.Geometry.width);
            pinchedDiameter = squeeze(states.Geometry.pinchedDiameter);
            t = squeeze(states.Time.values);
            septumLength = (width - pinchedDiameter)/2;
            
            cylindricalLength = zeros(size(t));
            maxTotalLengths = zeros(size(t));
            
            for i=1:length(cylindricalLength)
                if(pinchedDiameter(i) > 0)
                    [cylindricalLength(i), ~, maxTotalLengths(i)] = CellGeometry.calculateGeometry(...
                        width(i), ...
                        pinchedDiameter(i), ...
                        mass(i)/density);
                end
            end
            
            tidxs = [round(1:t(end)/4:t(end)) t(end)-1];
            
            scale = 1e-9;
            
            for i = 1:numel(tidxs)
                
                svg_width = maxTotalLengths(tidxs(i)) / scale + 10;
                svg_height = max(width) / scale + 10;
                
                fileNameSVG = sprintf('%s%sCellShape_%d.svg', outDirectory, filesep, i);
                fileNamePDF = sprintf('%s%sCellShape_%d.pdf', outDirectory, filesep, i);
                
                %open svg
                fid = fopen(fileNameSVG, 'w');
                fprintf(fid, '<?xml version="1.0" standalone="no"?>\n');
                fprintf(fid, '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n');
                fprintf(fid, '<svg width="%d" height="%d" version="1.1" xmlns="http://www.w3.org/2000/svg" viewBox="%d %d %d %d">\n', svg_width, svg_height, 0, 0, svg_width, svg_height);
                fprintf(fid, '<defs>\n');
                fprintf(fid, '<linearGradient id="cellBg" x1="0%%" y1="0%%" x2="0%%" y2="100%%">\n');
                fprintf(fid, '<stop style="stop-color:#97b0fb;stop-opacity:0.05" offset="0"/>\n');
                fprintf(fid, '<stop style="stop-color:#89b8d6;stop-opacity:1" offset="1"/>\n');
                fprintf(fid, '</linearGradient>');
                fprintf(fid, '</defs>\n');
                
                % White background
                fprintf(fid, '<rect x="%d" y="%d" width="%d" height="%d" style="fill:white; stroke:none"/>\n', ...
                    0, 0, svg_width, svg_height);
                
                fprintf(fid, Figures34.drawCell(scale,width(tidxs(i)),cylindricalLength(tidxs(i)),septumLength(tidxs(i)),max(maxTotalLengths),max(width)));
                
                %close svg
                fprintf(fid, '</svg>\n');
                fclose(fid);
                
                [status, result] = system(sprintf('inkscape "%s" --export-pdf="%s"', fileNameSVG, fileNamePDF));
                if status ~= 0
                    throw(MException('Figures34:error', 'Failed to convert svg to pdf: %s', result));
                end
            end
        end
        
        function normCellMass(simBatchDir, selectedSim, outDirectory, figHandle, position)
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if ~exist('selectedSim', 'var') || isempty(selectedSim)
                selectedSim = 1:SimulationDiskUtil.getNumSimulations(simBatchDir);
            end
            
            %% get constants
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            comp = sim.compartment;
            
            %% get data
            if exist([outDirectory 'normCellMass.mat'], 'file')
                load([outDirectory 'normCellMass.mat'])
            else
                stateNames = {
                    'Time' 'values'       ':' ':'
                    'Mass' 'dnaWt'        ':' '-sum'
                    'Mass' 'rnaWt'        ':' '-sum'
                    'Mass' 'proteinWt'    ':' '-sum'
                    'Mass' 'metaboliteWt' ':' [comp.cytosolIndexs; comp.membraneIndexs]
                    'Mass' 'cellDry'         ':' '-sum'
                    };
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSim);
                
                proteinWt = permute(states.Mass.proteinWt, [1 3 2]);
                dnaWt = permute(states.Mass.dnaWt, [1 3 2]);
                membraneWt = permute(states.Mass.metaboliteWt(:, 2, :), [1 3 2]);
                rnaWt = permute(states.Mass.rnaWt, [1 3 2]);
                cellWt = permute(states.Mass.cellDry, [1 3 2]);
                time = permute(states.Time.values, [1 3 2]) / 3600;
                
                proteinWt = proteinWt / proteinWt(1);
                dnaWt = dnaWt / dnaWt(1);
                membraneWt = membraneWt / membraneWt(1);
                rnaWt = rnaWt / rnaWt(1);
                cellWt = cellWt / cellWt(1);
                
                save([outDirectory 'normCellMass.mat'], 'proteinWt', 'dnaWt', 'membraneWt', 'rnaWt', 'cellWt', 'time');
            end
            
            %% plot
            if nargin >= 4
                axesHandle = subplot('position', position, 'parent', figHandle);
            else
                axesHandle = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            end
            
            lh = plot(axesHandle, time, [cellWt; dnaWt; rnaWt; proteinWt; membraneWt]);
            set(lh(2), 'LineWidth', 1, 'Color', 'g');
            uistack(lh(4), 'top');
            uistack(lh(1), 'top');
            uistack(lh(2), 'bottom');
            uistack(lh(3), 'bottom');
            
            xlim(axesHandle, [0 time(end)]);
            ylim(axesHandle, [0.92 2.05]);
            
            set(axesHandle, 'Box', 'off');
            set(axesHandle, 'TickDir', 'out');
            set(axesHandle, 'YTick', 1:0.5:2);
            set(axesHandle, 'XTick', 0:2:8);
            set(axesHandle, 'FontSize', 5);
            
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', 0.05, 'ytickoffset', 0.015);
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            set(xTicks, 'FontSize', 5);
            set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'right');
            delete(yTicks(2));
            delete(xTicks(2:2:end));
            
            xlabel(axesHandle, 'Time (h)', 'FontSize', 7);
            ylabel(axesHandle, 'Mass (norm)', 'FontSize', 7);
            xlabelPos = get(get(axesHandle, 'xlabel'), 'position');
            ylabelPos = get(get(axesHandle, 'ylabel'), 'position');
            set(get(axesHandle, 'xlabel'), 'position', [xlabelPos(1) 0.82 xlabelPos(3:end)]);
            set(get(axesHandle, 'ylabel'), 'position', [-0.50 ylabelPos(2:end)]);
            
            legendHandle = legend(lh, {'Total', 'DNA', 'RNA', 'Protein', 'Membrane'}, 'Location', 'NorthWest', 'FontSize', 7);
            legend(axesHandle, 'boxoff');
            legendPos = get(legendHandle, 'Position');
            set(legendHandle, 'Position', [0.45 1.05 legendPos(3:end)])
            set(legendHandle, 'units', 'pixels');
            legendPos = get(legendHandle, 'Position');
            set(legendHandle, 'Position', [legendPos(1:3) 49]);
            legendChildren = get(legendHandle, 'children');
            set(legendChildren(1:3:end), 'visible', 'off');
            for i = 1:3:numel(legendChildren)
                set(legendChildren(i+1), 'xdata', [0.03 0.10], 'LineWidth', 1)
                
                pos = get(legendChildren(i+2), 'Position');
                set(legendChildren(i+2), 'Position', [0.13 pos(2) 0])
            end
        end
        
        function populationMass(simBatchDir, selectedSim, figHandle, position, outDirectory)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            
            if exist([outDirectory 'populationMassDistribution.mat'], 'file')
                load([outDirectory 'populationMassDistribution.mat']);
            else
                data = load([SimulationDiskUtil.getSimulationBatchDir([simBatchDir filesep num2str(selectedSim)]) filesep 'population-MassDistribution.mat'], 'mass', 'time');
                mass = data.mass;
                time = data.time;
                save([outDirectory 'populationMassDistribution.mat'], 'mass', 'time');
            end
            
            % Determine division times
            simEndTimes = zeros(1, size(mass, 2));
            for i = 1:size(mass, 2)
                et = find(isnan(mass(:, i)), 1);
                if isempty(et)
                    et = NaN;
                end
                simEndTimes(i) = et - 1;
            end
            
            finalMass = zeros(1, size(mass, 2));
            for i = 1:size(finalMass, 2)
                if isnan(simEndTimes(i))
                    finalMass(i) = mass(end, i);
                else
                    finalMass(i) = mass(simEndTimes(i), i);
                end
            end
            
            %markIdx = find(simEndTimes == floor(nanmedian(simEndTimes)));
            markIdx = 10; %52;
            
            %pick every other cell
            simEndTimes = simEndTimes(2:2:end);
            mass = mass(:, 2:2:end);
            markIdx = markIdx / 2;
            
            %% plot
            axesSizes = [50; 5];
            
            options = struct;
            if nargin >= 4
                options.position = position;
            else
                [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
                clf(figHandle);
            end
            [axesHandles, xAxisHandle] = PlotUtil.multiElementPlot(figHandle, axesSizes, time([1 end])', options);
            xlim(xAxisHandle, [0 max(time)]);
            set(get(xAxisHandle, 'xlabel'), 'FontSize', 7);
            set(xAxisHandle, 'FontSize', 5);
            set(xAxisHandle, 'XTick', 0:5:10);
            tick2text(xAxisHandle, 'axis', 'x', 'xtickoffset', 50);
            xTicks = getappdata(xAxisHandle, 'XTickText');
            set(xTicks, 'FontSize', 5);
            xlabelPos = get(get(xAxisHandle, 'XLabel'), 'Position');
            set(get(xAxisHandle, 'XLabel'), 'Position', [xlabelPos(1) -120 xlabelPos(3:end)]);
            
            % Downsample
            [~, simOrder] = sort(mass(1, :));
            timeds = time(1:100:end);
            massds = mass(1:100:end, :);
            
            axesHandle = axesHandles(1);
            cla(axesHandle)
            hold(axesHandle, 'on');
            n1 = floor(numel(simOrder) / 2);
            n2 = ceil(numel(simOrder) / 2);
            color1 = 0.50 * [1 1 1];
            color2 = 0.25 * [1 1 1];
            color3 = 0.00 * [1 1 1];
            colorOrder = zeros(numel(simOrder), 3);
            colorOrder(simOrder, :) = [
                repmat(linspace(1, 0, n1)', 1, 3) .* repmat(color1, n1, 1)  +  repmat(linspace(0, 1, n1)', 1, 3) .* repmat(color2, n1, 1)
                repmat(linspace(1, 0, n2)', 1, 3) .* repmat(color2, n2, 1)  +  repmat(linspace(0, 1, n2)', 1, 3) .* repmat(color3, n2, 1)
                ];
            set(axesHandle, 'ColorOrder', colorOrder);
            h = plot(axesHandle, timeds, massds, 'LineWidth', 0.25);
            set(h(markIdx), 'Color', 'r', 'LineWidth', 2);
            uistack(h(markIdx), 'top');
            set(axesHandle, 'XLim', [0 max(time)], 'FontSize', 5);
            set(axesHandle, 'YTick', 20:10:30);
            tick2text(axesHandle, 'axis', 'y', 'ytickoffset', 0.035);
            yTicks = getappdata(axesHandle, 'YTickText');
            set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'right');
            ylabel(axesHandle, 'Mass (fg)', 'FontSize', 7);
            ylabelPos = get(get(axesHandle, 'ylabel'), 'position');
            set(get(axesHandle, 'ylabel'), 'position', [-1.05 ylabelPos(2:end)]);
            ylim(axesHandle, [min(massds(:)), max(massds(:))]);
            
            axesHandle = axesHandles(2);
            hold(axesHandle, 'on');
            
            edges = linspace(time(1), time(end), 100);
            freq = histc(simEndTimes / 3600, edges) / numel(simEndTimes) * 100;
            freq = smooth(freq, 9, 'lowess');
            line(edges, freq, 'Color', 'k', 'Parent', axesHandle)
            
            idx = find(edges <= simEndTimes(markIdx) / 3600, 1, 'last');
            y = freq(idx) + diff(freq(idx:idx+1)) / diff(edges(idx:idx+1)) * (simEndTimes(markIdx) / 3600 - edges(idx));
            
            axesPos = get(axesHandle, 'Position');
            arrowX = axesPos(1) + simEndTimes(markIdx) / 3600 / time(end) * axesPos(3);
            arrowY1 = axesPos(2);
            axesPos = get(axesHandles(1), 'Position');
            arrowY2 = axesPos(2) + axesPos(4)*(massds(floor(simEndTimes(markIdx)/100), markIdx) - min(massds(:))) / range(massds(:));
            y = linspace(arrowY1, arrowY2-0.01, 81);
            for i = 1:3:numel(y)
                if i >= 30 && i <= 50
                    continue;
                end
                h = annotation(get(axesHandle, 'Parent'), 'line', arrowX * [1 1], y([i+1 i]), ...
                    'Color', 'r');
                uistack(h, 'bottom');
            end
            
            annotation('textbox', [arrowX-0.06 y(26)+0.001 0.12 0.05], 'String', 'Median cell', ...
                'HorizontalAlign', 'center', 'verticalAlign', 'middle', ....
                'FontSize', 7, 'EdgeColor', 'none', 'BackgroundColor', 'none', 'Color', 'r');
            annotation('textbox', [arrowX-0.06 y(26)-0.014 0.12 0.05], 'String', sprintf('\\tau = %0.1f h', 8.9), ...
                'HorizontalAlign', 'center', 'verticalAlign', 'middle', ....
                'FontSize', 7, 'EdgeColor', 'none', 'BackgroundColor', 'none', 'Color', 'r');
            
            ylim(axesHandle, [0 max(freq)]);
            xlim(axesHandle, [time(1) time(end)]);
            
            ylabel(axesHandle, {'% Cell div'}, 'FontSize', 7);
            ylabelPos = get(get(axesHandle, 'ylabel'), 'position');
            set(get(axesHandle, 'ylabel'), 'position', [-1.05 ylabelPos(2:end)]);
            
            set(axesHandle, 'YTick', 0:4:8);
            set(axesHandle, 'XTickLabel','');
            set(axesHandle, 'XTick', []);
            set(axesHandle, 'Box', 'off');
            set(axesHandle, 'XLim', [0 max(time)]);
            
            tick2text(axesHandle, 'axis', 'y', 'ytickoffset', 0.035);
            yTicks = getappdata(axesHandle, 'YTickText');
            set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'right');
            delete(yTicks(2));
            
            %% offset axes
            PlotUtil.offsetYAxes(axesHandles, 0.02);
        end
        
        function growthMeasurement(figHandle, position, outDirectory)
            %% parameters
            [r, host] = system('hostname');
            
            if r ~= 0
                warning('WholeCell:warning:network', 'Could not determine host name');
                value = '';
                return;
            end
            host = deblank(host);
            switch host
                % Edit as appropriate (add jayodita and markus)
                case 'hydra'
                    baseDir = '/covertlab/home/share/Mycoplasma genitalium/MGenGrowthCurves';
                case 'silico'
                    baseDir = '/covertlab/home/share/Mycoplasma genitalium/MGenGrowthCurves';
                case 'covertlab-jkarr'
                    baseDir = 'Z:\share\Mycoplasma genitalium\MGenGrowthCurves';
                otherwise
                    if strcmp(computer, 'GLNXA64')
                        % Probably the cluster
                        baseDir = '/home/share/Mycoplasma genitalium/MGenGrowthCurves';
                    else
                        % Assume a PC with the correct Z drive mapping
                        baseDir = 'Z:\share\Mycoplasma genitalium\MGenGrowthCurves';
                        warning('WholeCell:warning:network', 'Assuming M. genitalium growth curves are in %s', baseDir);
                    end
            end
            
            plate = 93;
            condition = 1;
            replicate = 2;
            blankCondition = 2;
            blankReplicate = 1;
            nDilutions = 6;
            
            nReplicates = 3;
            nColumns = 10;
            excelRows = 3:8;
            %4: 1X
            %5: 1/5X
            %6: 1/25X
            excelCols = 3:12;
            %3: Cond 1, Rep 1;
            %4: Cond 1, Rep 2;
            %5: Cond 1, Rep 3;
            %6: Cond 2, Rep 1;
            %7: Cond 2, Rep 2;
            %8: Cond 2, Rep 3;
            %9: Cond 3, Rep 1;
            %10: Cond 3, Rep 2;
            %11: Cond 3, Rep 3;
            %12: Blank
            
            %% get file name and times
            files = dir(sprintf('%s%sPlate%d', baseDir, filesep, plate));
            files = files(3:end);
            [~, order] = sort(cellfun(@(x) str2double(x(end-5:end-4)), {files.name}));
            files = files(order);
            times = cellfun(@(x) datenum(...
                2000 + str2double(x(end-9:end-8)), ...   %year
                str2double(x(end-15:end-14)), ... %month
                str2double(x(end-12:end-11)), ... %day
                str2double(x(end-21:end-20)), ... %hour
                str2double(x(end-18:end-17)), ... %minute
                0 ... %second
                ), {files.name});
            times = times - times(1);
            
            %% get data
            %from plate reader files
            if exist([outDirectory 'growthMeasurement.mat'], 'file')
                load([outDirectory 'growthMeasurement.mat'])
            else
                data = zeros(nDilutions, nColumns, numel(files));
                for i = 1:numel(files)
                    [~, ~, raw] = xlsread(sprintf('%s%sPlate%d%s%s', baseDir, filesep, plate, filesep, files(i).name));
                    data(:, :, i) = cell2mat(raw(excelRows, excelCols));
                end
                save([outDirectory 'growthMeasurement.mat'], 'data');
            end
            
            % from Jayodita's excel sheet
            [~, ~, raw] = xlsread([outDirectory 'growthMeasurement.xls']);
            tfRows = [false; false; ~isnan(cell2mat(raw(3:end, 1)))];
            times = cell2mat(raw(tfRows, 1))';
            times = (times - times(1)) / 24;
            data = zeros(nDilutions, nColumns, numel(times));
            data(3:5, (condition-1) * nReplicates + replicate, :) = permute(cell2mat(raw(tfRows, 2:4)), [2 3 1]);
            data(2:6, (blankCondition-1) * nReplicates + blankReplicate, :) = repmat(permute(cell2mat(raw(tfRows, 5)), [2 3 1]), [5 1 1]);
            
            %% fit data
            t = NaN(nDilutions, 1);
            for i = 1:nDilutions
                tmpData = permute(data(i, (condition-1) * nReplicates + replicate, :), [1 3 2]);
                idx1 = find(tmpData > 0.35, 1, 'last');
                idx2 = find(tmpData < 0.2, 1, 'first');
                if ~isempty(idx1) && ~isempty(idx2)
                    t(i) = interp1(tmpData(:, idx1:idx2), times(:, idx1:idx2) * 24, 0.28, 'linear');
                end
            end
            deltaT = diff(t);
            tau = log(2) ./ (log(5) ./ deltaT);
            
            %% plot
            if nargin >= 2
                axesHandle = subplot('position', position, 'Parent', figHandle);
            else
                axesHandle = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            end
            hold(axesHandle, 'on');
            
            colors = 1 / 255 * [
                25 25 112
                0 0 255
                30 144 255
                139 137 137
                ];
            
            h = plot(axesHandle, times, permute(data(3:5, (condition-1) * nReplicates + replicate, :), [1 3 2]), 'LineWidth', 1); %dilutions
            plot(axesHandle, times, permute(mean(data(2:6, (blankCondition-1) * nReplicates + blankReplicate, :), 1), [1 3 2]), 'Color', colors(4, :), 'LineWidth', 1); %blank
            set(h(1), 'Color', colors(1, :)); %1x
            set(h(2), 'Color', colors(2, :)); %1/5 X
            set(h(3), 'Color', colors(3, :)); %1/25 X
            xlabel(axesHandle, 'Time (d)', 'FontSize', 7);
            ylabel(axesHandle, 'OD550', 'FontSize', 7);
            
            set(axesHandle, 'FontSize', 5, 'TickDir', 'out', 'YDir', 'Reverse', 'XTick', 0:5:15, 'YTick', 0:.2:0.6);
            box(axesHandle, 'off');
            xlim(axesHandle, [0 17]);
            ylim(axesHandle, [0 0.65]);
            
            text(0.3, 0.048, '\tau =', 'FontSize', 6, 'Interpreter', 'tex')
            text(4.2, 0.0125, 'ln(2) {\Delta}{\it{t}}', 'FontSize', 6, 'Interpreter', 'tex', 'HorizontalAlign', 'center')
            text(4.2, 0.0785, 'ln(dilution factor)', 'FontSize', 6, 'Interpreter', 'tex', 'HorizontalAlign', 'center')
            line([1.4 7], 0.051 * [1 1], 'Parent', axesHandle, 'Color', 'k');
            
            text(0.3, .35, '1X dilution', 'Parent', axesHandle, 'FontSize', 6, 'VerticalAlign', 'middle', 'Color', colors(1, :))
            text(0.3, .40, '5X dilution', 'Parent', axesHandle, 'FontSize', 6, 'VerticalAlign', 'middle', 'Color', colors(2, :))
            text(0.3, .45, '25X dilution', 'Parent', axesHandle, 'FontSize', 6, 'VerticalAlign', 'middle', 'Color', colors(3, :))
            text(0.3, .50, 'Blank', 'Parent', axesHandle, 'FontSize', 6, 'VerticalAlign', 'middle', 'Color', colors(4, :))
            
            axesPos = get(axesHandle, 'Position');
            
            y = 0.28;
            x = t(3)/24-0.3;
            text(x, y - 0.025, sprintf('{\\Delta}{\\it{t}} = %0.1f h', deltaT(3)), 'Parent', axesHandle, 'FontSize', 6, 'HorizontalAlign', 'right', 'verticalAlign', 'middle', 'Interpreter', 'tex')
            text(x-0.42, y + 0.025, sprintf('\\tau = %0.1f h', tau(3)), 'Parent', axesHandle, 'FontSize', 6, 'HorizontalAlign', 'right', 'verticalAlign', 'middle', 'Interpreter', 'tex')
            annotation('doublearrow', axesPos(1) + axesPos(3) * t([3 4])/24 / max(xlim(axesHandle)), ...
                axesPos(2) + axesPos(4) * (max(ylim(axesHandle))-y) / range(ylim(axesHandle)) * [1 1], ...
                'HeadLength', 1.5, 'HeadWidth', 2, 'HeadStyle', 'vback1')
            
            y = 0.28;
            x = t(5)/24+0.3;
            text(x, y - 0.025, sprintf('{\\Delta}{\\it{t}} = %0.1f h', deltaT(4)), 'Parent', axesHandle, 'FontSize', 6, 'HorizontalAlign', 'left', 'verticalAlign', 'middle', 'Interpreter', 'tex')
            text(x + 0.39, y + 0.025, sprintf('\\tau = %0.1f h', tau(4)), 'Parent', axesHandle, 'FontSize', 6, 'HorizontalAlign', 'left', 'verticalAlign', 'middle', 'Interpreter', 'tex')
            annotation('doublearrow', axesPos(1) + axesPos(3) * t([4 5])/24 / max(xlim(axesHandle)), ...
                axesPos(2) + axesPos(4) * (max(ylim(axesHandle))-y) / range(ylim(axesHandle)) * [1 1], ...
                'HeadLength', 1.5, 'HeadWidth', 2, 'HeadStyle', 'vback1')
            
            annotation('textbox', [0.37 0.8525+0.001 0.12 0.05], 'String', 'Mean', ...
                'HorizontalAlign', 'center', 'verticalAlign', 'middle', ....
                'FontSize', 7, 'EdgeColor', 'none', 'BackgroundColor', 'none', 'Color', 'r');
            annotation('textbox', [0.37 0.8525-0.014 0.12 0.05], 'String', sprintf('\\tau = %0.1f h', 9.03), ... %Jayodita: where do this come from?
                'HorizontalAlign', 'center', 'verticalAlign', 'middle', ....
                'FontSize', 7, 'EdgeColor', 'none', 'BackgroundColor', 'none', 'Color', 'r');
            
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', -1.05, 'ytickoffset', 0.015);
            xlabelPos = get(get(axesHandle, 'xlabel'), 'position');
            ylabelPos = get(get(axesHandle, 'ylabel'), 'position');
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            set(xTicks, 'FontSize', 5);
            set(yTicks, 'FontSize', 5, 'horizontalAlign', 'right');
            set(get(axesHandle, 'xlabel'), 'position', [xlabelPos(1) 0.71 xlabelPos(3:end)]);
            set(get(axesHandle, 'ylabel'), 'position', [-1.15 ylabelPos(2:end)]);
        end
        
        function doublingTime(outDirectory)
            % From
            % /home/projects/WholeCell/simulation/output/runSimulation/singleGeneDeletions/summary.xls
            modelMassDTMean = 8.92200510884003; % Col W
            modelMassDTStddev = 0.807282508239662; % Col X
            expMassDTMean = 9.03411809133663; % Col AX
            expMassDTStderr = 0.889827929955399; % Col AY
            
            [axesHandle, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            hold('on');
            bar(1, modelMassDTMean);
            bar(2, expMassDTMean, 'k');
            eh1 = errorbar([modelMassDTMean nan], [modelMassDTStddev nan], '.g');
            eh2 = errorbar([nan expMassDTMean], [nan expMassDTStderr], '.r');
            hold('off');
            legend([eh1 eh2], 'Standard Deviation', 'Standard Error', 'Location', 'BestOutside');
            ylabel('Mass Doubling Time (h)');
            set(axesHandle, 'XTick', [1 2]);
            set(axesHandle, 'XTickLabel', {'Model' 'Experiment'});
            set(axesHandle, 'YTick', [0 2 4 6 8 10]);
            set(axesHandle, 'Box', 'off');
            set(axesHandle, 'TickDir', 'out');
            
            set(figHandle, 'Renderer', 'Painters')
            print(figHandle, [outDirectory filesep 'MassDoublingTime.eps'], '-depsc', '-rgb');
        end
        
        function dryBiomassComposition(simBatchDir, figHandle, position, outDirectory)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            selectedSimulations = 1:SimulationDiskUtil.getNumSimulations(simBatchDir);
            
            %% get data
            if exist([outDirectory 'dryBiomassComposition.mat'], 'file')
                load([outDirectory 'dryBiomassComposition.mat'])
            else
                [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
                comp = sim.compartment;
                
                stateNames = {
                    'Time' 'values'       ':' ':'
                    'Mass' 'dnaWt'        ':' '-sum'
                    'Mass' 'rnaWt'        ':' '-sum'
                    'Mass' 'proteinWt'    ':' '-sum'
                    'Mass' 'metaboliteWt' ':' [comp.membraneIndexs]
                    'Mass' 'cellDry'         ':' '-sum'
                    };
                states = SimulationEnsemble.load(simBatchDir, stateNames, 3600, 3600, 1, 'extract', selectedSimulations);
                
                dnaWt = permute(states.Mass.dnaWt, [4 3 2 1]) ./ permute(states.Mass.cellDry, [4 3 2 1]) * 100;
                rnaWt = permute(states.Mass.rnaWt, [4 3 2 1]) ./ permute(states.Mass.cellDry, [4 3 2 1]) * 100;
                proteinWt = permute(states.Mass.proteinWt, [4 3 2 1]) ./ permute(states.Mass.cellDry, [4 3 2 1]) * 100;
                lipidWt = permute(states.Mass.metaboliteWt, [4 3 2 1]) ./ permute(states.Mass.cellDry, [4 3 2 1]) * 100;
                
                modelDNA = mean(dnaWt);
                modelRNA = mean(rnaWt);
                modelProtein = mean(proteinWt);
                modelLipid = mean(lipidWt);
                
                stdDNA = std(dnaWt);
                stdRNA = std(rnaWt);
                stdProtein = std(proteinWt);
                stdLipid = std(lipidWt);
                
                save([outDirectory 'dryBiomassComposition.mat'], ...
                    'modelDNA', 'modelRNA', 'modelProtein', 'modelLipid', ...
                    'stdDNA', 'stdRNA', 'stdProtein', 'stdLipid');
            end
            
            % From Jonathan
            expDNA = 4;
            expLipid = 11;
            expProtein = 80;
            expRNA = 8;
            
            if nargin >= 4
                axesHandle = subplot('Position', position, 'Parent', figHandle);
            else
                axesHandle = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            end
            hold(axesHandle, 'on');
            
            h = [
                bar(axesHandle, (0:3)*2.5-0.5, [modelDNA modelLipid modelProtein modelRNA], 0.37, 'b')
                bar(axesHandle, (0:3)*2.5+0.5, [expDNA expLipid expProtein expRNA], 0.37, 'k')
                ];
            errorbar((0:3)*2.5-0.5, [modelDNA modelLipid modelProtein modelRNA], [stdDNA stdLipid stdProtein stdRNA], ...
                'r', 'Parent', axesHandle, 'LineStyle', 'none');
            set(h, 'EdgeColor', 'none')
            
            set(axesHandle, 'XTick', (0:3)*2.5, 'Xlim', [-1.3 3*2.5+1.3]);
            set(axesHandle, 'XTickLabel', {'DNA' 'Lipid' 'Protein' 'RNA'});
            set(axesHandle, 'YTick', 0:25:75);
            tick2text(axesHandle, 'axis', 'y', 'ytickoffset', 0.015);
            yTicks = getappdata(axesHandle, 'YTickText');
            set(yTicks, 'FontSize', 5, 'horizontalAlign', 'right');
            
            set(axesHandle, 'Box', 'off');
            set(axesHandle, 'TickDir', 'out');
            set(axesHandle, 'FontSize', 5);
            
            legendHandle = legend({'Model', 'Experiment'}, 'Location', 'NorthWest', 'FontSize', 7);
            legend(axesHandle, 'boxoff');
            legendPos = get(legendHandle, 'Position');
            set(legendHandle, 'Position', [0.99 1.025 legendPos(3:end)])
            set(legendHandle, 'Units', 'pixels');
            legendPos = get(legendHandle, 'Position');
            set(legendHandle, 'Position', [legendPos(1:3) 24])
            legendChildren = get(legendHandle, 'children');
            for i = 1:2:numel(legendChildren)
                ptch = get(legendChildren(i), 'children');
                set(ptch, 'xdata', [0.05 0.05 0.12 0.12 0.05])
                
                pos = get(legendChildren(i+1), 'Position');
                set(legendChildren(i+1), 'Position', [0.15 pos(2) 0])
            end
            
            ylabel(axesHandle, 'Percent dry mass', 'FontSize', 7);
            ylabelPos = get(get(axesHandle, 'ylabel'), 'position');
            set(get(axesHandle, 'ylabel'), 'position', [-1.95 ylabelPos(2:end)]);
            
            tick2text(axesHandle, 'axis', 'x', 'xtickoffset', 0.05);
            xTicks = getappdata(axesHandle, 'XTickText');
            set(xTicks, 'FontSize', 7);
            set(xTicks(1), 'String', 'DNA');
            set(xTicks(2), 'String', 'Lipid');
            set(xTicks(3), 'String', 'Protein');
            set(xTicks(4), 'String', 'RNA');
        end
        
        function cellCycleLength(simBatchDir, selectedSim, figHandle, position)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SummaryLogger;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            
            %% get data
            data = load([SimulationDiskUtil.getSimulationBatchDir([simBatchDir filesep num2str(selectedSim)]) filesep 'population-CellCyclePhases.mat']);
            stats = data.stats;
            maxTime = data.maxTime;
            
            edges = linspace(0, maxTime, 50);
            n1 = histc(stats(:, SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME), edges);
            n2 = histc(stats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME), edges);
            n3 = histc(stats(:, SummaryLogger.SIM_STATUS_INDEX_REP_TIME) - stats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME), edges);
            n4 = histc(stats(:, SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME) - stats(:, SummaryLogger.SIM_STATUS_INDEX_REP_TIME), edges);
            
            lowIdx = find(n1 > 0, 1, 'first') - 2;
            highIdx = find(n1 > 0, 1, 'last') + 2;
            n1([1:lowIdx highIdx:end]) = NaN;
            lowIdx = find(n2 > 0, 1, 'first') - 2;
            highIdx = find(n2 > 0, 1, 'last') + 2;
            n2([1:lowIdx highIdx:end]) = NaN;
            lowIdx = find(n3 > 0, 1, 'first') - 2;
            highIdx = find(n3 > 0, 1, 'last') + 2;
            n3([1:lowIdx highIdx:end]) = NaN;
            lowIdx = find(n4 > 0, 1, 'first') - 2;
            highIdx = find(n4 > 0, 1, 'last') + 2;
            n4([1:lowIdx highIdx:end]) = NaN;
            
            %% plot data
            if nargin >= 4
                axesHandle = subplot('position', position, 'parent', figHandle);
            else
                axesHandle = PlotUtil.newAxesHandle();
            end
            
            n1 = smooth(n1, 'lowess', 13) / size(stats, 1) * 100;
            n2 = smooth(n2, 'lowess', 13) / size(stats, 1) * 100;
            n3 = smooth(n3, 'lowess', 13) / size(stats, 1) * 100;
            n4 = smooth(n4, 'lowess', 13) / size(stats, 1) * 100;
            h = plot(axesHandle, edges / 3600, [n1 n2 n3 n4], 'LineWidth', 1);
            set(h(2), 'Color', 'g');
            
            ylim(axesHandle, [0 max(max([n1 n2 n3 n4]))]);
            set(axesHandle, 'XTick', 0:5:10);
            set(axesHandle, 'YTick', 0:10:40);
            set(axesHandle, 'Box', 'off');
            set(axesHandle, 'TickDir', 'out');
            set(axesHandle, 'FontSize', 5);
            
            xlabel(axesHandle, 'Duration (h)', 'FontSize', 7);
            ylabel(axesHandle, '% Cells', 'FontSize', 7);
            xlabelPos = get(get(axesHandle, 'xlabel'), 'position');
            ylabelPos = get(get(axesHandle, 'ylabel'), 'position');
            set(get(axesHandle, 'xlabel'), 'position', [xlabelPos(1) -3.7 xlabelPos(3:end)]);
            set(get(axesHandle, 'ylabel'), 'position', [-1.35 ylabelPos(2:end)]);
            
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', 0.04, 'ytickoffset', 0.015);
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            set(xTicks, 'FontSize', 5);
            set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'right');
            delete(xTicks(2));
            delete(yTicks(2:2:end));
            
            legendHandle = legend(h, {'Cell cycle', 'Replication initiation', 'Replication', 'Cytokinesis'}, 'Location', 'NorthEast', 'FontSize', 7);
            legend(axesHandle, 'boxoff');
            legendPos = get(legendHandle, 'Position');
            set(legendHandle, 'Position', [0.295 0.43 legendPos(3:end)])
            set(legendHandle, 'units', 'pixels');
            legendPos = get(legendHandle, 'Position');
            set(legendHandle, 'Position', [legendPos(1:3) 40]);
            legendChildren = get(legendHandle, 'children');
            set(legendChildren(1:3:end), 'visible', 'off');
            for i = 1:3:numel(legendChildren)
                set(legendChildren(i+1), 'xdata', [0.03 0.08])
                
                pos = get(legendChildren(i+2), 'Position');
                set(legendChildren(i+2), 'Position', [0.11 pos(2) 0])
            end
        end
        
        function replicationInitiationVsReplicationDuration(simBatchDir, selectedSim, figHandle, position)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SummaryLogger;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            
            %% get data
            data = load([SimulationDiskUtil.getSimulationBatchDir([simBatchDir filesep num2str(selectedSim)]) filesep 'population-CellCyclePhases.mat']);
            repInitDur = data.stats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME) / 3600;
            repDur = (data.stats(:, SummaryLogger.SIM_STATUS_INDEX_REP_TIME) - data.stats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME)) / 3600;
            
            %% plot data
            if nargin >= 4
                axesHandle = subplot('position', position, 'parent', figHandle);
            else
                axesHandle = PlotUtil.newAxesHandle();
            end
            
            plot(axesHandle, repInitDur, repDur, 'b.')
            
            set(axesHandle, 'XTick', 0:5:10);
            set(axesHandle, 'YTick', 0:4:8);
            set(axesHandle, 'Box', 'off');
            set(axesHandle, 'TickDir', 'out');
            set(axesHandle, 'FontSize', 5);
            
            xlabel(axesHandle, {'Replication initiation' 'duration (h)'}, 'FontSize', 7);
            ylabel(axesHandle, {'Replication' 'duration (h)'}, 'FontSize', 7);
            xlabelPos = get(get(axesHandle, 'xlabel'), 'position');
            ylabelPos = get(get(axesHandle, 'ylabel'), 'position');
            set(get(axesHandle, 'xlabel'), 'position', [xlabelPos(1)*1.25 -0.7 xlabelPos(3:end)]);
            set(get(axesHandle, 'ylabel'), 'position', [-0.4 4 ylabelPos(3:end)]);
            
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', 0.04, 'ytickoffset', 0.015);
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            set(xTicks, 'FontSize', 5);
            set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'right');
            delete(xTicks(2));
            delete(yTicks(2:2:end));
        end
        
        function replicationInitiationVsReplicationDurationFit(simBatchDir, selectedSim, figHandle, position)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SummaryLogger;
            import edu.stanford.covert.util.ComputationUtil;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            
            %% get data
            data = load([SimulationDiskUtil.getSimulationBatchDir([simBatchDir filesep num2str(selectedSim)]) filesep 'population-CellCyclePhases.mat']);
            repInitDur = data.stats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME) / 3600;
            repDur = (data.stats(:, SummaryLogger.SIM_STATUS_INDEX_REP_TIME) - data.stats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME)) / 3600;
            
            %% fit
            [a, b, cutX, gof] = ComputationUtil.bilinearFit(repInitDur, repDur);
            
            %% plot
            if nargin >= 4
                axesHandle = subplot('position', position, 'parent', figHandle);
            else
                axesHandle = PlotUtil.newAxesHandle();
            end
            hold(axesHandle, 'on');
            x1 = [min(repInitDur) cutX];
            x2 = [cutX max(repInitDur)];
            h = [
                plot(axesHandle, repInitDur, repDur, 'b.')
                plot(axesHandle, x1, a(1) + b(1) * x1, 'r')
                plot(axesHandle, x2, a(2) + b(2) * x2, 'r')
                ];
            legend(h(1:2), {'Data', 'Fit'}, 'Location', 'NorthEast');
            xlabel(axesHandle, 'Replication Initiation Duration (h)')
            ylabel(axesHandle, 'Replication Duration (h)')
            
            text(10, 6.5, sprintf('R^2_{adj} = %0.2f', gof.adjrsquare(1)), 'Parent', axesHandle);
            text(0.7, 4, sprintf('R^2_{adj} = %0.2f', gof.adjrsquare(2)), 'Parent', axesHandle);
            text(8, 3, sprintf('R^2_{adj} = %0.2f', gof.adjrsquare(3)), 'Parent', axesHandle);
        end
        
        function proteinBursting(simBatchDir, selectedSim, outDirectory, figHandle, position)
            import edu.stanford.covert.cell.sim.util.MacromoleculeUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            
            %% get data
            if exist([outDirectory 'proteinBursting.mat'], 'file')
                load([outDirectory 'proteinBursting.mat']);
            else
                [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep num2str(selectedSim)]);
                
                comp = sim.compartment;
                rna = sim.state('Rna');
                pm = sim.state('ProteinMonomer');
                
                [~, ~, ~, ~, ~, ...
                    ~, ~, ~, ~, ...
                    ~, ~, ~, ~, ...
                    matureRnaIdxs, ~, ~, ~, ...
                    monomerIdxs, ~, ~, ~, ...
                    complexIdxs, ~, ~, ~] = ...
                    MacromoleculeUtil.getMacromoleculeIndexsIDsNames('MG_218', sim);
                
                stateNames = {
                    'Rna'               'counts' rna.matureIndexs(matureRnaIdxs)    comp.cytosolIndexs
                    'ProteinMonomer'    'counts' pm.matureIndexs(monomerIdxs)       '-sum'
                    
                    };
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSim);
                
                
                enzymeRNA = permute(states.Rna.counts, [1 3 2]);
                enzymeMonomer = permute(states.ProteinMonomer.counts, [1 3 2]);
                time = (1:numel(enzymeRNA)) / 3600;
                
                save([outDirectory 'proteinBursting.mat'], ...
                    'time', 'enzymeRNA', 'enzymeMonomer');
            end
            
            
            %% plot
            options = struct;
            if nargin >= 4
                options.position = position;
            else
                [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
                clf(figHandle);
            end
            
            options.xdata = time;
            options.ydata = {
                enzymeMonomer
                enzymeRNA
                };
            options.ylabelStr = {
                {'Protein (cnt)'}
                {'mRNA (cnt)'}
                };
            options.offsetYAxes = -0.02;
            [axesHandles, xAxisHandle, ~, offsetAxesHandles] = PlotUtil.multiElementPlot(figHandle, 8 * ones(2, 1), [0 time(end)], options);
            set(xAxisHandle, 'XTick', 0:2:8, 'FontSize', 5);
            set(get(xAxisHandle, 'XLabel'), 'FontSize', 7);
            xlabelPos = get(get(xAxisHandle, 'xlabel'), 'position');
            set(get(xAxisHandle, 'xlabel'), 'position', [xlabelPos(1) -151 xlabelPos(3:end)]);
            
            tick2text(xAxisHandle, 'axis', 'x', 'xtickoffset', 50);
            xTicks = getappdata(xAxisHandle, 'XTickText');
            set(xTicks, 'FontSize', 5);
            delete(xTicks(2:4));
            
            for i = 1:numel(options.ydata)
                set(get(offsetAxesHandles(i), 'ylabel'), 'FontSize', 7)
                set(offsetAxesHandles(i), 'FontSize', 5)
                set(axesHandles(i), 'FontSize', 5)
                
                yTicks = getappdata(offsetAxesHandles(i), 'YTickText');
                if isempty(yTicks)
                    tick2text(offsetAxesHandles(i), 'axis', 'y', 'ytickoffset', 0.14)
                    yTicks = getappdata(offsetAxesHandles(i), 'YTickText');
                end
                
                set(yTicks, 'FontSize', 5);
                for j = 1:numel(yTicks)
                    pos = get(yTicks(j), 'position');
                    set(yTicks(j), 'position', [-0.18 pos(2:end)], 'horizontalalign', 'right');
                end
                
                ylabelPos = get(get(offsetAxesHandles(i), 'ylabel'), 'position');
                set(get(offsetAxesHandles(i), 'ylabel'), 'position', [-0.6 ylabelPos(2:end)], 'FontSize', 7, 'VerticalAlign', 'bottom');
            end
            
            burstIdxs = find(diff(enzymeRNA) > 0);
            axesPos1 = get(axesHandles(1), 'position');
            axesPos2 = get(axesHandles(2), 'position');
            monYLims = ylim(axesHandles(1));
            rnaYLims = ylim(axesHandles(2));
            for i = 1:numel(burstIdxs)
                y = (axesPos2(2):0.0025:(axesPos1(2) + axesPos1(4) * (enzymeMonomer(burstIdxs(i)) - monYLims(1)) / range(monYLims)));
                x = (axesPos2(1) + axesPos2(3) * time(burstIdxs(i))/time(end)) * [1 1];
                for j = 1:3:3 * floor(numel(y) / 3)
                    %if any(...
                    %        (y(j:j+1)-axesPos2(2))/axesPos2(4)*range(rnaYLims) + rnaYLims(1) >= enzymeRNA(burstIdxs(i)) & ...
                    %        (y(j:j+1)-axesPos2(2))/axesPos2(4)*range(rnaYLims) + rnaYLims(1) <= enzymeRNA(burstIdxs(i)+1))
                    %    continue;
                    %end
                    annotation(figHandle, 'line', x, y(j:j+1), 'Color', [1 0 0]);
                end
            end
        end
        
        function geneToPhenotype(simBatchDir, selectedSim, outDirectory, figHandle, position)
            import edu.stanford.covert.cell.sim.analysis.SingleCell;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            
            %% get data
            if exist([outDirectory 'geneToPhenotype.mat'], 'file')
                load([outDirectory 'geneToPhenotype.mat']);
            else
                [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep num2str(selectedSim)]);
                
                [time, ~, ~, ~, ~, ~, ~, growth, reactionFlux, ...
                    enzymeComplex, enzymeMonomer, enzymeRNA, enzymeGene] = ...
                    SingleCell.calcGrowthData(sim, simBatchDir, selectedSim, 'GalE');
                enzymeGene = permute(enzymeGene, [1 3 2]);
                enzymeRNA = permute(enzymeRNA, [1 3 2]);
                enzymeMonomer = permute(enzymeMonomer, [1 3 2]);
                enzymeComplex = permute(enzymeComplex, [1 3 2]);
                reactionFlux = permute(reactionFlux, [1 3 2]);
                growth = permute(growth, [1 3 2]);
                
                save([outDirectory 'geneToPhenotype.mat'], ...
                    'time', 'enzymeGene', 'enzymeRNA', 'enzymeMonomer', 'enzymeComplex', 'reactionFlux', 'growth');
            end
            reactionFlux = -reactionFlux * 1000;
            
            %% plot
            options = struct;
            if nargin >= 4
                options.position = position;
            else
                [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
                clf(figHandle);
            end
            
            options.xdata = time;
            options.ydata = {
                growth
                reactionFlux
                enzymeComplex
                enzymeMonomer
                enzymeRNA
                enzymeGene
                };
            options.ylabelStr = {
                {'Growth' '(fg/h)'}
                {'Flux' '(10^3 rxn/h)'}
                {'Cpx' '(cnt)'}
                {'Mon' '(cnt)'}
                {'mRNA' '(cnt)'}
                {'Gene' '(cnt)'}
                };
            options.offsetYAxes = false;
            [axesHandles, xAxisHandle, ~, offsetAxesHandles] = PlotUtil.multiElementPlot(figHandle, 2 * ones(6, 1), time([1 end])', options);
            set(xAxisHandle, 'XTick', 0:2:8, 'FontSize', 5);
            set(get(xAxisHandle, 'XLabel'), 'FontSize', 7);
            xlabelPos = get(get(xAxisHandle, 'xlabel'), 'position');
            set(get(xAxisHandle, 'xlabel'), 'position', [xlabelPos(1) -80 xlabelPos(3:end)]);
            
            tick2text(xAxisHandle, 'axis', 'x', 'xtickoffset', 60);
            xTicks = getappdata(xAxisHandle, 'XTickText');
            set(xTicks, 'FontSize', 5);
            delete(xTicks(2:4));
            
            for i = 1:numel(axesHandles)
                set(get(axesHandles(i), 'ylabel'), 'FontSize', 7)
                set(axesHandles(i), 'FontSize', 5)
                
                yTicks = getappdata(axesHandles(i), 'YTickText');
                if isempty(yTicks)
                    tick2text(axesHandles(i), 'axis', 'y', 'ytickoffset', 0.06)
                    yTicks = getappdata(axesHandles(i), 'YTickText');
                end
                
                set(yTicks, 'FontSize', 5);
                for j = 1:numel(yTicks)
                    pos = get(yTicks(j), 'position');
                    set(yTicks(j), 'position', [-0.15 pos(2:end)], 'horizontalalign', 'right');
                end
                
                ylabelPos = get(get(axesHandles(i), 'ylabel'), 'position');
                set(get(axesHandles(i), 'ylabel'), 'position', [-0.8 ylabelPos(2:end)], 'FontSize', 5, 'VerticalAlign', 'bottom');
            end
            
            ylabelPos = get(get(axesHandles(1), 'ylabel'), 'position');
            set(get(axesHandles(1), 'ylabel'), 'position', [ylabelPos(1) 1.65 ylabelPos(3)]);
        end
        
        function spaceTimePlot(simBatchDir, selectedSim, outDirectory, figHandle, position)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.analysis.ChromosomeSpaceTimePlot;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if ~exist('selectedSim', 'var') || isempty(selectedSim)
                selectedSim = 1;
            end
            
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep num2str(selectedSim)]);
            c = sim.state('Chromosome');
            
            %% get data
            if exist([outDirectory 'spaceTimePlot.mat'], 'file')
                load([outDirectory 'spaceTimePlot.mat'])
            else
                [axesHandle, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
                [ensemble, states, dnaABoxBound] = ChromosomeSpaceTimePlot.plotSpaceTime(...
                    axesHandle, simBatchDir, selectedSim, [], false, [], false, true); %#ok<*NASGU,*ASGLU>
                save([outDirectory 'spaceTimePlot.mat'], 'ensemble', 'states', 'dnaABoxBound')
                close(figHandle);
            end
            
            %% space-time density plots
            %create axes
            if nargin < 4
                [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
                position = [0.1 0.1 0.8 0.8];
            end
            
            x = position(1);
            y = position(2);
            w = position(3);
            h = position(4);
            spcg = 0.08;
            figW = 1.45 * (0.5 - x);
            
            axesHandles = [
                subplot('position', [x            y  figW  h], 'parent', figHandle)
                subplot('position', [x+figW+spcg  y  w-figW-spcg  h], 'parent', figHandle)
                ];
            
            zoomTime = [1.34 1.41];
            zoomPos = [267000 273000];
            time = ensemble.stateData.time / 3600;
            
            zoomBoxHandle = [
                annotation(figHandle, 'line', x + [figW*zoomTime(2)/time(end) figW+spcg], y + [h*zoomPos(end)/c.sequenceLen h])
                annotation(figHandle, 'line', x + [figW*zoomTime(2)/time(end) figW+spcg], y + [h*zoomPos(1)/c.sequenceLen 0])
                annotation(figHandle, 'rectangle', [
                x + figW*zoomTime(1)/time(end)
                y + h*zoomPos(1)/c.sequenceLen
                figW*diff(zoomTime)/time(end)
                h*diff(zoomPos)/c.sequenceLen
                ]')
                ];
            set(zoomBoxHandle, 'LineWidth', 1);
            
            %% plot 1
            axesHandle = axesHandles(1);
            
            %plot
            
            [~, ~, ~, legendHandle] = ChromosomeSpaceTimePlot.plotSpaceTime(...
                axesHandle, simBatchDir, selectedSim, [], false, [], ...
                false, true, ensemble, states, dnaABoxBound, 2, 0.25, false);
            
            %style
            xlim(axesHandle, [0 time(end)]);
            ylim(axesHandle, [1 c.sequenceLen]);
            set(axesHandle, 'TickDir', 'out', 'XTick', 0:2:8, 'FontSize', 5);
            
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', 0.035, 'ytickoffset', 0.008);
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            set(xTicks, 'FontSize', 5);
            set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'right');
            delete(xTicks(2:4));
            set(yTicks([1 3]), 'String', 'terC');
            set(yTicks(2), 'String', 'oriC');
            
            xlabel(axesHandle, 'Time (h)', 'FontSize', 7);
            ylabel(axesHandle, 'Position', 'FontSize', 7);
            xLabelPos = get(get(axesHandle, 'xlabel'), 'position');
            yLabelPos = get(get(axesHandle, 'ylabel'), 'position');
            set(get(axesHandle, 'xlabel'), 'position', [xLabelPos(1) -25000 xLabelPos(3)])
            set(get(axesHandle, 'ylabel'), 'position', [-0.5 yLabelPos(2:3)])
            
            %legend
            legendChildren = get(legendHandle, 'children');
            set(legendHandle, 'FontSize', 7);
            legendPos = get(legendHandle, 'Position');
            set(legendHandle, 'Position', [0.037 0.346 legendPos(3:4)]);
            set(legendHandle, 'units', 'pixels');
            legendPos = get(legendHandle, 'Position');
            set(legendHandle, 'Position', [legendPos(1:2) 180 30]);
            tmp = findobj(legendChildren, 'type', 'line');
            set(tmp(1:2:end), 'visible', 'off');
            set(tmp(2:2:end), 'visible', 'on', 'LineWidth', 4)
            set(tmp(2:2:end), 'XData', [0.07 0.18])
            tmp = findobj(legendChildren, 'type', 'text');
            for i = 1:numel(tmp)
                pos = get(tmp(i), 'Position');
                set(tmp(i), 'Position', [0.24 pos(2) 0])
            end
            tmp = findobj(legendChildren, 'type', 'patch');
            line([0.07 0.18], [min(get(tmp, 'ydata')) max(get(tmp, 'ydata'))], 'LineWidth', 4, 'Parent', legendHandle, 'Color', 'r');
            set(tmp, 'visible', 'off');
            
            set(legendHandle, 'units', 'pixels');
            legendPos = get(legendHandle, 'position');
            set(legendHandle, 'position', [[1.4 1.53] .* legendPos(1:2) 0.47*legendPos(3) legendPos(4)]);
            
            %% plot 2
            axesHandle = axesHandles(2);
            
            %plot
            ChromosomeSpaceTimePlot.plotSpaceTime(axesHandle, simBatchDir, selectedSim, [], ...
                false, [], false, false, ensemble, states, dnaABoxBound, 0.5, 0.5, true); %#ok<*NODEF>
            
            xlim(axesHandle, zoomTime);
            ylim(axesHandle, zoomPos);
            
            %style
            set(axesHandle, 'TickDir', 'out', ...
                'XTick', linspace(zoomTime(1), zoomTime(2), 8), ...
                'YTick', linspace(zoomPos(1), zoomPos(2), 7)-38, ...
                'FontSize', 5);
            
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', 0.035, 'ytickoffset', 0.015);
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            set(xTicks, 'FontSize', 5);
            set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'right');
            set(yTicks(2), 'String', '558000')
            set(yTicks(4), 'String', '560000')
            set(yTicks(6), 'String', '562000')
            delete(xTicks(setdiff(1:end, [2 7])));
            delete(yTicks(1:2:end));
            
            xlabel(axesHandle, 'Time (h)', 'FontSize', 7);
            ylabel(axesHandle, 'Position (nt)', 'FontSize', 7);
            xLabelPos = get(get(axesHandle, 'xlabel'), 'position');
            yLabelPos = get(get(axesHandle, 'ylabel'), 'position');
            set(get(axesHandle, 'xlabel'), 'position', [xLabelPos(1) 266750 xLabelPos(3)])
            set(get(axesHandle, 'yLabel'), 'position', [1.321 yLabelPos(2:end)])
            
            figPos = get(figHandle, 'PaperSize');
            figW = figPos(1);
            figH = figPos(2);
            axisPos = get(axesHandle, 'position');
            x = axisPos(1);
            y = axisPos(2);
            w = axisPos(3);
            h = axisPos(4);
            annotation(figHandle, 'ellipse', [x+0.595*w  y+0.645*h  0.01 0.01*figW/figH], 'Color', 'k')
            annotation(figHandle, 'ellipse', [x+0.775*w  y+0.245*h  0.01 0.01*figW/figH], 'Color', 'k')
            annotation(figHandle, 'ellipse', [x+0.48*w  y+0.875*h  0.01 0.01*figW/figH], 'Color', 'k')
            annotation(figHandle, 'ellipse', [x+0.645*w  y+0.84*h  0.01 0.01*figW/figH], 'Color', 'k')
            annotation(figHandle, 'ellipse', [x+0.655*w  y+0.86*h  0.01 0.01*figW/figH], 'Color', 'k')
            annotation(figHandle, 'arrow', [x+0.535*w  x+0.60*w], [y+0.57*h  y+0.63*h], 'Color', 'k', ...
                'HeadLength', 3, 'HeadWidth', 5, 'HeadStyle', 'vback1')
            text(1.364, 270100, {'DNA pol-RNA' 'pol collision'}, 'parent', axesHandle, ...
                'FontSize', 6, 'VerticalAlign', 'middle', 'HorizontalAlign', 'center', 'Color', 'k');
            text(1.3860, 267400, {'Leading' 'DNA pol'}, 'parent', axesHandle, ...
                'FontSize', 6, 'VerticalAlign', 'middle', 'HorizontalAlign', 'center', 'Color', 'g');
            text(1.4020, 272500, {'Lagging' 'DNA pol'}, 'parent', axesHandle, ...
                'FontSize', 6, 'VerticalAlign', 'middle', 'HorizontalAlign', 'center', 'Color', 'g');
            
            %% zoom box
            uistack(zoomBoxHandle, 'top');
        end
        
        function metaboliteConcentrations(simBatchDir, selectedSim, outDirectory, figHandle, position)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.PrintUtil;
            import edu.stanford.covert.cell.sim.util.SummaryLogger;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            
            %% get constants
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep num2str(selectedSim)]);
            m = sim.state('Metabolite');
            metIdxs = [m.ntpIndexs; m.ndpIndexs; m.nmpIndexs; m.dntpIndexs; m.aminoAcidIndexs; m.phosphateIndexs; m.diphosphateIndexs; m.hydrogenIndexs];
            
            metIDs = m.wholeCellModelIDs;
            metIDs(m.dntpIndexs) = {'dATP', 'dCTP', 'dGTP', 'dTTP'};
            metIDs(m.aminoAcidIndexs) = {
                'Ala'
                'Arg'
                'Asn'
                'Asp'
                'Cys'
                'Gln'
                'Glu'
                'Gly'
                'His'
                'Ile'
                'Leu'
                'Lys'
                'Met'
                'Phe'
                'Pro'
                'Ser'
                'Thr'
                'Trp'
                'Tyr'
                'Val'
                'fMet'
                };
            metIDs{m.phosphateIndexs} = 'Pi';
            metIDs{m.diphosphateIndexs} = 'PPi';
            metIDs{m.hydrogenIndexs} = 'H^+';
            
            %% get data
            tmp = load([SimulationDiskUtil.getSimulationBatchDir([simBatchDir filesep '1']) filesep 'population-MetaboliteConcentrations.mat']);
            tmp2 = [tmp.meanConcs tmp.stdConcs];
            modelConcs = zeros(numel(m.wholeCellModelIDs), 2);
            modelConcs(m.aminoAcidIndexs, :) = tmp2(17:37, :);
            modelConcs(m.ntpIndexs, :) = tmp2(1:4, :);
            modelConcs(m.ndpIndexs, :) = tmp2(5:8, :);
            modelConcs(m.nmpIndexs, :) = tmp2(9:12, :);
            modelConcs(m.dntpIndexs, :) = tmp2(13:16, :);
            modelConcs(m.phosphateIndexs, :) = tmp2(38, :);
            modelConcs(m.diphosphateIndexs, :) = tmp2(39, :);
            modelConcs(m.hydrogenIndexs, :) = tmp2(40, :);
            
            [~, ~, data] = xlsread('documentation/reconstruction/Metabolite Concentrations.xls');
            tfs = ~cellfun(@(x) any(isnan(x)), data(:, 1));
            tfs(1:2) = false;
            data = data(tfs, :);
            
            [~, idxs] = ismember(data(:, 1), m.wholeCellModelIDs);
            expConcs = NaN(size(modelConcs, 1), 3);
            cyberCellConc = NaN(size(modelConcs, 1), 1);
            expConcs(idxs, :) = cell2mat(data(:, 2:4));
            cyberCellConc(idxs, :) = cell2mat(data(:, 5));
            
            sortedIdxs = [
                m.aminoAcidIndexs(1:end-1)
                m.ntpIndexs
                m.ndpIndexs
                m.nmpIndexs
                m.dntpIndexs
                m.phosphateIndexs
                m.diphosphateIndexs
                m.hydrogenIndexs
                ];
            allSortedIdxs = [
                m.dntpIndexs
                m.ndpIndexs
                m.nmpIndexs
                m.ntpIndexs
                m.aminoAcidIndexs
                m.phosphateIndexs
                m.diphosphateIndexs
                m.hydrogenIndexs
                ];
            
            %% plot
            if nargin >= 4
                axesHandle = subplot('position', position, 'parent', figHandle);
            else
                axesHandle = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            end
            
            hold(axesHandle, 'on');
            set(axesHandle, ...
                'YScale', 'log', ...
                'YMinorTick', 'on', ...
                'XTick', 1:numel(sortedIdxs), ...
                'XTickLabel', m.wholeCellModelIDs(sortedIdxs), ...
                'YTick', logspace(-2, 2, 5), ...
                'FontSize', 5, ...
                'TickDir', 'out');
            xlim(axesHandle, [0.5 numel(sortedIdxs)+0.5]);
            ylim(axesHandle, [1e-2/1.5 1.5e2]);
            
            plot(axesHandle, 1:numel(sortedIdxs), expConcs(sortedIdxs, 1) ./ modelConcs(sortedIdxs, 1), ...
                'rd', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
            plot(axesHandle, 1:numel(sortedIdxs), cyberCellConc(sortedIdxs) ./ modelConcs(sortedIdxs, 1), ...
                'gs', 'MarkerSize', 3, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
            nModel = numel(sortedIdxs);
            line([1:nModel; 1:nModel], max([ 1 + (modelConcs(sortedIdxs, 2) ./ modelConcs(sortedIdxs, 1))'; 1 - (modelConcs(sortedIdxs, 2) ./ modelConcs(sortedIdxs, 1))'], eps), 'Color', 'b');
            
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', 0.13, 'ytickoffset', 0.013);
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            set(xTicks, 'FontSize', 5, 'HorizontalAlign', 'left', 'VerticalAlign', 'middle', 'Interpreter', 'tex', 'rotation', 270);
            set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'right', 'Interpreter', 'tex');
            for i = 1:numel(xTicks)
                set(xTicks(i), 'String', metIDs{sortedIdxs(i)}, 'position', [i 5e-3 -1], 'Margin', 1e-6);
            end
            for i = 1:numel(yTicks)
                set(yTicks(i), 'String', sprintf('10^{%d}', log10(str2double(get(yTicks(i), 'string')))));
            end
            
            ylabel(axesHandle, 'Expt/model concentration', 'FontSize', 7);
            yLabelPos = get(get(axesHandle, 'ylabel'), 'position');
            set(get(axesHandle, 'ylabel'), 'position', [-1.6 0.8 yLabelPos(3)])
            
            legendHandle = legend({'Bennett {\it{et al.}}, 2009', 'Literature (CCDB)', 'Model s.d.'}, 'Location', 'SouthEast', 'FontSize', 6);
            legend(axesHandle, 'boxoff');
            legendPos = get(legendHandle, 'Position');
            set(legendHandle, 'Position', [0.088 0.092 legendPos(3:end)])
            set(legendHandle, 'units', 'pixels');
            legendPos = get(legendHandle, 'Position');
            set(legendHandle, 'Position', [legendPos(1:2) 82 28]);
            legendChildren = get(legendHandle, 'children');
            for i = 1:3:numel(legendChildren)
                if strcmp(get(legendChildren(i), 'marker'), 'none')
                    set(legendChildren(i), 'visible', 'off');
                    set(legendChildren(i+1), 'xdata', [0.03 0.10])
                else
                    set(legendChildren(i+1), 'visible', 'off');
                    set(legendChildren(i), 'xdata', 0.06)
                end
                
                pos = get(legendChildren(i+2), 'Position');
                set(legendChildren(i+2), 'Position', [0.13 pos(2) 0])
            end
            
            axisPos = get(axesHandle, 'Position');
            x = axisPos(1);
            w = axisPos(3);
            y = 0.021;
            annotation(figHandle, 'line', x + ([ 1 20]-0.5 + [-0.2 0.2])/39*w, y * [1 1])
            annotation(figHandle, 'line', x + ([21 24]-0.5 + [-0.2 0.2])/39*w, y * [1 1])
            annotation(figHandle, 'line', x + ([25 28]-0.5 + [-0.2 0.2])/39*w, y * [1 1])
            annotation(figHandle, 'line', x + ([29 32]-0.5 + [-0.2 0.2])/39*w, y * [1 1])
            annotation(figHandle, 'line', x + ([33 36]-0.5 + [-0.2 0.2])/39*w, y * [1 1])
            annotation(figHandle, 'line', x + ([37 39]-0.5 + [-0.2 0.2])/39*w, y * [1 1])
            
            y = 0.5e-3;
            text(mean([ 1 20]), y, 'Amino acid', 'FontSize', 7, 'HorizontalAlign', 'center', 'Parent', axesHandle)
            text(mean([21 24]), y, 'NTP',  'FontSize', 7, 'HorizontalAlign', 'center', 'Parent', axesHandle)
            text(mean([25 28]), y, 'NDP',  'FontSize', 7, 'HorizontalAlign', 'center', 'Parent', axesHandle)
            text(mean([29 32]), y, 'NMP',  'FontSize', 7, 'HorizontalAlign', 'center', 'Parent', axesHandle)
            text(mean([33 36]), y, 'dNTP', 'FontSize', 7, 'HorizontalAlign', 'center', 'Parent', axesHandle)
            text(mean([37 39]), y, 'Ion',  'FontSize', 7, 'HorizontalAlign', 'center', 'Parent', axesHandle)
            
            ptch = patch([0.5 numel(sortedIdxs)+0.5 numel(sortedIdxs)+0.5 0.5], ...
                [1e1 1e1 1e-1 1e-1], 1, ...
                'Parent', axesHandle, 'EdgeColor', 'none');
            set(ptch, 'FaceColor', 0.925 * [1 1 1], 'FaceAlpha', 1, 'EdgeAlpha', 1);
            h = line([0.5 0.5], [1/2*1e-1 2*1e1], 'LineWidth', 0.5, 'Color', 'k', 'Parent', axesHandle);
            uistack(ptch, 'bottom');
            uistack(h, 'top');
            
            %% table
            colLabels = {'Metabolite ID' 'Name', 'Whole Cell Model Mean' 'Whole Cell Model Std' 'Bennet et al. 2009 Mean' 'Bennet et al. 2009 Lower Bound' 'Bennet et al. 2009 Upper Bound' 'CyberCell'};
            content = [m.wholeCellModelIDs  m.names  ...
                num2cell([modelConcs expConcs cyberCellConc])];
            content = content(allSortedIdxs, :);
            
            PrintUtil.printToFile(content, colLabels, [outDirectory 'metaboliteConcentrations.xls'], 'Metabolite Concentrations');
        end
        
        function [exptGeneExpr, modelGeneExp] = geneExpression(simBatchDir, selectedSim, figHandle, position, logScale)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if ~exist('selectedSim', 'var') || isempty(selectedSim)
                selectedSim = 1;
            end
            
            %% get constants
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep num2str(selectedSim)]);
            rna = sim.state('Rna');
            g = sim.gene;
            
            %% get data
            exptGeneExpr = rna.expectedGeneExpression(:, 1);
            exptGeneExpr(exptGeneExpr == 0) = NaN;
            exptGeneExpr(g.mRNAIndexs) = exptGeneExpr(g.mRNAIndexs) / nansum(exptGeneExpr(g.mRNAIndexs));
            exptGeneExpr(setdiff(1:end, g.mRNAIndexs)) = exptGeneExpr(setdiff(1:end, g.mRNAIndexs)) / nansum(exptGeneExpr(setdiff(1:end, g.mRNAIndexs)));
            
            modelGeneExp = rna.matureRNAGeneComposition * rna.expression(rna.matureIndexs);
            modelGeneExp(g.mRNAIndexs) = modelGeneExp(g.mRNAIndexs) / sum(modelGeneExp(g.mRNAIndexs));
            modelGeneExp(setdiff(1:end, g.mRNAIndexs)) = modelGeneExp(setdiff(1:end, g.mRNAIndexs)) / sum(modelGeneExp(setdiff(1:end, g.mRNAIndexs)));
            
            %% plot
            if nargin >= 4
                axesHandle = subplot('Position', position, 'Parent', figHandle);
            else
                axesHandle = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            end
            
            if logScale
                set(axesHandle, 'XScale', 'log');
                set(axesHandle, 'YScale', 'log');
            end
            set(axesHandle, 'TickDir', 'out');
            set(axesHandle, 'FontSize', 5);
            hold(axesHandle, 'on');
            
            plot(axesHandle, exptGeneExpr(g.mRNAIndexs), modelGeneExp(g.mRNAIndexs, 1), '.', 'Color', [0 1 1])
            plot(axesHandle, exptGeneExpr(g.rRNAIndexs), modelGeneExp(g.rRNAIndexs, 1), '.', 'Color', [0 0 1])
            plot(axesHandle, exptGeneExpr(g.sRNAIndexs), modelGeneExp(g.sRNAIndexs, 1), '.', 'Color', [1 0 0])
            plot(axesHandle, exptGeneExpr(g.tRNAIndexs), modelGeneExp(g.tRNAIndexs, 1), '.', 'Color', [0 1 0])
            
            if logScale
                xlim(axesHandle, [6.5e-5 8e-2]);
            else
                xlim(axesHandle, [0 7.2e-2])
                set(axesHandle, 'XTick', 0:0.02:0.06, 'YTick', 0:0.02:0.06);
            end
            ylim(axesHandle, xlim(axesHandle));
            
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', 0.02, 'ytickoffset', 0.03);
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            set(xTicks, 'FontSize', 5);
            set(yTicks, 'FontSize', 5);
            if logScale
                for i = 1:numel(xTicks)
                    xTickPos = get(xTicks(i), 'position');
                    set(xTicks(i), 'position', [xTickPos(1) 5.5e-5 xTickPos(3)], ...
                        'String', sprintf('10^{%d}', log10(str2double(get(xTicks(i), 'string')))), ...
                        'Interpreter', 'tex');
                    
                    yTickPos = get(yTicks(i), 'position');
                    set(yTicks(i), 'position', [5e-5 yTickPos(2:3)], ...
                        'String', sprintf('10^{%d}', log10(str2double(get(yTicks(i), 'string')))), ...
                        'Interpreter', 'tex');
                end
            end
            
            xlabel(axesHandle, {'Observed gene expression (norm)'}, 'FontSize', 7);
            ylabel(axesHandle, {'Predicted gene expression (norm)'}, 'FontSize', 7);
            xLabelPos = get(get(axesHandle, 'xlabel'), 'position');
            yLabelPos = get(get(axesHandle, 'ylabel'), 'position');
            if logScale
                set(get(axesHandle, 'xlabel'), 'position', [xLabelPos(1) 4.75e-5 xLabelPos(3)])
                set(get(axesHandle, 'ylabel'), 'position', [4.1e-5 yLabelPos(2:3)])
            else
                set(get(axesHandle, 'xlabel'), 'position', [xLabelPos(1) -0.0032 xLabelPos(3)])
                set(get(axesHandle, 'ylabel'), 'position', [-0.0047 yLabelPos(2:3)])
            end
            
            if logScale
                text(3e-3, 3e-4, 'mRNA', 'Color', [0 1 1], 'Parent', axesHandle, 'FontSize', 7)
                text(7.8e-2, 2.5e-2, 'rRNA', 'Color', [0 0 1], 'Parent', axesHandle, 'FontSize', 7, 'HorizontalAlign', 'right')
                text(1.35e-2, 0.77e-2, 'sRNA', 'Color', [1 0 0], 'Parent', axesHandle, 'FontSize', 7)
                text(1.7e-2, 6e-2, 'tRNA', 'Color', [0 1 0], 'Parent', axesHandle, 'FontSize', 7, 'HorizontalAlign', 'right')
            else
                text(0.012, 0.003, 'mRNA', 'Color', [0 1 1], 'Parent', axesHandle, 'FontSize', 7)
                text(7e-2, 3.1e-2, 'rRNA', 'Color', [0 0 1], 'Parent', axesHandle, 'FontSize', 7, 'HorizontalAlign', 'right')
                text(1.35e-2, 0.77e-2, 'sRNA', 'Color', [1 0 0], 'Parent', axesHandle, 'FontSize', 7)
                text(4.5e-2, 2.7e-2, 'tRNA', 'Color', [0 1 0], 'Parent', axesHandle, 'FontSize', 7, 'HorizontalAlign', 'right')
            end
        end
        
        function geneExpressionRatios(exptGeneExpr, modelGeneExpr, axesHandle)
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            %% get constants
            sim = CachedSimulationObjectUtil.load();
            g = sim.gene;
            
            x1 = modelGeneExpr;
            x2 = exptGeneExpr;
            y1 = modelGeneExpr(g.mRNAIndexs);
            y2 = exptGeneExpr(g.mRNAIndexs);
            z1 = modelGeneExpr(setdiff(1:end, g.mRNAIndexs));
            z2 = exptGeneExpr(setdiff(1:end, g.mRNAIndexs));
            
            corr(x1(~isnan(x1) & ~isnan(x2)), x2(~isnan(x1) & ~isnan(x2)))
            corr(y1(~isnan(y1) & ~isnan(y2)), y2(~isnan(y1) & ~isnan(y2)))
            corr(z1(~isnan(z1) & ~isnan(z2)), z2(~isnan(z1) & ~isnan(z2)))
            
            %% plot
            ratio = modelGeneExpr ./ exptGeneExpr;
            
            spcg = 50;
            xPos = zeros(size(ratio));
            xPos(g.mRNAIndexs) = 1:numel(g.mRNAIndexs);
            xPos(g.rRNAIndexs) = (1:numel(g.rRNAIndexs)) + spcg + numel(g.mRNAIndexs);
            xPos(g.sRNAIndexs) = (1:numel(g.sRNAIndexs)) + 2*spcg + numel(g.mRNAIndexs) + numel(g.rRNAIndexs);
            xPos(g.tRNAIndexs) = (1:numel(g.tRNAIndexs)) + 3*spcg + numel(g.mRNAIndexs) + numel(g.rRNAIndexs) + numel(g.sRNAIndexs);
            
            hold(axesHandle, 'on');
            plot(axesHandle, xPos(g.mRNAIndexs), ratio(g.mRNAIndexs), 'r.');
            plot(axesHandle, xPos(g.rRNAIndexs), ratio(g.rRNAIndexs), 'g.');
            plot(axesHandle, xPos(g.sRNAIndexs), ratio(g.sRNAIndexs), 'b.');
            plot(axesHandle, xPos(g.tRNAIndexs), ratio(g.tRNAIndexs), 'c.');
            
            xlim(axesHandle, [-spcg numel(ratio)+4*spcg]);
            ylim(axesHandle, [1/2*1e-1 2*1e1]);
            ylabel(axesHandle, 'Model / Expt Gene Expression');
            xlabel(axesHandle, 'Gene');
            set(axesHandle, 'XTick', [], 'YScale', 'log', 'TickDir', 'out');
            box(axesHandle, 'off');
            
            lineY = 0.88 * max(ylim(axesHandle));
            textY = 0.99 * max(ylim(axesHandle));
            
            line(xPos(g.mRNAIndexs([1 end])) + [-10; 10], lineY * [1 1], 'Parent', axesHandle, 'Color', 'r')
            line(xPos(g.rRNAIndexs([1 end])) + [-10; 10], lineY * [1 1], 'Parent', axesHandle, 'Color', 'g')
            line(xPos(g.sRNAIndexs([1 end])) + [-10; 10], lineY * [1 1], 'Parent', axesHandle, 'Color', 'b')
            line(xPos(g.tRNAIndexs([1 end])) + [-10; 10], lineY * [1 1], 'Parent', axesHandle, 'Color', 'c')
            text(mean(xPos(g.mRNAIndexs([1 end]))), textY, 'mRNA', 'HorizontalAlign', 'center', 'Color', 'r');
            text(mean(xPos(g.rRNAIndexs([1 end]))), textY, 'rRNA', 'HorizontalAlign', 'center', 'Color', 'g');
            text(mean(xPos(g.sRNAIndexs([1 end]))), textY, 'sRNA', 'HorizontalAlign', 'center', 'Color', 'b');
            text(mean(xPos(g.tRNAIndexs([1 end]))), textY, 'tRNA', 'HorizontalAlign', 'center', 'Color', 'c');
        end
        
        function circosPlot(simBatchDir, selectedSim, outDirectory, figHandle, position)
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.analysis.ChromosomePositionHistogram;
            import edu.stanford.covert.cell.sim.analysis.ChromosomeSpaceTimePlot;
            import edu.stanford.covert.cell.sim.analysis.Figures34;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if ~exist('selectedSim', 'var') || isempty(selectedSim)
                selectedSim = 1;
            end
            
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep num2str(selectedSim)]);
            ch = sim.state('Chromosome');
            pc = sim.state('ProteinComplex');
            
            %% get data
            if exist([outDirectory 'circosPlot.mat'], 'file')
                load([outDirectory 'circosPlot.mat']);
            else
                stateNames = {
                    'Time'       'values'                   ':' ':'
                    'Chromosome' 'complexBoundSites'        ':' ':'
                    };
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSim);
                
                chromProteins = {...
                    %whole cell ids                strands  number mers
                    {'MG_094_HEXAMER'},            1:4      1              %DNA Pol
                    {'MG_469_1MER_ADP' ...                                 %DnaA
                    'MG_469_1MER_ATP' ...
                    'MG_469_2MER_ATP' ...
                    'MG_469_3MER_ATP' ...
                    'MG_469_4MER_ATP' ...
                    'MG_469_5MER_ATP' ...
                    'MG_469_6MER_ATP' ...
                    'MG_469_7MER_ATP'},           1:2      [1 1:7]
                    {'RNA_POLYMERASE'                                      %RNA Pol
                    'RNA_POLYMERASE_HOLOENZYME'}, 1:2      [1 1]
                    };
                
                densityBEMatrix = zeros(1000, size(chromProteins, 1));
                for i = 1:size(chromProteins, 1)
                    for j = 1:numel(chromProteins{i, 1})
                        protIdx = pc.getIndexs(chromProteins{i, 1}{j});
                        
                        %calculate frequency binding each position
                        tmp = full(sum(sum(states.Chromosome.complexBoundSites(:, chromProteins{i, 2}, :) == protIdx, 3), 2));
                        
                        %circularly convolve with footprint
                        tmp = cconv(tmp, chromProteins{i, 3}(j) * ones(ch.complexDNAFootprints(protIdx), 1), ch.sequenceLen);
                        
                        %bin to calculate density
                        densityBEMatrix(:, i) = ...
                            + densityBEMatrix(:, i) ...
                            + sum(reshape([tmp; tmp(1:1000 - mod(ch.sequenceLen, 1000), :)], [], 1000), 1)';
                    end
                end
                densityBEMatrix = densityBEMatrix / ceil(ch.sequenceLen/1000) / numel(states.Time.values); %proteins/(nt h)
                
                nBins = ch.sequenceLen - 1;
                densityBEMatrixOS = zeros(nBins, size(chromProteins, 1));
                for i = 1:size(chromProteins, 1)
                    for j = 1:numel(chromProteins{i, 1})
                        protIdx = pc.getIndexs(chromProteins{i, 1}{j});
                        
                        %calculate frequency binding each position
                        tmp = full(sum(sum(states.Chromosome.complexBoundSites(:, chromProteins{i, 2}, :) == protIdx, 3), 2));
                        
                        %circularly convolve with footprint
                        tmp = cconv(tmp, chromProteins{i, 3}(j) * ones(ch.complexDNAFootprints(protIdx), 1), ch.sequenceLen);
                        
                        %bin to calculate density
                        densityBEMatrixOS(:, i) = ...
                            + densityBEMatrixOS(:, i) ...
                            + sum(reshape([tmp; tmp(1:nBins - mod(ch.sequenceLen, nBins), :)], [], nBins), 1)';
                    end
                end
                
                densityBEMatrixOS = densityBEMatrixOS / ceil(ch.sequenceLen/nBins) / numel(states.Time.values); %proteins/(nt h)
                
                
                %save
                save([outDirectory 'circosPlot.mat'], 'densityBEMatrix', 'densityBEMatrixOS');
            end
            
            totalDensityMatrix = Figures34.calculateDNABoundProteinDensity(outDirectory, simBatchDir, selectedSim); %1/(nt s)
            totalDensityMatrix = sum(reshape([totalDensityMatrix; totalDensityMatrix(1:1000 - mod(numel(totalDensityMatrix), 1000), :)], [], 1000), 1)' ...
                / ceil(ch.sequenceLen/1000); %1/(nt s)
            
            densityBEMatrix = [densityBEMatrix totalDensityMatrix];
            %             densityBEMatrix = [densityBEMatrix [densityBEMatrixOS(579491:579560, 2); zeros(860, 1); densityBEMatrixOS(579421:579490, 2)]];
            %             clear('densityBEMatrixOS');
            
            
            %% plot
            if nargin < 5
                [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
                clf(figHandle);
                position = [0.1 0.1 0.8 0.8];
            end
            x = position(1);
            y = position(2);
            w = position(3);
            h = position(4);
            
            figPos = get(figHandle, 'papersize');
            figW = figPos(1);
            figH = figPos(2);
            
            y2 = y+(h-0.75*w*figW/figH)/2 + 0.81*w*figW/figH/2 - 0.5*h;
            
            colorbarHandles = [
                subplot('position', [x+0.727*w  y2            0.03*w  0.2*h], 'parent', figHandle)
                subplot('position', [x+0.727*w  y2+0.8*1/3*h  0.03*w  0.2*h], 'parent', figHandle)
                subplot('position', [x+0.727*w  y2+0.8*2/3*h  0.03*w  0.2*h], 'parent', figHandle)
                subplot('position', [x+0.727*w  y2+0.8*h      0.03*w  0.2*h], 'parent', figHandle)
                ];
            axesHandle = subplot('position', [x  y+(h-0.77*w*figW/figH)/2  0.70*w  0.81*w*figW/figH], 'parent', figHandle);
            
            cla(axesHandle);
            hold(axesHandle, 'on');
            set(axesHandle, 'visible', 'off');
            
            cmap = cell(4, 1);
            cmapLims = zeros(4, 2);
            cmapTicks = cell(4, 1);
            cmapFuncs = {
                [0 1 0]
                [1 0 0]
                [0 0 1]
                [0 1 1]
                };
            totcmap = zeros(0, 3);
            
            nTheta = size(densityBEMatrix, 1);
            for i = 1:size(densityBEMatrix, 2)
                x = zeros(10, size(densityBEMatrix, 1));
                y = zeros(10, size(densityBEMatrix, 1));
                c = zeros(size(densityBEMatrix, 1), 1);
                
                data = densityBEMatrix(:, i);
                
                cmapLims(i, :) = [
                    max(-6, floor(log10(min(data(data > 0)))))
                    ceil(log10(max(data)))
                    ];
                
                if i == 2
                    cmapLims(i, :) = [-5 0];
                end
                
                nBins = 2 + 10 * diff(cmapLims(i, :));
                cmap{i} = ones(nBins, 3) - repmat([1 1 1] - cmapFuncs{i}, nBins, 1) .* repmat(linspace(0, 1, nBins)', 1, 3);
                cmapTicks{i} = [0 logspace(cmapLims(i, 1), cmapLims(i, 2), 1 + 10 * diff(cmapLims(i, :)))];
                
                for j = 1:nTheta
                    x(:, j) = [
                        (i+1)*cos(linspace(2*pi/nTheta*j, 2*pi/nTheta*(j+1), 5)+pi/2) ...
                        (i+2)*cos(linspace(2*pi/nTheta*(j+1), 2*pi/nTheta*j, 5)+pi/2) ...
                        ];
                    y(:, j) = [
                        (i+1)*sin(linspace(2*pi/nTheta*j, 2*pi/nTheta*(j+1), 5)+pi/2) ...
                        (i+2)*sin(linspace(2*pi/nTheta*(j+1), 2*pi/nTheta*j, 5)+pi/2) ...
                        ];
                    c(j) = find(cmapTicks{i} >= data(j), 1, 'first');
                end
                
                c = c + size(totcmap, 1);
                totcmap = [totcmap; cmap{i}]; %#ok<AGROW>
                
                ptch = patch(x, y, ones(1, nTheta), 'EdgeColor', 'none', 'Parent', axesHandle);
                set(ptch, 'FaceColor', 'flat', 'FaceVertexCData', c, 'CDataMapping', 'direct', 'FaceAlpha', 1, 'EdgeAlpha', 1);
            end
            colormap(axesHandle, totcmap);
            
            xlims = [-7.8 6.2];
            xlim(axesHandle, xlims);
            ylim(axesHandle, range(xlims)*0.81/0.7/2 * [-1 1]);
            
            %oriC label
            line([0 0], [2 6.2], 'Color', 'k', 'Parent', axesHandle);
            text(0, 6.18, 'oriC', ...
                'Parent', axesHandle, ...
                'VerticalAlign', 'bottom', ...
                'HorizontalAlign', 'center', ...
                'FontSize', 5);
            
            % rRNA label
            
            % "Box"
            s1 = 0.85;
            s2 = 1 + (1 - s1);
            theta1 = (170007)/580076 * (2*pi);
            theta2 = (170007+4786)/580076 * (2*pi);
            nDivs = 500;
            x = zeros(nDivs, 1);
            y = zeros(nDivs, 1);
            z = ones(nDivs, 1);
            for i = 1:nDivs
                x(:, i) =  [
                    6.25 * cos(linspace(s1 * theta1, s2 * theta2, nDivs/2) + pi/2)...
                    7.65 * cos(linspace(s2 * theta2, s1 * theta1, nDivs/2) + pi/2)...
                    ];
                y(:, i) =  [
                    6.25 * sin(linspace(s1 * theta1, s2 * theta2, nDivs/2) + pi/2)...
                    7.65 * sin(linspace(s2 * theta2, s1 * theta1, nDivs/2) + pi/2)...
                    ];
            end
            p = patch(x, y, z, 'EdgeColor', 'k', 'Parent', axesHandle);
            set(p, 'FaceColor', 'flat', 'CDataMapping', 'direct', 'FaceAlpha', 1);
            
            line([(4.9 * (cos(theta1 + pi/2))) (6.25 * (cos(s1 * theta1 + pi/2)))], ...
                [(4.9 * (sin(theta1 + pi/2))) (6.25 * (sin(s1 * theta1 + pi/2)))], ...
                'Color', 'k', 'Parent', axesHandle);
            line([(4.9 * (cos(theta2 + pi/2))) (6.25 * (cos(s2 * theta2 + pi/2)))], ...
                [(4.9 * (sin(theta2 + pi/2))) (6.25 * (sin(s2 * theta2 + pi/2)))], ...
                'Color', 'k', 'Parent', axesHandle);
            
            theta = mean([theta1 theta2]);
            text((6.3 * (cos(theta + pi/2))), (6.25 * (sin(theta + pi/2))), ...
                {'rRNA' '5S, 16S, 23S'}, ...
                'Parent', axesHandle, ...
                'VerticalAlign', 'bottom', ...
                'HorizontalAlign', 'center', ...
                'FontSize', 5, ...
                'Rotation', theta * 180/pi);
            
            % DnaA label
            
            % "Box"
            s1 = 0.97;
            s2 = 1 + (1 - s1);
            theta1 = 997/1000 * (2*pi);
            theta2 = 997/1000 * (2*pi);
            nDivs = 500;
            off = pi / 10;
            x = zeros(nDivs, 1);
            y = zeros(nDivs, 1);
            z = ones(nDivs, 1);
            for i = 1:nDivs
                x(:, i) =  [
                    6.25 * cos(linspace(s1 * theta1, s2 * theta2, nDivs/2) + pi/2 - off)...
                    7.65 * cos(linspace(s2 * theta2, s1 * theta1, nDivs/2) + pi/2 - off)...
                    ];
                y(:, i) =  [
                    6.25 * sin(linspace(s1 * theta1, s2 * theta2, nDivs/2) + pi/2 - off)...
                    7.65 * sin(linspace(s2 * theta2, s1 * theta1, nDivs/2) + pi/2 - off)...
                    ];
            end
            p = patch(x, y, z, 'EdgeColor', 'k', 'Parent', axesHandle);
            set(p, 'FaceColor', 'flat', 'CDataMapping', 'direct', 'FaceAlpha', 1);
            
            line([(3.7 * (cos(theta1 + pi/2))) (6.25 * (cos(s1 * theta1 + pi/2 - off)))], ...
                [(3.7 * (sin(theta1 + pi/2))) (6.25 * (sin(s1 * theta1 + pi/2 - off)))], ...
                'Color', 'k', 'Parent', axesHandle);
            line([(3.7 * (cos(theta2 + pi/2))) (6.25 * (cos(s2 * theta2 + pi/2 - off)))], ...
                [(3.7 * (sin(theta2 + pi/2))) (6.25 * (sin(s2 * theta2 + pi/2 - off)))], ...
                'Color', 'k', 'Parent', axesHandle);
            
            data = densityBEMatrixOS(579434:579550, 2);
            nTheta = numel(data) + 8;
            x = zeros(10, numel(data));
            y = zeros(10, numel(data));
            c = zeros(numel(data), 1);
            
            arcLen = s2 * theta2 - s1 * theta1;
            off = off + arcLen / 2;
            for j = 1:numel(data)
                x(:, j) = [
                    (6.25)*cos(linspace(arcLen/nTheta*j, arcLen/nTheta*(j+1), 5)+pi/2-off) ...
                    (7.65)*cos(linspace(arcLen/nTheta*(j+1), arcLen/nTheta*j, 5)+pi/2-off) ...
                    ];
                y(:, j) = [
                    (6.25)*sin(linspace(arcLen/nTheta*j, arcLen/nTheta*(j+1), 5)+pi/2-off) ...
                    (7.65)*sin(linspace(arcLen/nTheta*(j+1), arcLen/nTheta*j, 5)+pi/2-off) ...
                    ];
                c(j) = find(cmapTicks{2} >= data(j), 1, 'first');
            end
            
            c = c + size(cmap{1},1);
            
            ptch = patch(x, y, ones(1, numel(data)), 'EdgeColor', 'none', 'Parent', axesHandle);
            set(ptch, 'FaceColor', 'flat', 'FaceVertexCData', c, 'CDataMapping', 'direct', 'FaceAlpha', 1, 'EdgeAlpha', 1);
            
            %             text(5.4, 6.6, ...
            %                 {'DnaA complex'}, ...
            %                 'Parent', axesHandle, ...
            %                 'VerticalAlign', 'bottom', ...
            %                 'HorizontalAlign', 'center', ...
            %                 'FontSize', 5, ...
            %                 'Rotation', 0);
            
            theta = mean([(s1 * theta1) (s2 * theta2)]) - off + arcLen/2;
            text(2.9, 7.80, ...
                {'DnaA complex'}, ...
                'Parent', axesHandle, ...
                'VerticalAlign', 'bottom', ...
                'HorizontalAlign', 'center', ...
                'FontSize', 5, ...
                'Rotation', theta * 180/pi + 2);
            
            text(1.34, 7.42, '*', 'Parent', axesHandle, 'VerticalAlign', 'bottom', ...
                'HorizontalAlign', 'center', 'FontSize', 5, 'Rotation', 0);
            text(2.30, 7.19, '*', 'Parent', axesHandle, 'VerticalAlign', 'bottom', ...
                'HorizontalAlign', 'center', 'FontSize', 5, 'Rotation', 0);
            text(2.64, 7.07, '*', 'Parent', axesHandle, 'VerticalAlign', 'bottom', ...
                'HorizontalAlign', 'center', 'FontSize', 5, 'Rotation', 0);
            text(3.75, 6.58, '*', 'Parent', axesHandle, 'VerticalAlign', 'bottom', ...
                'HorizontalAlign', 'center', 'FontSize', 5, 'Rotation', 0);
            
            %color bars
            ylabels = {
                {'DNA pol' 'prob. bound'}
                {'DnaA' 'prob. bound'}
                {'RNA pol' 'prob. bound'}
                {'All proteins' 'prob. bound'}
                };
            for i = 1:numel(colorbarHandles)
                image(permute(cmap{i}, [1 3 2]), 'Parent', colorbarHandles(i));
                
                set(colorbarHandles(i), ...
                    'TickDir', 'out', ...
                    'box', 'off', ...
                    'FontSize', 5, ...
                    'YAxisLocation', 'right', ...
                    'XTick', [], ...
                    'YTick', 2:10:size(cmap{i}, 1), ...
                    'XLim', [0.5 1.5], ...
                    'visible', 'on', ...
                    'YDir', 'normal');
                ylim(colorbarHandles(i), [0 size(cmap{i}, 1)]+0.5);
                
                line([min(xlim(colorbarHandles(i))) * [1 1]   max(xlim(colorbarHandles(i)))], ...
                    [ylim(colorbarHandles(i)) max(ylim(colorbarHandles(i)))], ...
                    'parent', colorbarHandles(i), 'color', 'k');
                
                tick2text(colorbarHandles(i), 'axis', 'y', 'ytickoffset', -1.25);
                yTicks = getappdata(colorbarHandles(i), 'YTickText');
                set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'left', 'Interpreter', 'tex');
                
                for j = 1:numel(yTicks)
                    set(yTicks(j), 'String', sprintf('10^{%d}', cmapLims(i, 1)+j-1));
                end
                if numel(yTicks) == 6
                    delete(yTicks(1:2:5));
                end
                
                ylabel(colorbarHandles(i), ylabels{i}, 'FontSize', 5);
                ylabelpos = get(get(colorbarHandles(i), 'ylabel'), 'position');
                set(get(colorbarHandles(i), 'ylabel'), 'position', [3.75 ylabelpos(2:end)]);
            end
        end
        
        function dnaBoundProteinDensityVsDisplacement(simBatchDir, selectedSim, outDirectory, figHandle, position)
            import edu.stanford.covert.cell.sim.analysis.Figures34;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if ~exist('selectedSim', 'var') || isempty(selectedSim)
                selectedSim = 1;
            end
            
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep num2str(selectedSim)]);
            c = sim.state('Chromosome');
            
            %% get data
            dnaBoundProteinDensity = Figures34.calculateDNABoundProteinDensity(outDirectory, simBatchDir, selectedSim);
            
            dnaBoundProteinDisplacements = histc(...
                Figures34.calculateDNABoundProteinDisplacements(outDirectory, simBatchDir, selectedSim), ...
                1:c.sequenceLen);
            dnaBoundProteinDisplacements = conv(dnaBoundProteinDisplacements, ones(1000, 1), 'same');
            
            %% plot
            if nargin >= 5
                x = position(1);
                y = position(2);
                w = position(3);
                h = position(4);
                axesHandle = subplot('position', [x y 0.82*w h], 'parent', figHandle);
                colorBarHandle = subplot('position', [x+0.85*w y 0.05*w h], 'parent', figHandle);
            else
                axesHandle = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
                colorBarHandle = colorbar();
            end
            
            % data
            [cnts, centers] = hist3([dnaBoundProteinDensity dnaBoundProteinDisplacements], [100 100]);
            cnts = cnts / c.sequenceLen;
            edges = logspace(-6, -1, 51);
            colormap = jet(numel(edges) + 1);
            colors = zeros(100, 100, 3);
            
            [i, j] = find(cnts < edges(1));
            colors(sub2ind(size(colors), i, j, ones(size(i)))) = colormap(1, 1);
            colors(sub2ind(size(colors), i, j, 2*ones(size(i)))) = colormap(1, 2);
            colors(sub2ind(size(colors), i, j, 3*ones(size(i)))) = colormap(1, 3);
            
            [i, j] = find(cnts > edges(end-1));
            colors(sub2ind(size(colors), i, j, ones(size(i)))) = colormap(end, 1);
            colors(sub2ind(size(colors), i, j, 2*ones(size(i)))) = colormap(end, 2);
            colors(sub2ind(size(colors), i, j, 3*ones(size(i)))) = colormap(end, 3);
            
            for idx = 1:numel(edges)-1
                [i, j] = find(cnts >= edges(idx) & cnts <= edges(idx+1));
                colors(sub2ind(size(colors), i, j, ones(size(i)))) = colormap(idx+1, 1);
                colors(sub2ind(size(colors), i, j, 2*ones(size(i)))) = colormap(idx+1, 2);
                colors(sub2ind(size(colors), i, j, 3*ones(size(i)))) = colormap(idx+1, 3);
            end
            
            image(permute(colors, [2 1 3]), 'Parent', axesHandle);
            set(axesHandle, 'XDir', 'normal');
            set(axesHandle, 'YDir', 'normal');
            
            set(axesHandle, 'Box', 'on', 'FontSize', 5, 'TickDir', 'out');
            xlim(axesHandle, [0.5 0.4/centers{1}(end)*100])
            ylim(axesHandle, [0.5 1500/centers{2}(end)*100])
            
            set(axesHandle, 'XTick', 0.5 + (0:0.1:0.4)/centers{1}(end)*100);
            set(axesHandle, 'YTick', 0.5 + (0:500:1500)/centers{2}(end)*100);
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', 0.06, 'ytickoffset', 0.015);
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            set(xTicks, 'FontSize', 5);
            set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'right');
            set(xTicks(1), 'String', '0')
            set(xTicks(end), 'String', '0.4')
            set(yTicks(1), 'String', '0')
            set(yTicks(end), 'String', '1500')
            delete(yTicks(2:3));
            delete(xTicks(2:4));
            
            xlabel(axesHandle, 'Density', 'FontSize', 7);
            ylabel(axesHandle, 'Displacement', 'FontSize', 7);
            xlabelPos = get(get(axesHandle, 'xlabel'), 'position');
            ylabelPos = get(get(axesHandle, 'ylabel'), 'position');
            set(get(axesHandle, 'xlabel'), 'position', [xlabelPos(1) -3 xlabelPos(3:end)]);
            set(get(axesHandle, 'ylabel'), 'position', [-3 ylabelPos(2:end)]);
            
            % color bar
            image(permute(colormap, [1 3 2]), 'Parent', colorBarHandle);
            
            box(colorBarHandle, 'on');
            set(colorBarHandle, 'XTick', []);
            set(colorBarHandle, 'TickDir', 'out');
            set(colorBarHandle, 'YAxisLocation', 'right');
            set(colorBarHandle, 'YDir', 'normal')
            set(colorBarHandle, 'FontSize', 5);
            
            set(colorBarHandle, 'YTick', 2:10:52);
            tick2text(colorBarHandle, 'axis', 'y', 'ytickoffset', -1.15);
            yTicks = getappdata(colorBarHandle, 'YTickText');
            set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'left', 'Interpreter', 'tex');
            for i = 1:6
                set(yTicks(i), 'String', sprintf('10^{%d}', -7+i));
            end
            
            ylabel(colorBarHandle, '% Bases', 'FontSize', 7);
            ylabelPos = get(get(colorBarHandle, 'ylabel'), 'position');
            set(get(colorBarHandle, 'ylabel'), 'position', [3.0 ylabelPos(2:end)]);
        end
        
        function [posFreq, protDensity] = analyzeDNABoundProteinDisplacements(simBatchDir, outDirectory, figHandle, position)
            import edu.stanford.covert.cell.sim.analysis.Figures34;
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.PrintUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SummaryLogger;
            
            %% constants
            sim = CachedSimulationObjectUtil.load();
            
            g = sim.gene;
            c = sim.state('Chromosome');
            transcript = sim.state('Transcript');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            rna = sim.state('Rna');
            transcription = sim.process('Transcription');
            
            %% get data
            nMon = numel(pm.matureIndexs);
            nCpx = numel(pc.matureIndexs);
            nProt = nMon + nCpx;            
            
            if ~exist([outDirectory filesep 'dnaBoundProteinDisplacement.mat'], 'file')
                nSims = SimulationDiskUtil.getNumSimulations(simBatchDir);
                
                ensemble = SimulationEnsemble(simBatchDir, cell(0, 2), [], 1:nSims);
                simEndTimes = ensemble.stateData.simulationEndTimes;
                
                nDisplacements = 0;
                allDisplacements = zeros(0, 5);
                protByProtFreq = zeros(nProt);
                timeFreq = zeros(1000, 1);
                posFreq = zeros(c.sequenceLen, 2);
                for i = 1:nSims
                    dnaBoundProteinDisplacements = Figures34.calculateDNABoundProteinDisplacements(outDirectory, simBatchDir, i);
                    
                    nDisplacements = nDisplacements + size(dnaBoundProteinDisplacements, 1);
                    allDisplacements = [allDisplacements; dnaBoundProteinDisplacements]; %#ok<AGROW>
                    
                    tmp = dnaBoundProteinDisplacements(:, 1:2);
                    tmp(tmp < 0) = -tmp(tmp < 0) + nMon;
                    tmpCnt = hist3(tmp, 'Edges', {1:nProt 1:nProt});
                    assert(sum(tmpCnt(:)) == size(dnaBoundProteinDisplacements, 1));
                    protByProtFreq = protByProtFreq + tmpCnt;
                    
                    timeFreq = timeFreq + histc(dnaBoundProteinDisplacements(:, 5), linspace(1, simEndTimes(i), 1000));
                    posFreq(:, 1) = posFreq(:, 1) + histc(dnaBoundProteinDisplacements(:, 3), 1:c.sequenceLen);
                    posFreq(:, 2) = posFreq(:, 2) + histc(dnaBoundProteinDisplacements(:, 4), 1:c.sequenceLen);
                end
                
                simStats = SummaryLogger.getSimulationStatistics(simBatchDir);
                repInitDuration = simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME);
                repDuration = diff(simStats(:, [SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME SummaryLogger.SIM_STATUS_INDEX_REP_TIME]), 1, 2);
                cytokinesisDuration = diff(simStats(:, [SummaryLogger.SIM_STATUS_INDEX_REP_TIME SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME]), 1, 2);
                
                save([outDirectory filesep 'dnaBoundProteinDisplacement.mat'], ...
                    'nSims', 'simEndTimes', 'allDisplacements', 'nDisplacements', ...
                    'protByProtFreq', 'timeFreq', 'posFreq', ...
                    'repInitDuration', 'repDuration', 'cytokinesisDuration')
            else
                load([outDirectory filesep 'dnaBoundProteinDisplacement.mat']);
            end
            
            [mat, allIds, ids, names, shortNames] = Figures34.groupProteinDisplacements(protByProtFreq, simEndTimes, sim);
            
            monIds = intersect(allIds, pm.wholeCellModelIDs);
            cpxIds = intersect(allIds, pm.wholeCellModelIDs);
            
            [monTfs, monIdxs] = ismember(allIds, pm.wholeCellModelIDs(pm.matureIndexs));
            monIdxs = monIdxs(monTfs);
            
            [cpxTfs, cpxIdxs] = ismember(allIds, pc.wholeCellModelIDs(pc.matureIndexs));
            cpxIdxs = cpxIdxs(cpxTfs);
            
            if ~exist([outDirectory filesep 'proteinDensity.mat'], 'file')
                protDensity = zeros(numel(allIds), c.sequenceLen);
                for i = 1:nSims
                    states = SimulationEnsemble.load(simBatchDir, {
                        'Chromosome' 'monomerBoundSites'
                        'Chromosome' 'complexBoundSites'
                        }, [], [], 1, 'extract', i);
                    
                    [subs, vals] = find(states.Chromosome.monomerBoundSites);
                    tfs = ismember(vals, monIdxs);
                    protDensity(monTfs, :) = protDensity(monTfs, :) + hist3([vals(tfs) subs(tfs, 1)], {monIdxs 1:c.sequenceLen});
                    
                    [subs, vals] = find(states.Chromosome.complexBoundSites);
                    tfs = ismember(vals, cpxIdxs);
                    protDensity(cpxTfs, :) = protDensity(cpxTfs, :) + hist3([vals(tfs) subs(tfs, 1)], {cpxIdxs 1:c.sequenceLen});
                    
                    clear states;
                end
                
                save([outDirectory filesep 'proteinDensity.mat'], 'protDensity');
            else
                load([outDirectory filesep 'proteinDensity.mat']);
            end
            
            %%figure
            x = position(1);
            y = position(2);
            W = position(3);
            H = position(4);
            paperSize = get(figHandle, 'PaperSize');
            
            %A
            w = H*paperSize(2)/paperSize(1);
            Figures34.plotCollisionProteinDistribution(figHandle, [x y w H], [x+w+0.02 y 0.02 H], mat, shortNames);
            axesHandle = findobj(figHandle, 'Position', [x y w H]);
            xlabelPos = get(get(axesHandle, 'xlabel'), 'position');
            ylabelPos = get(get(axesHandle, 'ylabel'), 'position');
            set(get(axesHandle, 'xlabel'), 'FontSize', 7, 'position', [xlabelPos(1) 16.4 xlabelPos(3:end)]);
            set(get(axesHandle, 'ylabel'), 'FontSize', 7, 'position', [-2.3 ylabelPos(2:end)]);
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            for i = 1:numel(xTicks)
                set(xTicks(i), 'position', [i 13.75 0])
                set(yTicks(i), 'position', [0.25 i 0])
            end
            
            %B
            x = x+w+0.18;
            axesHandle = subplot('Position', [x y 1-x H], 'Parent', figHandle);
            tmpX = sum(reshape(sum(protDensity(:, 1:580000), 1), 1000, []) / sum(simEndTimes) / 1000 * 1e3, 1);
            tmpY = sum(reshape(posFreq(1:580000, 1)', 1000, []) / sum(simEndTimes) / 1000 * 3600 * 1e3, 1);
            plot(axesHandle, tmpX, tmpY, '.');
            xlim(axesHandle, [min(tmpX) 2.05] + range(tmpX)*0.025*[-1 0]);
            ylim(axesHandle, [min(tmpY) 14.2] + range(tmpY)*0.025*[-1 1]);
            set(axesHandle, 'FontSize', 5, 'box', 'off', 'tickDir', 'out', 'XTick', 0.5:0.5:3.5);
            xlabel(axesHandle, 'Density (knt^{-1})', 'FontSize', 7);
            ylabel(axesHandle, 'Collisions (knt^{-1} h^{-1})', 'FontSize', 7);
            xlabelPos = get(get(axesHandle, 'xlabel'), 'position');
            ylabelPos = get(get(axesHandle, 'ylabel'), 'position');
            set(get(axesHandle, 'xlabel'), 'FontSize', 7, 'position', [xlabelPos(1) -1.2 xlabelPos(3:end)]);
            set(get(axesHandle, 'ylabel'), 'FontSize', 7, 'position', [0.03 ylabelPos(2:end)]);
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', 0.05, 'ytickoffset', 0.03);
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');            
            set(xTicks, 'FontSize', 5);
            set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'right');
            delete(xTicks(1:2:end));
            
            %% supplement
            fprintf('Displacement rate: %0.2f (1/s)\n', nDisplacements / sum(simEndTimes));
            
            tmp = num2cell(mat);
            tmp(mat == 0) = {''};
            content = [
                [cell(3, 3) [ids names shortNames]']
                [ids names shortNames tmp]
                ];
            PrintUtil.printToFile(content, {}, [outDirectory filesep 'dnaBoundProteinDisplacement.xls'], 'DNA Bound Protein Displacement');
            
            monIdxs = setdiff(unique([
                find(any(c.reactionMonomerCatalysisMatrix, 1))'
                c.reactionBoundMonomer
                ]), 0);
            cpxIdxs = setdiff(unique([
                find(any(c.reactionComplexCatalysisMatrix, 1))'
                c.reactionBoundComplex
                ]), 0);
            content = [[
                pm.wholeCellModelIDs(pm.matureIndexs(monIdxs))
                pc.wholeCellModelIDs(pc.matureIndexs(cpxIdxs))
                ] [
                pm.names(pm.matureIndexs(monIdxs))
                pc.names(pc.matureIndexs(cpxIdxs))
                ] cell(numel(monIdxs) + numel(cpxIdxs), 1)];
            tfs = ismember(content(:, 1), allIds);
            content(tfs, 3) = {'Y'};
            PrintUtil.printToFile(content, {'ID' 'Name' 'Displaced/Displacing'}, [outDirectory filesep 'dnaBoundProteinDisplacement.xls'], 'Displac(ed|ing) Proteins');
                        
            %protein distribution
            w = 11.4; h = 9.7; % Cell 1.5 column
            [~, figHandle] = PlotUtil.newAxesHandle([w h]);
            clf(figHandle);
            
            Figures34.plotCollisionProteinDistribution(figHandle, [0.12 0.125 0.74 0.74*w/h], [0.90 0.125 0.03 0.74*w/h], mat, shortNames);
            
            saveas(figHandle, [outDirectory filesep 'dnaBoundProteinDisplacement-proteins.pdf']);
            close(figHandle);
            
            %time, spatial distribution
            w = 11.4; h = 5; % Cell 1.5 column
            [~, figHandle] = PlotUtil.newAxesHandle([w h]);
            clf(figHandle);
            
            
            axesHandle = subplot('Position', [0.07 0.12 0.40 0.835], 'Parent', figHandle);
            x = linspace(1, 100, 1000);
            plot(axesHandle, x(2:end-1), timeFreq(2:end-1) / nSims / (sim.state('Time').cellCycleLength / 1000));
            xlabel(axesHandle, 'Time (% Cell Cycle)', 'FontSize', 7);
            ylabel(axesHandle, {'Displacements (s^{-1})'}, 'FontSize', 7);
            xlim(axesHandle, [0 100])
            ylim(axesHandle, [0.65 1.3]);
            set(axesHandle, 'box', 'off', 'tickdir', 'out', 'FontSize', 5, 'XTick', 0:25:100, 'YTick', 0.75:0.25:1.25);
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', 0.035, 'ytickoffset', 0.02);
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            set(xTicks, 'FontSize', 5);
            set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'right');
            ylabelPos = get(get(axesHandle, 'ylabel'), 'position');
            set(get(axesHandle, 'ylabel'), 'position', [-10 ylabelPos(2:end)]);
            
            cellCycleLength = nanmean(repInitDuration + repDuration + cytokinesisDuration);
            y = [ylim(axesHandle) fliplr(ylim(axesHandle))];
            x = [0 0 nanmean(repInitDuration) * [1 1]] / cellCycleLength * 100;
            ptch = patch(x, y, 1, 'Parent', axesHandle, 'EdgeColor', 'none', 'FaceColor', 1 - [0.1 0.1 0.1], 'FaceAlpha', 1, 'EdgeAlpha', 1);
            uistack(ptch, 'bottom');
            x = [nanmean(repInitDuration) * [1 1]  nanmean(repInitDuration+repDuration) * [1 1]] / cellCycleLength * 100;
            ptch = patch(x, y, 1, 'Parent', axesHandle, 'EdgeColor', 'none', 'FaceColor', 1 - [0.15 0.15 0.15], 'FaceAlpha', 1, 'EdgeAlpha', 1);
            uistack(ptch, 'bottom');
            x = [nanmean(repInitDuration+repDuration) / cellCycleLength * [1 1] 1 1] * 100;
            ptch = patch(x, y, 1, 'Parent', axesHandle, 'EdgeColor', 'none', 'FaceColor', 1 - [0.1 0.1 0.1], 'FaceAlpha', 1, 'EdgeAlpha', 1);
            uistack(ptch, 'bottom');
            line(xlim(axesHandle), min(ylim(axesHandle)) * [1 1], 'Color', 'k', 'Parent', axesHandle);
            line(min(xlim(axesHandle)) * [1 1], ylim(axesHandle), 'Color', 'k', 'Parent', axesHandle);
            text(nanmean(0.5 * repInitDuration)/ cellCycleLength * 100, min(ylim(axesHandle)) + 1.03 * range(ylim(axesHandle)), ...
                'Replication Initiation', ...
                'FontSize', 6, 'HorizontalAlign', 'center');
            text(nanmean(repInitDuration + 0.5 * repDuration)/ cellCycleLength * 100, min(ylim(axesHandle)) + 1.03 * range(ylim(axesHandle)), ...
                'Replication', ...
                'FontSize', 6, 'HorizontalAlign', 'center');
            text(nanmean(repInitDuration + repDuration + 0.5 * cytokinesisDuration)/ cellCycleLength * 100, min(ylim(axesHandle)) + 1.03 * range(ylim(axesHandle)), ...
                'Cytokinesis', ...
                'FontSize', 6, 'HorizontalAlign', 'center');
            
            
            axesHandle = subplot('Position', [0.60 0.12 0.40 0.845], 'Parent', figHandle);
            tmp = [
                cconv(posFreq(:, 1) / nSims, ones(1000, 1) / 1000, c.sequenceLen) ...
                cconv(posFreq(:, 2) / nSims, ones(1000, 1) / 1000, c.sequenceLen)
                ] / sum(simEndTimes) * 1e9;
            h = plot(axesHandle, (1:c.sequenceLen)' / 1e3, tmp);
            xlabel(axesHandle, 'Position (kb)', 'FontSize', 7);
            ylabel(axesHandle, {'Displacements' '(10^{-9} s^{-1})'}, 'FontSize', 7);
            xlim(axesHandle, [1 c.sequenceLen] / 1e3)
            ylim(axesHandle, [5 45]);
            set(axesHandle, 'box', 'off', 'tickdir', 'out', 'FontSize', 5, 'YTick', 10:10:40);            
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', 0.035, 'ytickoffset', 0.02);
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            set(xTicks, 'FontSize', 5);
            set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'right');            
            legend(h, {'Unbinding'; 'Binding'}, 'Location', 'NorthEast')
            legend(axesHandle, 'boxoff');
            
            [peaks, locs] = findpeaks(tmp(:, 1), 'NPEAKS', 5, 'SORTSTR', 'descend', 'MINPEAKDISTANCE', 5000);
            
            g = sim.gene;
            for i = 1:numel(locs)
                gIdx = find(g.startCoordinates <= locs(i) & g.startCoordinates + g.lengths - 1 >= locs(i));
                text(locs(i) / 1e3, peaks(i), sprintf('%s\n(%s)', g.names{gIdx}, g.wholeCellModelIDs{gIdx}), ...
                    'FontSize', 5, 'HorizontalAlign', 'center', 'VerticalAlign', 'bottom', 'Interpreter', 'none');
            end
            
            saveas(figHandle, [outDirectory filesep 'dnaBoundProteinDisplacement-time-space.pdf']);
            close(figHandle);
            
            %% get data
            dnaBoundProteinDisplacements = histc(allDisplacements(:, 3), 1:c.sequenceLen);
            dnaBoundProteinDisplacements = cconv(dnaBoundProteinDisplacements / nSims, ones(1000, 1) / 1000, c.sequenceLen) / sum(simEndTimes) * 1e9;
            
            %% find meaning of peaks
            proteinIdxs = zeros(size(ids));
            [~, tmp] = ismember(ids, pm.wholeCellModelIDs(pm.matureIndexs));
            proteinIdxs(tmp ~= 0) = tmp(tmp ~= 0);
            [~, tmp] = ismember(ids, pc.wholeCellModelIDs(pc.matureIndexs));
            proteinIdxs(tmp ~= 0) = -tmp(tmp ~= 0);
                        
            [~, peaks] = findpeaks(dnaBoundProteinDisplacements, 'MINPEAKDISTANCE', 5000, 'NPEAKS', 10, 'SORTSTR', 'descend');
            tuStarts = transcript.transcriptionUnitFivePrimeCoordinates - (transcript.transcriptionUnitLengths - 1) .* (transcript.transcriptionUnitDirections == 0);
            tuEnds = transcript.transcriptionUnitFivePrimeCoordinates + (transcript.transcriptionUnitLengths - 1) .* (transcript.transcriptionUnitDirections == 1);
            peakGenes = NaN(size(peaks));
            peakTUs = NaN(size(peaks));
            for i = 1:numel(peaks)
                tmp = find(g.startCoordinates <= peaks(i) & g.startCoordinates + g.lengths - 1 >= peaks(i));
                if ~isempty(tmp)
                    peakGenes(i) = tmp;
                end
                tmp = find(tuStarts <= peaks(i) & tuEnds >= peaks(i));
                if ~isempty(tmp)
                    peakTUs(i) = tmp;
                end
                
                if isnan(peakGenes(i))
                    continue;
                end
                
                tfs = ...
                    allDisplacements(:, 3) >= g.startCoordinates(peakGenes(i)) - 500  & ...
                    allDisplacements(:, 3) <= g.startCoordinates(peakGenes(i)) + g.lengths(peakGenes(i)) - 1 + 500;
                idxs = allDisplacements(tfs, 1:2);
                idxs(idxs < 0) = -idxs(idxs < 0) + nMon;
                protByProtFreq = hist3(idxs, 'Edges', {1:nMon+nCpx  1:nMon+nCpx});
                [mat, allIds, ids, names, shortNames] = Figures34.groupProteinDisplacements(protByProtFreq, simEndTimes, sim);
                
                w = 11.4; h = 9.7; % Cell 1.5 column
                [~, figHandle] = PlotUtil.newAxesHandle([w h]);
                clf(figHandle);
                Figures34.plotCollisionProteinDistribution(figHandle, [0.12 0.125 0.74 0.74*w/h], [0.90 0.125 0.03 0.74*w/h], mat, shortNames);
                saveas(figHandle, [outDirectory filesep sprintf('dnaBoundProteinDisplacement-proteins-%s.pdf', g.wholeCellModelIDs{peakGenes(i)})]);
                close(figHandle);
            end
            [g.wholeCellModelIDs(peakGenes(~isnan(peakGenes))) ...  %most displaced genes
                g.names(peakGenes(~isnan(peakGenes))) ...
                num2cell(dnaBoundProteinDisplacements(peaks(~isnan(peakGenes)))) ...
                num2cell(rna.nascentRNAGeneComposition(peakGenes(~isnan(peakGenes)), :) * transcription.transcriptionUnitBindingProbabilities)
                ] %#ok<NOPRT>
            [transcript.transcriptionUnitWholeCellModelIDs(peakTUs(~isnan(peakTUs))) ... %most displaced transcription units
                num2cell(dnaBoundProteinDisplacements(peaks(~isnan(peakTUs)))) ...
                num2cell(transcription.transcriptionUnitBindingProbabilities(peakTUs(~isnan(peakTUs))))] %#ok<NOPRT>
            
            %% correlate displacement spatial frequency with SMC density            
            [monTfs, monIdxs] = ismember(allIds, pm.wholeCellModelIDs(pm.matureIndexs));
            [cpxTfs, cpxIdxs] = ismember(allIds, pc.wholeCellModelIDs(pc.matureIndexs));
            tmp = zeros(size(protDensity));
            for i = 1:size(protDensity, 1)
                if monTfs(i)
                    tmp(i, :) = cconv(protDensity(i, :)', ones(c.monomerDNAFootprints(monIdxs(i)), 1), c.sequenceLen);
                else
                    tmp(i, :) = cconv(protDensity(i, :)', ones(c.complexDNAFootprints(cpxIdxs(i)), 1), c.sequenceLen);
                end
            end
            
            figHandle = figure();
            plot(sum(reshape(sum(protDensity(:, 1:580000), 1), 1000, []) / sum(simEndTimes) / 1000, 1), ...
                sum(reshape(posFreq(1:580000, 1)', 1000, []) / sum(simEndTimes) / 1000, 1), ...
                '.');
            xlabel('Protein Density (nt^{-1})');
            ylabel('Collisions Density (nt^{-1} s^{-1})');
            saveas(figHandle, [outDirectory filesep 'proteinVsCollisionDensity.pdf']);
            close(figHandle);
            
            figHandle = figure();
            plot(sum(reshape(protDensity(strcmp('MG_213_214_298_6MER_ADP', allIds), 1:580000), 1000, []) / sum(simEndTimes) / 1000, 1), ...
                sum(reshape(posFreq(1:580000, 1)', 1000, []) / sum(simEndTimes) / 1000, 1), ...
                '.');
            xlabel('SMC Density');
            ylabel('Collisions Density');
            saveas(figHandle, [outDirectory filesep 'smcVsCollisionDensity.pdf']);
            close(figHandle);
        end
        
        function [mat, allIds, ids, names, shortNames] = groupProteinDisplacements(protByProtFreq, simEndTimes, sim)
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            nMon = numel(pm.matureIndexs);
            
            [~, idxs] = ismember({
                'MG_101_MONOMER'
                'MG_236_MONOMER'
                'DNA_GYRASE'
                'DNA_POLYMERASE_2CORE_BETA_CLAMP_GAMMA_COMPLEX_PRIMASE'
                'DNA_POLYMERASE_CORE_BETA_CLAMP_GAMMA_COMPLEX'
                'DNA_POLYMERASE_CORE_BETA_CLAMP_PRIMASE'
                'MG_001_DIMER'
                'MG_091_OCTAMER'
                'MG_094_HEXAMER'
                'MG_203_204_TETRAMER'
                'MG_205_DIMER'
                'MG_213_214_298_6MER_ADP'
                'MG_428_DIMER'
                'MG_469_1MER_ATP'
                'MG_469_1MER_ADP'
                'MG_469_7MER_ATP'
                'RNA_POLYMERASE'
                'RNA_POLYMERASE_HOLOENZYME'
                }, [pm.wholeCellModelIDs(pm.matureIndexs); pc.wholeCellModelIDs(pc.matureIndexs)]);            
            idxs = unique([idxs; find(any(protByProtFreq, 2) | any(protByProtFreq, 1)')]);
            
            mat = protByProtFreq(idxs, idxs) / sum(simEndTimes) * 3600;
            allIds = [
                pm.wholeCellModelIDs(pm.matureIndexs(idxs(idxs <= nMon)))
                pc.wholeCellModelIDs(pc.matureIndexs(idxs(idxs > nMon) - nMon))
                ];
            ids = allIds;
            names = [
                pm.names(pm.matureIndexs(idxs(idxs <= nMon)))
                pc.names(pc.matureIndexs(idxs(idxs > nMon) - nMon))
                ];
            groups = {
                {'RNA_POLYMERASE', 'RNA_POLYMERASE_HOLOENZYME'}
                {'DNA_POLYMERASE_2CORE_BETA_CLAMP_GAMMA_COMPLEX_PRIMASE', 'DNA_POLYMERASE_CORE_BETA_CLAMP_GAMMA_COMPLEX', 'DNA_POLYMERASE_CORE_BETA_CLAMP_PRIMASE'}
                {'MG_469_1MER_ATP', 'MG_469_1MER_ADP' 'MG_469_7MER_ATP'}
                };
            for i = 1:numel(groups)
                [~, idxs] = ismember(groups{i}, ids);
                idxs = idxs(idxs > 0);
                ids = ids(setdiff(1:end, idxs(2:end)));
                names = names(setdiff(1:end, idxs(2:end)));
                mat(:, idxs(1)) = mat(:, idxs(1)) + sum(mat(:, idxs(2:end)), 2);
                mat(idxs(1), :) = mat(idxs(1), :) + sum(mat(idxs(2:end), :), 1);
                mat = mat(setdiff(1:end, idxs(2:end)), setdiff(1:end, idxs(2:end)));
            end
            
            shortNames = cell(size(ids));
            shortNames(strcmp(ids, 'MG_101_MONOMER')) = {'GntR'};
            shortNames(strcmp(ids, 'MG_236_MONOMER')) = {'Fur'};
            shortNames(strcmp(ids, 'DNA_GYRASE')) = {'GyrAB'};
            shortNames(strcmp(ids, 'DNA_POLYMERASE_2CORE_BETA_CLAMP_GAMMA_COMPLEX_PRIMASE')) = {'DNA Pol'};
            shortNames(strcmp(ids, 'DNA_POLYMERASE_CORE_BETA_CLAMP_GAMMA_COMPLEX')) = {'DNA Pol'};
            shortNames(strcmp(ids, 'DNA_POLYMERASE_CORE_BETA_CLAMP_PRIMASE')) = {'DNA Pol'};
            shortNames(strcmp(ids, 'MG_001_DIMER')) = {'DnaN'};
            shortNames(strcmp(ids, 'MG_091_OCTAMER')) = {'SSB'};
            shortNames(strcmp(ids, 'MG_094_HEXAMER')) = {'DnaB'};
            shortNames(strcmp(ids, 'MG_203_204_TETRAMER')) = {'Topo IV'};
            shortNames(strcmp(ids, 'MG_205_DIMER')) = {'HrcA'};
            shortNames(strcmp(ids, 'MG_213_214_298_6MER_ADP')) = {'SMC'};
            shortNames(strcmp(ids, 'MG_428_DIMER')) = {'LuxR'};
            shortNames(strcmp(ids, 'MG_469_1MER_ATP')) = {'DnaA'};
            shortNames(strcmp(ids, 'MG_469_1MER_ADP')) = {'DnaA'};
            shortNames(strcmp(ids, 'MG_469_7MER_ATP')) = {'DnaA'};
            shortNames(strcmp(ids, 'RNA_POLYMERASE')) = {'RNA Pol'};
            shortNames(strcmp(ids, 'RNA_POLYMERASE_HOLOENZYME')) = {'RNA Pol'};
                        
            [~, order] = sort(shortNames);
            ids = ids(order);
            names = names(order);
            shortNames = shortNames(order);
            mat = mat(order, order);
        end
        
        function plotCollisionProteinDistribution(figHandle, position, position2, mat, shortNames)
            paperSize = get(figHandle, 'PaperSize');
            w = paperSize(1);
            h = paperSize(2);
            
            axesHandle = subplot('Position', position, 'Parent', figHandle);
            
            tmp = (log10(mat) - -3.4) / (3.4 - -3.4);
            colorsB = 0.85 * (1 - tmp);
            colorsG = 0.85 * (1 - tmp);
            colorsR = ones(size(tmp));
            colorsR(mat == 0) = 1;
            colorsG(mat == 0) = 1;
            colorsB(mat == 0) = 1;
            colors = cat(3, colorsR, colorsG, colorsB);
            for i = 1:size(colors, 1)
                for j = 1:size(colors, 2)
                    if mat(i, j)
                        patch([j-0.5 j-0.5 j+0.5 j+0.5], [i-0.5 i+0.5 i+0.5 i-0.5], 1, ...
                            'FaceColor', permute(colors(i, j, :), [1 3 2]), ...
                            'EdgeColor', 'none', 'FaceAlpha', 1, 'EdgeAlpha', 1, ...
                            'Parent', axesHandle);
                    end
                end
            end
            set(axesHandle, ...
                'FontSize', 5, 'TickDir', 'out', 'box', 'off', ...
                'GridLineStyle', 'none', ...
                'MinorGridLineStyle', '-', ...
                'XTick', 1:size(mat, 1), 'XAxisLocation', 'bottom', ...
                'YTick', 1:size(mat, 1), 'YDir', 'reverse', 'visible', 'on');
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', -1.02, 'ytickoffset', 0.02);
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            set(xTicks, 'FontSize', 5, 'HorizontalAlign', 'right', 'interpreter', 'none', 'rotation', 90);
            set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'right', 'interpreter', 'none');
            for i = 1:numel(shortNames)
                set(xTicks(i), 'String', shortNames{i}(1:min(end, 24)))
                set(yTicks(i), 'String', shortNames{i}(1:min(end, 24)))
            end
            xlabel(axesHandle, 'Binding', 'FontSize', 7);
            ylabel(axesHandle, 'Unbinding', 'FontSize', 7);
            xlabelPos = get(get(axesHandle, 'xlabel'), 'position');
            ylabelPos = get(get(axesHandle, 'ylabel'), 'position');
            set(get(axesHandle, 'xlabel'), 'position', [xlabelPos(1) 14.8 xlabelPos(3:end)]);
            set(get(axesHandle, 'ylabel'), 'position', [-1 ylabelPos(2:end)]);
            xlim(axesHandle, [0.5 size(mat, 1)+0.5]);
            ylim(axesHandle, [0.5 size(mat, 1)+0.5]);
            for i = 0.5:size(mat, 1)+0.5
                line([0.5 size(mat, 1)+0.5], [i i], 'Parent', axesHandle, 'Color', 'k');
                line([i i], [0.5 size(mat, 1)+0.5], 'Parent', axesHandle, 'Color', 'k');
            end
            axis(axesHandle, 'square');
            
            axesHandle = subplot('Position', position2, 'Parent', figHandle);            
            image(cat(3, ...
                ones(69, 1), ...
                0.85 * (1 - linspace(0, 1, 69))', ...
                0.85 * (1 - linspace(0, 1, 69))'), ...
                'Parent', axesHandle);
            set(axesHandle, ...
                'FontSize', 5, 'TickDir', 'out', 'box', 'off', ...
                'XTick', [], ...
                'YTick', 15:20:55, ...
                'YAxisLocation', 'right', 'YDir', 'normal');
            tick2text(axesHandle, 'axis', 'y', 'ytickoffset', -1.2);
            yTicks = getappdata(axesHandle, 'YTickText');
            set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'left', 'interpreter', 'tex');
            yTickLabels = {'10^{-2}', '10^{0}', '10^{2}'};
            for i = 1:numel(yTicks)
                set(yTicks(i), 'String', yTickLabels{i});
            end
            xlim(axesHandle, [0.5 1.5]);
            ylim(axesHandle, [0.5 69.5]);
            ylabel(axesHandle, 'Freq (h^{-1})', 'FontSize', 7);
            line([0.5 0.5], [0.5 69.5], 'Parent', axesHandle, 'Color', 'k')
            line([0.5 1.5], [69.5 69.5], 'Parent', axesHandle, 'Color', 'k')
        end
        
        function totalDensityMatrix = calculateDNABoundProteinDensity(outDirectory, simBatchDir, selectedSim)
            import edu.stanford.covert.cell.sim.analysis.ChromosomePositionHistogram;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            
            if exist([outDirectory 'DNABoundProteinDensity.mat'], 'file')
                load([outDirectory 'DNABoundProteinDensity.mat'])
            else
                %% get constants
                [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep num2str(selectedSim)]);
                
                stateNames = {
                    'Chromosome' 'complexBoundSites'        ':' ':'
                    'Chromosome' 'monomerBoundSites'        ':' ':'
                    };
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSim);
                
                complexIdxs = full(unique(states.Chromosome.complexBoundSites));
                complexIdxs = complexIdxs(complexIdxs > 0);
                
                monomerIdxs = full(unique(states.Chromosome.monomerBoundSites));
                monomerIdxs = monomerIdxs(monomerIdxs > 0);
                
                boundIndices = struct('complex', {num2cell(complexIdxs')}, 'monomer', {num2cell(monomerIdxs')});
                
                [cDensityMatrix, mDensityMatrix] = ChromosomePositionHistogram.makeDensityMatrix(states, boundIndices, sim);
                totalDensityMatrix = sum(cDensityMatrix, 2) + sum(mDensityMatrix, 2);
                
                save([outDirectory 'DNABoundProteinDensity.mat'], 'totalDensityMatrix')
            end
        end
        
        function displacements = calculateDNABoundProteinDisplacements(outDirectory, simBatchDir, selectedSim, iJob, nJobs)
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            if nargin < 5
                nJobs = 1;
            end
            
            relSimBatchDir = simBatchDir;
            % Make simBatchDir a relative path
            while relSimBatchDir(end) == filesep
                relSimBatchDir(end) = [];
            end
            relSimBatchDir = [relSimBatchDir filesep];
            p = find(relSimBatchDir == filesep, 2, 'last');
            relSimBatchDir = relSimBatchDir(p(1) + 1:p(2) - 1);
            
            outDirectory = [outDirectory filesep];
            
            %% name of file where results will be cached
            if ~exist(sprintf('%sdnaBoundProteinDisplacement', outDirectory), 'dir')
                mkdir(sprintf('%sdnaBoundProteinDisplacement', outDirectory));
            end
            if nJobs == 1 || nargin < 4 || isempty(iJob)
                fileName = sprintf('%sdnaBoundProteinDisplacement%s%s-%d.mat', outDirectory, filesep, relSimBatchDir, selectedSim);
            else
                fileName = sprintf('%sdnaBoundProteinDisplacement%s%s-%d-%d-%d.mat', outDirectory, filesep, relSimBatchDir, selectedSim, iJob, nJobs);
            end
            
            %% get parameters
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep num2str(selectedSim)]);
            c = sim.state('Chromosome');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            %% if data already computed, just load data
            if exist(fileName, 'file')
                load(fileName);
                return;
            end
            
            %% stitch together jobs
            if nargin >= 5 && isempty(iJob)
                displacements = zeros(0, 5);
                for iJob = 1:nJobs
                    tmp = load(sprintf('%sdnaBoundProteinDisplacement-%s-%d-%d-%d.mat', outDirectory, simBatchDir, selectedSim, iJob, nJobs));
                    displacements = [displacements; tmp.displacements]; %#ok<AGROW>
                end
                save(fileName, 'displacements');
                return;
            end
            
            %% get data
            stateNames = {
                'Chromosome' 'monomerBoundSites'
                'Chromosome' 'complexBoundSites'
                };
            displacements = zeros(0, 5);
            
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSim);
            
            monBndSites = permute(max(states.Chromosome.monomerBoundSites(:, 1:2, :), [], 2), [1 3 2]);
            cpxBndSites = permute(max(states.Chromosome.complexBoundSites(:, 1:2, :), [], 2), [1 3 2]);
            
            monMonDisplacement = zeros(numel(pm.matureIndexs), numel(pm.matureIndexs));
            monCpxDisplacement = zeros(numel(pm.matureIndexs), numel(pc.matureIndexs));
            cpxMonDisplacement = zeros(numel(pc.matureIndexs), numel(pm.matureIndexs));
            cpxCpxDisplacement = zeros(numel(pc.matureIndexs), numel(pc.matureIndexs));
            
            [i, j] = find(c.reactionMonomerCatalysisMatrix);
            tfs = c.reactionBoundMonomer(i) ~= 0;
            monMonDisplacement(sub2ind(size(monMonDisplacement), c.reactionBoundMonomer(i(tfs)), j(tfs))) = 1;
            cpxMonDisplacement(sub2ind(size(cpxMonDisplacement), c.reactionBoundComplex(i(~tfs)), j(~tfs))) = 1;
            
            [i, j] = find(c.reactionComplexCatalysisMatrix);
            tfs = c.reactionBoundMonomer(i) ~= 0;
            monCpxDisplacement(sub2ind(size(monCpxDisplacement), c.reactionBoundMonomer(i(tfs)), j(tfs))) = 1;
            cpxCpxDisplacement(sub2ind(size(cpxCpxDisplacement), c.reactionBoundComplex(i(~tfs)), j(~tfs))) = 1;
            
            monDNAFootprints = c.monomerDNAFootprints;
            cpxDNAFootprints = c.complexDNAFootprints;
            cpxDNAFootprints(pc.rnaPolymeraseIndexs) = max(cpxDNAFootprints(pc.rnaPolymeraseIndexs), sim.process('Transcription').rnaPolymeraseElongationRate);
            cpxDNAFootprints(pc.dnaPolymeraseIndexs) = max(cpxDNAFootprints(pc.dnaPolymeraseIndexs), sim.process('Replication').dnaPolymeraseElongationRate);
            
            monDNANegFootprints = zeros(size(c.monomerDNAFootprints));
            cpxDNANegFootprints = zeros(size(c.complexDNAFootprints));
            cpxDNANegFootprints(pc.rnaPolymeraseIndexs) = sim.process('Transcription').rnaPolymeraseElongationRate - c.complexDNAFootprints(pc.rnaPolymeraseIndexs);
            cpxDNANegFootprints(pc.dnaPolymeraseIndexs) = sim.process('Replication').dnaPolymeraseElongationRate - c.complexDNAFootprints(pc.dnaPolymeraseIndexs);
            cpxDNANegFootprints = max(0, cpxDNANegFootprints);
            
            for i = iJob:nJobs:size(states.Chromosome.monomerBoundSites, 3) - 1
                monBndSites1 = monBndSites(:, i);
                monBndSites2 = monBndSites(:, i+1);
                [subs1, vals1] = find(monBndSites1);
                [subs2, vals2] = find(monBndSites2);
                monBndSites1(subs2(monBndSites1(subs2) == vals2, :)) = 0;
                monBndSites2(subs1(monBndSites2(subs1) == vals1, :)) = 0;
                [monPos1, monIdxs1] = find(monBndSites1);
                [monPos2, monIdxs2] = find(monBndSites2);
                monStarts1 = monPos1(:, 1) - monDNANegFootprints(monIdxs1);
                monStarts2 = monPos2(:, 1) - monDNANegFootprints(monIdxs2);
                monEnds1 = monPos1(:, 1) + monDNAFootprints(monIdxs1);
                monEnds2 = monPos2(:, 1) + monDNAFootprints(monIdxs2);
                
                cpxBndSites1 = cpxBndSites(:, i);
                cpxBndSites2 = cpxBndSites(:, i+1);
                [subs1, vals1] = find(cpxBndSites1);
                [subs2, vals2] = find(cpxBndSites2);
                cpxBndSites1(subs2(cpxBndSites1(subs2) == vals2, :)) = 0;
                cpxBndSites2(subs1(cpxBndSites2(subs1) == vals1, :)) = 0;
                [cpxPos1, cpxIdxs1] = find(cpxBndSites1);
                [cpxPos2, cpxIdxs2] = find(cpxBndSites2);
                cpxStarts1 = cpxPos1(:, 1) - cpxDNANegFootprints(cpxIdxs1);
                cpxStarts2 = cpxPos2(:, 1) - cpxDNANegFootprints(cpxIdxs2);
                cpxEnds1 = cpxPos1(:, 1) + cpxDNAFootprints(cpxIdxs1);
                cpxEnds2 = cpxPos2(:, 1) + cpxDNAFootprints(cpxIdxs2);
                
                if any(any(monMonDisplacement(monIdxs1, monIdxs2)))
                    for j = 1:size(monPos1, 1)
                        tfs = ...
                            (monStarts1(j) >= monStarts2 &  monStarts1(j) <= monEnds2) | ...
                            (monEnds1(j) >= monStarts2  &  monEnds1(j) <= monEnds2) | ...
                            (monStarts1(j) <= monStarts2  &  monEnds1(j) >= monEnds2);
                        if any(monMonDisplacement(monIdxs1(j), monIdxs2(tfs)))
                            displacements = [displacements;
                                monIdxs1(j(ones(sum(tfs), 1))) monIdxs2(tfs) monStarts1(j(ones(sum(tfs), 1))) monStarts2(tfs) i(ones(sum(tfs), 1))]; %#ok<AGROW>
                        end
                    end
                end
                
                if any(any(monCpxDisplacement(monIdxs1, cpxIdxs2)))
                    for j = 1:size(monPos1, 1)
                        tfs = ...
                            (monStarts1(j) >= cpxStarts2 &  monStarts1(j) <= cpxEnds2) | ...
                            (monEnds1(j) >= cpxStarts2  &  monEnds1(j) <= cpxEnds2) | ...
                            (monStarts1(j) <= cpxStarts2  &  monEnds1(j) >= cpxEnds2);
                        if any(monCpxDisplacement(monIdxs1(j), cpxIdxs2(tfs)))
                            displacements = [displacements;
                                monIdxs1(j(ones(sum(tfs), 1))) -cpxIdxs2(tfs) monStarts1(j(ones(sum(tfs), 1))) cpxStarts2(tfs) i(ones(sum(tfs), 1))]; %#ok<AGROW>
                        end
                    end
                end
                
                if any(any(cpxMonDisplacement(cpxIdxs1, monIdxs2)))
                    for j = 1:size(monPos2, 1)
                        tfs = ...
                            (cpxStarts1 >= monStarts2(j) &  cpxStarts1 <= monEnds2(j)) | ...
                            (cpxEnds1 >= monStarts2(j)  &  cpxEnds1 <= monEnds2(j)) | ...
                            (cpxStarts1 <= monStarts2(j)  &  cpxEnds1 >= monEnds2(j));
                        if any(cpxMonDisplacement(cpxIdxs1(tfs), monIdxs2(j)))
                            displacements = [displacements;
                                -cpxIdxs1(tfs) monIdxs2(j(ones(sum(tfs), 1))) cpxStarts1(tfs) monStarts2(j(ones(sum(tfs), 1))) i(ones(sum(tfs), 1))]; %#ok<AGROW>
                        end
                    end
                end
                
                if any(any(cpxCpxDisplacement(cpxIdxs1, cpxIdxs2)))
                    for j = 1:size(cpxPos1, 1)
                        tfs = ...
                            (cpxStarts1(j) >= cpxStarts2 &  cpxStarts1(j) <= cpxEnds2) | ...
                            (cpxEnds1(j) >= cpxStarts2  &  cpxEnds1(j) <= cpxEnds2) | ...
                            (cpxStarts1(j) <= cpxStarts2  &  cpxEnds1(j) >= cpxEnds2);
                        if any(cpxCpxDisplacement(cpxIdxs1(j), cpxIdxs2(tfs)))
                            displacements = [displacements;
                                -cpxIdxs1(j(ones(sum(tfs), 1))) -cpxIdxs2(tfs) cpxStarts1(j(ones(sum(tfs), 1))) cpxStarts2(tfs) i(ones(sum(tfs), 1))]; %#ok<AGROW>
                        end
                    end
                end
            end
            
            %% save
            save(fileName, 'displacements');
        end
        
        function metabolicMap(simBatchDir, selectedSim, outDirectory, figHandle, position, showAllLabels)
            import edu.stanford.covert.cell.kb.KnowledgeBaseUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.PrintUtil;
            import edu.stanford.covert.cell.sim.analysis.Figures34;
            import edu.stanford.covert.db.MySQLDatabase;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if ~exist('selectedSim', 'var') || isempty(selectedSim)
                selectedSim = 1;
            end
            if nargin < 6
                showAllLabels = false;
            end
            
            %% get constants
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep num2str(selectedSim)]);
            comp = sim.compartment;
            met = sim.state('Metabolite');
            mr = sim.state('MetabolicReaction');
            metabolism = sim.process('Metabolism');
            
            %% get data
            if exist([outDirectory 'metabolicMap.mat'], 'file')
                load([outDirectory 'metabolicMap.mat'])
            else
                dbConnectionParameters = config();
                database = MySQLDatabase(dbConnectionParameters);
                database.setNullValue(0);
                
                kbWID = KnowledgeBaseUtil.selectLatestKnowledgeBase(database);
                
                database.prepareStatement('CALL get_metabolicmapmetabolites("{Si}", null)', kbWID);
                data = database.query();
                mapMets = struct;
                [~, mapMets.idxs] = ismember(data.Metabolite, met.wholeCellModelIDs);
                [~, mapMets.compIdxs] = ismember(data.Compartment, comp.wholeCellModelIDs);
                [~, mapMets.uidxs] = ismember(mapMets.idxs, unique(mapMets.idxs));
                mapMets.x = data.X;
                mapMets.y = cellfun(@str2double, data.Y);
                
                database.prepareStatement('CALL get_metabolicmapreactions("{Si}", null)', kbWID);
                data = database.query();
                mapRxns = struct;
                [~, mapRxns.idxs] = ismember(data.Reaction, mr.reactionWholeCellModelIDs);
                [~, mapRxns.uidxs] = ismember(mapRxns.idxs, unique(mapRxns.idxs));
                mapRxns.path = data.Path;
                mapRxns.labelX = data.LabelX;
                mapRxns.labelY = data.LabelY;
                mapRxns.valueX = data.ValueX;
                mapRxns.valueY = data.ValueY;
                
                %% get data
                stateNames = {
                    'MetabolicReaction'  'fluxs'  ':'  ':'
                    };
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSim);
                
                states.MetabolicReaction.fluxs = full(states.MetabolicReaction.fluxs);
                
                %% table
                colLabels = {'Reaction ID', 'Reaction Name', 'Flux Mean (rxn/s)', 'Flux Std (rxn/s)'};
                content = [
                    metabolism.reactionWholeCellModelIDs  ...
                    metabolism.reactionNames  ...
                    num2cell(nanmean(states.MetabolicReaction.fluxs, 3)) ...
                    num2cell(nanstd(states.MetabolicReaction.fluxs, [], 3))
                    ];
                PrintUtil.printToFile(content, colLabels, [outDirectory 'reactionFluxes.xls'], 'Reaction Fluxes');
                
                %% iso fluxs
                isoFluxIDs = {
                    {'Adk1', 'Adk2', 'Adk3', 'Adk4'}
                    {'PfkA1', 'PfkA2', 'PfkA3', 'PfkA4', 'PfkA5'}
                    {'Pgk1', 'Pgk2', 'Pgk3', 'Pgk4', 'Pgk5'}
                    {'Pyk_ADP', 'Pyk_CDP', 'Pyk_DADP', 'Pyk_DCDP', 'Pyk_DGDP', 'Pyk_DTDP', 'Pyk_DUDP', 'Pyk_GDP', 'Pyk_IDP', 'Pyk_TDP', 'Pyk_UDP'}
                    };
                isoFluxs = false(size(mr.reactionWholeCellModelIDs));
                for i = 1:numel(isoFluxIDs)
                    tfs2 = ismember(mr.reactionWholeCellModelIDs, isoFluxIDs{i}(2:end));
                    isoFluxs = isoFluxs | tfs2;
                    
                    [~, idxs1] = ismember(isoFluxIDs{i}(1), mr.reactionWholeCellModelIDs);
                    [~, idxs2] = ismember(isoFluxIDs{i}(2:end), mr.reactionWholeCellModelIDs);
                    states.MetabolicReaction.fluxs(idxs1, :, :) = ...
                        + states.MetabolicReaction.fluxs(idxs1, :, :) ...
                        + sum(states.MetabolicReaction.fluxs(idxs2, :, :), 1);
                end
                
                fluxMeans = nanmean(states.MetabolicReaction.fluxs(mapRxns.idxs, :, :), 3);
                fluxStds = nanstd(states.MetabolicReaction.fluxs(mapRxns.idxs, :, :), [], 3) ./ abs(fluxMeans); % (Standard Deviation)/|Mean|
                
                %% save
                save([outDirectory 'metabolicMap.mat'], 'mapMets', 'mapRxns', 'fluxMeans', 'fluxStds', 'isoFluxs');
            end
            
            %% plot map
            H = 600;
            W = 800;
            metRadius = 3;
            arrowRadius = 14;
            barW = 20;
            reactionLims = [0 2 4];
            strokes = [0.25 1.125 1];
            
            %align
            mapRxns.path{mapRxns.idxs == metabolism.reactionIndexs('PfkA1')} = 'M 230,145 L 230,178';
            mapRxns.path{mapRxns.idxs == metabolism.reactionIndexs('Pgk1')} = 'M 270,265 L 270,293';
            mapRxns.path{mapRxns.idxs == metabolism.reactionIndexs('GpmA')} = 'M 270,305 L 270,331';
            mapRxns.path{mapRxns.idxs == metabolism.reactionIndexs('Pyk_ADP')} = 'M 270,385 L 270,418';
            mapRxns.path{mapRxns.idxs == metabolism.reactionIndexs('Adk1')} = 'M 761,80 L 815,80';
            
            %scale
            minX = min(mapMets.x)-metRadius;
            minY = min(mapMets.y)+metRadius;
            scaleX = W / (range(mapMets.x) + 2*metRadius);
            scaleY = H / (range(mapMets.y) + 2*metRadius);
            mapMets.x = scaleX * (mapMets.x - minX);
            mapMets.y = scaleY * (mapMets.y - minY);
            for i = 1:numel(mapRxns.path)
                [starts, ends, tokens] = regexp(mapRxns.path{i}, '(-*\d*\.*\d*),(-*\d*\.*\d*)', 'start', 'end', 'tokens');
                for j = 1:numel(starts)
                    valX = num2str(scaleX * (str2double(tokens{j}{1}) - minX));
                    valY = num2str(scaleY * (str2double(tokens{j}{2}) - minY));
                    mapRxns.path{i} = [...
                        mapRxns.path{i}(1:starts(j)-1) ...
                        valX ...
                        ',' ...
                        valY ...
                        mapRxns.path{i}(ends(j)+1:end) ...
                        ];
                    starts = starts + numel(valX) + numel(valY) - numel(tokens{j}{1}) - numel(tokens{j}{2});
                    ends = ends + numel(valX) + numel(valY) - numel(tokens{j}{1}) - numel(tokens{j}{2});
                end
            end
            mapRxns.labelX = scaleX * (mapRxns.labelX - minX);
            mapRxns.labelY = scaleY * (mapRxns.labelY - minY);
            mapRxns.valueX = scaleX * (mapRxns.valueX - minX);
            mapRxns.valueY = scaleY * (mapRxns.valueY - minY);
            
            if nargin >= 5
                axesHandle = subplot('position', position, 'parent', figHandle);
            else
                axesHandle = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            end
            hold(axesHandle, 'on');
            set(axesHandle, 'YDir', 'reverse');
            xlim(axesHandle, [0 W+10]);
            ylim(axesHandle, [-8 H-10]);
            set(axesHandle, 'visible', 'off')
            
            %metabolites
            tfs = true(size(mapMets.x));
            tfs([59:63 14 39 85 53 67 18 69 35 40 33 44 37 70 83]) = false;
            plot(mapMets.x(tfs), mapMets.y(tfs), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 0.2 * [1 1 1], 'MarkerSize', metRadius);
            if showAllLabels
                for i = 1:numel(mapMets.x)
                    if ~tfs(i)
                        continue;
                    end
                    text(mapMets.x(i)+7, mapMets.y(i)+7, met.wholeCellModelIDs{mapMets.idxs(i)}, ...
                        'FontSize', 5, 'Color', 'k', 'Interpreter', 'none');
                end
            end
            
            %reactions
            relMean = log10(abs(fluxMeans));
            relMean(isinf(relMean)) = NaN;
            relMean = max(0, min(4, relMean)) / 4;
            strokeWidth = 0.25 * ones(size(relMean));
            color1 = 0.2 * [0 0 1] + 0.8 * (0.9 * [1 1 1]);
            color2 = [0 0 1];
            colors = ...
                + repmat(color1, numel(relMean), 1) .* repmat(1-relMean, 1, 3) ...
                + repmat(color2, numel(relMean), 1) .* repmat(relMean, 1, 3);
            
            [~, order] = sort(relMean);
            relMean = relMean(order, :);
            fluxMeans = fluxMeans(order, :);
            colors = colors(order, :);
            strokeWidth = strokeWidth(order, :);
            mapRxns.idxs = mapRxns.idxs(order);
            mapRxns.path = mapRxns.path(order);
            mapRxns.labelX = mapRxns.labelX(order);
            mapRxns.labelY = mapRxns.labelY(order);
            mapRxns.valueX = mapRxns.valueX(order);
            mapRxns.valueY = mapRxns.valueY(order);
            
            arrowRadii = arrowRadius * ones(size(mapRxns.idxs));
            arrowRadii(mapRxns.idxs == metabolism.reactionIndexs('AspC1')) = 0.6 * arrowRadius;
            
            [~, hideFluxIdxs] = ismember({'GlpQ1', 'GlpQ2', 'GlpQ3', 'GlpQ4', 'GlpQ5', 'DeoD', 'Rpe', 'RpiB'}, mr.reactionWholeCellModelIDs);
            for i = 1:numel(mapRxns.idxs)
                if isoFluxs(mapRxns.idxs(i)) || any(mapRxns.idxs(i) == hideFluxIdxs) || fluxMeans(i) == 0
                    continue;
                end
                
                endX = 0;
                endY = 0;
                endAngle = NaN;
                tokens = regexp(mapRxns.path{i}, '([A-Za-z])( (-*\d*\.*\d*),(-*\d*\.*\d*))+', 'tokens');
                for j = 1:numel(tokens)
                    tokens2 = regexp(tokens{j}{2}, '(-*\d*\.*\d*),(-*\d*\.*\d*)', 'tokens');
                    switch tokens{j}{1}
                        case 'M'
                            endX = str2double(tokens2{1}{1});
                            endY = str2double(tokens2{1}{2});
                        case 'L'
                            startX = endX;
                            startY = endY;
                            endX2 = str2double(tokens2{1}{1});
                            endY2 = str2double(tokens2{1}{2});
                            line([endX endX2], [endY endY2], ...
                                'Parent', axesHandle, ...
                                'Color', colors(i, :) * (~showAllLabels), ...
                                'lineWidth', strokeWidth(i));
                            startAngle = atan2(endY-endY2, endX-endX2);
                            endAngle = atan2(endY2 - endY, endX2 - endX);
                            endX = endX2;
                            endY = endY2;
                        case 'C'
                            startX = endX;
                            startY = endY;
                            endX2 = str2double(tokens2{1}{1});
                            endY2 = str2double(tokens2{1}{2});
                            endX3 = str2double(tokens2{2}{1});
                            endY3 = str2double(tokens2{2}{2});
                            endX4 = str2double(tokens2{3}{1});
                            endY4 = str2double(tokens2{3}{2});
                            [h, qx, qy] = Funct_Bezier([endX endX2 endX3 endX4], [endY endY2 endY3 endY4], 50, axesHandle);
                            set(h, 'color', colors(i, :) * (~showAllLabels), 'lineWidth', strokeWidth(i));
                            startAngle = atan2(qy(1)-qy(2), qx(1)-qx(2));
                            endAngle = atan2(qy(end)-qy(end-1), qx(end)-qx(end-1));
                            endX = endX4;
                            endY = endY4;
                        otherwise
                            throw(MException('Figures34:error', 'undefined command %s', tokens{j}{1}))
                    end
                    if j == 2 && fluxMeans(i) < 0
                        x = startX + 2*cos(startAngle) + [0  arrowRadii(i)*cos(startAngle+pi/3+pi/2)  arrowRadii(i)*cos(startAngle+pi*2/3+pi/2)];
                        y = startY + 2*sin(startAngle) + [0  arrowRadii(i)*sin(startAngle+pi/3+pi/2)  arrowRadii(i)*sin(startAngle+pi*2/3+pi/2)];
                        
                        arrowPtch = patch(x, y, 1, 'Parent', axesHandle, 'EdgeColor', 'none');
                        set(arrowPtch, 'FaceColor', colors(i, :) * (~showAllLabels), 'FaceAlpha', 1, 'EdgeAlpha', 1);
                    end
                end
                if ~isnan(endAngle) && fluxMeans(i) > 0
                    x = endX + 2*cos(endAngle) + [0  arrowRadii(i)*cos(endAngle+pi/3+pi/2)  arrowRadii(i)*cos(endAngle+pi*2/3+pi/2)];
                    y = endY + 2*sin(endAngle) + [0  arrowRadii(i)*sin(endAngle+pi/3+pi/2)  arrowRadii(i)*sin(endAngle+pi*2/3+pi/2)];
                    
                    arrowPtch = patch(x, y, 1, 'Parent', axesHandle, 'EdgeColor', 'none');
                    set(arrowPtch, 'FaceColor', colors(i, :) * (~showAllLabels), 'FaceAlpha', 1, 'EdgeAlpha', 1);
                end
                
                if showAllLabels
                    text(mapRxns.labelX(i), mapRxns.labelY(i), ...
                        metabolism.reactionWholeCellModelIDs(mapRxns.idxs(i)), ...
                        'FontSize', 5, 'Interpreter', 'none');
                end
            end
            
            %labels
            labels = {
                [90 90 210 210]            [60 400 400 60]          [0 1 0]  [105 335] 'Glycolysis'                 0.2  0
                [230 230 370 370]          [55 230 230 55]          [1 0 0]  [205 190] {'Pentose' 'phosphate'}      0.2  90
                [265 265 425 425]          [230 340 340 230]        [0 0 1]  [375 310] {'ATP' 'synthesis'}          0.2  0
                [450 680 680 820 820 450]  [20 20 320 320 580 580]  [0 1 1]  [635 305] 'Nucleotide metabolism'     0.1  0
                [170 390 390 170]          [420 420 810 810]        [1 0.5 0]  [280 470] {'Pyruvate' 'metabolism'} 0.3  0
                };
            h = zeros(size(labels, 1), 1);
            for i = 1:size(labels, 1)
                if showAllLabels
                    ptch = patch(labels{i, 1}, labels{i, 2}, 1, 'Parent', axesHandle, 'EdgeColor', 'none');
                    set(ptch, 'FaceColor', labels{i,6}*labels{i, 3} + (1-labels{i,6})*[1 1 1], 'FaceAlpha', 1, 'EdgeAlpha', 1);
                    uistack(ptch, 'bottom');
                    
                    h(i) = text(labels{i, 4}(1), labels{i, 4}(2), labels{i, 5}, ...
                        'HorizontalAlign', 'center',  ...
                        'VerticalAlign', 'middle',  ...
                        'FontSize', 5, ...
                        'Color', labels{i, 3});
                else
                    h(i) = text(labels{i, 4}(1), labels{i, 4}(2), labels{i, 5}, ...
                        'HorizontalAlign', 'center',  ...
                        'VerticalAlign', 'middle',  ...
                        'FontSize', 5, ...
                        'Color', 'k', 'rotation', labels{i, 7});
                end
            end
            set(h(2:3), 'HorizontalAlign', 'left')
            
            %colormap
            if ~showAllLabels
                x1 = W-100;
                x2 = x1 + 20;
                h = 6.5;
                y = 0 + h*40;
                colorMap = ...
                    + repmat(color1, 40, 1) .* repmat(1 - linspace(0, 1, 40)', 1, 3) ...
                    + repmat(color2, 40, 1) .* repmat(    linspace(0, 1, 40)', 1, 3);
                for i = 1:40
                    ptch = patch([x1 x2 x2 x1], y - h*[i-1 i-1 i i], 1, 'EdgeColor', 'none', 'Parent', axesHandle);
                    set(ptch, 'FaceColor', colorMap(i, :), 'FaceAlpha', 1, 'EdgeAlpha', 1);
                end
                patch([x1 x2 x2 x1], y - [0 0 40*h 40*h], 1, 'EdgeColor', 'k', 'FaceColor', 'none');
                
                line([x2 x2+5], y + [0 0], 'Color', 'k', 'parent', axesHandle)
                line([x2 x2+5], y-20*h + [0 0], 'Color', 'k', 'parent', axesHandle)
                line([x2 x2+5], y-40*h + [0 0], 'Color', 'k', 'parent', axesHandle)
                
                text(x2+75, y-20*h, 'Flux (rxn s^{-1})', ...
                    'Color', 'k', 'parent', axesHandle, 'Interpreter', 'tex', ...
                    'HorizontalAlign', 'center', 'VerticalAlign', 'middle', 'FontSize', 7, 'rotation', 270);
                text(x2+7, y, sprintf('10^{%d}', 0), ...
                    'Color', 'k', 'parent', axesHandle, 'Interpreter', 'tex', ...
                    'HorizontalAlign', 'left', 'VerticalAlign', 'middle', 'FontSize', 5);
                text(x2+7, y-20*h, sprintf('10^{%d}', 2), ...
                    'Color', 'k', 'parent', axesHandle, 'Interpreter', 'tex', ...
                    'HorizontalAlign', 'left', 'VerticalAlign', 'middle', 'FontSize', 5);
                text(x2+7, y-40*h, sprintf('10^{%d}', 4), ...
                    'Color', 'k', 'parent', axesHandle, 'Interpreter', 'tex', ...
                    'HorizontalAlign', 'left', 'VerticalAlign', 'middle', 'FontSize', 5);
            end
            
            %% labels
            if ~showAllLabels
                idxGpsA = find(mapRxns.idxs == sim.process('Metabolism').reactionIndexs('GpsA'), 1, 'first');
                idxTalA = find(mapRxns.idxs == sim.process('Metabolism').reactionIndexs('TalA'), 1, 'first');
                idxPgi = find(mapRxns.idxs == sim.process('Metabolism').reactionIndexs('Pgi'), 1, 'first');
                
                ptch = patch([335 335 450 450], [90 190 190 90], 1, 'Parent', axesHandle, 'EdgeColor', 'none');
                set(ptch, 'FaceColor', 0.2*[1 0 0] + (1-0.2)*[1 1 1], 'FaceAlpha', 1, 'EdgeAlpha', 1);
                uistack(ptch, 'bottom');
                
                ptch = patch([0 0 107 107], [190 285 285 190], 1, 'Parent', axesHandle, 'EdgeColor', 'none');
                set(ptch, 'FaceColor', 0.2*[1 0 0] + (1-0.2)*[1 1 1], 'FaceAlpha', 1, 'EdgeAlpha', 1);
                uistack(ptch, 'bottom');
                
                text(95, 220, 'GpsA', ...
                    'FontSize', 7, ...
                    'HorizontalAlign', 'right')
                text(95, 260, sprintf('(%.2f%%)', abs(fluxMeans(idxGpsA)/fluxMeans(idxPgi)*100)), ...
                    'FontSize', 5, ...
                    'HorizontalAlign', 'right')
                text(390, 120, 'TalA', ...
                    'FontSize', 7, ...
                    'HorizontalAlign', 'center')
                text(390, 160, sprintf('(%.2f%%)', abs(fluxMeans(idxTalA)/fluxMeans(idxPgi)*100)), ...
                    'FontSize', 5, ...
                    'HorizontalAlign', 'center')
            end
        end
        
        function chromosomeSampling(simBatchDir, selectedSim, outDirectory, figHandle, position)
            import edu.stanford.covert.cell.sim.analysis.Figures34;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if ~exist('selectedSim', 'var') || isempty(selectedSim)
                selectedSim = 1;
            end
            
            %% get data
            if exist([outDirectory 'chromosomeSampling.mat'], 'file')
                load([outDirectory 'chromosomeSampling.mat']);
            else
                [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep num2str(selectedSim)]);
                ch = sim.state('Chromosome');
                pc = sim.state('ProteinComplex');
                rnaPol = sim.state('RNAPolymerase');
                
                stateNames = {
                    'Time'        'values'             ':'  ':'
                    'Chromosome'  'monomerBoundSites'  ':'  ':'
                    'Chromosome'  'complexBoundSites'  ':'  ':'
                    'RNAPolymerase' 'positionStrands'  ':'  ':'
                    'RNAPolymerase' 'states'           ':'  ':'
                    };
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSim);
                time = permute(states.Time.values, [3 1 2]) / 3600;
                
                %monomers
                [subs, vals] = find(states.Chromosome.monomerBoundSites);
                times = subs(:, 3);
                startPos = subs(:, 1);
                endPos = subs(:, 1) + ch.monomerDNAFootprints(vals) - 1;
                idxs = find(endPos > ch.sequenceLen);
                vals = [vals; vals(idxs)];
                times = [times; times(idxs)];
                startPos = [startPos; ones(size(idxs))];
                endPos = [endPos; endPos(idxs) - ch.sequenceLen];
                endPos(idxs) = ch.sequenceLen;
                monTimeBound = Figures34.calcChromosomeSampling(startPos, endPos, times, 1, ch.sequenceLen, ch.sequenceLen);
                
                %complexes
                complexDNAFootprints = ch.complexDNAFootprints;
                complexDNAFootprints(pc.dnaPolymeraseIndexs) = max(complexDNAFootprints(pc.dnaPolymeraseIndexs), sim.process('Replication').dnaPolymeraseElongationRate);
                complexDNAFootprints(pc.getIndexs('MG_094_HEXAMER')) = max(complexDNAFootprints(pc.getIndexs('MG_094_HEXAMER')), sim.process('Replication').dnaPolymeraseElongationRate);
                
                [subs, vals] = find(states.Chromosome.complexBoundSites);
                subs(ismember(vals, pc.rnaPolymeraseIndexs), :) = [];
                vals(ismember(vals, pc.rnaPolymeraseIndexs), :) = [];
                times = subs(:, 3);
                startPos = subs(:, 1);
                endPos = subs(:, 1) + complexDNAFootprints(vals) - 1;
                
                tmp = [reshape(states.RNAPolymerase.positionStrands(:, 1, :), [], 1)  ...
                    reshape(repmat(permute(1:size(states.RNAPolymerase.positionStrands, 3), [1 3 2]), [size(states.RNAPolymerase.positionStrands, 1), 1 1]), [], 1) ...
                    reshape(states.RNAPolymerase.states(:, 1, :), [], 1)];
                tmp(tmp(:, 3) == rnaPol.freeValue | tmp(:, 3) == rnaPol.notExistValue, :) = [];
                times = [times; tmp(:, 2)];
                vals = [vals; pc.rnaPolymeraseIndexs(ones(size(tmp, 1), 1))];
                startPos = [startPos; tmp(:, 1)];
                endPos = [endPos;
                    tmp(:, 1) ...
                    + ch.complexDNAFootprints(pc.rnaPolymeraseIndexs(1)) * (tmp(:, 1) < rnaPol.activelyTranscribingValue) ...
                    + sim.process('Transcription').rnaPolymeraseElongationRate * (tmp(:, 1) >= rnaPol.activelyTranscribingValue)
                    ];
                
                idxs = find(endPos > ch.sequenceLen);
                vals = [vals; vals(idxs)];
                times = [times; times(idxs)];
                startPos = [startPos; ones(size(idxs))];
                endPos = [endPos; endPos(idxs) - ch.sequenceLen];
                endPos(idxs) = ch.sequenceLen;
                cpxTimeBound = Figures34.calcChromosomeSampling(startPos, endPos, times, 1, ch.sequenceLen, ch.sequenceLen);
                
                %fraction sampled
                timeBound = min(monTimeBound, cpxTimeBound);
                fractionSampled = cumsum(histc(timeBound, (1:size(states.Time.values, 3))')) / ch.sequenceLen;
                
                %fraction sampled - individual proteins
                cpxIdxs = pc.getIndexs({'MG_213_214_298_6MER_ADP', 'RNA_POLYMERASE', 'MG_203_204_TETRAMER', 'MG_094_HEXAMER'});
                fractionSampled_IndivProtein = zeros(size(states.Time.values, 3), numel(cpxIdxs));
                for i = 1:numel(cpxIdxs)
                    timeBound = Figures34.calcChromosomeSampling(startPos(vals == cpxIdxs(i)), endPos(vals == cpxIdxs(i)), times(vals == cpxIdxs(i)), 1, ch.sequenceLen, ch.sequenceLen);
                    fractionSampled_IndivProtein(:, i) = cumsum(histc(timeBound, (1:size(states.Time.values, 3))')) / ch.sequenceLen;
                end
                
                %fraction sampled - RNA Pol
                nTime = size(states.RNAPolymerase.positionStrands, 3);
                nPols = size(states.RNAPolymerase.positionStrands, 1);
                rnaPolRate = sim.process('Transcription').rnaPolymeraseElongationRate;
                rnaPolFtpt = ch.complexDNAFootprints(pc.rnaPolymeraseIndexs(1));
                fractionSampled_RNAPol = zeros(nTime, 3);
                
                tmpPosTime = [reshape(permute(states.RNAPolymerase.positionStrands(:, 1, :), [1 3 2]), [], 1)  reshape(repmat(1:nTime, nPols, 1), [], 1)];
                tmpStates = reshape(permute(states.RNAPolymerase.states, [1 3 2]), [], 1);
                
                tfs = tmpStates >= rnaPol.activelyTranscribingValue;
                timeBound = Figures34.calcChromosomeSampling(tmpPosTime(tfs, 1), tmpPosTime(tfs, 1)+rnaPolRate, tmpPosTime(tfs, 2), 1, ch.sequenceLen, ch.sequenceLen);
                fractionSampled_RNAPol(:, 1) = cumsum(histc(timeBound, (1:nTime)')) / ch.sequenceLen;
                
                tfs = false(ch.sequenceLen, 1);
                for j = 1:numel(rnaPol.transcripts.transcriptionUnitFivePrimeCoordinates)
                    if rnaPol.transcripts.transcriptionUnitDirections(j) == 1
                        tfs(rnaPol.transcripts.transcriptionUnitFivePrimeCoordinates(j)+(1:rnaPol.transcripts.transcriptionUnitLengths(j))-1) = true;
                    else
                        tfs(rnaPol.transcripts.transcriptionUnitFivePrimeCoordinates(j)-(1:rnaPol.transcripts.transcriptionUnitLengths(j))+1) = true;
                    end
                end
                timeBound = timeBound(tfs);
                fractionSampled_RNAPol(:, 3) = cumsum(histc(timeBound, (1:nTime)')) / numel(timeBound);
                
                tfs = tmpStates == rnaPol.nonSpecificallyBoundValue | tmpStates == rnaPol.specificallyBoundValue;
                timeBound = Figures34.calcChromosomeSampling(tmpPosTime(tfs, 1), tmpPosTime(tfs, 1)+rnaPolFtpt, tmpPosTime(tfs, 2), 1, ch.sequenceLen, ch.sequenceLen);
                fractionSampled_RNAPol(:, 2) = cumsum(histc(timeBound, (1:nTime)')) / ch.sequenceLen;
                
                %fit
                expFittype = fittype('1 - a * exp(-log(2) * t / tau)', ...
                    'independent', 't', ...
                    'coefficients', {'a', 'tau'});
                expFitOptions = fitoptions(...
                    'Method', 'NonlinearLeastSquares', ...
                    'Startpoint', [.25 2], ...
                    'Lower', [0 0], ...
                    'Upper', [1 10]);
                
                [fitResult, gof] = fit(time, fractionSampled, expFittype, expFitOptions);
                
                tau = zeros(size(fractionSampled_IndivProtein, 2), 1);
                for i = 1:size(fractionSampled_IndivProtein, 2)
                    tmp = fit(time, fractionSampled_IndivProtein(:, i), expFittype, expFitOptions);
                    tau(i) = tmp.tau;
                end
                
                tau_RNAPol = zeros(size(fractionSampled_RNAPol, 2), 1);
                for i = 1:size(fractionSampled_RNAPol, 2)
                    tmp = fit(time, fractionSampled_RNAPol(:, i), expFittype, expFitOptions);
                    tau_RNAPol(i) = tmp.tau;
                end
                
                save([outDirectory 'chromosomeSampling.mat'], 'time', 'fractionSampled', ...
                    'fractionSampled_IndivProtein', 'fractionSampled_RNAPol', ...
                    'fitResult', 'gof', 'tau', 'tau_RNAPol');
            end
            
            %% plot
            if nargin >= 5
                axesHandle = subplot('position', position, 'parent', figHandle);
            else
                axesHandle = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            end
            
            labels = {'All', 'DNA pol', 'Gyrase', 'RNA pol', 'SMC'};
            colors = [
                0 1 1
                0 1 0
                1 0 0
                0 0 1
                1 0 1
                ];
            data = [fractionSampled  fractionSampled_IndivProtein(:, end:-1:1)];
            tau = [fitResult.tau; tau(end:-1:1, 1)] * 60 * log(2);
            
            show = true(size(labels));
            %show([3 5]) = false;
            labels = labels(show);
            colors = colors(show, :);
            data = data(:, show);
            tau = tau(show, :);
            
            h = plot(axesHandle, time, 100 * data, 'LineWidth', 1);
            
            [~, order] = sort(tau);
            labels = labels(order);
            tau = tau(order);
            h = h(order);
            colors = colors(order, :);
            data = data(:, order);
            
            tableFont = 6;
            rowSep = 12;
            text(7.4, 88, 'Protein', 'FontSize', tableFont, 'HorizontalAlign', 'left', 'VerticalAlign', 'middle', 'Parent', axesHandle);
            line(7.4 + [0 1.4], 82 + [0 0], 'Parent', axesHandle, 'Color', 'k')
            for i = 1:numel(labels)
                set(h(i), 'Color', colors(i, :));
                line(6.9 + [0 0.3], 87-i*rowSep + [0 0], 'Parent', axesHandle, 'Color', colors(i, :), 'LineWidth', 2)
                text(7.4, 87-i*rowSep, labels{i}, 'FontSize', tableFont, 'HorizontalAlign', 'left', 'VerticalAlign', 'middle', 'Parent', axesHandle);
            end
            
            xlim(axesHandle, [0 time(end)]);
            ylim(axesHandle, [0 100]);
            
            set(axesHandle, 'FontSize', 5, 'TickDir', 'out', 'XTick', 0:2:8, 'YTick', 0:25:100, 'box', 'off');
            
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', 0.08, 'ytickoffset', 0.025);
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            set(xTicks, 'FontSize', 5);
            set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'right');
            delete(xTicks(2:4));
            delete(yTicks(2:4));
            
            xlabel(axesHandle, 'Time (h)', 'FontSize', 7);
            ylabel(axesHandle, {'% Chromosome' 'explored'}, 'FontSize', 7);
            xlabelPos = get(get(axesHandle, 'xlabel'), 'position');
            ylabelPos = get(get(axesHandle, 'ylabel'), 'position');
            set(get(axesHandle, 'xlabel'), 'position', [xlabelPos(1) -11 xlabelPos(3:end)]);
            set(get(axesHandle, 'ylabel'), 'position', [-0.55 ylabelPos(2:end)]);
        end
        
        function plotRNAPolSampling(simBatchDir, selectedSim, outDirectory, figHandle, position)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            %% lood constants
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep num2str(selectedSim)]);
            ch = sim.state('Chromosome');
            pc = sim.state('ProteinComplex');
            rnaPol = sim.state('RNAPolymerase');
            
            %% load data
            load([outDirectory 'chromosomeSampling.mat']);
            load([outDirectory 'rnaSampling.mat']);
            
            %% plot
            if nargin >= 5
                axesHandle = subplot('position', position, 'parent', figHandle);
            else
                axesHandle = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            end
            
            h = plot(axesHandle, time, 100 * [fractionSampled_IndivProtein(:, 2) fractionSampled_RNAPol(:, [1 3 2]) fracRNASampled]); %#ok<NODEF>
            legend(h, {
                sprintf('All RNA Polymerase (\\tau = %0.1f h)', tau(2))
                sprintf('Active/SB - of whole Chromosome (\\tau = %0.1f h)', tau_RNAPol(1))
                sprintf('Active/SB - of TUs (\\tau = %0.1f h)', tau_RNAPol(3))
                sprintf('Non-SB (\\tau = %0.1f h)', tau_RNAPol(2))
                sprintf('Expressed RNA (\\tau = %0.1f h)', rnaFitResult.tau)
                }, ...
                'Location', 'SouthEast', 'FontSize', 8);
            xlabel(axesHandle, 'Time (h)')
            ylabel(axesHandle, '% Explored / expressed')
            xlim(axesHandle, [0 max(time)])
            ylim(axesHandle, [0 100])
        end
        
        function timeBound = calcChromosomeSampling(startPos, endPos, times, pos1, pos2, L)
            import edu.stanford.covert.cell.sim.analysis.Figures34;
            
            if pos2-pos1 > 500
                pos21 = floor((pos1+pos2)/2);
                pos12 = pos21+1;
                tfs1 = startPos <= pos21;
                tfs2 = endPos >= pos12;
                timeBound = min(...
                    Figures34.calcChromosomeSampling(startPos(tfs1), endPos(tfs1), times(tfs1), pos1, pos21, L), ...
                    Figures34.calcChromosomeSampling(startPos(tfs2), endPos(tfs2), times(tfs2), pos12, pos2, L));
            else
                timeBound = NaN(L, 1);
                for i = pos1:pos2
                    tf = startPos <= i & endPos >= i;
                    if ~any(tf)
                        continue;
                    end
                    timeBound(i) = min(times(tf));
                end
            end
        end
        
        function rnaSampling(simBatchDir, selectedSim, outDirectory, figHandle, position)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if ~exist('selectedSim', 'var') || isempty(selectedSim)
                selectedSim = 1;
            end
            
            %% data
            if exist([outDirectory 'rnaSampling.mat'], 'file')
                load([outDirectory 'rnaSampling.mat']);
            else
                [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep num2str(selectedSim)]);
                comp = sim.compartment;
                rna = sim.state('Rna');
                
                stateNames = {
                    'Time'       'values'        ':'  ':'
                    'Rna'        'counts'        rna.matureIndexs  comp.cytosolIndexs
                    };
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSim);
                time = permute(states.Time.values, [3 1 2]) / 3600;
                
                [i, j] = find(permute(states.Rna.counts, [1 3 2]));
                timeFirstExpressed = NaN(size(rna.matureIndexs));
                timeFirstExpressed(i(end:-1:1)) = j(end:-1:1);
                fracRNASampled = cumsum(histc(timeFirstExpressed, (1:numel(time)))) / numel(rna.matureIndexs);
                
                expFittype = fittype('1 - a * exp(-log(2) * t / tau)', ...
                    'independent', 't', ...
                    'coefficients', {'a', 'tau'});
                expFitOptions = fitoptions(...
                    'Method', 'NonlinearLeastSquares', ...
                    'Startpoint', [.25 2], ...
                    'Lower', [0 0], ...
                    'Upper', [1 10]);
                [rnaFitResult, rnaGof] = fit(time, fracRNASampled, expFittype, expFitOptions);
                
                save([outDirectory 'rnaSampling.mat'], 'time', 'fracRNASampled', 'rnaFitResult', 'rnaGof')
            end
            
            %% plot
            if nargin >= 5
                axesHandle = subplot('position', position, 'parent', figHandle);
            else
                axesHandle = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            end
            
            plot(axesHandle, time, 100 * fracRNASampled, 'Color', 'b', 'LineWidth', 1);
            
            xlim(axesHandle, [0 time(end)]);
            ylim(axesHandle, [0 100]);
            
            set(axesHandle, 'FontSize', 5, 'TickDir', 'out', 'XTick', 0:2:8, 'YTick', 0:25:100, 'box', 'off');
            
            t_50 = min(time(fracRNASampled >= 0.50)) * 60;
            text(8.84, 98, sprintf('t_{50} = %.0f min', t_50), ...
                'Parent', axesHandle, ...
                'FontSize', 6, ...
                'horizontalalign', 'right', ...
                'verticalalign', 'top');
            
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', 0.08, 'ytickoffset', 0.025);
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            set(xTicks, 'FontSize', 5);
            set(yTicks, 'FontSize', 5, 'HorizontalAlign', 'right');
            delete(xTicks(2:4));
            delete(yTicks(2:4));
            
            xlabel(axesHandle, 'Time (h)', 'FontSize', 7);
            ylabel(axesHandle, {'% RNA' 'expressed'}, 'FontSize', 7);
            xlabelPos = get(get(axesHandle, 'xlabel'), 'position');
            ylabelPos = get(get(axesHandle, 'ylabel'), 'position');
            set(get(axesHandle, 'xlabel'), 'position', [xlabelPos(1) -11 xlabelPos(3:end)]);
            set(get(axesHandle, 'ylabel'), 'position', [-0.55 ylabelPos(2:end)]);
        end
        
        function mrnaSampling(simBatchDir, selectedSim, outDirectory, figHandle, position)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if ~exist('selectedSim', 'var') || isempty(selectedSim)
                selectedSim = 1;
            end
            
            %% data
            if exist([outDirectory 'mrnaSampling.mat'], 'file')
                load([outDirectory 'mrnaSampling.mat']);
            else
                [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep num2str(selectedSim)]);
                comp = sim.compartment;
                rna = sim.state('Rna');
                
                stateNames = {
                    'Time'       'values'        ':'  ':'
                    'Rna'        'counts'        rna.matureIndexs  comp.cytosolIndexs
                    };
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSim);
                time = permute(states.Time.values, [3 1 2]) / 3600;
                
                [i, j] = find(permute(states.Rna.counts, [1 3 2]));
                timeFirstExpressed = NaN(size(rna.matureIndexs));
                timeFirstExpressed(i(end:-1:1)) = j(end:-1:1);
                fracMRNASampled = cumsum(histc(timeFirstExpressed(rna.matureMRNAIndexs), (1:numel(time)))) / numel(rna.matureMRNAIndexs);
                
                expFittype = fittype('1 - a * exp(-log(2) * t / tau)', ...
                    'independent', 't', ...
                    'coefficients', {'a', 'tau'});
                expFitOptions = fitoptions(...
                    'Method', 'NonlinearLeastSquares', ...
                    'Startpoint', [.25 2], ...
                    'Lower', [0 0], ...
                    'Upper', [1 10]);
                [mrnaFitResult, mrnaGof] = fit(time, fracMRNASampled, expFittype, expFitOptions);
                
                save([outDirectory 'mrnaSampling.mat'], 'time', 'fracMRNASampled', 'mrnaFitResult', 'mrnaGof')
            end
            
            %% plot mRNA
            if nargin >= 5
                axesHandle = subplot('position', position, 'parent', figHandle);
            else
                axesHandle = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            end
            
            h = plot(axesHandle, time, 100 * [fracMRNASampled   1-mrnaFitResult.a*exp(-log(2)*time/mrnaFitResult.tau)]);
            set(h(1), 'Color', 'b');
            set(h(2), 'Color', 'k');
            uistack(h(1), 'top');
            
            xlim(axesHandle, time([1 end]));
            ylim(axesHandle, [0 100]);
            
            xlabel(axesHandle, 'Time (h)', 'FontSize', 7);
            ylabel(axesHandle, '% mRNA expressed', 'FontSize', 7);
            set(axesHandle, 'FontSize', 5, 'TickDir', 'out', 'XTick', [0 5], 'YTick', [0 50 100], 'box', 'off');
            
            text(time(end), 98.5, sprintf('\\gamma = %.0f min', mrnaFitResult.tau * 60 * log(2)), ...
                'Parent', axesHandle, ...
                'FontSize', 7, ...
                'horizontalalign', 'right', ...
                'verticalalign', 'top');
            text(time(end), 93.5, sprintf('R^2_{adj} = %.2f', mrnaGof.adjrsquare), ...
                'Parent', axesHandle, ...
                'FontSize', 7, ...
                'horizontalalign', 'right', ...
                'verticalalign', 'top');
        end
    end
    
    methods (Static = true);
        function rgbStr = colorFunc(normVal)
            r = 0; g = 0; b = 0;
            
            if normVal <= 0.5
                b = 255 * (1 - 2 * normVal);
                g = 255 * (2 * normVal);
            else
                r = 255 * (2 * normVal - 1);
                g = 255 * (2 - 2 * normVal);
            end
            
            rgbStr = sprintf('rgb(%d,%d,%d)', uint8(r), uint8(g), uint8(b));
        end
        
        function svg = drawCell(scale, width, cylindricalLength, septumLength, maxTotalLength, maxWidth) %#ok<INUSD>
            width = width / scale;
            cylindricalLength = cylindricalLength / scale;
            septumLength = septumLength / scale;
            H = maxTotalLength / scale;
            W = width;
            
            svg = [...
                sprintf('  <path d="\n') ...
                sprintf('    M%0.4f,%0.4f a%0.4f,%0.4f 0 1 0 0,%0.4f\n', W/2+5, 5, width/2, width/2, width), ...
                sprintf('    l%0.4f,0\n', cylindricalLength/2) ...
                sprintf('    a%0.4f,%0.4f 0 0 0 %0.4f,-%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('    a%0.4f,%0.4f 0 0 0 %0.4f,%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('    l%0.4f,0\n', cylindricalLength/2) ...
                sprintf('    a%0.4f,%0.4f 0 1 0 0,-%0.4f\n', width/2,width/2, width) ...
                sprintf('    l-%0.4f,0\n', cylindricalLength/2) ...
                sprintf('    a%0.4f,%0.4f 0 0 0 -%0.4f,%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('    a%0.4f,%0.4f 0 0 0 -%0.4f,-%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('    z\n') ...
                '    " style="stroke:darkblue;stroke-width:2;stroke-linecap:round;stroke-linejoin:miter;fill:url(#cellBg)"/>\n'];
        end
    end
end