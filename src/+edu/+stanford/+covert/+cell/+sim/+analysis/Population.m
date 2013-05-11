%Population
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Derek Macklin, macklin@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 7/6/2011
classdef Population
    methods (Static = true)
        function run(simBatchDir, fileName, selectedSimulations, methodsToRun)
            import edu.stanford.covert.cell.sim.analysis.Population;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.PrintUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            if nargin < 1 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            %evaluate each method and optionally save at the result as
            %pdf/xls files
            metaData = meta.class.fromName('edu.stanford.covert.cell.sim.analysis.Population');
            for i = 1:numel(metaData.Methods)
                %only plot methods matching the signature
                %  plotFcn(figHandle, simBatchDir, selectedSimulations)
                if ...
                        numel(metaData.Methods{i}.InputNames) < 3 || ~isequal(metaData.Methods{i}.InputNames(1:3), {'figHandle'; 'simBatchDir'; 'selectedSimulations'}) || ...
                        numel(metaData.Methods{i}.OutputNames) < 1 || ~isequal(metaData.Methods{i}.OutputNames(1), {'figData'})
                    continue;
                end
                
                %get method name
                methodName = metaData.Methods{i}.Name;
                MethodName = [upper(methodName(1)) methodName(2:end)];
                if nargin >= 4 && ~ismember(metaData.Methods{i}.Name, methodsToRun)
                    continue;
                end
                
                %run plot method
                [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
                try
                    argout = cell(size(metaData.Methods{i}.OutputNames));
                    [argout{:}] = Population.(methodName)(figHandle, simBatchDir, selectedSimulations);
                    
                    % save plot
                    if nargin >= 2
                        figData = argout{1}; %#ok<NASGU>
                        save([fileName '-' MethodName '.mat'], '-struct', 'figData');
                        
                        orient(figHandle, 'portrait');
                        saveas(figHandle, [fileName '-' MethodName '.pdf']);
                    end
                    
                    % save table
                    if numel(argout) > 1
                        if nargin >= 2
                            PrintUtil.printToFile(argout{2}, argout{3}, [fileName '-' MethodName '.xls'], MethodName);
                        else
                            PrintUtil.printToStdIO(argout{2}, argout{3});
                        end
                    end
                catch exception
                    warning('WholeCell:warning', 'Unable to make %s plot:\n%s', methodName, exception.getReport());
                end
                
                if nargin >= 2
                    close(figHandle);
                    clear figHandle;
                end
                clear argout figData;
            end
        end
        
        function figData = growth(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.MacromoleculeUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.analysis.CellOverview;
            import edu.stanford.covert.util.ConstantUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            comp = sim.compartment;
            g = sim.gene;
            r = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            %% get data
            stateNames = {
                'mass'              'Mass'
                'ploidy'            {'Chromosome' 'Copy Number'}
                'rnas'              'RNA'
                'proteins'          'Proteins'
                'lipids'            'Lipids'
                'growth_rate'       'Growth'
                'ribosomes'         'Ribosomes'
                };
            ensemble = SimulationEnsemble(simBatchDir, stateNames, [], selectedSimulations);
            simEndTimes = ensemble.stateData.simulationEndTimes;
            
            time = ensemble.stateData.time / 3600;
            mass = ensemble.stateData.values(ensemble.getPropertyIndices('mass'), :, :, :) * sim.state('Mass').cellInitialDryWeight / (1-sim.state('Mass').fractionWetWeight) * 1e15;
            growth = ensemble.stateData.values(ensemble.getPropertyIndices('growth_rate'), :, :, :) * 3600 * sim.state('Mass').cellInitialDryWeight / (1 - sim.state('Mass').fractionWetWeight) * 1e15;
            ploidy = ensemble.stateData.values(ensemble.getPropertyIndices('ploidy'), :, :, :);
            rnas = ensemble.stateData.values(ensemble.getPropertyIndices('rnas'), :, :, :);
            proteins = ensemble.stateData.values(ensemble.getPropertyIndices('proteins'), :, :, :);
            lipids = ensemble.stateData.values(ensemble.getPropertyIndices('lipids'), :, :, :);
            ribosomes = ensemble.stateData.values(ensemble.getPropertyIndices('ribosomes'), :, :, :);
            
            [~, simOrder] = sort(mass(:, 1, min(ensemble.stateData.simulationEndTimes), :));
            simOrder = simOrder(:);
            
            clear ensemble;
            
            ribMonomerIdxs = find(any(pc.proteinComplexComposition(g.mRNAIndexs, pc.ribosome70SIndexs, :), 3));
            [~, tmp] = min(r.matureRNAGeneComposition(g.mRNAIndexs(ribMonomerIdxs), :) * r.expression(r.matureIndexs));
            [~, ~, ~, ~, ~, ...
                ~, ~, ~, ~, ...
                ~, ~, ~, ~, ...
                ribMatureRnaIdxs, ~, ~, ~, ...
                ribMonomerIdxs, ~, ~, ribMonomerNames] = ...
                MacromoleculeUtil.getMacromoleculeIndexsIDsNames(pm.wholeCellModelIDs{pm.matureIndexs(g.mRNAIndexs(ribMonomerIdxs(tmp)))}, sim);
            ribMonomerNames = strrep(ribMonomerNames, 'ribosomal protein ', '');
            stateNames = {
                'Rna'               'counts'  r.matureIndexs(ribMatureRnaIdxs)   comp.cytosolIndexs
                'ProteinMonomer'    'counts'  pm.matureIndexs(ribMonomerIdxs)    comp.cytosolIndexs
                'Mass'              'proteinWt' ':' ':'
                };
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations);
            ribRNAs = states.Rna.counts;
            ribMonomers = states.ProteinMonomer.counts;
            proteinWt = states.Mass.proteinWt;
            
            for i = 1:numel(simEndTimes)
                ribRNAs(:, :, simEndTimes(i)+1:end, i) = NaN;
                ribMonomers(:, :, simEndTimes(i)+1:end, i) = NaN;
                proteinWt(:, :, simEndTimes(i)+1:end, i) = NaN;
            end
            
            clear states;
            
            %% clear figure
            clf(figHandle);
            
            %% left panel
            
            %layout plots
            axesHandlesL = PlotUtil.multiElementPlot(figHandle, 3.5 * ones(5, 1), [0 time(end)], struct(...
                'titleStr', 'Composition', ...
                'position', [0.12 0.52 0.21 0.4], ...
                'colorOrder', PlotUtil.getRedGreenColorOrder(simOrder)));
            
            %plot data
            PlotUtil.plotLine(axesHandlesL(1), time, permute(mass, [4 3 2 1]), false, true, false);
            ylabel(axesHandlesL(1), {'Mass' '(fg)'});
            
            PlotUtil.plotLine(axesHandlesL(2), time, permute(ploidy, [4 3 2 1]), false, true, false);
            ylabel(axesHandlesL(2), {'Chrosome' 'Copy' 'Number'});
            
            PlotUtil.plotLine(axesHandlesL(3), time, permute(rnas, [4 3 2 1]), false, true, false);
            ylabel(axesHandlesL(3), 'RNAs');
            
            PlotUtil.plotLine(axesHandlesL(4), time, permute(proteins, [4 3 2 1]), false, true, false);
            ylabel(axesHandlesL(4), 'Proteins');
            
            PlotUtil.plotLine(axesHandlesL(5), time, permute(lipids, [4 3 2 1]), false, true, false);
            ylabel(axesHandlesL(5), 'Lipids');
            
            % Format Y-axis
            PlotUtil.alignYAxesLabels(axesHandlesL);
            PlotUtil.offsetYAxes(axesHandlesL, 0.04);
            PlotUtil.labelSubplots(axesHandlesL(1), 'A', -0.42, 1.2);
            
            %% right panel
            
            %layout plots
            axesHandlesR = PlotUtil.multiElementPlot(figHandle, 3.5 * ones(5, 1), [0 time(end)], struct(...
                'titleStr', 'Growth', ...
                'position', [0.45 0.52 0.21 0.4], ...
                'colorOrder', PlotUtil.getRedGreenColorOrder(simOrder)));
            
            %plot data
            PlotUtil.plotLine(axesHandlesR(1), time, permute(growth, [4 3 2 1]), false, true, false);
            ylabel(axesHandlesR(1), {'Growth' '(fg h^{-1})'});
            
            PlotUtil.plotLine(axesHandlesR(2), time(3:end), diff(permute(sum(proteinWt, 2), [4 3 2 1]), 1, 2), false, true, false);
            ylabel(axesHandlesR(2), {'Protein' 'Synthesis'});
            
            PlotUtil.plotLine(axesHandlesR(3), time, permute(ribosomes, [4 3 2 1]), false, true, false);
            ylabel(axesHandlesR(3), 'Ribosomes');
            
            PlotUtil.plotLine(axesHandlesR(4), time(2:end), permute(ribRNAs, [4 3 2 1]), false, true, false);
            ylabel(axesHandlesR(4), [ribMonomerNames 'Protein']);
            
            PlotUtil.plotLine(axesHandlesR(5), time(2:end), permute(ribMonomers, [4 3 2 1]), false, true, false);
            ylabel(axesHandlesR(5), [ribMonomerNames 'mRNA']);
            
            % Format Y-axis
            PlotUtil.alignYAxesLabels(axesHandlesR);
            PlotUtil.offsetYAxes(axesHandlesR, 0.04);
            PlotUtil.labelSubplots(axesHandlesR(1), 'B', -0.34, 1.2);
            
            %% organize figure data
            figData = struct;
            figData.simOrder = simOrder;
            figData.time = time;
            figData.mass = mass;
            figData.ploidy = ploidy;
            figData.rnas = rnas;
            figData.proteins = proteins;
            figData.lipids = lipids;
            figData.growth = growth;
            figData.proteinWt = proteinWt;
            figData.ribosomes = ribosomes;
            figData.ribRNAs = ribRNAs;
            figData.ribMonomers = ribMonomers;
            figData.ribMonomerNames = ribMonomerNames;
        end
        
        function figData = energyProduction(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            %% constants
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            comp = sim.compartment;
            m = sim.state('Metabolite');
            met = sim.process('Metabolism');
            mass = sim.state('Mass');
            
            metIdxs = sub2ind([numel(m.wholeCellModelIDs) comp.count], ...
                [m.ntpIndexs; m.ndpIndexs; m.nmpIndexs], repmat(comp.cytosolIndexs, 12, 1));
            metProcessIdx = sim.processIndex(met.wholeCellModelID);
            
            %% get data
            stateNames = {
                'Time'               'values'         ':'     ':'
                'MetabolicReaction'  'growth'         ':'     ':'
                'Mass'               'cell'           ':'     '-sum'
                'Metabolite'         'processUsages'  metIdxs  metProcessIdx
                };
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations);
            simEndTimes = permute(max(states.Time.values, [], 3), [4 1 2 3]);
            
            timeIdx = min(simEndTimes);
            [~, simOrder] = sort(states.Mass.cell(:, :, timeIdx, :));
            simOrder = simOrder(:);
            
            time = permute(max(states.Time.values, [], 4), [3 4 1 2]) / 3600;
            states.Metabolite.processUsages = -full(states.Metabolite.processUsages);
            for i = 1:numel(simEndTimes)
                states.MetabolicReaction.growth(:, :, simEndTimes(i)+1:end, i) = NaN;
                states.Metabolite.processUsages(:, :, simEndTimes(i)+1:end, i) = NaN;
            end
            
            %% plot data
            clf(figHandle);
            options = struct;
            options.position = [0.02 0.32 0.95 0.4];
            options.titleStr = {'Growth'; 'A'; 'C'; 'G'; 'U'};
            options.colorOrder = PlotUtil.getRedGreenColorOrder(simOrder);
            options.xdata = time;
            options.ydata = {
                {permute(states.MetabolicReaction.growth * 3600 * mass.cellInitialDryWeight / (1 - mass.fractionWetWeight) * 1e15, [3 4 1 2])}
                permute(mat2cell(permute(states.Metabolite.processUsages(1:4:end, :, :, :), [3 4 1 2]), numel(time), numel(simEndTimes), ones(3, 1)), [3 1 2])
                permute(mat2cell(permute(states.Metabolite.processUsages(2:4:end, :, :, :), [3 4 1 2]), numel(time), numel(simEndTimes), ones(3, 1)), [3 1 2])
                permute(mat2cell(permute(states.Metabolite.processUsages(3:4:end, :, :, :), [3 4 1 2]), numel(time), numel(simEndTimes), ones(3, 1)), [3 1 2])
                permute(mat2cell(permute(states.Metabolite.processUsages(4:4:end, :, :, :), [3 4 1 2]), numel(time), numel(simEndTimes), ones(3, 1)), [3 1 2])
                };
            options.ylabelStr = {
                {'Growth (fg h^{-1})'}
                {'ATP'; 'ADP'; 'AMP'}
                {'CTP'; 'CDP'; 'CMP'}
                {'GTP'; 'GDP'; 'GMP'}
                {'UTP'; 'UDP'; 'UMP'}
                };
            PlotUtil.multiElementPlot(figHandle, {30 [10; 10; 10] [10; 10; 10] [10; 10; 10] [10; 10; 10]}, [0 time(end)], options);
            
            clear options;
            
            %% organize figure data
            figData = struct;
            figData.simOrder = simOrder;
            figData.time = time;
            figData.growth = states.MetabolicReaction.growth;
            figData.processUsages = states.Metabolite.processUsages;
        end
        
        function figData = cellMassDistribution(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            %% get data
            stateNames = {
                'mass' 'Mass'
                };
            ensemble = SimulationEnsemble(simBatchDir, stateNames, [], selectedSimulations);
            mass = ensemble.stateData.values(ensemble.getPropertyIndices('mass'), :, :, :);
            
            simEndTimes = ensemble.stateData.simulationEndTimes;
            
            initMass = permute(mass(:, :, 1, :), [4 3 1 2]);
            finMass = mass(sub2ind(size(mass), ...
                ones(size(simEndTimes)), ...
                ones(size(simEndTimes)), ...
                simEndTimes, ...
                (1:numel(simEndTimes))'));
            
            %% plot data
            clf(figHandle);
            
            mass = linspace(min([initMass; finMass - 1]), max([initMass; finMass - 1]), 20);
            initMassFreq = histc(initMass, mass);
            finMassFreq = histc(finMass - 1, mass);
            
            axesHandle = subplot(2, 1, 1, 'Parent', figHandle);
            bar(axesHandle, mass, initMassFreq);
            xlabel(axesHandle, 'Initial Mass (cell dry weight)');
            ylabel(axesHandle, 'Frequency');
            
            axesHandle = subplot(2, 1, 2, 'Parent', figHandle);
            bar(axesHandle, mass + 1, finMassFreq);
            xlabel(axesHandle, 'Final Mass (cell dry weight)');
            ylabel(axesHandle, 'Frequency');
            
            %% figure data
            figData = struct;
            figData.initMass = initMass;
            figData.finMass = finMass;
        end
        
        function figData = massDistribution(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            comp = sim.compartment;
            
            %% get data
            stateNames = {
                'Time'      'values'             ':' ':'
                'Mass'      'cell'               ':' '-sum'
                'Mass'      'dnaWt'              ':' comp.cytosolIndexs
                'Mass'      'rnaWt'              ':' comp.cytosolIndexs
                'Mass'      'proteinWt'          ':' '-sum'
                'Mass'      'metaboliteWt'       ':' [comp.cytosolIndexs; comp.membraneIndexs]
                };
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations);
            
            time = permute(max(states.Time.values, [], 4), [3 4 1 2]) / 3600;
            mass = permute(states.Mass.cell, [3 4 1 2]) * 1e15;
            dnaWt = permute(states.Mass.dnaWt, [3 4 1 2]) * 1e15;
            rnaWt = permute(states.Mass.rnaWt, [3 4 1 2]) * 1e15;
            proteinWt = permute(states.Mass.proteinWt, [3 4 1 2]) * 1e15;
            membraneWt = permute(states.Mass.metaboliteWt(:, 2, :, :), [3 4 1 2]) * 1e15;
            metaboliteWt = permute(states.Mass.metaboliteWt(:, 1, :, :), [3 4 1 2]) * 1e15;
            
            dnaWt(mass == 0) = NaN;
            rnaWt(mass == 0) = NaN;
            proteinWt(mass == 0) = NaN;
            membraneWt(mass == 0) = NaN;
            metaboliteWt(mass == 0) = NaN;
            mass(mass == 0) = NaN;
            
            minSimEndTime = min(max(states.Time.values, [], 4));
            [~, simOrder] = sort(mass(minSimEndTime, :));
            
            clear states;
            
            %% layout figure
            clf(figHandle);
            
            options = struct();
            options.titleStr = 'Cell Mass';
            options.colorOrder = PlotUtil.getRedGreenColorOrder(simOrder);
            axesHandles = PlotUtil.multiElementPlot(figHandle, 3 * ones(6, 1), [0 time(end)], options);
            
            %% plot data
            
            PlotUtil.plotLine(axesHandles(1), time, mass, false, true, false);
            ylabel(axesHandles(1), 'Mass (fg)');
            
            PlotUtil.plotLine(axesHandles(2), time, dnaWt, false, true, false);
            ylabel(axesHandles(2), 'DNA (fg)');
            
            PlotUtil.plotLine(axesHandles(3), time, rnaWt, false, true, false);
            ylabel(axesHandles(3), 'RNA (fg)');
            
            PlotUtil.plotLine(axesHandles(4), time, proteinWt, false, true, false);
            ylabel(axesHandles(4), 'Protein (fg)');
            
            PlotUtil.plotLine(axesHandles(5), time, membraneWt, false, true, false);
            ylabel(axesHandles(5), 'Membrane (fg)');
            
            PlotUtil.plotLine(axesHandles(6), time, metaboliteWt, false, true, false);
            ylabel(axesHandles(6), 'Metabolites (fg)');
            
            %% Format Y-axis
            PlotUtil.alignYAxesLabels(axesHandles);
            PlotUtil.offsetYAxes(axesHandles, 0.03);
            
            %% organize figure data
            figData = struct;
            figData.simOrder = simOrder;
            figData.time = time;
            figData.mass = mass;
            figData.dnaWt = dnaWt;
            figData.rnaWt = rnaWt;
            figData.proteinWt = proteinWt;
            figData.membraneWt = membraneWt;
            figData.metaboliteWt = metaboliteWt;
        end
        
        function figData = cellCyclePhases(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.SummaryLogger;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            %% get data
            [simDir, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            simBatchDir = simDir(1:find(simDir == filesep, 1, 'last') - 1);
            mass = sim.state('Mass');
            stats = SummaryLogger.getSimulationStatistics(simBatchDir);
            nSims = size(stats, 1);
            
            ensemble = SimulationEnsemble(simBatchDir, cell(0, 2), [], selectedSimulations);
            maxTime = ensemble.stateData.time(end);
            
            clear ensemble;
            
            %% clear figure
            clf(figHandle);
            
            %% top panel
            
            %layout plots
            edges = 15 * 60 * (0 : ceil(maxTime / (15 * 60)));
            n1 = histc(stats(:, SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME), edges);
            n2 = histc(stats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME), edges);
            n3 = histc(stats(:, SummaryLogger.SIM_STATUS_INDEX_REP_TIME) - stats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME), edges);
            n4 = histc(stats(:, SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME) - stats(:, SummaryLogger.SIM_STATUS_INDEX_REP_TIME), edges);
            
            [axesHandlesT, xAxesHandleT] = PlotUtil.multiElementPlot(figHandle, 3 * ones(4, 1), edges([1 end]) / 3600', ...
                struct('position', [0.12 0.68 0.21 0.24], 'xlabelStr', 'Cell Cycle Phase Duration (h)')); %#ok<NASGU>
            
            %plot data
            bar(axesHandlesT(1), edges / 3600, n1 / nSims * 100, 'histc');
            ylabel(axesHandlesT(1), {'Cell' 'Cycle'});
            ytick = get(axesHandlesT(1), 'YTick');
            set(axesHandlesT(1), 'YTick', ytick([1 end]));
            ylim = get(axesHandlesT(1), 'YLim');
            set(axesHandlesT(1), 'YLim', ylim);
            
            bar(axesHandlesT(2), edges / 3600, n2 / nSims * 100, 'histc');
            ylabel(axesHandlesT(2), {'Rep' 'Init'});
            ytick = get(axesHandlesT(2), 'YTick');
            set(axesHandlesT(2), 'YTick', ytick([1 end]));
            ylim = get(axesHandlesT(2), 'YLim');
            set(axesHandlesT(2), 'YLim', ylim);
            
            bar(axesHandlesT(3), edges / 3600, n3 / nSims * 100, 'histc');
            ylabel(axesHandlesT(3), {'Repli-' 'cation'});
            ytick = get(axesHandlesT(3), 'YTick');
            set(axesHandlesT(3), 'YTick', ytick([1 end]));
            ylim = get(axesHandlesT(3), 'YLim');
            set(axesHandlesT(3), 'YLim', ylim);
            
            bar(axesHandlesT(4), edges / 3600, n4 / nSims * 100, 'histc');
            ylabel(axesHandlesT(4), {'Cyto-' 'kinesis'});
            ytick = get(axesHandlesT(4), 'YTick');
            set(axesHandlesT(4), 'YTick', ytick([1 end]));
            ylim = get(axesHandlesT(4), 'YLim');
            set(axesHandlesT(4), 'YLim', ylim);
            
            %xlabelPos = get(get(xAxesHandleT, 'XLabel'), 'Position');
            %set(get(xAxesHandleT, 'XLabel'), 'Position', xlabelPos+[0 17 0]);
            
            %% bottom panel
            edges = linspace(min(stats(:, SummaryLogger.SIM_STATUS_INDEX_FIN_MASS)), max(stats(:, SummaryLogger.SIM_STATUS_INDEX_FIN_MASS)), 50);
            n = histc(stats(:, SummaryLogger.SIM_STATUS_INDEX_FIN_MASS), edges);
            
            [axesHandlesB, xAxesHandleB] = PlotUtil.multiElementPlot(figHandle, 3, ...
                edges([1 end]) * mass.cellInitialDryWeight / (1-mass.fractionWetWeight) * 1e15, ...
                struct('position', [0.12 0.52 0.21 0.06], 'xlabelStr', 'Cell Mass (fg)')); %#ok<NASGU>
            
            bar(axesHandlesB(1), edges * mass.cellInitialDryWeight / (1-mass.fractionWetWeight) * 1e15, n / nSims * 100, 'histc');
            ylabel(axesHandlesB(1), 'Freq');
            ytick = get(axesHandlesB(1), 'YTick');
            set(axesHandlesB(1), 'YTick', ytick([1 end]));
            ylim = get(axesHandlesB(1), 'YLim');
            set(axesHandlesB(1), 'YLim', ylim);
            
            %xlabelPos = get(get(xAxesHandleB, 'XLabel'), 'Position');
            %set(get(xAxesHandleB, 'XLabel'), 'Position', xlabelPos+[0 15 0]);
            
            %% Format Y-axis
            PlotUtil.alignYAxesLabels(axesHandlesT);
            PlotUtil.offsetYAxes([axesHandlesT; axesHandlesB], 0.04);
            PlotUtil.labelSubplots(axesHandlesT(1), 'A', -0.25, 1.2);
            PlotUtil.labelSubplots(axesHandlesB(1), 'B', -0.25, 1.2);
            
            %% organize figure data
            figData = struct;
            figData.nSims = nSims;
            figData.maxTime = maxTime;
            figData.mass = mass;
            figData.stats = stats;
        end
        
        function [figData, tabContent, tabColLabels] = nucleotideDistribution(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.MacromoleculeUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.analysis.CellOverview;
            import edu.stanford.covert.util.ConstantUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            comp = sim.compartment;
            m = sim.state('Metabolite');
            
            %% get data
            stateNames = {
                'mass' 'Mass'
                };
            ensemble = SimulationEnsemble(simBatchDir, stateNames, [], selectedSimulations);            
            mass = ensemble.stateData.values(ensemble.getPropertyIndices('mass'), :, :, :) * sim.state('Mass').cellInitialDryWeight / (1-sim.state('Mass').fractionWetWeight) * 1e15;
            
            simEndTimes = ensemble.stateData.simulationEndTimes;
            timeIdx = ceil(min(simEndTimes) / ensemble.stateData.downsampleStepSec);
            [~, simOrder] = sort(mass(:, 1, timeIdx, :));
            simOrder = simOrder(:);
            
            clear ensemble mass;
            
            stateNames = {
                'Time'       'values' ':'  ':'
                'Metabolite' 'counts' [m.ntpIndexs; m.ndpIndexs; m.nmpIndexs; m.aminoAcidIndexs] comp.cytosolIndexs
                };
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations);
            time = max(full(states.Time.values), [], 4) / 3600;
            ntps = full(states.Metabolite.counts(1:4, :, :, :));
            ndps = full(states.Metabolite.counts(5:8, :, :, :));
            nmps = full(states.Metabolite.counts(9:12, :, :, :));
            aas = full(states.Metabolite.counts(13:end, :, :, :));
            
            for i = 1:numel(simEndTimes)
                ntps(:, :, simEndTimes(i)+1:end, i) = NaN;
                ndps(:, :, simEndTimes(i)+1:end, i) = NaN;
                nmps(:, :, simEndTimes(i)+1:end, i) = NaN;
                aas(:, :, simEndTimes(i)+1:end, i) = NaN;
            end
            
            clear states;
            
            %% clear figure
            clf(figHandle);
            
            %% left panel
            
            %layout plots
            options = struct();
            options.titleStr = {'Nucleotides'; 'AXP'; 'CXP'; 'GXP'; 'UXP'};
            options.xdata = permute(time, [4 3 2 1]);
            options.ydata = {{
                permute(sum(aas, 1), [4 3 2 1])
                permute(sum(ntps, 1), [4 3 2 1])
                permute(sum(ndps, 1), [4 3 2 1])
                permute(sum(nmps, 1), [4 3 2 1])
                };{
                []
                permute(ntps(1, :, :, :), [4 3 2 1])
                permute(ndps(1, :, :, :), [4 3 2 1])
                permute(nmps(1, :, :, :), [4 3 2 1])
                };{
                []
                permute(ntps(2, :, :, :), [4 3 2 1])
                []
                permute(nmps(2, :, :, :), [4 3 2 1])
                };{
                []
                permute(ntps(3, :, :, :), [4 3 2 1])
                permute(ndps(3, :, :, :), [4 3 2 1])
                permute(nmps(3, :, :, :), [4 3 2 1])
                }; {
                []
                permute(ntps(4, :, :, :), [4 3 2 1])
                []
                permute(nmps(4, :, :, :), [4 3 2 1])
                }};
            options.ylabelStr = {
                {'Amino Acid'; 'NTP'; 'NDP'; 'NMP'};
                {}; {}; {}; {}
                };
            options.colorOrder = PlotUtil.getRedGreenColorOrder(simOrder);
            [axesHandles, ~, ~, offsetAxesHandles] = PlotUtil.multiElementPlot(figHandle, repmat({3.5 * ones(4, 1)}, 5, 1), [0 time(end)], options);
            
            set(cellfun(@(x) x(1), axesHandles(2:end)), 'visible', 'off');
            set(cellfun(@(x) x(1), offsetAxesHandles(2:end)), 'visible', 'off');
            set(cellfun(@(x) x(3), axesHandles([3 5])), 'visible', 'off');
            set(cellfun(@(x) x(3), offsetAxesHandles([3 5])), 'visible', 'off');
            
            clear options;
            
            %% table
            aaSums = sum(aas, 1);
            medianNtps = [
                median(ntps(sub2ind(size(ntps), 1 * ones(size(simEndTimes)), ones(size(simEndTimes)), simEndTimes, (1:numel(simEndTimes))')))
                median(ntps(sub2ind(size(ntps), 2 * ones(size(simEndTimes)), ones(size(simEndTimes)), simEndTimes, (1:numel(simEndTimes))')))
                median(ntps(sub2ind(size(ntps), 3 * ones(size(simEndTimes)), ones(size(simEndTimes)), simEndTimes, (1:numel(simEndTimes))')))
                median(ntps(sub2ind(size(ntps), 4 * ones(size(simEndTimes)), ones(size(simEndTimes)), simEndTimes, (1:numel(simEndTimes))')))
                ];
            medianNdps = [
                median(ndps(sub2ind(size(ndps), 1 * ones(size(simEndTimes)), ones(size(simEndTimes)), simEndTimes, (1:numel(simEndTimes))')))
                median(ndps(sub2ind(size(ndps), 2 * ones(size(simEndTimes)), ones(size(simEndTimes)), simEndTimes, (1:numel(simEndTimes))')))
                median(ndps(sub2ind(size(ndps), 3 * ones(size(simEndTimes)), ones(size(simEndTimes)), simEndTimes, (1:numel(simEndTimes))')))
                median(ndps(sub2ind(size(ndps), 4 * ones(size(simEndTimes)), ones(size(simEndTimes)), simEndTimes, (1:numel(simEndTimes))')))
                ];
            medianNmps = [
                median(nmps(sub2ind(size(nmps), 1 * ones(size(simEndTimes)), ones(size(simEndTimes)), simEndTimes, (1:numel(simEndTimes))')))
                median(nmps(sub2ind(size(nmps), 2 * ones(size(simEndTimes)), ones(size(simEndTimes)), simEndTimes, (1:numel(simEndTimes))')))
                median(nmps(sub2ind(size(nmps), 3 * ones(size(simEndTimes)), ones(size(simEndTimes)), simEndTimes, (1:numel(simEndTimes))')))
                median(nmps(sub2ind(size(nmps), 4 * ones(size(simEndTimes)), ones(size(simEndTimes)), simEndTimes, (1:numel(simEndTimes))')))
                ];
            medianAAs = median(aaSums(sub2ind(size(aaSums), ones(size(simEndTimes)), ones(size(simEndTimes)), simEndTimes, (1:numel(simEndTimes))')));
            nxps = [
                median(ntps(:, :, 1, :), 4)' medianNtps'
                median(ndps(:, :, 1, :), 4)' medianNdps'
                median(nmps(:, :, 1, :), 4)' medianNmps'
                median(aaSums(:, :, 1, :)) zeros(1, 3) medianAAs zeros(1, 3)
                ];
            tabColLabels = {[] 't = 0' [] [] [] ['t = ' num2str(median(simEndTimes))] [] [] []};
            tabContent = [{[] 'A' 'C' 'G' 'U' 'A' 'C' 'G' 'U';}
                {'NTP'; 'NDP'; 'NMP'; 'AAs'} num2cell(nxps)];
            tabContent(end, [3:5 7:9]) = {[]};
            
            %% organize figure data
            figData = struct;
            figData.simEndTimes = simEndTimes;
            figData.simOrder = simOrder;
            figData.timeIdx = timeIdx;
            figData.time = time;
            figData.aas = aas;
            figData.ntps = ntps;
            figData.ndps = ndps;
            figData.nmps = nmps;
        end
        
        function figData = aminoAcidCounts(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            comp = sim.compartment;
            m = sim.state('Metabolite');
            
            %% get data
            stateNames = {
                'Time'        'values'  ':'                ':'
                'Mass'        'cell'    ':'                '-sum'
                'Metabolite'  'counts'  m.aminoAcidIndexs  comp.cytosolIndexs
                };
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations);
            
            time = permute(max(states.Time.values, [], 4), [4 3 1 2]) / 3600;
            simEndTimes = squeeze(max(states.Time.values, [], 3));
            [~, simOrder] = sort(states.Mass.cell(sub2ind(size(states.Mass.cell), ones(size(simEndTimes)), ones(size(simEndTimes)), simEndTimes, (1:numel(simEndTimes))')));
            
            aminoAcids = full(permute(states.Metabolite.counts, [4 3 1 2]));
            for i = 1:numel(simEndTimes)
                aminoAcids(i, simEndTimes(i) + 1:end, :, :) = NaN;
            end
            
            %% plot data
            clf(figHandle);
            
            options = struct;
            options.colorOrder = PlotUtil.getRedGreenColorOrder(simOrder);
            options.xdata = time;
            options.ylabelStr = {
                m.wholeCellModelIDs(m.aminoAcidIndexs(1:6));
                cell(6, 1);
                m.wholeCellModelIDs(m.aminoAcidIndexs(7:12));
                cell(6, 1);
                m.wholeCellModelIDs(m.aminoAcidIndexs(13:18));
                cell(6, 1);
                [m.wholeCellModelIDs(m.aminoAcidIndexs(19:21)); {[]; []; []}]
                cell(6, 1);
                };
            options.ydata = {
                permute(mat2cell(aminoAcids(:, :, 1:6), size(aminoAcids, 1), size(aminoAcids, 2), ones(1, 6)), [3 1 2])
                repmat({[]}, 6, 1)
                permute(mat2cell(aminoAcids(:, :, 7:12), size(aminoAcids, 1), size(aminoAcids, 2), ones(1, 6)), [3 1 2])
                repmat({[]}, 6, 1)
                permute(mat2cell(aminoAcids(:, :, 13:18), size(aminoAcids, 1), size(aminoAcids, 2), ones(1, 6)), [3 1 2])
                repmat({[]}, 6, 1)
                [permute(mat2cell(aminoAcids(:, :, 19:21), size(aminoAcids, 1), size(aminoAcids, 2), ones(1, 3)), [3 1 2]); {[]; []; []}]
                repmat({[]}, 6, 1)
                };
            options.colWidths = repmat([4; 1], 4, 1);
            options.yAxesLabelWidths = repmat([0.1; 0.01], 4, 1);
            [axesHandles, xAxesHandles, ~, offsetAxesHandles] = PlotUtil.multiElementPlot(figHandle, repmat({10 * ones(6, 1)}, 8, 1), [0 time(end)], options);
            
            clear options;
            
            delete(axesHandles{7}(4:6));
            delete(axesHandles{8}(4:6));
            delete(offsetAxesHandles{7}(4:6));
            delete(offsetAxesHandles{8}(4:6));
            delete(cell2mat(xAxesHandles(2:2:end)));
            
            for i = 1:numel(m.aminoAcidIndexs)
                row = mod(i - 1, 6) + 1;
                col = 2 * ceil(i / 6);
                
                cnts = aminoAcids(sub2ind(size(aminoAcids), (1:numel(simEndTimes))', simEndTimes, i * ones(size(simEndTimes))));
                edges = linspace(min(min(aminoAcids(:, :, i))), max(max(aminoAcids(:, :, i))), 10);
                edges = edges(:);
                freq = histc(cnts, edges);
                
                axesHandle = axesHandles{col}(row);
                offsetAxesHandle = offsetAxesHandles{col}(row);
                
                cla(axesHandle)
                cla(offsetAxesHandle)
                
                barh(axesHandle, edges, freq)
                xlim(axesHandle, [0 max(freq(:))]);
                ylim(axesHandle, [min(min(aminoAcids(:, :, i))) max(max(aminoAcids(:, :, i)))]);
                set(axesHandle, 'YTick', [])
                set(axesHandle, 'XTick', [])
                set(offsetAxesHandle, 'YTick', [])
                set(offsetAxesHandle, 'XTick', [])
            end
            
            %% data
            figData = struct;
            figData.simOrder = simOrder;
            figData.time = time;
            figData.aminoAcids = aminoAcids;
        end
        
        function [figData, tabContent, tabColLabels] = metabolism(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            met = sim.process('Metabolism');
            
            %% get data
            corrEnzymeFlux = zeros(numel(met.reactionWholeCellModelIDs), numel(selectedSimulations));
            corrEnzymeGrowth = zeros(numel(met.enzymeWholeCellModelIDs), numel(selectedSimulations));
            corrFluxGrowth = zeros(numel(met.reactionWholeCellModelIDs), numel(selectedSimulations));
            stateNames = {
                'MetabolicReaction'  'fluxs'   ':'                                            ':'
                'MetabolicReaction'  'growth'  ':'                                            ':'
                'ProteinMonomer'     'counts'  pm.matureIndexs(met.enzymeMonomerGlobalIndexs) '-sum'
                'ProteinComplex'     'counts'  pc.matureIndexs(met.enzymeComplexGlobalIndexs) '-sum'
                };
            for i = 1:numel(selectedSimulations)
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations(i));
                growth = permute(states.MetabolicReaction.growth, [1 3 2]);
                enzymes = zeros(numel(met.enzymeWholeCellModelIDs), size(states.ProteinMonomer.counts, 3));
                enzymes(met.enzymeMonomerLocalIndexs, :, :) = full(permute(states.ProteinMonomer.counts, [1 3 2]));
                enzymes(met.enzymeComplexLocalIndexs, :, :) = full(permute(states.ProteinComplex.counts, [1 3 2]));
                fluxs = permute(states.MetabolicReaction.fluxs, [1 3 2]);
                
                for j = 1:numel(met.reactionIndexs_fba)
                    k = met.reactionIndexs_fba(j);
                    corrEnzymeFlux(k, i) = corr(fluxs(k, :)', (met.reactionCatalysisMatrix(k, :) * enzymes)');
                end
                corrEnzymeGrowth(:, i) = corr(enzymes', growth');
                corrFluxGrowth(:, i) = corr(fluxs', growth');
                
                clear states;
            end
            
            corrEnzymeFlux(isnan(corrEnzymeFlux)) = 0;
            corrEnzymeGrowth(isnan(corrEnzymeGrowth)) = 0;
            corrFluxGrowth(isnan(corrFluxGrowth)) = 0;
            
            corrEnzymeFluxMean = nanmean(corrEnzymeFlux, 2);
            corrEnzymeGrowthMean = nanmean(corrEnzymeGrowth, 2);
            corrFluxGrowthMean = nanmean(corrFluxGrowth, 2);
            corrEnzymeFluxStd = nanstd(corrEnzymeFlux, [], 2);
            corrEnzymeGrowthStd = nanstd(corrEnzymeGrowth, [], 2);
            corrFluxGrowthStd = nanstd(corrFluxGrowth, [], 2);
            
            [~, enzIdx] = max(corrEnzymeGrowthMean);
            [~, simOrder] = sort(corrEnzymeGrowth(enzIdx, :));
            
            clear corrEnzymeFlux corrEnzymeGrowth corrFluxGrowth;
            
            %% get temporal data
            rxnIdx = find(met.reactionCatalysisMatrix(:, enzIdx));
            monIdx = pm.matureIndexs(met.enzymeMonomerGlobalIndexs(met.enzymeMonomerLocalIndexs == enzIdx));
            cpxIdx = pc.matureIndexs(met.enzymeComplexGlobalIndexs(met.enzymeComplexLocalIndexs == enzIdx));
            
            stateNames = {
                'Time'               'values'  ':'     ':'
                'MetabolicReaction'  'fluxs'   rxnIdx  ':'
                };
            if ~isempty(monIdx)
                stateNames = [stateNames; {'ProteinMonomer' 'counts' monIdx '-sum'}];
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations);
                enzymes = full(permute(states.ProteinMonomer.counts, [3 4 1 2]));
            else
                stateNames = [stateNames; {'ProteinComplex' 'counts' cpxIdx '-sum'}];
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations);
                enzymes = full(permute(states.ProteinComplex.counts, [3 4 1 2]));
            end
            time = permute(max(states.Time.values, [], 4), [3 4 1 2]) / 3600;
            simEndTimes = squeeze(max(states.Time.values, [], 3));
            fluxs = full(permute(states.MetabolicReaction.fluxs, [3 4 1 2]));
            
            clear states;
            
            for i = 1:numel(simEndTimes)
                enzymes(simEndTimes(i)+1:end, i) = NaN;
                fluxs(simEndTimes(i)+1:end, i) = NaN;
            end
            
            colorOrder = PlotUtil.getRedGreenColorOrder(simOrder);
            
            %% plot
            clf(figHandle);
            
            axesHandles = repmat({zeros(4, 1)}, 2, 1);
            
            %enzyme-fluxs
            axesHandle = subplot(4, 2, 1);
            axesHandles{1}(1) = axesHandle;
            [~, order] = sort(corrEnzymeFluxMean);
            errorbar((1:numel(corrEnzymeFluxMean))', corrEnzymeFluxMean(order), corrEnzymeFluxStd(order), ...
                'Parent', axesHandle, 'LineStyle', 'none', 'Marker', '.');
            xlabel(axesHandle, 'Flux', 'FontSize', 8);
            ylabel(axesHandle, {'Corr' '(Flux, Enzyme)'}, 'FontSize', 8);
            xlim(axesHandle, [0.5 numel(corrEnzymeFluxMean) + 0.5]);
            ylim(axesHandle, [min(corrEnzymeFluxMean - corrEnzymeFluxStd) max(corrEnzymeFluxMean - corrEnzymeFluxStd)]);
            
            axesHandle = subplot(4, 2, 2);
            axesHandles{2}(1) = axesHandle;
            hist(axesHandle, corrEnzymeFluxMean);
            xlabel(axesHandle, 'Corr (Flux, Enzyme)', 'FontSize', 8);
            ylabel(axesHandle, 'Freq', 'FontSize', 8);
            
            %fluxes-growth
            axesHandle = subplot(4, 2, 3);
            axesHandles{1}(2) = axesHandle;
            [~, order] = sort(corrFluxGrowthMean);
            errorbar((1:numel(corrFluxGrowthMean))', corrFluxGrowthMean(order), corrFluxGrowthStd(order), ...
                'Parent', axesHandle, 'LineStyle', 'none', 'Marker', '.');
            xlabel(axesHandle, 'Flux', 'FontSize', 8);
            ylabel(axesHandle, {'Corr' '(Flux, Growth)'}, 'FontSize', 8);
            xlim(axesHandle, [0.5 numel(corrFluxGrowthMean) + 0.5]);
            ylim(axesHandle, [min(corrFluxGrowthMean - corrFluxGrowthStd) max(corrFluxGrowthMean - corrFluxGrowthStd)]);
            
            axesHandle = subplot(4, 2, 4);
            axesHandles{2}(2) = axesHandle;
            hist(axesHandle, corrFluxGrowthMean);
            xlabel(axesHandle, 'Corr (Flux, Growth)', 'FontSize', 8);
            ylabel(axesHandle, 'Freq', 'FontSize', 8);
            
            %enzymes-growth
            axesHandle = subplot(4, 2, 5);
            axesHandles{1}(3) = axesHandle;
            [~, order] = sort(corrEnzymeGrowthMean);
            errorbar((1:numel(corrEnzymeGrowthMean))', corrEnzymeGrowthMean(order), corrEnzymeGrowthStd(order), ...
                'Parent', axesHandle, 'LineStyle', 'none', 'Marker', '.');
            xlabel(axesHandle, 'Flux', 'FontSize', 8);
            ylabel(axesHandle, {'Corr' '(Growth, Enzyme)'}, 'FontSize', 8);
            xlim(axesHandle, [0.5 numel(corrEnzymeGrowthMean) + 0.5]);
            ylim(axesHandle, [min(corrEnzymeGrowthMean - corrEnzymeGrowthStd) max(corrEnzymeGrowthMean - corrEnzymeGrowthStd)]);
            
            axesHandle = subplot(4, 2, 6);
            axesHandles{2}(3) = axesHandle;
            hist(axesHandle, corrEnzymeGrowthMean);
            xlabel(axesHandle, 'Corr (Growth, Enzyme)', 'FontSize', 8);
            ylabel(axesHandle, 'Freq', 'FontSize', 8);
            
            %flux-time
            axesHandle = subplot(4, 2, 7);
            axesHandles{1}(4) = axesHandle;
            hold(axesHandle, 'on');
            set(axesHandle, 'ColorOrder', colorOrder);
            PlotUtil.plotLine(axesHandle, time, fluxs, true, true);
            xlabel(axesHandle, 'Time (h)', 'FontSize', 8);
            ylabel(axesHandle, [met.reactionWholeCellModelIDs(rxnIdx) 'Flux'], ...
                'FontSize', 8, 'Interpreter', 'none');
            
            %enzyme-time
            axesHandle = subplot(4, 2, 8);
            axesHandles{2}(4) = axesHandle;
            hold(axesHandle, 'on');
            set(axesHandle, 'ColorOrder', colorOrder);
            PlotUtil.plotLine(axesHandle, time, enzymes, true, true);
            xlabel(axesHandle, 'Time (h)', 'FontSize', 8);
            ylabel(axesHandle, [met.enzymeWholeCellModelIDs(enzIdx) 'Enzyme'], ...
                'FontSize', 8, 'Interpreter', 'none');
            
            %% format plots
            for i = 1:numel(axesHandles)
                for j = 1:numel(axesHandles{i})
                    axesHandle = axesHandles{i}(j);
                    set(axesHandle, 'FontSize', 6);
                    set(get(axesHandle, 'YLabel'), 'Units', 'normalized')
                end
            end
            PlotUtil.alignYAxesLabels(axesHandles);
            
            %% layout table
            tabColLabels = {
                'Reaction' 'Mean Enzyme Correlation' 'Std Enzyme Correlation' ...
                [] ...
                'Reaction' 'Mean Growth Correlation' 'Std Growth Correlation' ...
                [] ...
                'Enzyme'   'Mean Growth Correlation' 'Std Growth Correlation' ...
                };
            tabContent = [
                met.reactionWholeCellModelIDs num2cell(corrEnzymeFluxMean) num2cell(corrEnzymeFluxStd) ...
                cell(size(met.reactionWholeCellModelIDs)) ...
                met.reactionWholeCellModelIDs num2cell(corrFluxGrowthMean) num2cell(corrFluxGrowthStd) ...
                cell(size(met.reactionWholeCellModelIDs)) ...
                [met.enzymeWholeCellModelIDs num2cell(corrEnzymeGrowthMean) num2cell(corrEnzymeGrowthStd); cell(numel(met.reactionWholeCellModelIDs) - numel(met.enzymeWholeCellModelIDs), 3)] ...
                ];
            
            %% organize figure data
            figData = struct;
            figData.simOrder = simOrder;
            figData.time = time;
            figData.enzymes = enzymes;
            figData.fluxs = fluxs;
        end
        
        function figData = metaboliteConcentrations(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            comp = sim.compartment;
            m = sim.state('Metabolite');
            
            %% get data
            metIdxs = [m.ntpIndexs; m.ndpIndexs; m.nmpIndexs; m.dntpIndexs; m.aminoAcidIndexs; m.phosphateIndexs; m.diphosphateIndexs; m.hydrogenIndexs];
            stateNames = {
                'Geometry'    'volume'  ':'      ':'
                'Metabolite'  'counts'  metIdxs  comp.cytosolIndexs
                };
            meanMetConcs = zeros(numel(metIdxs), 1, 1, numel(selectedSimulations));
            for i = 1:numel(selectedSimulations)
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations(i));
                meanMetConcs(:, :, :, i) = mean(states.Metabolite.counts ./ states.Geometry.volume(ones(size(metIdxs)), :, :), 3) / ConstantUtil.nAvogadro * 1000;
                
                clear states;
            end
            
            meanConcs = mean(meanMetConcs, 4);
            stdConcs = std(meanMetConcs, [], 4);
            
            clear meanMetConcs;
            
            %% plot
            clf(figHandle);
            
            axesHandle = subplot(1, 1, 1, 'Position', [0.1300    0.1500    0.7750    0.8150]);
            errorbar((1:numel(metIdxs))', meanConcs, meanConcs - max(1e-3, meanConcs - stdConcs), stdConcs, ...
                'Parent', axesHandle, 'LineStyle', 'none', 'Marker', '.');
            set(axesHandle, 'YScale', 'log');
            xlim(axesHandle, [0.5  numel(metIdxs) + 0.5]);
            set(axesHandle, ...
                'XTick', 1:numel(metIdxs), ...
                'XTickLabel', [], ...
                'FontSize', 8, ...
                'TickDir', 'out');
            tick2text(axesHandle, 'axis', 'x');
            for i = 1:numel(metIdxs)
                text((i-0.5) / numel(metIdxs), -0.02, m.wholeCellModelIDs{metIdxs(i)}, ...
                    'Parent', axesHandle, ...
                    'FontSize', 6, ...
                    'Rotation', 270, ...
                    'VerticalAlignment', 'middle', ...
                    'HorizontalAlignment', 'left', ...
                    'Units', 'normalized');
            end
            axesPos = get(axesHandle, 'Position');
            
            annotation(figHandle, 'line', axesPos(1) + axesPos(3) * ((-0.5 + [1 4] + [-0.25 0.35]) / numel(metIdxs)), axesPos(2) + axesPos(4) * (-0.09 * [1 1]), ...
                'Units', 'normalized');
            annotation(figHandle, 'line', axesPos(1) + axesPos(3) * ((-0.5 + [5 8] + [-0.25 0.35]) / numel(metIdxs)), axesPos(2) + axesPos(4) * (-0.09 * [1 1]), ...
                'Units', 'normalized');
            annotation(figHandle, 'line', axesPos(1) + axesPos(3) * ((-0.5 + [9 12] + [-0.25 0.35]) / numel(metIdxs)), axesPos(2) + axesPos(4) * (-0.09 * [1 1]), ...
                'Units', 'normalized');
            annotation(figHandle, 'line', axesPos(1) + axesPos(3) * ((-0.5 + [13 16] + [-0.25 0.35]) / numel(metIdxs)), axesPos(2) + axesPos(4) * (-0.09 * [1 1]), ...
                'Units', 'normalized');
            annotation(figHandle, 'line', axesPos(1) + axesPos(3) * ((-0.5 + [17 37] + [-0.25 0.35]) / numel(metIdxs)), axesPos(2) + axesPos(4) * (-0.09 * [1 1]), ...
                'Units', 'normalized');
            annotation(figHandle, 'line', axesPos(1) + axesPos(3) * ((-0.5 + [38 40] + [-0.25 0.35]) / numel(metIdxs)), axesPos(2) + axesPos(4) * (-0.09 * [1 1]), ...
                'Units', 'normalized');
            
            text(2 / numel(metIdxs), -0.092, 'NTPs', ...
                'Parent', axesHandle, ...
                'FontSize', 6, ...
                'VerticalAlignment', 'top', ...
                'HorizontalAlignment', 'center', ...
                'Units', 'normalized');
            text(6 / numel(metIdxs), -0.092, 'NDPs', ...
                'Parent', axesHandle, ...
                'FontSize', 6, ...
                'VerticalAlignment', 'top', ...
                'HorizontalAlignment', 'center', ...
                'Units', 'normalized');
            text(10 / numel(metIdxs), -0.092, 'NMPs', ...
                'Parent', axesHandle, ...
                'FontSize', 6, ...
                'VerticalAlignment', 'top', ...
                'HorizontalAlignment', 'center', ...
                'Units', 'normalized');
            text(14 / numel(metIdxs), -0.092, 'dNTPs', ...
                'Parent', axesHandle, ...
                'FontSize', 6, ...
                'VerticalAlignment', 'top', ...
                'HorizontalAlignment', 'center', ...
                'Units', 'normalized');
            text(26.5 / numel(metIdxs), -0.092, 'AAs', ...
                'Parent', axesHandle, ...
                'FontSize', 6, ...
                'VerticalAlignment', 'top', ...
                'HorizontalAlignment', 'center', ...
                'Units', 'normalized');
            text(38.5 / numel(metIdxs), -0.092, 'Other', ...
                'Parent', axesHandle, ...
                'FontSize', 6, ...
                'VerticalAlignment', 'top', ...
                'HorizontalAlignment', 'center', ...
                'Units', 'normalized');
            
            xlabel(axesHandle, 'Metabolite', 'FontSize', 12);
            set(get(axesHandle, 'xlabel'), 'Position', [0.5 -0.12 1], 'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
            ylabel(axesHandle, 'Concentration (mM)', 'FontSize', 12);
            
            %% organize figure data
            figData = struct;
            figData.meanConcs = meanConcs;
            figData.stdConcs = stdConcs;
        end
        
        function figData = translation(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            comp = sim.compartment;
            g = sim.gene;
            
            t = sim.state('Time');
            mass = sim.state('Mass');
            m = sim.state('Metabolite');
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            rib = sim.state('Ribosome');
            pol = sim.state('Polypeptide');
            
            tl = sim.process('Translation');
            ta = sim.process('tRNAAminoacylation');
            
            expTotRnas = mass.initialFractionNTPsInRNAs * mass.dryWeightFractionRNA * mass.cellInitialDryWeight * rna.expression(rna.matureIndexs) / ...
                (rna.expression(rna.matureIndexs)' * rna.molecularWeights(rna.matureIndexs) / ConstantUtil.nAvogadro);
            monExp = (rna.matureRNAGeneComposition(g.mRNAIndexs, :) * rna.expression(rna.matureIndexs)) ./ ...
                (log(2) / t.cellCycleLength + pm.decayRates(pm.matureIndexs));
            expTotMons = mass.initialFractionAAsInMonomers * mass.dryWeightFractionProtein * mass.cellInitialDryWeight * ...
                monExp / (monExp' * pm.molecularWeights(pm.matureIndexs) / ConstantUtil.nAvogadro);
            expTotCpxs = min(repmat(expTotMons, [1 numel(pc.matureIndexs)]) ./ sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3), [], 1);
            
            %% get data
            %mass, energy
            ensemble = SimulationEnsemble(simBatchDir, {'mass' 'mass'; 'atp' 'atp'}, [], selectedSimulations);
            time = ensemble.stateData.time / 3600;
            mass = ensemble.stateData.values(ensemble.getPropertyIndices('mass'), :, :, :) * sim.state('Mass').cellInitialDryWeight / (1-sim.state('Mass').fractionWetWeight) * 1e15;
            simEndTimes = ensemble.stateData.simulationEndTimes;
            
            [~, simOrder] = sort(mass(:, 1, min(simEndTimes), :));
            simOrder = simOrder(:);
            
            %cleanup
            clear ensemble mass;
            
            %metabolites, rna, proteins
            monIdxs = [tl.enzymeMonomerGlobalIndexs; ta.enzymeMonomerGlobalIndexs];
            cpxIdxs = [tl.enzymeComplexGlobalIndexs; ta.enzymeComplexGlobalIndexs];
            stateNames = {
                'Metabolite'     'counts' [m.ntpIndexs([1 3]); m.ndpIndexs([1 3]); m.nmpIndexs(1); m.aminoAcidIndexs; m.getIndexs('FTHF10')] comp.cytosolIndexs
                'Rna'            'counts' [rna.matureIndexs(rna.matureTRNAIndexs); rna.aminoacylatedIndexs(rna.matureTRNAIndexs)] comp.cytosolIndexs
                'ProteinMonomer' 'counts' [pm.matureIndexs(monIdxs) pm.boundIndexs(monIdxs)] comp.cytosolIndexs
                'ProteinComplex' 'counts' [pc.matureIndexs(cpxIdxs) pc.boundIndexs(cpxIdxs)] comp.cytosolIndexs
                };
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations);
            
            tmp = full(states.Metabolite.counts);
            atp = permute(tmp(1, :, :, :), [4 3 1 2]);
            gtp = permute(tmp(2, :, :, :), [4 3 1 2]);
            adp = permute(tmp(3, :, :, :), [4 3 1 2]);
            gdp = permute(tmp(4, :, :, :), [4 3 1 2]);
            amp = permute(tmp(5, :, :, :), [4 3 1 2]);
            aas = permute(sum(tmp(6:26, :, :, :), 1), [4 3 1 2]);
            fthf10 = permute(sum(tmp(27:end, :, :, :), 1), [4 3 1 2]);
            
            monCnts = full(states.ProteinMonomer.counts);
            cpxCnts = full(states.ProteinComplex.counts);
            
            enzProcesses = [
                repmat({'TL'}, numel(tl.enzymeWholeCellModelIDs), 1)
                repmat({'TA'}, numel(ta.enzymeWholeCellModelIDs)-numel(ta.enzymeRNALocalIndexs), 1)
                ];
            enzProcessIdxs = [
                (1:numel(tl.enzymeWholeCellModelIDs))'
                (1:numel(ta.enzymeWholeCellModelIDs)-numel(ta.enzymeRNALocalIndexs))'
                ];
            
            enzymes = zeros(0, 1, size(monCnts, 3), size(monCnts, 4));
            expEnzymes = zeros(0, 1);
            
            tmp = zeros(numel(tl.enzymes), 1, size(monCnts, 3), size(monCnts, 4));
            tmp(tl.enzymeMonomerLocalIndexs, :, :, :) = monCnts(1:numel(tl.enzymeMonomerGlobalIndexs), :, :, :);
            tmp(tl.enzymeComplexLocalIndexs, :, :, :) = cpxCnts(1:numel(tl.enzymeComplexGlobalIndexs), :, :, :);
            tmp2 = zeros(numel(tl.enzymes), 1);
            tmp2(tl.enzymeMonomerLocalIndexs) = expTotMons(tl.enzymeMonomerGlobalIndexs);
            tmp2(tl.enzymeComplexLocalIndexs) = expTotCpxs(tl.enzymeComplexGlobalIndexs);
            enzymes = [enzymes; tmp];
            expEnzymes = [expEnzymes; tmp2];
            
            tmp = zeros(numel(ta.enzymes), 1, size(monCnts, 3), size(monCnts, 4));
            tmp(ta.enzymeMonomerLocalIndexs, :, :, :) = monCnts(end-numel(ta.enzymeMonomerGlobalIndexs)+1:end, :, :, :);
            tmp(ta.enzymeComplexLocalIndexs, :, :, :) = cpxCnts(end-numel(ta.enzymeComplexGlobalIndexs)+1:end, :, :, :);
            tmp2 = zeros(numel(ta.enzymes), 1);
            tmp2(ta.enzymeMonomerLocalIndexs) = expTotMons(ta.enzymeMonomerGlobalIndexs);
            tmp2(ta.enzymeComplexLocalIndexs) = expTotCpxs(ta.enzymeComplexGlobalIndexs);
            tmp(ta.enzymeRNALocalIndexs, :, :, :) = [];
            tmp2(ta.enzymeRNALocalIndexs) = [];
            enzymes = [enzymes; tmp];
            expEnzymes = [expEnzymes; tmp2];
            
            expEnzymes = repmat(expEnzymes, [1 1 numel(time)-1 1]) .* exp(log(2) * repmat(permute(time(2:end),[1 3 2]), size(expEnzymes, 1), 1) / (t.cellCycleLength/3600));
            
            trnas = states.Rna.counts(1:end/2, :, :, :);
            aatrnas = states.Rna.counts(end/2+1:end, :, :, :);
            expTRNAs = repmat(expTotRnas(rna.matureTRNAIndexs), [1 1 numel(time)-1 1]) .* exp(log(2) * repmat(permute(time(2:end),[1 3 2]), size(rna.matureTRNAIndexs, 1), 1) / (t.cellCycleLength/3600));
            
            for i = 1:numel(simEndTimes)
                atp(i, simEndTimes(i)+1:end) = NaN;
                gtp(i, simEndTimes(i)+1:end) = NaN;
                adp(i, simEndTimes(i)+1:end) = NaN;
                gdp(i, simEndTimes(i)+1:end) = NaN;
                amp(i, simEndTimes(i)+1:end) = NaN;
                aas(i, simEndTimes(i)+1:end) = NaN;
                fthf10(i, simEndTimes(i)+1:end) = NaN;
                enzymes(:, :, simEndTimes(i)+1:end, i) = NaN;
                trnas(:, :, simEndTimes(i)+1:end, i) = NaN;
                aatrnas(:, :, simEndTimes(i)+1:end, i) = NaN;
            end
            
            %cleanup
            clear states rnaCnts monCnts cpxCnts;
            
            %ribosomes
            pols = NaN(4, 1, numel(time)-1, numel(selectedSimulations));
            stateNames = {
                'Ribosome'  'states'          ':' ':'
                'Ribosome'  'boundMRNAs'      ':' ':'
                'Ribosome'  'mRNAPositions'   ':' ':'
                'Ribosome'  'tmRNAPositions'  ':' ':'
                };
            for j = 1:numel(selectedSimulations)
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations(j));
                
                for i = 1:size(states.Ribosome.states, 3)
                    tfs = states.Ribosome.states(:, :, i) == rib.activeValue | states.Ribosome.states(:, :, i) == rib.stalledValue;
                    boundMRNAs = states.Ribosome.boundMRNAs(tfs, :, i);
                    lengths1 = states.Ribosome.mRNAPositions(tfs, :, i);
                    lengths2 = states.Ribosome.tmRNAPositions(tfs, :, i);
                    
                    pols(1, 1, i, j) = sum(pm.lengths(pm.nascentIndexs(boundMRNAs)) == lengths1 | lengths2 == pol.proteolysisTagLength);
                    pols(2, 1, i, j) = sum(lengths1 == 0 & lengths2 == 0);
                    pols(4, 1, i, j) = sum(states.Ribosome.states(tfs, :, i) == rib.stalledValue & lengths2 ~= pol.proteolysisTagLength);
                    pols(3, 1, i, j) = sum(tfs) - sum(pols([1 2 4], 1, i, j));
                end
                
                clear states tfs boundMRNAs lengths1 lengths2;
            end
            
            %% clear figure
            clf(figHandle);
            C = 70;
            
            %% left panel
            
            %layout plots
            h = (C - 10 + 1) / 10;
            axesHandlesL = PlotUtil.multiElementPlot(figHandle, h * ones(11, 1), [0 time(end)], struct(...
                'titleStr', 'Metabolites', ...
                'position', [0.10 0.10 0.20 0.8], ...
                'colorOrder', PlotUtil.getRedGreenColorOrder(simOrder)));
            
            %plot data
            PlotUtil.plotLine(axesHandlesL(1), time(2:end), atp, false, true, false);
            ylabel(axesHandlesL(1), 'ATP');
            
            PlotUtil.plotLine(axesHandlesL(2), time(2:end), gtp, false, true, false);
            ylabel(axesHandlesL(2), 'GTP');
            
            PlotUtil.plotLine(axesHandlesL(3), time(2:end), adp, false, true, false);
            ylabel(axesHandlesL(3), 'ADP');
            
            PlotUtil.plotLine(axesHandlesL(4), time(2:end), gdp, false, true, false);
            ylabel(axesHandlesL(4), 'GDP');
            
            PlotUtil.plotLine(axesHandlesL(5), time(2:end), amp, false, true, false);
            ylabel(axesHandlesL(5), 'AMP');
            
            PlotUtil.plotLine(axesHandlesL(6), time(2:end), aas, false, true, false);
            ylabel(axesHandlesL(6), 'AAs');
            
            PlotUtil.plotLine(axesHandlesL(7), time(2:end), fthf10, false, true, false);
            ylabel(axesHandlesL(7), 'FTHF10');
            
            PlotUtil.plotLine(axesHandlesL(8), time(2:end), permute(pols(2, :, :, :), [4 3 1 2]), false, true, false);
            ylabel(axesHandlesL(8), {'Init' 'Ribs'});
            
            PlotUtil.plotLine(axesHandlesL(9), time(2:end), permute(pols(3, :, :, :), [4 3 1 2]), false, true, false);
            ylabel(axesHandlesL(9), {'Elng' 'Ribs'});
            
            PlotUtil.plotLine(axesHandlesL(10), time(2:end), permute(pols(1, :, :, :), [4 3 1 2]), false, true, false);
            ylabel(axesHandlesL(10), {'Term' 'Ribs'});
            
            PlotUtil.plotLine(axesHandlesL(11), time(2:end), permute(pols(4, :, :, :), [4 3 1 2]), false, true, false);
            ylabel(axesHandlesL(11), {'Stalled' 'Ribs'});
            
            % Format Y-axis
            PlotUtil.alignYAxesLabels(axesHandlesL);
            PlotUtil.offsetYAxes(axesHandlesL, 0.04);
            PlotUtil.labelSubplots(axesHandlesL(1), 'A', -0.42, 1.2);
            
            %% middle panels
            
            nPlots1 = ceil(size(enzymes, 1) / 2);
            nPlots2 = floor(size(enzymes, 1) / 2);
            
            %layout plots
            colorOrder = [PlotUtil.getRedGreenColorOrder(simOrder); 0.5 0.5 0.5];
            h = (C - ceil(size(enzymes, 1)/2) + 1) / ceil(size(enzymes, 1) / 2);
            axesHandlesM1 = PlotUtil.multiElementPlot(figHandle, h * ones(nPlots1, 1), [0 time(end)], struct(...
                'titleStr', 'Enzymes', ...
                'position', [0.36 0.10 0.10 0.8], ...
                'colorOrder', colorOrder));
            axesHandlesM2 = PlotUtil.multiElementPlot(figHandle, h * ones(nPlots1, 1), [0 time(end)], struct(...
                'position', [0.53 0.10 0.10 0.8], ...
                'colorOrder', colorOrder));
            
            %plot data
            for i = 1:nPlots1
                PlotUtil.plotLine(axesHandlesM1(i), time(2:end), permute(cat(4, enzymes(i, :, :, :), expEnzymes(i, :, :, :)), [4 3 1 2]), false, true, false);
                ylabel(axesHandlesM1(i), [enzProcesses(i) num2str(enzProcessIdxs(i))]);
            end
            for i = 1:nPlots2
                PlotUtil.plotLine(axesHandlesM2(i), time(2:end), permute(cat(4, enzymes(nPlots1+i, :, :, :), expEnzymes(nPlots1+i, :, :, :)), [4 3 1 2]), false, true, false);
                ylabel(axesHandlesM2(i), [enzProcesses(nPlots1+i) num2str(enzProcessIdxs(nPlots1+i))]);
            end
            
            % Format Y-axis
            PlotUtil.alignYAxesLabels(axesHandlesM1);
            PlotUtil.alignYAxesLabels(axesHandlesM2);
            axesHandlesM12 = PlotUtil.offsetYAxes(axesHandlesM1, 0.04);
            axesHandlesM22 = PlotUtil.offsetYAxes(axesHandlesM2, 0.04);
            PlotUtil.labelSubplots(axesHandlesM1(1), 'B', -0.42, 1.2);
            
            set(axesHandlesM1(nPlots1+1:end), 'visible', 'off')
            set(axesHandlesM2(nPlots2+1:end), 'visible', 'off')
            set(axesHandlesM12(nPlots1+1:end), 'visible', 'off')
            set(axesHandlesM22(nPlots2+1:end), 'visible', 'off')
            
            %% right panels
            nPlots1 = ceil(size(trnas, 1) / 2);
            nPlots2 = floor(size(trnas, 1) / 2);
            
            %layout plots
            colorOrder = [PlotUtil.getRedGreenColorOrder(simOrder); PlotUtil.getRedGreenColorOrder(simOrder); 0.5 0.5 0.5];
            h = (C - ceil(size(trnas, 1)/2) + 1) / ceil(size(trnas, 1) / 2);
            axesHandlesR1 = PlotUtil.multiElementPlot(figHandle, h * ones(nPlots1, 1), [0 time(end)], struct(...
                'titleStr', 'tRNA', ...
                'position', [0.70 0.10 0.10 0.8], ...
                'colorOrder', colorOrder));
            axesHandlesR2 = PlotUtil.multiElementPlot(figHandle, h * ones(nPlots1, 1), [0 time(end)], struct(...
                'position', [0.87 0.10 0.10 0.8], ...
                'colorOrder', colorOrder));
            
            %plot data
            for i = 1:nPlots1
                PlotUtil.plotLine(axesHandlesR1(i), time(2:end), permute(cat(4, trnas(i, :, :, :), aatrnas(i, :, :, :), expTRNAs(i, :, :, :)), [4 3 1 2]), false, true, false);
                ylabel(axesHandlesR1(i), num2str(i));
            end
            for i = 1:nPlots2
                PlotUtil.plotLine(axesHandlesR2(i), time(2:end), permute(cat(4, trnas(nPlots1+i, :, :, :), aatrnas(nPlots1+i, :, :, :), expTRNAs(nPlots1+i, :, :, :)), [4 3 1 2]), false, true, false);
                ylabel(axesHandlesR2(i), num2str(nPlots1+i));
            end
            
            % Format Y-axis
            PlotUtil.alignYAxesLabels(axesHandlesR1);
            PlotUtil.alignYAxesLabels(axesHandlesR2);
            axesHandlesR12 = PlotUtil.offsetYAxes(axesHandlesR1, 0.04);
            axesHandlesR22 = PlotUtil.offsetYAxes(axesHandlesR2, 0.04);
            PlotUtil.labelSubplots(axesHandlesR1(1), 'C', -0.42, 1.2);
            
            set(axesHandlesR1(nPlots1+1:end), 'visible', 'off')
            set(axesHandlesR2(nPlots2+1:end), 'visible', 'off')
            set(axesHandlesR12(nPlots1+1:end), 'visible', 'off')
            set(axesHandlesR22(nPlots2+1:end), 'visible', 'off')
            
            %% organize figure data
            figData = struct;
            figData.selectedSimulations = selectedSimulations;
            figData.simOrder = simOrder;
            figData.time = time;
            figData.atp = atp;
            figData.gtp = gtp;
            figData.adp = adp;
            figData.gdp = gdp;
            figData.amp = amp;
            figData.aas = aas;
            figData.fthf10 = fthf10;
            figData.enzProcesses = enzProcesses;
            figData.enzProcessIdxs = enzProcessIdxs;
            figData.enzymes = enzymes;
            figData.expEnzymes = expEnzymes;
            figData.trnas = trnas;
            figData.aatrnas = aatrnas;
            figData.expTRNAs = expTRNAs;
            figData.pols = pols;
        end
        
        function figData = proteinSynthesis(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            comp = sim.compartment;
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            ppI = sim.process('ProteinProcessingI');
            ppII = sim.process('ProteinProcessingII');
            pt = sim.process('ProteinTranslocation');
            pf = sim.process('ProteinFolding');
            pmod = sim.process('ProteinModification');
            pd = sim.process('ProteinDecay');
            
            enzProcesses = [
                repmat({'PP-I'}, size(ppI.enzymeWholeCellModelIDs))
                repmat({'PT'}, size(pt.enzymeWholeCellModelIDs))
                repmat({'PP-II'}, size(ppII.enzymeWholeCellModelIDs))
                repmat({'PF'}, size(pf.enzymeWholeCellModelIDs))
                repmat({'PMod'}, size(pmod.enzymeWholeCellModelIDs))
                repmat({'PD'}, size(pd.enzymeWholeCellModelIDs))
                ];
            enzProcessIdxs = [
                (1:size(ppI.enzymeWholeCellModelIDs))'
                (1:size(pt.enzymeWholeCellModelIDs))'
                (1:size(ppII.enzymeWholeCellModelIDs))'
                (1:size(pf.enzymeWholeCellModelIDs))'
                (1:size(pmod.enzymeWholeCellModelIDs))'
                (1:size(pd.enzymeWholeCellModelIDs))'
                ];
            
            %% get data
            ensemble = SimulationEnsemble(simBatchDir, {'mass' 'mass'; 'immatureMonomers' 'immatureMonomers'}, [], selectedSimulations);
            time = ensemble.stateData.time / 3600;
            mass = ensemble.stateData.values(ensemble.getPropertyIndices('mass'), :, :, :) * sim.state('Mass').cellInitialDryWeight / (1-sim.state('Mass').fractionWetWeight) * 1e15;
            simEndTimes = ensemble.stateData.simulationEndTimes;
            
            [~, simOrder] = sort(mass(:, 1, min(simEndTimes), :));
            simOrder = simOrder(:);
            
            clear ensemble;
            
            stateNames = {
                'Polypeptide'    'abortedPolypeptides'  '-nnz'           1
                'ProteinMonomer' 'counts'               ':'              ':'
                'ProteinComplex' 'counts'               pc.matureIndexs  '-sum'
                };
            nscnt = NaN(numel(selectedSimulations), max(simEndTimes));
            procsI = NaN(numel(selectedSimulations), max(simEndTimes));
            trans = NaN(numel(selectedSimulations), max(simEndTimes));
            procsII = NaN(numel(selectedSimulations), max(simEndTimes));
            folds = NaN(numel(selectedSimulations), max(simEndTimes));
            signals = NaN(numel(selectedSimulations), max(simEndTimes));
            inacts = NaN(numel(selectedSimulations), max(simEndTimes));
            misfolds = NaN(numel(selectedSimulations), max(simEndTimes));
            damages = NaN(numel(selectedSimulations), max(simEndTimes));
            aborts = NaN(numel(selectedSimulations), max(simEndTimes));
            enzymes = NaN(numel(enzProcesses), 1, max(simEndTimes), numel(selectedSimulations));
            for i = 1:numel(selectedSimulations)
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations(i));
                
                monCnts = full(sum(states.ProteinMonomer.counts, 2));
                nscnt(i, 1:simEndTimes(i)) = permute(sum(monCnts(pm.nascentIndexs, :, :, :), 1), [4 3 1 2]);
                procsI(i, 1:simEndTimes(i)) = full(permute(sum(states.ProteinMonomer.counts(...
                    pm.processedIIndexs(pm.compartments(pm.processedIIIndexs) ~= comp.cytosolIndexs), ...
                    comp.cytosolIndexs, :, :), 1), [4 3 1 2]));
                trans(i, 1:simEndTimes(i)) = permute(sum(monCnts(pm.processedIIndexs, :, :, :), 1), [4 3 1 2]) - procsI(i, 1:simEndTimes(i));
                procsII(i, 1:simEndTimes(i)) = permute(sum(monCnts(pm.processedIIIndexs, :, :, :), 1), [4 3 1 2]);
                folds(i, 1:simEndTimes(i)) = permute(sum(monCnts(pm.foldedIndexs, :, :, :), 1), [4 3 1 2]);
                signals(i, 1:simEndTimes(i)) = permute(sum(monCnts(pm.signalSequenceIndexs, :, :, :), 1), [4 3 1 2]);
                inacts(i, 1:simEndTimes(i)) = permute(sum(monCnts(pm.inactivatedIndexs, :, :, :), 1), [4 3 1 2]);
                misfolds(i, 1:simEndTimes(i)) = permute(sum(monCnts(pm.misfoldedIndexs, :, :, :), 1), [4 3 1 2]);
                damages(i, 1:simEndTimes(i)) = permute(sum(monCnts(pm.damagedIndexs, :, :, :), 1), [4 3 1 2]);
                aborts(i, 1:simEndTimes(i)) = permute(full(states.Polypeptide.abortedPolypeptides), [4 3 1 2]);
                
                monCnts = monCnts(pm.matureIndexs, :, :, :);
                cpxCnts = full(states.ProteinComplex.counts);
                
                tmpEnzymes = zeros(0, 1, size(monCnts, 3), size(monCnts, 4));
                
                tmp = zeros(numel(ppI.enzymes), 1, size(monCnts, 3), size(monCnts, 4));
                tmp(ppI.enzymeMonomerLocalIndexs, :, :, :) = monCnts(ppI.enzymeMonomerGlobalIndexs, :, :, :);
                tmp(ppI.enzymeComplexLocalIndexs, :, :, :) = cpxCnts(ppI.enzymeComplexGlobalIndexs, :, :, :);
                tmpEnzymes = [tmpEnzymes; tmp]; %#ok<AGROW>
                
                tmp = zeros(numel(pt.enzymes), 1, size(monCnts, 3), size(monCnts, 4));
                tmp(pt.enzymeMonomerLocalIndexs, :, :, :) = monCnts(pt.enzymeMonomerGlobalIndexs, :, :, :);
                tmp(pt.enzymeComplexLocalIndexs, :, :, :) = cpxCnts(pt.enzymeComplexGlobalIndexs, :, :, :);
                tmpEnzymes = [tmpEnzymes; tmp]; %#ok<AGROW>
                
                tmp = zeros(numel(ppII.enzymes), 1, size(monCnts, 3), size(monCnts, 4));
                tmp(ppII.enzymeMonomerLocalIndexs, :, :, :) = monCnts(ppII.enzymeMonomerGlobalIndexs, :, :, :);
                tmp(ppII.enzymeComplexLocalIndexs, :, :, :) = cpxCnts(ppII.enzymeComplexGlobalIndexs, :, :, :);
                tmpEnzymes = [tmpEnzymes; tmp]; %#ok<AGROW>
                
                tmp = zeros(numel(pf.enzymes), 1, size(monCnts, 3), size(monCnts, 4));
                tmp(pf.enzymeMonomerLocalIndexs, :, :, :) = monCnts(pf.enzymeMonomerGlobalIndexs, :, :, :);
                tmp(pf.enzymeComplexLocalIndexs, :, :, :) = cpxCnts(pf.enzymeComplexGlobalIndexs, :, :, :);
                tmpEnzymes = [tmpEnzymes; tmp]; %#ok<AGROW>
                
                tmp = zeros(numel(pmod.enzymes), 1, size(monCnts, 3), size(monCnts, 4));
                tmp(pmod.enzymeMonomerLocalIndexs, :, :, :) = monCnts(pmod.enzymeMonomerGlobalIndexs, :, :, :);
                tmp(pmod.enzymeComplexLocalIndexs, :, :, :) = cpxCnts(pmod.enzymeComplexGlobalIndexs, :, :, :);
                tmpEnzymes = [tmpEnzymes; tmp]; %#ok<AGROW>
                
                tmp = zeros(numel(pd.enzymes), 1, size(monCnts, 3), size(monCnts, 4));
                tmp(pd.enzymeMonomerLocalIndexs, :, :, :) = monCnts(pd.enzymeMonomerGlobalIndexs, :, :, :);
                tmp(pd.enzymeComplexLocalIndexs, :, :, :) = cpxCnts(pd.enzymeComplexGlobalIndexs, :, :, :);
                tmpEnzymes = [tmpEnzymes; tmp]; %#ok<AGROW>
                
                enzymes(:, :, 1:simEndTimes(i), i) = tmpEnzymes;
                
                clear states monCnts cpxCnts tmpEnzymes tmp;
            end
            
            %% clear figure
            clf(figHandle);
            
            %% left panel
            
            %layout plots
            axesHandlesL = PlotUtil.multiElementPlot(figHandle, 3.5 * ones(10, 1), [0 time(end)], struct(...
                'titleStr', 'Monomers', ...
                'position', [0.12 0.10 0.21 0.8], ...
                'colorOrder', PlotUtil.getRedGreenColorOrder(simOrder)));
            
            %plot data
            PlotUtil.plotLine(axesHandlesL(1), time(2:end), nscnt, false, true, false);
            ylabel(axesHandlesL(1), 'Nascent');
            
            PlotUtil.plotLine(axesHandlesL(2), time(2:end), procsI, false, true, false);
            ylabel(axesHandlesL(2), 'Processed-I');
            
            PlotUtil.plotLine(axesHandlesL(3), time(2:end), trans, false, true, false);
            ylabel(axesHandlesL(3), 'Translocated');
            
            PlotUtil.plotLine(axesHandlesL(4), time(2:end), procsII, false, true, false);
            ylabel(axesHandlesL(4), 'Processed-II');
            
            PlotUtil.plotLine(axesHandlesL(5), time(2:end), folds, false, true, false);
            ylabel(axesHandlesL(5), 'Folded');
            
            PlotUtil.plotLine(axesHandlesL(6), time(2:end), signals, false, true, false);
            ylabel(axesHandlesL(6), {'Signal' 'Sequences'});
            
            PlotUtil.plotLine(axesHandlesL(7), time(2:end), inacts, false, true, false);
            ylabel(axesHandlesL(7), 'Inactivated');
            
            PlotUtil.plotLine(axesHandlesL(8), time(2:end), misfolds, false, true, false);
            ylabel(axesHandlesL(8), 'Misfolded');
            
            PlotUtil.plotLine(axesHandlesL(9), time(2:end), damages, false, true, false);
            ylabel(axesHandlesL(9), 'Damaged');
            
            PlotUtil.plotLine(axesHandlesL(10), time(2:end), aborts, false, true, false);
            ylabel(axesHandlesL(10), 'Aborted');
            
            % Format Y-axis
            PlotUtil.alignYAxesLabels(axesHandlesL);
            PlotUtil.offsetYAxes(axesHandlesL, 0.04);
            PlotUtil.labelSubplots(axesHandlesL(1), 'A', -0.42, 1.2);
            
            %% right panels
            
            nPlots1 = ceil(size(enzymes, 1) / 2);
            nPlots2 = floor(size(enzymes, 1) / 2);
            
            %layout plots
            axesHandlesM = PlotUtil.multiElementPlot(figHandle, 3.5 * ones(nPlots1, 1), [0 time(end)], struct(...
                'titleStr', 'Enzymes', ...
                'position', [0.45 0.10 0.21 0.8], ...
                'colorOrder', PlotUtil.getRedGreenColorOrder(simOrder)));
            axesHandlesR = PlotUtil.multiElementPlot(figHandle, 3.5 * ones(nPlots2, 1), [0 time(end)], struct(...
                'position', [0.78 0.10 0.21 0.8], ...
                'colorOrder', PlotUtil.getRedGreenColorOrder(simOrder)));
            
            %plot data
            for i = 1:nPlots1
                PlotUtil.plotLine(axesHandlesM(i), time(2:end), permute(enzymes(i, :, :, :), [4 3 1 2]), false, true, false);
                ylabel(axesHandlesM(i), [enzProcesses(i) num2str(enzProcessIdxs(i))]);
            end
            for i = 1:nPlots2
                PlotUtil.plotLine(axesHandlesR(i), time(2:end), permute(enzymes(nPlots1+i, :, :, :), [4 3 1 2]), false, true, false);
                ylabel(axesHandlesR(i), [enzProcesses(nPlots1+i) num2str(enzProcessIdxs(nPlots1+i))]);
            end
            
            % Format Y-axis
            PlotUtil.alignYAxesLabels(axesHandlesM);
            PlotUtil.alignYAxesLabels(axesHandlesR);
            PlotUtil.offsetYAxes(axesHandlesM, 0.04);
            PlotUtil.offsetYAxes(axesHandlesR, 0.04);
            PlotUtil.labelSubplots(axesHandlesM(1), 'B', -0.42, 1.2);
            
            %% organize figure data
            figData = struct;
            figData.simOrder = simOrder;
            figData.time = time;
            figData.nscnt = nscnt;
            figData.procsI = procsI;
            figData.trans = trans;
            figData.procsII = procsII;
            figData.folds = folds;
            figData.signals = signals;
            figData.inacts = inacts;
            figData.misfolds = misfolds;
            figData.damages = damages;
            figData.aborts = aborts;
            figData.enzymes = enzymes;
            figData.enzProcesses = enzProcesses;
        end
        
        function figData = rnaSynthesis(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            g = sim.gene;
            
            t = sim.state('Time');
            mass = sim.state('Mass');
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            rp = sim.process('RNAProcessing');
            rmod = sim.process('RNAModification');
            rd = sim.process('RNADecay');
            
            monExp = (rna.matureRNAGeneComposition(g.mRNAIndexs, :) * rna.expression(rna.matureIndexs)) ./ ...
                (log(2) / t.cellCycleLength + pm.decayRates(pm.matureIndexs));
            expTotMons = mass.initialFractionAAsInMonomers * mass.dryWeightFractionProtein * mass.cellInitialDryWeight * ...
                monExp / (monExp' * pm.molecularWeights(pm.matureIndexs) / ConstantUtil.nAvogadro);
            expTotCpxs = min(repmat(expTotMons, [1 numel(pc.matureIndexs)]) ./ sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3), [], 1);
            
            enzProcesses = [
                repmat({'RP'}, size(rp.enzymeWholeCellModelIDs))
                repmat({'RMod'}, size(rmod.enzymeWholeCellModelIDs))
                repmat({'RD'}, size(rd.enzymeWholeCellModelIDs))
                ];
            enzProcessIdxs = [
                (1:size(rp.enzymeWholeCellModelIDs))'
                (1:size(rmod.enzymeWholeCellModelIDs))'
                (1:size(rd.enzymeWholeCellModelIDs))'
                ];
            
            %% get data
            ensemble = SimulationEnsemble(simBatchDir, {'mass' 'mass'; 'immatureRnas' 'immatureRnas'}, [], selectedSimulations);
            time = ensemble.stateData.time / 3600;
            mass = ensemble.stateData.values(ensemble.getPropertyIndices('mass'), :, :, :) * sim.state('Mass').cellInitialDryWeight / (1-sim.state('Mass').fractionWetWeight) * 1e15;
            simEndTimes = ensemble.stateData.simulationEndTimes();
            
            [~, simOrder] = sort(mass(:, 1, min(simEndTimes), :));
            simOrder = simOrder(:);
            
            clear ensemble;
            
            stateNames = {
                'Rna'            'counts'  ':'               '-sum'
                'Transcript'     'abortedTranscripts' '-nnz'  1
                'ProteinMonomer' 'counts'  pm.matureIndexs  '-sum'
                'ProteinComplex' 'counts'  pc.matureIndexs  '-sum'
                };
            nscnt = NaN(numel(selectedSimulations), max(simEndTimes));
            procs = NaN(numel(selectedSimulations), max(simEndTimes));
            misfolds = NaN(numel(selectedSimulations), max(simEndTimes));
            damages = NaN(numel(selectedSimulations), max(simEndTimes));
            aborts = NaN(numel(selectedSimulations), max(simEndTimes));
            enzymes = NaN(numel(enzProcessIdxs), 1, max(simEndTimes), numel(selectedSimulations));
            for i = 1:numel(selectedSimulations)
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations(i));
                
                rnaCnts = full(states.Rna.counts);
                monCnts = full(states.ProteinMonomer.counts);
                cpxCnts = full(states.ProteinComplex.counts);
                nscnt(i, 1:simEndTimes(i)) = permute(sum(rnaCnts(rna.nascentIndexs, :, :, :), 1), [4 3 1 2]);
                procs(i, 1:simEndTimes(i)) = permute(sum(rnaCnts(rna.processedIndexs, :, :, :), 1), [4 3 1 2]);
                misfolds (i, 1:simEndTimes(i))= permute(sum(rnaCnts(rna.misfoldedIndexs, :, :, :), 1), [4 3 1 2]);
                damages(i, 1:simEndTimes(i)) = permute(sum(rnaCnts(rna.damagedIndexs, :, :, :), 1), [4 3 1 2]);
                aborts(i, 1:simEndTimes(i)) = permute(full(states.Transcript.abortedTranscripts), [4 3 1 2]);
                
                tmpEnzymes = zeros(0, 1, size(monCnts, 3), size(monCnts, 4));
                
                tmp = zeros(numel(rp.enzymes), 1, size(monCnts, 3), size(monCnts, 4));
                tmp(rp.enzymeRNALocalIndexs, :, :, :) = rnaCnts(numel(rna.nascentIndexs)+numel(rna.processedIndexs) + rp.enzymeRNAGlobalIndexs, :, :, :);
                tmp(rp.enzymeMonomerLocalIndexs, :, :, :) = monCnts(rp.enzymeMonomerGlobalIndexs, :, :, :);
                tmp(rp.enzymeComplexLocalIndexs, :, :, :) = cpxCnts(rp.enzymeComplexGlobalIndexs, :, :, :);
                tmpEnzymes = [tmpEnzymes; tmp]; %#ok<AGROW>
                
                tmp = zeros(numel(rmod.enzymes), 1, size(monCnts, 3), size(monCnts, 4));
                tmp(rmod.enzymeRNALocalIndexs, :, :, :) = rnaCnts(numel(rna.nascentIndexs)+numel(rna.processedIndexs) + rmod.enzymeRNAGlobalIndexs, :, :, :);
                tmp(rmod.enzymeMonomerLocalIndexs, :, :, :) = monCnts(rmod.enzymeMonomerGlobalIndexs, :, :, :);
                tmp(rmod.enzymeComplexLocalIndexs, :, :, :) = cpxCnts(rmod.enzymeComplexGlobalIndexs, :, :, :);
                tmpEnzymes = [tmpEnzymes; tmp]; %#ok<AGROW>
                
                tmp = zeros(numel(rd.enzymes), 1, size(monCnts, 3), size(monCnts, 4));
                tmp(rd.enzymeRNALocalIndexs, :, :, :) = rnaCnts(numel(rna.nascentIndexs)+numel(rna.processedIndexs) + rd.enzymeRNAGlobalIndexs, :, :, :);
                tmp(rd.enzymeMonomerLocalIndexs, :, :, :) = monCnts(rd.enzymeMonomerGlobalIndexs, :, :, :);
                tmp(rd.enzymeComplexLocalIndexs, :, :, :) = cpxCnts(rd.enzymeComplexGlobalIndexs, :, :, :);
                tmpEnzymes = [tmpEnzymes; tmp]; %#ok<AGROW>
                
                enzymes(:, :, 1:simEndTimes(i), i) = tmpEnzymes;
                
                clear states rnaCnts monCnts cpxCnts tmpEnzymes tmp;
            end
            
            expEnzymes = zeros(0, 1);
            
            tmp2 = zeros(numel(rp.enzymes), 1);
            tmp2(rp.enzymeMonomerLocalIndexs) = expTotMons(rp.enzymeMonomerGlobalIndexs);
            tmp2(rp.enzymeComplexLocalIndexs) = expTotCpxs(rp.enzymeComplexGlobalIndexs);
            expEnzymes = [expEnzymes; tmp2];
            
            tmp2 = zeros(numel(rmod.enzymes), 1);
            tmp2(rmod.enzymeMonomerLocalIndexs) = expTotMons(rmod.enzymeMonomerGlobalIndexs);
            tmp2(rmod.enzymeComplexLocalIndexs) = expTotCpxs(rmod.enzymeComplexGlobalIndexs);
            expEnzymes = [expEnzymes; tmp2];
            
            tmp2 = zeros(numel(rd.enzymes), 1);
            tmp2(rd.enzymeMonomerLocalIndexs) = expTotMons(rd.enzymeMonomerGlobalIndexs);
            tmp2(rd.enzymeComplexLocalIndexs) = expTotCpxs(rd.enzymeComplexGlobalIndexs);
            expEnzymes = [expEnzymes; tmp2];
            
            expEnzymes = repmat(expEnzymes, [1 1 numel(time)-1 1]) .* exp(log(2) * repmat(permute(time(2:end),[1 3 2]), size(expEnzymes, 1), 1) / (t.cellCycleLength/3600));
            
            %% clear figure
            clf(figHandle);
            C = 70;
            
            %% left panel
            
            %layout plots
            h = (C - 2 + 1) / 2;
            axesHandlesL = PlotUtil.multiElementPlot(figHandle, h* ones(5, 1), [0 time(end)], struct(...
                'titleStr', 'RNAs', ...
                'position', [0.12 0.10 0.21 0.8], ...
                'colorOrder', PlotUtil.getRedGreenColorOrder(simOrder)));
            
            %plot data
            PlotUtil.plotLine(axesHandlesL(1), time(2:end), nscnt, false, true, false);
            ylabel(axesHandlesL(1), 'Nascent');
            
            PlotUtil.plotLine(axesHandlesL(2), time(2:end), procs, false, true, false);
            ylabel(axesHandlesL(2), 'Processed');
            
            PlotUtil.plotLine(axesHandlesL(3), time(2:end), misfolds, false, true, false);
            ylabel(axesHandlesL(3), 'Misfolded');
            
            PlotUtil.plotLine(axesHandlesL(4), time(2:end), damages, false, true, false);
            ylabel(axesHandlesL(4), 'Damaged');
            
            PlotUtil.plotLine(axesHandlesL(5), time(2:end), aborts, false, true, false);
            ylabel(axesHandlesL(5), 'Aborted');
            
            % Format Y-axis
            PlotUtil.alignYAxesLabels(axesHandlesL);
            PlotUtil.offsetYAxes(axesHandlesL, 0.04);
            PlotUtil.labelSubplots(axesHandlesL(1), 'A', -0.42, 1.2);
            
            %% right panels
            
            nPlots1 = ceil(size(enzymes, 1) / 2);
            nPlots2 = floor(size(enzymes, 1) / 2);
            
            %layout plots
            colorOrder = [PlotUtil.getRedGreenColorOrder(simOrder); 0.5 0.5 0.5];
            h = (C - ceil(size(enzymes, 1)/2) + 1) / ceil(size(enzymes, 1) / 2);
            axesHandlesM = PlotUtil.multiElementPlot(figHandle, h * ones(nPlots1, 1), [0 time(end)], struct(...
                'titleStr', 'Enzymes', ...
                'position', [0.45 0.10 0.21 0.8], ...
                'colorOrder', colorOrder));
            axesHandlesR = PlotUtil.multiElementPlot(figHandle, h * ones(nPlots2, 1), [0 time(end)], struct(...
                'position', [0.78 0.10 0.21 0.8], ...
                'colorOrder', colorOrder));
            
            %plot data
            for i = 1:nPlots1
                PlotUtil.plotLine(axesHandlesM(i), time(2:end), permute(cat(4, enzymes(i, :, :, :), expEnzymes(i, :, :, :)), [4 3 1 2]), false, true, false);
                ylabel(axesHandlesM(i), [enzProcesses(i) num2str(enzProcessIdxs(i))]);
            end
            for i = 1:nPlots2
                PlotUtil.plotLine(axesHandlesR(i), time(2:end), permute(cat(4, enzymes(nPlots1+i, :, :, :), expEnzymes(nPlots1+i, :, :, :)), [4 3 1 2]), false, true, false);
                ylabel(axesHandlesR(i), [enzProcesses(nPlots1+i) num2str(enzProcessIdxs(nPlots1+i))]);
            end
            
            % Format Y-axis
            PlotUtil.alignYAxesLabels(axesHandlesM);
            PlotUtil.alignYAxesLabels(axesHandlesR);
            PlotUtil.offsetYAxes(axesHandlesM, 0.04);
            PlotUtil.offsetYAxes(axesHandlesR, 0.04);
            PlotUtil.labelSubplots(axesHandlesM(1), 'B', -0.42, 1.2);
            
            %% organize figure data
            figData = struct;
            figData.simOrder = simOrder;
            figData.time = time;
            figData.nscnt = nscnt;
            figData.procs = procs;
            figData.misfolds = misfolds;
            figData.damages = damages;
            figData.aborts = aborts;
            figData.enzymes = enzymes;
            figData.expEnzymes = expEnzymes;
            figData.enzProcessIdxs = enzProcessIdxs;
        end
        
        function figData = cellShape(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            %% get data
            stateNames = {
                'Time'      'values'             ':' ':'
                'Mass'      'cell'               ':' '-sum'
                'Geometry'  'volume'             ':' ':'
                'Geometry'  'width'              ':' ':'
                'Geometry'  'totalLength'        ':' ':'
                'Geometry'  'surfaceArea'        ':' ':'
                'Geometry'  'pinchedDiameter'    ':' ':'
                };
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations);
            simEndTimes = permute(max(states.Time.values, [], 3), [4 3 1 2]);
            
            states.Geometry.width(states.Geometry.width == -1) = NaN;
            states.Geometry.totalLength(states.Geometry.totalLength == -1) = NaN;
            states.Geometry.surfaceArea(states.Geometry.surfaceArea == -1) = NaN;
            
            time = permute(max(states.Time.values, [], 4), [3 4 1 2]) / 3600;
            mass = permute(states.Mass.cell, [3 4 1 2]) * 1e15;
            volume = permute(states.Geometry.volume, [3 4 1 2]) * 1e15;
            width = permute(states.Geometry.width, [3 4 1 2]) * 1e6;
            totalLength = permute(states.Geometry.totalLength, [3 4 1 2]) * 1e6;
            surfaceArea = permute(states.Geometry.surfaceArea, [3 4 1 2]) * 1e12;
            pinchedDiameter = permute(states.Geometry.pinchedDiameter, [3 4 1 2]) * 1e6;
            
            minSimEndTime = min(max(states.Time.values, [], 4));
            [~, simOrder] = sort(mass(minSimEndTime, :));
            
            for i = 1:numel(simEndTimes)
                mass(simEndTimes(i)+1:end, i) = NaN;
                volume(simEndTimes(i)+1:end, i) = NaN;
                width(simEndTimes(i)+1:end, i) = NaN;
                totalLength(simEndTimes(i)+1:end, i) = NaN;
                surfaceArea(simEndTimes(i)+1:end, i) = NaN;
                pinchedDiameter(simEndTimes(i)+1:end, i) = NaN;
            end
            
            clear states;
            
            %% plot data
            clf(figHandle);
            
            options = struct();
            options.titleStr = 'Cell Shape';
            options.colorOrder = PlotUtil.getRedGreenColorOrder(simOrder);
            options.xdata = time;
            options.ylabelStr = {
                'Mass (fg)'
                'Volume (fL)'
                {'Surface Area' '({\mu}m^2)'}
                'Length ({\mu}m)'
                'Width ({\mu}m)'
                {'Pinched' 'Diameter ({\mu}m)'}
                };
            options.ydata = {
                mass
                volume
                surfaceArea
                totalLength
                width
                pinchedDiameter
                };
            PlotUtil.multiElementPlot(figHandle, 3 * ones(6, 1), [0 time(end)], options);
            
            clear options;
            
            %% organize figure data
            figData = struct;
            figData.simOrder = simOrder;
            figData.time = time;
            figData.mass = mass;
            figData.volume = volume;
            figData.surfaceArea = surfaceArea;
            figData.totalLength = totalLength;
            figData.width = width;
            figData.pinchedDiameter = pinchedDiameter;
        end
        
        function figData = rnaPolymerases(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            %% constants
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            g = sim.gene;
            time = sim.state('Time');
            rna = sim.state('Rna');
            pc = sim.state('ProteinComplex');
            rnaPol = sim.state('RNAPolymerase');
            tr = sim.process('TranscriptionalRegulation');
            sc = sim.process('DNASupercoiling');
            pcComp = sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3);
            
            fcRNAIdxs = unique([tr.tuIndexs(:); sc.tuIndexs(:)]);
            
            expectedRNASyn = rna.expression(rna.matureIndexs);
            expectedRNASyn = expectedRNASyn / sum(expectedRNASyn);
            
            expectedRNAExp = expectedRNASyn ./ (log(2) / time.cellCycleLength + rna.decayRates(rna.matureIndexs));
            expectedRNAExp = expectedRNAExp / sum(expectedRNAExp);
            
            %% get data
            %RNA polymerases
            stateNames = {
                'Time'           'values'                                    ':' ':'
                'Mass'           'cell'                                      ':' '-sum'
                'RNAPolymerase'  'nActive'                                   ':' ':'
                'RNAPolymerase'  'nSpecificallyBound'                        ':' ':'
                'RNAPolymerase'  'nNonSpecificallyBound'                     ':' ':'
                'RNAPolymerase'  'nFree'                                     ':' ':'
                'RNAPolymerase'  'transcriptionFactorBindingProbFoldChange'  fcRNAIdxs 1
                'RNAPolymerase'  'supercoilingBindingProbFoldChange'         fcRNAIdxs 1
                };
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations);
            time = permute(max(states.Time.values, [], 4), [3 4 1 2]) / 3600;
            mass = permute(states.Mass.cell, [3 4 1 2]);
            
            simEndTimes = permute(max(states.Time.values, [], 3), [4 3 1 2]);
            minSimEndTime = min(simEndTimes);
            [~, simOrder] = sort(mass(minSimEndTime, :));
            
            foldChanges = ...
                states.RNAPolymerase.transcriptionFactorBindingProbFoldChange .* ...
                states.RNAPolymerase.supercoilingBindingProbFoldChange;
            for i = 1:numel(simEndTimes)
                foldChanges(:, :, simEndTimes(i)+1:end, i) = NaN;
            end
            
            actPols = states.RNAPolymerase.nActive;
            sbPols = states.RNAPolymerase.nSpecificallyBound;
            nsbPols = states.RNAPolymerase.nNonSpecificallyBound;
            freePols = states.RNAPolymerase.nFree;
            tmp = (actPols + sbPols + nsbPols + freePols);
            actPols = actPols ./ tmp * 100;
            sbPols = sbPols ./ tmp * 100;
            nsbPols = nsbPols ./ tmp * 100;
            freePols = freePols ./ tmp * 100;
            
            for i = 1:numel(simEndTimes)
                actPols(:, :, simEndTimes(i)+1:end, i) = NaN;
                sbPols(:, :, simEndTimes(i)+1:end, i) = NaN;
                nsbPols(:, :, simEndTimes(i)+1:end, i) = NaN;
                freePols(:, :, simEndTimes(i)+1:end, i) = NaN;
            end
            
            %cleanup
            clear clear states mass minSimEndTime tmp;
            
            %transcription units
            stateNames = {
                'RNAPolymerase'  'states'                   ':'  ':'
                'Transcript'     'boundTranscriptionUnits'  ':'  ':'
                'Transcript'     'boundTranscriptProgress'  ':'  ':'
                'Rna'            'counts'                   ':'  '-sum'
                'ProteinComplex' 'counts'                   ':'  '-sum'
                };
            rnaSyn = zeros(numel(rna.matureIndexs), 1, 1, numel(selectedSimulations));
            rnaExp = zeros(numel(rna.matureIndexs), 1, 1, numel(selectedSimulations));
            
            rnaPolAbsPosDensity = zeros(numel(rna.nascentIndexs), 1 + max(rna.lengths(rna.nascentIndexs)));
            
            for i = 1:numel(selectedSimulations)
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations(i));
                
                tfs = cat(3, diff(states.Transcript.boundTranscriptProgress, [], 3) < 0, states.Transcript.boundTranscriptionUnits(:, :, end) ~= rnaPol.notExistValue);
                rnaSyn(:, :, :, i) = histc(states.Transcript.boundTranscriptionUnits(tfs), (1:numel(rna.matureIndexs)));
                
                tmp = sum(states.Rna.counts, 3);
                rnaExp(:, :, :, i) = full(...
                    + tmp(rna.processedIndexs, :) ...
                    + tmp(rna.matureIndexs, :) ...
                    + tmp(rna.aminoacylatedIndexs, :) ...
                    + tmp(rna.boundIndexs, :) ...
                    + tmp(rna.misfoldedIndexs, :) ...
                    + tmp(rna.damagedIndexs, :));
                rnaExp(setdiff(1:end, rna.matureMRNAIndexs), :, :, i) = ...
                    + rnaExp(setdiff(1:end, rna.matureMRNAIndexs), :, :, i) ...
                    + pcComp * full(sum(reshape(sum(states.ProteinComplex.counts, 3), ...
                    [numel(pc.matureIndexs)  numel(pc.wholeCellModelIDs)/numel(pc.matureIndexs)]), 2));
                
                states.RNAPolymerase.states = full(states.RNAPolymerase.states);
                states.Transcript.boundTranscriptionUnits = full(states.Transcript.boundTranscriptionUnits);
                states.Transcript.boundTranscriptProgress = full(states.Transcript.boundTranscriptProgress);
                states.RNAPolymerase.states = states.RNAPolymerase.states(:);
                states.Transcript.boundTranscriptionUnits = states.Transcript.boundTranscriptionUnits(:);
                states.Transcript.boundTranscriptProgress = states.Transcript.boundTranscriptProgress(:);
                tfs = states.RNAPolymerase.states >= rnaPol.activelyTranscribingValue;
                rnaPolAbsPosDensity = ...
                    + rnaPolAbsPosDensity ...
                    + accumarray([states.Transcript.boundTranscriptionUnits(tfs)  states.Transcript.boundTranscriptProgress(tfs) + 1], ones(sum(tfs), 1), size(rnaPolAbsPosDensity));
                
                clear states tfs tmp;
            end
            
            totRNASyn = sum(rnaSyn, 1);
            rnaSyn = rnaSyn ./ totRNASyn(ones(size(rnaSyn, 1), 1), :, :, :);
            meanRNASyn = mean(rnaSyn, 4);
            stdRNASyn = std(rnaSyn, [], 4);
            
            totRNAExp = sum(rnaExp, 1);
            rnaExp = rnaExp ./ totRNAExp(ones(size(rnaExp, 1), 1), :, :, :);
            meanRNAExp = mean(rnaExp, 4);
            stdRNAExp = std(rnaExp, [], 4);
            
            rnaPolRelPosDensity = zeros(numel(rna.nascentIndexs), 100);
            for i = 1:numel(rna.nascentIndexs)
                rnaPolRelPosDensity(i, :) = accumarray(...
                    [ones(rna.lengths(rna.nascentIndexs(i)) + 1, 1)  ceil((1:rna.lengths(rna.nascentIndexs(i)) + 1)' / (rna.lengths(rna.nascentIndexs(i)) + 1) * 100)], ...
                    rnaPolAbsPosDensity(i, 1:rna.lengths(rna.nascentIndexs(i)) + 1)', [1 100]);
            end
            
            %cleanup
            clear rnaSyn rnaExp totRNASyn totRNAExp;
            
            %% plot data
            clf(figHandle);
            
            options = struct();
            options.titleStr = {'Transcriptional Regulation'; 'RNA Polymerase'; 'RNA'};
            options.colorOrder = PlotUtil.getRedGreenColorOrder(simOrder);
            axesHandles = PlotUtil.multiElementPlot(figHandle, {10 * ones(size(fcRNAIdxs)) [10 10 10 10] [10 10 10 10]}, [0 time(end)], options);
            
            %left
            for i = 1:numel(fcRNAIdxs)
                axesHandle = axesHandles{1}(i);
                PlotUtil.plotLine(axesHandle, time, permute(foldChanges(i, :, :, :), [3 4 1 2]));
                xlim(axesHandle, [0 time(end)]);
                ylabel(axesHandle, rna.wholeCellModelIDs{rna.nascentIndexs(fcRNAIdxs(i))});
            end
            
            %middle
            axesHandle = axesHandles{2}(1);
            plot(axesHandle, time, permute(actPols, [3 4 1 2]));
            xlim(axesHandle, [0 time(end)])
            ylim(axesHandle, [0 100])
            set(axesHandle, 'YTick', [0 50 100]);
            ylabel(axesHandle, 'Active');
            
            axesHandle = axesHandles{2}(2);
            plot(axesHandle, time, permute(sbPols, [3 4 1 2]));
            xlim(axesHandle, [0 time(end)])
            ylim(axesHandle, [0 100])
            set(axesHandle, 'YTick', [0 50 100]);
            ylabel(axesHandle, {'Specifically' 'bound'});
            
            axesHandle = axesHandles{2}(3);
            plot(axesHandle, time, permute(nsbPols, [3 4 1 2]));
            xlim(axesHandle, [0 time(end)])
            ylim(axesHandle, [0 100])
            set(axesHandle, 'YTick', [0 50 100]);
            ylabel(axesHandle, {'Non-specifically' 'bound'});
            
            axesHandle = axesHandles{2}(4);
            plot(axesHandle, time, permute(freePols, [3 4 1 2]));
            xlim(axesHandle, [0 time(end)])
            ylim(axesHandle, [0 100])
            set(axesHandle, 'YTick', [0 50 100]);
            ylabel(axesHandle, 'Free');
            
            %right
            axesHandle = axesHandles{3}(1);
            cla(axesHandle);
            density = smooth(sum(rnaPolAbsPosDensity, 1) / sum(simEndTimes) * 3600, 11, 'moving');
            plot(axesHandle, 0:max(rna.lengths(rna.nascentIndexs)), density, 'Color', 'r');
            xlim(axesHandle, [-0.5 max(rna.lengths(rna.nascentIndexs)) + 0.5]);
            ylim(axesHandle, [0 max(density)]);
            set(axesHandle, 'XTick', [0 max(rna.lengths(rna.nascentIndexs))]);
            set(axesHandle, 'YTickMode', 'auto');
            set(axesHandle, 'XColor', 'k');
            set(axesHandle, 'YScale', 'log')
            xlabel(axesHandle, 'Absolute Position');
            ylabel(axesHandle, {'Density'  '(h^{-1})'});
            
            axesHandle = axesHandles{3}(2);
            density = smooth(sum(rnaPolRelPosDensity, 1) / sum(simEndTimes) * 3600, 11, 'moving');
            plot(axesHandle, (1:100) - 0.5, density, 'Color', 'r');
            xlim(axesHandle, [0 100]);
            ylim(axesHandle, [0 max(density)]);
            set(axesHandle, 'XTick', [0 100]);
            set(axesHandle, 'YTickMode', 'auto');
            set(axesHandle, 'XColor', 'k');
            set(axesHandle, 'YScale', 'log')
            xlabel(axesHandle, 'Relative Position (%)');
            ylabel(axesHandle, {'Density'  '(h^{-1})'});
            
            axesHandle = axesHandles{3}(3);
            [~, order] = sort(expectedRNASyn);
            errorbar((1:numel(rna.matureIndexs))', meanRNASyn(order), stdRNASyn(order), 'Parent', axesHandle, 'Color', 'g');
            plot(axesHandle, (1:numel(rna.matureIndexs))', expectedRNASyn(order), 'Color', 'b');
            xlabel(axesHandle, 'Transcription Units');
            xlim(axesHandle, [0.5 numel(rna.matureIndexs)+0.5])
            ylabel(axesHandle, 'Synthesis');
            ylim(axesHandle, [0 0.05]);
            set(axesHandle, 'YTick', 0:0.05:0.05);
            
            axesHandle = axesHandles{3}(4);
            [~, order] = sort(expectedRNAExp);
            errorbar((1:numel(rna.matureIndexs))', meanRNAExp(order), stdRNAExp(order), 'Parent', axesHandle, 'Color', 'g');
            plot(axesHandle, (1:numel(rna.matureIndexs))', expectedRNAExp(order), 'Color', 'b');
            xlabel(axesHandle, 'Transcription Units');
            xlim(axesHandle, [0.5 numel(rna.matureIndexs)+0.5])
            ylabel(axesHandle, 'Expression');
            ylim(axesHandle, [0 0.05]);
            set(axesHandle, 'YTick', 0:0.05:0.05);
            
            %y axes
            PlotUtil.alignYAxesLabels(axesHandles{1});
            PlotUtil.alignYAxesLabels(axesHandles{2});
            PlotUtil.alignYAxesLabels(axesHandles{3});
            PlotUtil.offsetYAxes(axesHandles{1}, 0.03);
            PlotUtil.offsetYAxes(axesHandles{2}, 0.03);
            PlotUtil.offsetYAxes(axesHandles{3}, 0.03);
            
            %% organize figure data
            figData = struct;
            figData.simEndTimes = simEndTimes;
            figData.simOrder = simOrder;
            figData.time = time;
            figData.foldChanges = foldChanges;
            figData.actPols = actPols;
            figData.sbPols = sbPols;
            figData.nsbPols = nsbPols;
            figData.freePols = freePols;
            figData.meanRNASyn = meanRNASyn;
            figData.stdRNASyn = stdRNASyn;
            figData.meanRNAExp = meanRNAExp;
            figData.stdRNAExp = stdRNAExp;
            figData.expectedRNASyn = expectedRNASyn;
            figData.expectedRNAExp = expectedRNAExp;
            figData.rnaPolAbsPosDensity = rnaPolAbsPosDensity;
            figData.rnaPolRelPosDensity = rnaPolRelPosDensity;
        end
        
        function figData = ribosomes(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            %% constants
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            g = sim.gene;
            time = sim.state('Time');
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            rib = sim.state('Ribosome');
            pcComp = sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3);
            
            expectedMonSyn = rna.matureRNAGeneComposition(g.mRNAIndexs, :) * rna.expression(rna.matureIndexs);
            expectedMonSyn = expectedMonSyn / sum(expectedMonSyn);
            
            expectedMonExp = expectedMonSyn ./ (log(2) / time.cellCycleLength + pm.decayRates(pm.matureIndexs));
            expectedMonExp = expectedMonExp / sum(expectedMonExp);
            
            %% get data
            %ribosome states
            stateNames = {
                'Time'           'values'          ':' ':'
                'Mass'           'cell'            ':' '-sum'
                'Ribosome'       'nActive'         ':' ':'
                'Ribosome'       'nStalled'        ':' ':'
                };
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations);
            time = permute(max(states.Time.values, [], 4), [3 4 1 2]) / 3600;
            mass = permute(states.Mass.cell, [3 4 1 2]);
            
            simEndTimes = permute(max(states.Time.values, [], 3), [4 3 1 2]);
            minSimEndTime = min(simEndTimes);
            [~, simOrder] = sort(mass(minSimEndTime, :));
            
            actRibs = states.Ribosome.nActive;
            stalledRibs = states.Ribosome.nStalled;
            tmp = (actRibs + stalledRibs);
            actRibs = actRibs ./ tmp * 100;
            stalledRibs = stalledRibs ./ tmp * 100;
            
            for i = 1:numel(simEndTimes)
                actRibs(:, :, simEndTimes(i)+1:end, i) = NaN;
                stalledRibs(:, :, simEndTimes(i)+1:end, i) = NaN;
            end
            
            %cleanup
            clear states mass minSimEndTime tmp;
            
            %synthesis, expression
            stateNames = {
                'Ribosome'       'states'          ':'  ':'
                'Ribosome'       'boundMRNAs'      ':'  ':'
                'Ribosome'       'mRNAPositions'   ':'  ':'
                'ProteinMonomer' 'counts'          ':'  '-sum'
                'ProteinComplex' 'counts'          ':'  '-sum'
                };
            monSyn = zeros(numel(pm.matureIndexs), 1, 1, numel(selectedSimulations));
            monExp = zeros(numel(pm.matureIndexs), 1, 1, numel(selectedSimulations));
            cpxExp = zeros(numel(pc.matureIndexs), 1, 1, numel(selectedSimulations));
            ribAbsPosDensity = zeros(numel(pm.nascentIndexs), 1 + max(pm.lengths(pm.nascentIndexs)));
            for i = 1:numel(selectedSimulations)
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations(i));
                
                boundMRNAs = states.Ribosome.boundMRNAs;
                mRNAPositions = states.Ribosome.mRNAPositions;
                tfs = cat(3, diff(mRNAPositions, [], 3) < 0, boundMRNAs(:, :, end) ~= rib.notExistValue);
                monSyn(:, :, :, i) = histc(boundMRNAs(tfs), (1:numel(pm.matureIndexs)));
                
                monExp(:, :, :, i) = full(sum(reshape(sum(states.ProteinMonomer.counts, 3), ...
                    [numel(pm.matureIndexs)  numel(pm.wholeCellModelIDs)/numel(pm.matureIndexs)]), 2));
                cpxExp(:, :, :, i) = full(sum(reshape(sum(states.ProteinComplex.counts, 3), ...
                    [numel(pc.matureIndexs)  numel(pc.wholeCellModelIDs)/numel(pc.matureIndexs)]), 2));
                
                states.Ribosome.states = full(states.Ribosome.states);
                states.Ribosome.boundMRNAs = full(states.Ribosome.boundMRNAs);
                states.Ribosome.mRNAPositions = full(states.Ribosome.mRNAPositions);
                states.Ribosome.states = states.Ribosome.states(:);
                states.Ribosome.boundMRNAs = states.Ribosome.boundMRNAs(:);
                states.Ribosome.mRNAPositions = states.Ribosome.mRNAPositions(:);
                tfs = states.Ribosome.states == rib.activeValue;
                ribAbsPosDensity = ...
                    + ribAbsPosDensity ...
                    + accumarray([states.Ribosome.boundMRNAs(tfs)  states.Ribosome.mRNAPositions(tfs) + 1], ones(sum(tfs), 1), size(ribAbsPosDensity));
                
                clear states boundMRNAs mRNAPositions tfs;
            end
            
            totMonSyn = sum(monSyn, 1);
            monSyn = monSyn ./ totMonSyn(ones(size(monSyn, 1), 1), :, :, :);
            meanMonSyn = mean(monSyn, 4);
            stdMonSyn = std(monSyn, [], 4);
            
            for j = 1:size(monExp, 4)
                monExp(:, :, :, j) = ...
                    + monExp(:, :, :, j) ...
                    + pcComp * cpxExp(:, :, :, j);
            end
            totMonExp = sum(monExp, 1);
            monExp = monExp ./ totMonExp(ones(size(monExp, 1), 1), :, :, :);
            meanMonExp = mean(monExp, 4);
            stdMonExp = std(monExp, [], 4);
            
            ribRelPosDensity = zeros(numel(pm.nascentIndexs), 100);
            for i = 1:numel(pm.nascentIndexs)
                ribRelPosDensity(i, :) = accumarray(...
                    [ones(pm.lengths(pm.nascentIndexs(i)) + 1, 1)  ceil((1:pm.lengths(pm.nascentIndexs(i)) + 1)' / (pm.lengths(pm.nascentIndexs(i)) + 1) * 100)], ...
                    ribAbsPosDensity(i, 1:pm.lengths(pm.nascentIndexs(i)) + 1)', [1 100]);
            end
            
            %cleanup
            clear monExp monSyn cpxExp totMonSyn totMonExp;
            
            %% plot data
            clf(figHandle);
            
            options = struct();
            options.titleStr = {'States'; 'Monomers'};
            options.colorOrder = PlotUtil.getRedGreenColorOrder(simOrder);
            axesHandles = PlotUtil.multiElementPlot(figHandle, {[10 10] [10 10 10 10]}, [0 time(end)], options);
            
            %left
            axesHandle = axesHandles{1}(1);
            plot(axesHandle, time, permute(actRibs, [3 4 1 2]));
            xlim(axesHandle, [0 time(end)])
            ylim(axesHandle, [0 100])
            set(axesHandle, 'YTick', [0 50 100]);
            ylabel(axesHandle, 'Active');
            
            axesHandle = axesHandles{1}(2);
            plot(axesHandle, time, permute(stalledRibs, [3 4 1 2]));
            xlim(axesHandle, [0 time(end)])
            ylim(axesHandle, [0 100])
            set(axesHandle, 'YTick', [0 50 100]);
            ylabel(axesHandle, 'Stalled');
            
            %right
            axesHandle = axesHandles{2}(1);
            density = smooth(sum(ribAbsPosDensity, 1) / sum(simEndTimes) * 3600, 11, 'moving');
            plot(axesHandle, 0:max(pm.lengths(pm.nascentIndexs)), density, 'Color', 'r');
            xlim(axesHandle, [-0.5 max(pm.lengths(pm.nascentIndexs)) + 0.5]);
            ylim(axesHandle, [0 max(density)]);
            set(axesHandle, 'XTick', [0 max(pm.lengths(pm.nascentIndexs))]);
            set(axesHandle, 'YTickMode', 'auto', 'YScale', 'log');
            set(axesHandle, 'XColor', 'k');
            xlabel(axesHandle, 'Absolute Position');
            ylabel(axesHandle, {'Density'  '(h^{-1})'});
            
            axesHandle = axesHandles{2}(2);
            density = smooth(sum(ribRelPosDensity, 1) / sum(simEndTimes) * 3600, 11, 'moving');
            plot(axesHandle, (1:100) - 0.5, density, 'Color', 'r');
            xlim(axesHandle, [0 100]);
            ylim(axesHandle, [0 max(density)]);
            set(axesHandle, 'XTick', [0 100]);
            set(axesHandle, 'YTickMode', 'auto', 'YScale', 'log');
            set(axesHandle, 'XColor', 'k');
            xlabel(axesHandle, 'Relative Position (%)');
            ylabel(axesHandle, {'Density'  '(h^{-1})'});
            
            axesHandle = axesHandles{2}(3);
            [~, order] = sort(expectedMonSyn);
            errorbar((1:numel(pm.matureIndexs))', meanMonSyn(order), stdMonSyn(order), 'Parent', axesHandle, 'Color', 'g');
            plot(axesHandle, (1:numel(pm.matureIndexs))', expectedMonSyn(order), 'Color', 'b');
            xlabel(axesHandle, 'Protein-Coding Genes');
            xlim(axesHandle, [0.5 numel(pm.matureIndexs)+0.5])
            ylabel(axesHandle, 'Synthesis');
            ylim(axesHandle, [min(meanMonSyn(order) - stdMonSyn(order)) max(meanMonSyn(order) + stdMonSyn(order))]);
            set(axesHandle, 'YTick', 0:0.01:0.02);
            
            axesHandle = axesHandles{2}(4);
            [~, order] = sort(expectedMonExp);
            errorbar((1:numel(pm.matureIndexs))', meanMonExp(order), stdMonExp(order), 'Parent', axesHandle, 'Color', 'g');
            plot(axesHandle, (1:numel(pm.matureIndexs))', expectedMonExp(order), 'Color', 'b');
            xlabel(axesHandle, 'Protein-Coding Genes');
            xlim(axesHandle, [0.5 numel(pm.matureIndexs)+0.5])
            ylabel(axesHandle, 'Expression');
            ylim(axesHandle, [min(meanMonExp(order) - stdMonExp(order)) max(meanMonExp(order) + stdMonExp(order))]);
            set(axesHandle, 'YTick', 0:0.01:0.02);
            
            %y axes
            PlotUtil.alignYAxesLabels(axesHandles{1});
            PlotUtil.alignYAxesLabels(axesHandles{2});
            PlotUtil.offsetYAxes(axesHandles{1}, 0.03);
            PlotUtil.offsetYAxes(axesHandles{2}, 0.03);
            
            %% organize figure data
            figData = struct;
            figData.simEndTimes = simEndTimes;
            figData.simOrder = simOrder;
            figData.time = time;
            figData.actRibs = actRibs;
            figData.stalledRibs = stalledRibs;
            figData.expectedMonSyn = expectedMonSyn;
            figData.expectedMonExp = expectedMonExp;
            figData.meanMonSyn = meanMonSyn;
            figData.stdMonSyn = stdMonSyn;
            figData.meanMonExp = meanMonExp;
            figData.stdMonExp = stdMonExp;
            figData.ribAbsPosDensity = ribAbsPosDensity;
            figData.ribRelPosDensity = ribRelPosDensity;
        end
        
        function figData = macromoleculeExpression(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            %% get data
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            %% constants
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            g = sim.gene;
            m = sim.state('Metabolite');
            rna = sim.state('Rna');
            mon = sim.state('ProteinMonomer');
            cpx = sim.state('ProteinComplex');
            
            rnaPcComp = sum(cpx.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3);
            monPcComp = sum(cpx.proteinComplexComposition(g.mRNAIndexs, :, :), 3);
            
            essGene = strcmp(g.essential, 'Y');
            
            %% get data
            stateNames = {
                'Rna'             'counts'
                'ProteinMonomer'  'counts'
                'ProteinComplex'  'counts'
                };
            rnaFreq = zeros(1000, numel(rna.matureIndexs));
            monFreq = zeros(10000, numel(mon.matureIndexs));
            for i = 1:numel(selectedSimulations)
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations(i));
                
                tmp = permute(sum(states.Rna.counts, 2), [1 3 2]);
                rnaCnts = ...
                    + rna.nascentRNAMatureRNAComposition * full(tmp(rna.nascentIndexs, :)) ...
                    + full(...
                    + tmp(rna.processedIndexs, :) ...
                    + tmp(rna.matureIndexs, :) ...
                    + tmp(rna.boundIndexs, :) ...
                    + tmp(rna.misfoldedIndexs, :) ...
                    + tmp(rna.damagedIndexs, :) ...
                    + tmp(rna.aminoacylatedIndexs, :) ...
                    );
                clear tmp;
                
                monCnts = full(permute(sum(reshape(sum(states.ProteinMonomer.counts, 2), [numel(mon.matureIndexs)  size(states.ProteinMonomer.counts, 1) / numel(mon.matureIndexs) size(states.ProteinMonomer.counts, 3)]), 2), [1 3 2]));
                cpxCnts = full(permute(sum(reshape(sum(states.ProteinComplex.counts, 2), [numel(cpx.matureIndexs)  size(states.ProteinComplex.counts, 1) / numel(cpx.matureIndexs) size(states.ProteinComplex.counts, 3)]), 2), [1 3 2]));
                rnaCnts(setdiff(1:end, rna.matureMRNAIndexs), :) = ...
                    + rnaCnts(setdiff(1:end, rna.matureMRNAIndexs), :) ...
                    + rnaPcComp * cpxCnts;
                monCnts = ...
                    + monCnts ...
                    + monPcComp * cpxCnts;
                
                idxs = (1:numel(rna.matureIndexs))';
                idxs = idxs(:, ones(size(rnaCnts, 2), 1));
                rnaFreq = ...
                    + rnaFreq ...
                    + accumarray([rnaCnts(:)+1 idxs(:)], ones(numel(rnaCnts), 1), size(rnaFreq));
                
                idxs = (1:numel(mon.matureIndexs))';
                idxs = idxs(:, ones(size(monCnts, 2), 1));
                monFreq = ...
                    + monFreq ...
                    + accumarray([monCnts(:)+1 idxs(:)], ones(numel(monCnts), 1), size(monFreq));
                
                clear states rnaCnts monCnts cpxCnts idxs;
            end
            
            rnaFreq = rnaFreq(1:find(any(rnaFreq, 2), 1, 'last'), :);
            monFreq = monFreq(1:find(any(monFreq, 2), 1, 'last'), :);
            
            meanRnaCnt = ((0:size(rnaFreq, 1)-1) * rnaFreq) ./ sum(rnaFreq, 1);
            stdRnaCnt = sqrt((((1:size(rnaFreq, 1)).^2) * rnaFreq) ./ sum(rnaFreq, 1) - meanRnaCnt.^2);
            meanMRnaCnt = meanRnaCnt * rna.matureRNAGeneComposition(g.mRNAIndexs, :)';
            stdMRnaCnt = stdRnaCnt * rna.matureRNAGeneComposition(g.mRNAIndexs, :)';
            meanMonCnt = ((0:size(monFreq, 1)-1) * monFreq) ./ sum(monFreq, 1);
            stdMonCnt = sqrt((((1:size(monFreq, 1)).^2) * monFreq) ./ sum(monFreq, 1) - meanMonCnt.^2);
            
            rnaLognParam = NaN(numel(rna.matureIndexs), 2);
            rnaLognParamCI = NaN(numel(rna.matureIndexs), 2, 2);
            rnaNonNormPValue = zeros(numel(rna.matureIndexs), 1);
            for i = 1:numel(rna.matureIndexs)
                cnts = rude(rnaFreq(:, i), 0:size(rnaFreq, 1) - 1);
                if ~rnaFreq(1, i)
                    [rnaLognParam(i, :), rnaLognParamCI(i, :, :)] = lognfit(cnts);
                end
                if cnts > 5e3
                    cnts = cnts(randsample(numel(cnts), 5e3, false));
                end
                [~, rnaNonNormPValue(i)] = swtest(log10(cnts), 0.05, 0);
            end
            
            monGamParam = zeros(numel(mon.matureIndexs), 2);
            monGamParamCI = zeros(numel(mon.matureIndexs), 2, 2);
            monNonNormPValue = zeros(numel(mon.matureIndexs), 1);
            for i = 1:numel(mon.matureIndexs)
                cnts = rude(monFreq(:, i), 0:size(monFreq, 1) - 1);
                [monGamParam(i, :), monGamParamCI(i, :, :)] = gamfit(cnts);
                if cnts > 5e3
                    cnts = cnts(randsample(numel(cnts), 5e3, false));
                end
                [~, monNonNormPValue(i)] = swtest(cnts, 0.05, 0);
            end
            
            %% plot data
            clf(figHandle);
            
            [axesHandles, xAxesHandles, ~, offsetAxesHandles] = PlotUtil.multiElementPlot(figHandle, {1 * ones(3, 1); 1 * ones(2,1); 3 * ones(3, 1); 3 * ones(4, 1)}, [0 1]);
            delete(cell2mat(xAxesHandles));
            
            set(cell2mat(offsetAxesHandles), 'visible', 'on');
            set(cell2mat(axesHandles), 'xcolor', 'k', 'xtickmode', 'auto', 'xlimmode', 'auto', 'ylimmode', 'auto');
            
            %mRNA count distribution for particular gene
            [~, idx] = max(meanMRnaCnt);
            tmp = rnaFreq(:, rna.matureMRNAIndexs(idx)) / sum(rnaFreq(:, rna.matureMRNAIndexs(idx)));
            axesHandle = axesHandles{1}(1);
            cla(axesHandle);
            bar(axesHandle, (1:size(tmp, 1))'-1, tmp, 'EdgeColor', 'none', 'BarWidth', 0.9);
            minX = find(tmp, 1, 'first') - 1.5;
            maxX = find(tmp, 1, 'last') - 0.5;
            xlim(axesHandle, [minX maxX]);
            ylim(axesHandle, [0 max(tmp)]);
            xlabel(axesHandle, 'mRNA Count');
            ylabel(axesHandle, 'Freq');
            
            %rRNA count distribution for particular gene
            x = rnaFreq(:, rna.matureRRNAIndexs(1));
            cnts = [];
            for j = 0:numel(x) - 1
                cnts = [cnts; j(ones(x(j+1), 1), 1)]; %#ok<AGROW>
            end
            axesHandle = axesHandles{1}(2);
            cla(axesHandle);
            hold(axesHandle, 'on');
            
            edges = logspace(log10(min(cnts)), log10(max(cnts)), 200);
            freq = histc(cnts, edges);
            freq = freq / sum(freq);
            bar(axesHandle, log10(edges), freq, 'EdgeColor', 'none', 'BarWidth', 2);
            set(axesHandle, 'XScale', 'linear');
            xlim(axesHandle, log10(edges([1 end])))
            set(axesHandle, 'XTickMode', 'auto');
            set(axesHandle, 'XTick', 10.^get(axesHandle, 'XTick'));
            ylim(axesHandle, [0 max(freq)])
            
            y = lognpdf(edges, rnaLognParam(rna.matureRRNAIndexs(2), 1), rnaLognParam(rna.matureRRNAIndexs(2), 2));
            plot(axesHandle, log10(edges), y, 'r');
            
            xlabel(axesHandle, 'rRNA Count');
            ylabel(axesHandle, 'Freq');
            
            %monomer count distribution for particular gene
            tmp = monFreq(:, 1) ./ sum(monFreq(:, 1));
            axesHandle = axesHandles{1}(3);
            cla(axesHandle);
            bar(axesHandle, (1:size(tmp, 1))'-1, tmp, 'EdgeColor', 'none', 'BarWidth', 2);
            minX = find(tmp, 1, 'first') - 1.5;
            maxX = find(tmp, 1, 'last') - 0.5;
            xlim(axesHandle, [minX maxX])
            ylim([0 max(tmp)]);
            
            x = linspace(max(0, minX), maxX, 100);
            y = gampdf(x, monGamParam(1, 1), monGamParam(1, 2));
            plot(axesHandle, x, y, 'r');
            
            xlabel(axesHandle, 'Protein Count');
            ylabel(axesHandle, 'Freq');
            
            %freq of strains with mean protein number
            axesHandle = axesHandles{2}(1);
            cla(axesHandle);
            edges = logspace(log10(min(meanMonCnt)), log10(max(meanMonCnt)), 10);
            cnts = [
                histc(meanMonCnt, edges)
                histc(meanMonCnt(essGene(g.mRNAIndexs)), edges)
                ];
            h = bar(axesHandle, log10(edges'), cnts', 'grouped');
            set(h, 'BarWidth', 2);
            set(h(1), 'EdgeColor', 'none')
            set(h(2), 'EdgeColor', 'none')
            set(axesHandle, 'XScale', 'linear');
            xlim(axesHandle, log10(edges([1 end])))
            set(axesHandle, 'XTickMode', 'auto');
            set(axesHandle, 'XTick', 10.^get(axesHandle, 'XTick'));
            xlabel(axesHandle, 'Prot Exp');
            ylabel(axesHandle, 'Freq');
            
            %noise vs mean protein number
            axesHandle = axesHandles{2}(2);
            cla(axesHandle);
            hold(axesHandle, 'on');
            noise = (stdMonCnt ./ meanMonCnt).^2;
            x = quantile(meanMonCnt, [0.01 0.99])';
            y = 5./quantile(meanMonCnt, [0.01 0.99])';
            plot(axesHandle, meanMonCnt, noise, 'k.');
            plot(axesHandle, x, y, 'r-');
            plot(axesHandle, [min(meanMonCnt) max(meanMonCnt)], quantile(noise, 0.01) * [1 1], 'b-');
            ylim(axesHandle, [min(noise) max(noise)])
            xlim(axesHandle, [min(meanMonCnt) max(meanMonCnt)])
            set(axesHandle, 'XScale', 'log');
            set(axesHandle, 'YScale', 'log');
            xlabel(axesHandle, 'Mean Protein Number (\mu)');
            ylabel(axesHandle, 'Protein Noise (\sigma^2 / \mu^2)');
            
            %Protein vs RNA expression
            axesHandle = axesHandles{3}(1);
            cla(axesHandle);
            plot(axesHandle, meanMRnaCnt, meanMonCnt, '.');
            xlim(axesHandle, [min(meanMRnaCnt) max(meanMRnaCnt)])
            ylim(axesHandle, [min(meanMonCnt) max(meanMonCnt)])
            xlabel(axesHandle, 'Mean mRNA Level');
            ylabel(axesHandle, 'Mean Protein Number');
            set(axesHandle, 'XScale', 'log');
            set(axesHandle, 'YScale', 'log');
            
            %RNA expression vs noise
            axesHandle = axesHandles{3}(2);
            cla(axesHandle);
            noise = (stdMRnaCnt ./ meanMRnaCnt) .^2;
            plot(axesHandle, meanMRnaCnt, noise, '.');
            xlim(axesHandle, [min(meanMRnaCnt) max(meanMRnaCnt)])
            ylim(axesHandle, [min(noise) max(noise)])
            xlabel(axesHandle, 'Mean MRNA Number');
            ylabel(axesHandle, 'mRNA Noise');
            set(axesHandle, 'XScale', 'log');
            set(axesHandle, 'YScale', 'log');
            
            %mRNA fano factor frequence
            axesHandle = axesHandles{3}(3);
            cla(axesHandle);
            tfs = ~isnan(stdMRnaCnt) & meanMRnaCnt > 0;
            hist(axesHandle, (stdMRnaCnt(tfs) .^2) ./ meanMRnaCnt(tfs), 200);
            h = findobj(axesHandle, 'Type', 'patch');
            set(h, 'EdgeColor', 'none');
            xlim(axesHandle, [2 20]);
            xlabel(axesHandle, 'mRNA Fano Factor');
            ylabel(axesHandle, 'Number of Strains');
            
            %mRNA gamma a vs b
            axesHandle = axesHandles{4}(1);
            cla(axesHandle)
            plot(axesHandle, monGamParam(:, 1), monGamParam(:, 2), '.');
            xlim(axesHandle, [min(monGamParam(:, 1)) max(monGamParam(:, 1))])
            ylim(axesHandle, [min(monGamParam(:, 2)) max(monGamParam(:, 2))])
            xlabel(axesHandle, 'a Value');
            ylabel(axesHandle, 'b Value');
            set(axesHandle, 'XScale', 'log');
            set(axesHandle, 'XScale', 'log');
            
            axesHandle = axesHandles{4}(2);
            cla(axesHandle)
            plot(axesHandle, meanMonCnt, monGamParam(:, 2), '.');
            xlim(axesHandle, [min(meanMonCnt) max(meanMonCnt)])
            ylim(axesHandle, [min(monGamParam(:, 2)) max(monGamParam(:, 2))])
            xlabel(axesHandle, 'Mean Protein Number');
            ylabel(axesHandle, 'b Value');
            set(axesHandle, 'XScale', 'log');
            set(axesHandle, 'XScale', 'log');
            
            %protein expression vs gc content
            axesHandle = axesHandles{4}(3);
            cla(axesHandle)
            tmp = rna.baseCounts(rna.matureIndexs, m.nmpIndexs);
            gcContent = rna.matureRNAGeneComposition(g.mRNAIndexs, :) * (sum(tmp(:, [2 3]), 2) ./ sum(tmp, 2));
            plot(axesHandle, meanMonCnt, gcContent, '.');
            xlim(axesHandle, [min(meanMonCnt) max(meanMonCnt)])
            ylim(axesHandle, [min(gcContent) max(gcContent)])
            xlabel('Mean Protein Number');
            ylabel('GC Content');
            set(axesHandle, 'XScale', 'log');
            
            %protein expression vs mRNA half life
            axesHandle = axesHandles{4}(4);
            cla(axesHandle)
            plot(axesHandle, meanMonCnt, rna.matureRNAGeneComposition(g.mRNAIndexs, :) * rna.halfLives(rna.matureIndexs), '.');
            xlim(axesHandle, [min(meanMonCnt) max(meanMonCnt)])
            ylim(axesHandle, [min(rna.halfLives(rna.matureIndexs(rna.matureMRNAIndexs))) max(rna.halfLives(rna.matureIndexs(rna.matureMRNAIndexs)))])
            xlabel('Mean Protein Number');
            ylabel('mRNA Life Time (min)');
            set(axesHandle, 'XScale', 'log');
            set(axesHandle, 'YScale', 'log');
            
            %% data
            figData = struct;
            figData.rnaFreq = rnaFreq;
            figData.monFreq = monFreq;
            figData.meanRnaCnt = meanRnaCnt;
            figData.stdRnaCnt = stdRnaCnt;
            figData.meanMRnaCnt = meanMRnaCnt;
            figData.stdMRnaCnt = stdMRnaCnt;
            figData.meanMonCnt = meanMonCnt;
            figData.stdMonCnt = stdMonCnt;
            figData.rnaLognParam = rnaLognParam;
            figData.rnaLognParamCI = rnaLognParamCI;
            figData.rnaNonNormPValue = rnaNonNormPValue;
            figData.monGamParam = monGamParam;
            figData.monGamParamCI = monGamParamCI;
            figData.monNonNormPValue = monNonNormPValue;
        end
        
        function figData = dnaBoundProteinDisplacement(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            %% get data
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            %% constants
            nSims = numel(selectedSimulations);
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            c = sim.state('Chromosome');
            pc = sim.state('ProteinComplex');
            
            %% get data
            stateNames = {
                'Chromosome' 'monomerBoundSites'
                'Chromosome' 'complexBoundSites'
                };
            displacements = zeros(numel(c.reactionBoundMonomer), nSims);
            simEndTimes = zeros(1, nSims);
            for i = 1:nSims
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations(i));
                
                simEndTimes(i) = size(states.Chromosome.monomerBoundSites, 3);
                [monPosStrndTimes, monIdxs] = find(states.Chromosome.monomerBoundSites);
                [cpxPosStrndTimes, cpxIdxs] = find(states.Chromosome.complexBoundSites);
                
                rxnIdxs = find(...
                    (ismember(c.reactionBoundComplex, [pc.dnaPolymeraseIndexs; pc.rnaPolymeraseIndexs]) | ...
                    any(c.reactionComplexCatalysisMatrix(:, [pc.dnaPolymeraseIndexs; pc.rnaPolymeraseIndexs]), 2)) & ...
                    (ismember(c.reactionBoundMonomer, monIdxs) | ismember(c.reactionBoundComplex, cpxIdxs)) & ...
                    (any(c.reactionMonomerCatalysisMatrix(:, unique(monIdxs)), 2) | any(c.reactionComplexCatalysisMatrix(:, unique(cpxIdxs)), 2)) ...
                    );
                for j = 1:numel(rxnIdxs)
                    if c.reactionBoundMonomer(rxnIdxs(j))
                        tfs = monIdxs == c.reactionBoundMonomer(rxnIdxs(j));
                        boundPosStrdTimes = monPosStrndTimes(tfs, :);
                        boundFtpt = c.monomerDNAFootprints(c.reactionBoundMonomer(rxnIdxs(j)));
                    else
                        tfs = cpxIdxs == c.reactionBoundComplex(rxnIdxs(j));
                        boundPosStrdTimes = cpxPosStrndTimes(tfs, :);
                        boundFtpt = c.complexDNAFootprints(c.reactionBoundComplex(rxnIdxs(j)));
                    end
                    if isempty(boundPosStrdTimes)
                        continue;
                    end
                    
                    if any(c.reactionMonomerCatalysisMatrix(rxnIdxs(j), :))
                        tfs = monIdxs == find(c.reactionMonomerCatalysisMatrix(rxnIdxs(j), :));
                        bindingPosStrdTimes = monPosStrndTimes(tfs, :);
                        bindingFtpt = c.monomerDNAFootprints(c.reactionMonomerCatalysisMatrix(rxnIdxs(j), :) ~= 0);
                    else
                        tfs = cpxIdxs == find(c.reactionComplexCatalysisMatrix(rxnIdxs(j), :));
                        bindingPosStrdTimes = cpxPosStrndTimes(tfs, :);
                        bindingFtpt = c.complexDNAFootprints(c.reactionComplexCatalysisMatrix(rxnIdxs(j), :) ~= 0);
                    end
                    if isempty(bindingPosStrdTimes)
                        continue;
                    end
                    
                    times = intersect(unique(boundPosStrdTimes(:, 3)), unique(bindingPosStrdTimes(:, 3) - 1));
                    for k = 1:numel(times)
                        tmpboundPosStrdTimes = boundPosStrdTimes(boundPosStrdTimes(:, 3) == times(k), 1:2);
                        tmpbindingPosStrdTimes = bindingPosStrdTimes(bindingPosStrdTimes(:, 3) == times(k) + 1, 1:2);
                        
                        if numel(tmpboundPosStrdTimes) <= numel(tmpbindingPosStrdTimes)
                            for l = 1:size(tmpboundPosStrdTimes, 1)
                                if ...
                                        any(tmpbindingPosStrdTimes(tmpbindingPosStrdTimes(:, 1) >= tmpboundPosStrdTimes(l, 1), 1) - tmpboundPosStrdTimes(l, 1) < boundFtpt) || ...
                                        any(tmpbindingPosStrdTimes(:, 1) + c.sequenceLen - tmpboundPosStrdTimes(l, 1) < boundFtpt) || ...
                                        any(tmpboundPosStrdTimes(l, 1) - tmpbindingPosStrdTimes(tmpbindingPosStrdTimes(:, 1) < tmpboundPosStrdTimes(l, 1), 1) < bindingFtpt) || ...
                                        any(tmpboundPosStrdTimes(l, 1) + c.sequenceLen - tmpbindingPosStrdTimes(:, 1) < boundFtpt)
                                    displacements(rxnIdxs(j), i) = ...
                                        displacements(rxnIdxs(j), i) + 1;
                                end
                            end
                        else
                            for l = 1:size(tmpbindingPosStrdTimes, 1)
                                if ...
                                        any(tmpboundPosStrdTimes(tmpboundPosStrdTimes(:, 1) >= tmpbindingPosStrdTimes(l, 1), 1) - tmpbindingPosStrdTimes(l, 1) < bindingFtpt) || ...
                                        any(tmpboundPosStrdTimes(:, 1) + c.sequenceLen - tmpbindingPosStrdTimes(l, 1) < bindingFtpt) || ...
                                        any(tmpbindingPosStrdTimes(l, 1) - tmpboundPosStrdTimes(tmpboundPosStrdTimes(:, 1) < tmpbindingPosStrdTimes(l, 1), 1) < boundFtpt) || ...
                                        any(tmpbindingPosStrdTimes(l, 1) + c.sequenceLen - tmpboundPosStrdTimes(:, 1) < boundFtpt)
                                    displacements(rxnIdxs(j), i) = ...
                                        displacements(rxnIdxs(j), i) + 1;
                                end
                            end
                        end
                    end
                    
                    clear tfs boundPosStrdTimes boundFtpt bindingPosStrdTimes bindingFtpt tmpboundPosStrdTimes tmpbindingPosStrdTimes;
                end
                
                clear states monPosStrndTimes cpxPosStrndTimes monIdxs cpxIdxs;
            end
            
            %% plot data
            clf(figHandle);
            
            meanDisplacements = mean(displacements ./ simEndTimes(ones(size(displacements, 1), 1), :), 2) * 3600;
            stdDisplacements = std(displacements ./ simEndTimes(ones(size(displacements, 1), 1), :), [], 2) * 3600;
            axesHandles = zeros(2, 2);
            
            axesHandle = subplot(2, 2, 1);
            axesHandles(1, 1) = axesHandle;
            hold(axesHandle, 'on');
            means = meanDisplacements(any(c.reactionComplexCatalysisMatrix(:, pc.dnaPolymeraseIndexs), 2));
            stds = stdDisplacements(any(c.reactionComplexCatalysisMatrix(:, pc.dnaPolymeraseIndexs), 2));
            bar(axesHandle, (1:numel(means))', means);
            errorbar((1:numel(means))', means, stds, 'Parent', axesHandle, 'LineStyle', 'none', 'Marker', 'none');
            title(axesHandle, 'DNA Pol-displaced Proteins', 'FontSize', 8);
            xlabel(axesHandle, 'Displaced Protein', 'FontSize', 8);
            ylabel(axesHandle, 'Collisions h^{-1}', 'FontSize', 8);
            set(axesHandle, 'XTick', []);
            ylim(axesHandle, [0 max(means + stds)])
            
            axesHandle = subplot(2, 2, 2);
            axesHandles(1, 2) = axesHandle;
            hold(axesHandle, 'on');
            means = meanDisplacements(any(c.reactionComplexCatalysisMatrix(:, pc.rnaPolymeraseIndexs), 2));
            stds = stdDisplacements(any(c.reactionComplexCatalysisMatrix(:, pc.rnaPolymeraseIndexs), 2));
            bar(axesHandle, (1:numel(means))', means);
            errorbar((1:numel(means))', means, stds, 'Parent', axesHandle, 'LineStyle', 'none', 'Marker', 'none');
            title(axesHandle, 'RNA Pol-displaced Proteins', 'FontSize', 8);
            xlabel(axesHandle, 'Displaced Protein', 'FontSize', 8);
            ylabel(axesHandle, 'Collisions h^{-1}', 'FontSize', 8);
            set(axesHandle, 'XTick', []);
            ylim(axesHandle, [0 max(means + stds)])
            
            axesHandle = subplot(2, 2, 3);
            axesHandles(2, 1) = axesHandle;
            hold(axesHandle, 'on');
            means = meanDisplacements(ismember(c.reactionBoundComplex, pc.rnaPolymeraseIndexs));
            stds = stdDisplacements(ismember(c.reactionBoundComplex, pc.rnaPolymeraseIndexs));
            bar(axesHandle, (1:numel(means))', means);
            errorbar((1:numel(means))', means, stds, 'Parent', axesHandle, 'LineStyle', 'none', 'Marker', 'none');
            title(axesHandle, 'RNA Pol Displacements', 'FontSize', 8);
            xlabel(axesHandle, 'Binding Protein', 'FontSize', 8);
            ylabel(axesHandle, 'Collisions h^{-1}', 'FontSize', 8);
            set(axesHandle, 'XTick', []);
            ylim(axesHandle, [0 max(means + stds)])
            
            axesHandle = subplot(2, 2, 4);
            axesHandles(2, 2) = axesHandle;
            hold(axesHandle, 'on');
            tmp = sum(displacements(ismember(c.reactionBoundComplex, pc.rnaPolymeraseIndexs) & any(c.reactionComplexCatalysisMatrix(:, pc.dnaPolymeraseIndexs), 2), :), 1) ./ simEndTimes * 3600;
            bar(axesHandle, tmp);
            title(axesHandle, 'RNA Pol Displacements by DNA Pol', 'FontSize', 8);
            xlabel(axesHandle, 'Collisions h^{-1}', 'FontSize', 8);
            ylabel(axesHandle, 'Frequency', 'FontSize', 8);
            
            %y-axes
            set(axesHandles, 'FontSize', 6);
            set(cell2mat(get(axesHandles, 'Ylabel')), 'Units', 'normalized')
            PlotUtil.alignYAxesLabels(axesHandles(:, 1));
            PlotUtil.alignYAxesLabels(axesHandles(:, 2));
            
            %% organize figure data
            figData = struct;
            figData.simEndTimes = simEndTimes;
            figData.displacements = displacements;
        end
        
        function figData = macromoleculeHalfLives(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            %% constants
            nSims = numel(selectedSimulations);
            
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            g = sim.gene;
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            pcComp = sum(pc.proteinComplexComposition, 3);
            
            expRNAHalfLives = rna.halfLives(rna.matureIndexs) / 60;
            expMonHalfLives = pm.halfLives(pm.matureIndexs) / 3600;
            expCpxHalfLives = pc.halfLives(pc.matureIndexs) / 3600;
            
            expRNAHalfLives(isinf(expRNAHalfLives)) = NaN;
            expMonHalfLives(isinf(expMonHalfLives)) = NaN;
            expCpxHalfLives(isinf(expCpxHalfLives)) = NaN;
            
            %% get data
            rnaCounts = zeros(numel(rna.matureIndexs), nSims);
            monCounts = zeros(numel(pm.matureIndexs),  nSims);
            cpxCounts = zeros(numel(pc.matureIndexs),  nSims);
            rnaDecays = zeros(numel(rna.matureIndexs), nSims);
            monDecays = zeros(numel(pm.matureIndexs),  nSims);
            cpxDecays = zeros(numel(pc.matureIndexs),  nSims);
            
            stateNames = {
                'Rna'             'counts'  ':'  '-sum'
                'ProteinMonomer'  'counts'  ':'  '-sum'
                'ProteinComplex'  'counts'  ':'  '-sum'
                };
            for i = 1:nSims
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations(i));
                
                if isa(states.Rna.counts, 'edu.stanford.covert.util.SparseMat')
                    [subs, vals] = find(states.Rna.counts);
                    tfs = ...
                        (subs(:, 1) >= rna.nascentIndexs(1) & subs(:, 1) <= rna.nascentIndexs(end)) | ...
                        (subs(:, 1) >= rna.intergenicIndexs(1) & subs(:, 1) <= rna.intergenicIndexs(end));
                    subs = subs(~tfs, :);
                    vals = vals(~tfs, :);
                    
                    subs(subs(:, 1) > rna.intergenicIndexs(end)) = subs(subs(:, 1) > rna.intergenicIndexs(end)) - numel(rna.intergenicIndexs);
                    subs(subs(:, 1) > rna.nascentIndexs(end)) = subs(subs(:, 1) > rna.nascentIndexs(end)) - numel(rna.nascentIndexs);
                    
                    tmp = edu.stanford.covert.util.SparseMat(subs, vals, [6 * numel(rna.matureIndexs) 1 size(states.Rna.counts, 3)]);
                    states.Rna.counts = full(permute(sum(reshape(tmp, [numel(rna.matureIndexs) 6 size(states.Rna.counts, 3)]), 2), [1 3 2 4]));
                    
                    clear tmp tfs subs vals;
                else
                    states.Rna.counts = permute(...
                        + states.Rna.counts(rna.processedIndexs, :, :) ...
                        + states.Rna.counts(rna.matureIndexs, :, :) ...
                        + states.Rna.counts(rna.aminoacylatedIndexs, :, :) ...
                        + states.Rna.counts(rna.boundIndexs, :, :) ...
                        + states.Rna.counts(rna.misfoldedIndexs, :, :) ...
                        + states.Rna.counts(rna.damagedIndexs, :, :)...
                        , [1 3 2]);
                end
                
                states.ProteinMonomer.counts = full(permute(sum(reshape(states.ProteinMonomer.counts, ...
                    [numel(pm.matureIndexs)  numel(pm.wholeCellModelIDs)/numel(pm.matureIndexs)  size(states.ProteinMonomer.counts, 3)]), 2), [1 3 2]));
                states.ProteinComplex.counts = full(permute(sum(reshape(states.ProteinComplex.counts, ...
                    [numel(pc.matureIndexs)  numel(pc.wholeCellModelIDs)/numel(pc.matureIndexs)  size(states.ProteinComplex.counts, 3)]), 2), [1 3 2]));
                
                tmp = states.Rna.counts;
                tmp(setdiff(1:end, rna.matureMRNAIndexs), :) = ...
                    + tmp(setdiff(1:end, rna.matureMRNAIndexs), :) ...
                    + pcComp(setdiff(1:end, g.mRNAIndexs), :) * states.ProteinComplex.counts;
                rnaCounts(:, i) = sum(tmp(:, 1:end-1), 2);
                rnaDecays(:, i) = -sum(min(0, diff(tmp, 1, 2)), 2);
                
                tmp = ...
                    + states.ProteinMonomer.counts ...
                    + pcComp(g.mRNAIndexs, :) * states.ProteinComplex.counts;
                monCounts(:, i) = sum(tmp(:, 1:end-1), 2);
                monDecays(:, i) = -sum(min(0, diff(tmp, 1, 2)), 2);
                
                tmp = states.ProteinComplex.counts;
                cpxCounts(:, i) = sum(tmp(:, 1:end-1), 2);
                cpxDecays(:, i) = -sum(min(0, diff(tmp, 1, 2)), 2);
                
                clear states simEndTimes tmp;
            end
            
            avgRnaHalfLives = 1 / 60   * log(2) ./ (sum(rnaDecays, 2) ./ sum(rnaCounts, 2));
            avgMonHalfLives = 1 / 3600 * log(2) ./ (sum(monDecays, 2) ./ sum(monCounts, 2));
            avgCpxHalfLives = 1 / 3600 * log(2) ./ (sum(cpxDecays, 2) ./ sum(cpxCounts, 2));
            stdRnaHalfLives = std(1 / 60   * log(2) ./ (rnaDecays ./ rnaCounts), [], 2);
            stdMonHalfLives = std(1 / 3600 * log(2) ./ (monDecays ./ monCounts), [], 2);
            stdCpxHalfLives = std(1 / 3600 * log(2) ./ (cpxDecays ./ cpxCounts), [], 2);
            
            avgRnaHalfLives(isinf(avgRnaHalfLives) | isnan(expRNAHalfLives)) = NaN;
            avgMonHalfLives(isinf(avgMonHalfLives) | isnan(expMonHalfLives)) = NaN;
            avgCpxHalfLives(isinf(avgCpxHalfLives) | isnan(expCpxHalfLives)) = NaN;
            stdRnaHalfLives(isinf(stdRnaHalfLives) | isnan(expRNAHalfLives)) = NaN;
            stdMonHalfLives(isinf(stdMonHalfLives) | isnan(expMonHalfLives)) = NaN;
            stdCpxHalfLives(isinf(stdCpxHalfLives) | isnan(expCpxHalfLives)) = NaN;
            
            avgMonHalfLive = nanmean(avgMonHalfLives);
            avgCpxHalfLive = nanmean(avgCpxHalfLives);
            
            clear rnaCounts monCounts cpxCounts rnaDecays monDecays cpxDecays;
            
            %% plot data
            clf(figHandle);
            [axesHandles, xAxesHandles] = PlotUtil.multiElementPlot(figHandle, {100 100 100}, [0 1], struct(...
                'position', [0.02 0.62 0.95 0.2417]));
            
            %RNA
            axesHandle = axesHandles{1};
            xAxesHandle = xAxesHandles{1};
            hold(axesHandle, 'on');
            h = zeros(4, 1);
            h(1) = errorbar(expRNAHalfLives(rna.matureMRNAIndexs), avgRnaHalfLives(rna.matureMRNAIndexs), stdRnaHalfLives(rna.matureMRNAIndexs), ...
                'Parent', axesHandle, 'LineStyle', 'none', 'Color', 'r', 'Marker', '.');
            h(2) = errorbar(expRNAHalfLives(rna.matureRRNAIndexs), avgRnaHalfLives(rna.matureRRNAIndexs), stdRnaHalfLives(rna.matureRRNAIndexs), ...
                'Parent', axesHandle, 'LineStyle', 'none', 'Color', 'g', 'Marker', '.');
            h(3) = errorbar(expRNAHalfLives(rna.matureSRNAIndexs), avgRnaHalfLives(rna.matureSRNAIndexs), stdRnaHalfLives(rna.matureSRNAIndexs), ...
                'Parent', axesHandle, 'LineStyle', 'none', 'Color', 'b', 'Marker', '.');
            h(4) = errorbar(expRNAHalfLives(rna.matureTRNAIndexs), avgRnaHalfLives(rna.matureTRNAIndexs), stdRnaHalfLives(rna.matureTRNAIndexs), ...
                'Parent', axesHandle, 'LineStyle', 'none', 'Color', 'c', 'Marker', '.');
            minVal = min(min(avgRnaHalfLives - stdRnaHalfLives, expRNAHalfLives - 1));
            maxVal = max(max(avgRnaHalfLives - stdRnaHalfLives, expRNAHalfLives + 1));
            xlim(xAxesHandle, [minVal maxVal]);
            xlim(axesHandle, [minVal maxVal]);
            ylim(axesHandle, [minVal maxVal]);
            title(axesHandle, 'RNAs', 'FontSize', 12);
            xlabel(xAxesHandle, 'ExpHalf Life (min)', 'FontSize', 8);
            ylabel(axesHandle, 'Sim Half Life (min)', 'FontSize', 8);
            set(axesHandle, 'FontSize', 6);
            set(axesHandle, 'YTick', get(xAxesHandle, 'XTick'));
            legend(h, {'m', 'r', 's', 't'}, 'Location', 'SouthEast');
            
            %Protein Monomers
            axesHandle = axesHandles{2};
            xAxesHandle = xAxesHandles{2};
            [~, order] = sort(avgMonHalfLives);
            avgMonHalfLives = avgMonHalfLives(order);
            stdMonHalfLives = stdMonHalfLives(order);
            expMonHalfLives = expMonHalfLives(order);
            hold(axesHandle, 'on');
            line([0.5 numel(avgMonHalfLives)+0.5], expMonHalfLives([1 1]), 'Parent', axesHandle, 'Color', 'r');
            line([0.5 numel(avgMonHalfLives)+0.5], avgMonHalfLive([1 1]), 'Parent', axesHandle, 'Color', 'g');
            idxs = find(~isnan(stdMonHalfLives));
            errorbar(idxs, avgMonHalfLives(idxs), stdMonHalfLives(idxs), ...
                'Parent', axesHandle, ...
                'LineStyle', 'none', ...
                'Marker', '.', ...
                'Color', 'b');
            idxs = find(isnan(stdMonHalfLives) & ~isnan(avgMonHalfLives));
            plot(axesHandle, idxs, avgMonHalfLives(idxs), 'b.');
            minVal = max(0, min(min(avgMonHalfLives - stdMonHalfLives, expMonHalfLives - 1)));
            maxVal = max(max(avgMonHalfLives - stdMonHalfLives, expMonHalfLives + 1));
            xlim(xAxesHandle, [0.5 numel(avgMonHalfLives)+0.5]);
            xlim(axesHandle, [0.5 numel(avgMonHalfLives)+0.5]);
            ylim(axesHandle, [minVal maxVal]);
            title(axesHandle, 'Monomers', 'FontSize', 12);
            xlabel(xAxesHandle, 'Monomer', 'FontSize', 8);
            ylabel(axesHandle, 'Sim Half Life (h)', 'FontSize', 8);
            set(axesHandle, 'FontSize', 6);
            
            %Protein Complexes
            axesHandle = axesHandles{3};
            xAxesHandle = xAxesHandles{3};
            [~, order] = sort(avgCpxHalfLives);
            avgCpxHalfLives = avgCpxHalfLives(order);
            stdCpxHalfLives = stdCpxHalfLives(order);
            expCpxHalfLives = expCpxHalfLives(order);
            hold(axesHandle, 'on');
            line([0.5 numel(avgCpxHalfLives)+0.5], expCpxHalfLives([1 1]), 'Parent', axesHandle, 'Color', 'r');
            line([0.5 numel(avgCpxHalfLives)+0.5], avgCpxHalfLive([1 1]), 'Parent', axesHandle, 'Color', 'g');
            idxs = find(~isnan(stdCpxHalfLives));
            errorbar(idxs, avgCpxHalfLives(idxs), stdCpxHalfLives(idxs), ...
                'Parent', axesHandle, ...
                'LineStyle', 'none', ...
                'Marker', '.', ...
                'Color', 'b');
            idxs = find(isnan(stdCpxHalfLives) & ~isnan(avgCpxHalfLives));
            plot(axesHandle, idxs, avgCpxHalfLives(idxs), 'b.');
            minVal = max(0, min(min(avgCpxHalfLives - stdCpxHalfLives, expCpxHalfLives - 1)));
            maxVal = max(max(avgCpxHalfLives - stdCpxHalfLives, expCpxHalfLives + 1));
            xlim(xAxesHandle, [0.5 numel(avgCpxHalfLives)+0.5]);
            xlim(axesHandle, [0.5 numel(avgCpxHalfLives)+0.5]);
            ylim(axesHandle, [minVal maxVal]);
            title(axesHandle, 'Complexes', 'FontSize', 12);
            xlabel(xAxesHandle, 'Complex', 'FontSize', 8);
            ylabel(axesHandle, 'Sim Half Life (h)', 'FontSize', 8);
            set(axesHandle, 'FontSize', 6);
            
            %align
            PlotUtil.offsetYAxes(axesHandles, 0.04);
            PlotUtil.labelSubplots(cell2mat(axesHandles));
            
            %% organize figure data
            figData = struct;
            figData.expRNAHalfLives = expRNAHalfLives;
            figData.avgRnaHalfLives = avgRnaHalfLives;
            figData.stdRnaHalfLives = stdRnaHalfLives;
            figData.avgMonHalfLives = avgMonHalfLives;
            figData.stdMonHalfLives = stdMonHalfLives;
            figData.expMonHalfLives = expMonHalfLives;
            figData.avgCpxHalfLives = avgCpxHalfLives;
            figData.stdCpxHalfLives = stdCpxHalfLives;
            figData.expCpxHalfLives = expCpxHalfLives;
        end
        
        function figData = rnaSynthesisDuration(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            %% constants
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            g = sim.gene;
            comp = sim.compartment;
            c = sim.state('Chromosome');
            r = sim.state('Rna');
            pc = sim.state('ProteinComplex');
            rnaPol = sim.state('RNAPolymerase');
            transcript = sim.state('Transcript');
            
            nSims = numel(selectedSimulations);
            pcComp = sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3);
            
            %% get data
            stateNames = {
                'Rna'             'counts'                   ':'  comp.cytosolIndexs
                'Transcript'      'boundTranscriptionUnits'  ':'  ':'
                'Transcript'      'boundTranscriptProgress'  ':'  ':'
                'RNAPolymerase'   'positionStrands'          ':'  ':'
                'RNAPolymerase'   'states'                   ':'  ':'
                'Ribosome'        'boundMRNAs'               ':'  ':'
                'ProteinComplex'  'counts'                   ':'  '-sum'
                };
            rnaSynthesisTimes = zeros(numel(r.matureIndexs), 3, nSims);
            mrnaFracBounds = zeros(numel(g.mRNAIndexs), 1, nSims);
            rnaFracAminoacylated = zeros(numel(r.matureIndexs), 1, nSims);
            rnaPauses = zeros(4, 1000);
            for i = 1:nSims
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations(i));
                
                states.RNAPolymerase.states = permute(states.RNAPolymerase.states, [1 3 2]);
                states.RNAPolymerase.positionStrands = permute(states.RNAPolymerase.positionStrands, [1 3 2]);
                states.Transcript.boundTranscriptionUnits = permute(states.Transcript.boundTranscriptionUnits, [1 3 2]);
                states.Transcript.boundTranscriptProgress = permute(states.Transcript.boundTranscriptProgress, [1 3 2]);
                states.Ribosome.boundMRNAs = permute(states.Ribosome.boundMRNAs, [1 3 2]);
                
                states.Rna.counts = permute(states.Rna.counts, [1 3 2]);
                
                nascs = full(states.Rna.counts(r.nascentIndexs, :));
                procs = full(states.Rna.counts(r.processedIndexs, :));
                matrs = full(...
                    + states.Rna.counts(r.matureIndexs, :) ...
                    + states.Rna.counts(r.boundIndexs, :) ...
                    + states.Rna.counts(r.misfoldedIndexs, :) ...
                    + states.Rna.counts(r.damagedIndexs, :) ...
                    + states.Rna.counts(r.aminoacylatedIndexs, :));
                matrs(setdiff(1:end, r.matureMRNAIndexs), :) = ...
                    + matrs(setdiff(1:end, r.matureMRNAIndexs), :) ...
                    + pcComp * full(permute(sum(reshape(states.ProteinComplex.counts, ...
                    [numel(pc.matureIndexs)  numel(pc.wholeCellModelIDs)/numel(pc.matureIndexs)  size(states.ProteinComplex.counts, 3)]), 2), [1 3 2]));
                amincs = full(sum(states.Rna.counts(r.aminoacylatedIndexs, :), 2));
                
                %transcription times
                [startTimes, startPols] = find(cat(2, ...
                    states.RNAPolymerase.states(:, 1) >= rnaPol.activelyTranscribingValue, ...
                    states.RNAPolymerase.states(:, 2:end)   >= rnaPol.activelyTranscribingValue & ...
                    states.RNAPolymerase.states(:, 1:end-1) < rnaPol.activelyTranscribingValue)');
                [endTimes, endPols] = find(cat(2, ...
                    states.RNAPolymerase.states(:, 1:end-1) >= rnaPol.activelyTranscribingValue & ...
                    states.RNAPolymerase.states(:, 2:end)   < rnaPol.activelyTranscribingValue, ...
                    states.RNAPolymerase.states(:, end) >= rnaPol.activelyTranscribingValue)');
                assert(all(startPols == endPols));
                assert(all(endTimes >= startTimes));
                boundTUs = states.Transcript.boundTranscriptionUnits(sub2ind(...
                    size(states.Transcript.boundTranscriptionUnits), ...
                    startPols, startTimes));
                rnaSynthesisTimes(1:numel(r.nascentIndexs), 1, i) = ...
                    accumarray(boundTUs, endTimes - startTimes + 1, [numel(r.nascentIndexs) 1]) ./ histc(boundTUs, (1:numel(r.nascentIndexs))');
                
                %transcription pauses
                pause = states.RNAPolymerase.states(:, 1:end-1) >= rnaPol.activelyTranscribingValue & diff(states.Transcript.boundTranscriptProgress, 1, 2) == 0;
                [startTimes, startPols] = find([pause(:, 1)  diff(pause, 1, 2) > 0]');
                [endTimes, endPols] = find([diff(pause, 1, 2) < 0  pause(:, end)]');
                assert(all(startPols == endPols));
                assert(all(endTimes >= startTimes));
                
                [~, bases] = ismember(c.sequence.subsequence(states.RNAPolymerase.positionStrands(sub2ind(...
                    size(states.RNAPolymerase.positionStrands), ...
                    [startPols startPols], ...
                    [startTimes startTimes], ...
                    [ones(size(startPols, 1), 1) 2 * ones(size(startPols, 1), 1)]))), 'ACGT');
                
                if max(endTimes - startTimes + 1) > size(rnaPauses, 2)
                    rnaPauses = [rnaPauses zeros(4, max(endTimes - startTimes + 1) - size(rnaPauses, 2))]; %#ok<AGROW>
                end
                rnaPauses = ...
                    + rnaPauses...
                    + accumarray([bases  endTimes - startTimes + 1], ones(size(startTimes)), [4 size(rnaPauses, 2)]);
                
                %processing time
                rnaSynthesisTimes(1:numel(r.nascentIndexs), 2, i) = log(2) ./ (sum(max(0, -diff(nascs, 1, 2)), 2) ./ sum(nascs(:, 1:end-1), 2));
                
                %modification time
                rnaSynthesisTimes(:, 3, i) = log(2) ./ (sum(max(0, -diff(procs, 1, 2)), 2) ./ sum(procs(:, 1:end-1), 2));
                
                %fraction time bound
                [~, j, val] = find(states.Ribosome.boundMRNAs);
                tmp = zeros(numel(g.mRNAIndexs), size(states.Ribosome.boundMRNAs, 2));
                tmp(sub2ind(size(tmp), val, j)) = 1;
                mrnaFracBounds(:, 1, i) = sum(tmp, 2) ./ sum(r.matureRNAGeneComposition(g.mRNAIndexs, :) * matrs > 0, 2);
                
                %fraction aminoacylated bound
                rnaFracAminoacylated(:, 1, i) = amincs ./ sum(matrs, 2);
                
                clear states nascs procs matrs amincs startTimes startPols endTimes endPols boundTUs pause j val tmp;
            end
            
            rnaSynthesisTimes(isinf(rnaSynthesisTimes)) = NaN;
            mrnaFracBounds(isinf(mrnaFracBounds)) = NaN;
            rnaFracAminoacylated(isinf(rnaFracAminoacylated)) = NaN;
            
            %% plot
            clf(figHandle);
            
            %synthesis times
            axesHandle = subplot(3, 2, 1, 'Parent', figHandle);
            axesHandles(1, 1) = axesHandle;
            tmp = [zeros(numel(r.matureIndexs), 2)  nanmean(rnaSynthesisTimes(:, 3, :), 3)];
            for i = 1:numel(r.nascentIndexs)
                tmpTfs = r.nascentRNAMatureRNAComposition(:, i) == 1;
                tmp(tmpTfs, 1) = nanmean(rnaSynthesisTimes(i, 1, :), 3);
                tmp(tmpTfs, 2) = nanmean(rnaSynthesisTimes(i, 2, :), 3);
            end
            [~, order] = sort(nansum(tmp, 2));
            h = bar(axesHandle, 1:numel(r.matureIndexs),  tmp(order, :), 'stacked');
            xlim(axesHandle, [0.5 numel(r.matureIndexs) + 0.5]);
            ylim(axesHandle, [0 max(nansum(tmp, 2))]);
            legend(h, {'Transcription', 'Processing', 'Modification'}, 'Location', 'NorthWest');
            xlabel(axesHandle, 'RNA', 'FontSize', 8);
            ylabel(axesHandle, {'Synthesis' 'Duration (s)'}, 'FontSize', 8);
            
            %transcription time vs length
            axesHandle = subplot(3, 2, 3, 'Parent', figHandle);
            axesHandles(2, 1) = axesHandle;
            meanVal = nanmean(rnaSynthesisTimes(1:numel(r.nascentIndexs), 1, :), 3);
            stdVal = nanstd(rnaSynthesisTimes(1:numel(r.nascentIndexs), 1, :), [], 3);
            errorbar(r.lengths(r.nascentIndexs), meanVal, stdVal, ...
                'LineStyle', 'none', 'Parent', axesHandle);
            xlim(axesHandle, [min(r.lengths(r.nascentIndexs)) max(r.lengths(r.nascentIndexs))]);
            ylim(axesHandle, [min(meanVal - stdVal) max(meanVal + stdVal)]);
            xlabel(axesHandle, 'Length', 'FontSize', 8);
            ylabel(axesHandle, {'Transcription' 'Duration (s)'}, 'FontSize', 8)
            
            %transcription time vs G/C content
            axesHandle = subplot(3, 2, 5, 'Parent', figHandle);
            axesHandles(3, 1) = axesHandle;
            gcContent = cellfun(@(seq) sum(ismember(seq, 'GC')) / numel(seq), transcript.transcriptionUnitSequences(r.nascentIndexs));
            meanVal = nanmean(rnaSynthesisTimes(1:numel(r.nascentIndexs), 1, :), 3);
            stdVal = nanstd(rnaSynthesisTimes(1:numel(r.nascentIndexs), 1, :), [], 3);
            errorbar(gcContent, ...
                nanmean(rnaSynthesisTimes(1:numel(r.nascentIndexs), 1, :), 3), ...
                nanstd(rnaSynthesisTimes(1:numel(r.nascentIndexs), 1, :), [], 3), ...
                'LineStyle', 'none', 'Parent', axesHandle);
            xlim(axesHandle, [min(gcContent) max(gcContent)]);
            ylim(axesHandle, [min(meanVal - stdVal) max(meanVal + stdVal)]);
            xlabel(axesHandle, 'G/C Content', 'FontSize', 8);
            ylabel(axesHandle, {'Transcription' 'Duration'}, 'FontSize', 8)
            
            %pauses
            axesHandle = subplot(3, 2, 2, 'Parent', figHandle);
            axesHandles(1, 2) = axesHandle;
            tmp = rnaPauses ./ repmat(sum(rnaPauses, 2), 1, size(rnaPauses, 2));
            h = bar(axesHandle, 1:size(rnaPauses, 2), tmp', 'group');
            xlim(axesHandle, [0.5 find(any(tmp > 0.01, 1), 1, 'last') + 0.5]);
            legend(h, {'A', 'C', 'G', 'T'});
            ylim(axesHandle, [0 max(tmp(:))]);
            xlabel(axesHandle, 'Pause Duration (s)', 'FontSize', 8);
            ylabel(axesHandle, 'Freq', 'FontSize', 8)
            
            %fraction bound
            axesHandle = subplot(3, 2, 4, 'Parent', figHandle);
            axesHandles(2, 2) = axesHandle;
            [~, order] = sort(nanmean(mrnaFracBounds, 3));
            errorbar((1:numel(g.mRNAIndexs))', ...
                nanmean(mrnaFracBounds(order, :, :), 3), ...
                nanstd(mrnaFracBounds(order, :, :), [], 3), ...
                'LineStyle', 'none', 'Parent', axesHandle);
            xlim(axesHandle, [0.5 numel(g.mRNAIndexs) + 0.5])
            ylim(axesHandle, [
                min(nanmean(mrnaFracBounds, 3) - nanstd(mrnaFracBounds, [], 3))  ...
                max(nanmean(mrnaFracBounds, 3) + nanstd(mrnaFracBounds, [], 3))
                ]);
            xlabel(axesHandle, 'mRNA', 'FontSize', 8)
            ylabel(axesHandle, {'Frac' 'Bound'}, 'FontSize', 8);
            
            %fraction aminoacylated
            axesHandle = subplot(3, 2, 6, 'Parent', figHandle);
            axesHandles(3, 2) = axesHandle;
            expr = r.expression(r.matureIndexs(r.matureTRNAIndexs));
            expr = expr / sum(expr);
            errorbar(expr, ...
                nanmean(rnaFracAminoacylated(r.matureTRNAIndexs, :, :), 3), ...
                nanstd(rnaFracAminoacylated(r.matureTRNAIndexs, :, :), [], 3), ...
                'LineStyle', 'none', 'Parent', axesHandle);
            xlim(axesHandle, [min(expr) max(expr)]);
            ylim(axesHandle, [
                min(nanmean(rnaFracAminoacylated(r.matureTRNAIndexs, :, :), 3) - nanstd(rnaFracAminoacylated(r.matureTRNAIndexs, :, :), [], 3)) ...
                max(nanmean(rnaFracAminoacylated(r.matureTRNAIndexs, :, :), 3) + nanstd(rnaFracAminoacylated(r.matureTRNAIndexs, :, :), [], 3))
                ]);
            xlabel(axesHandle, 'tRNA Exp', 'FontSize', 8);
            ylabel(axesHandle, {'Frac' 'Aminoacylated'}, 'FontSize', 8);
            
            %align axes labels
            set(axesHandles, 'FontSize', 6);
            PlotUtil.alignYAxesLabels(axesHandles(:, 1));
            PlotUtil.alignYAxesLabels(axesHandles(:, 2));
            
            %% organize figure data
            figData = struct;
            figData.rnaSynthesisTimes = rnaSynthesisTimes;
            figData.mrnaFracBounds = mrnaFracBounds;
            figData.rnaFracAminoacylated = rnaFracAminoacylated;
            figData.rnaPauses = rnaPauses;
        end
        
        function figData = proteinSynthesisDuration(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            %% constants
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            g = sim.gene;
            comp = sim.compartment;
            m = sim.state('Metabolite');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            rib = sim.state('Ribosome');
            pol = sim.state('Polypeptide');
            
            nSims = numel(selectedSimulations);
            pcComp = sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3);
            
            %% get data
            stateNames = {
                'ProteinMonomer'  'counts'                   ':'  ':'
                'ProteinComplex'  'counts'                   ':'  '-sum'
                'Ribosome'        'states'                   ':'  ':'
                'Ribosome'        'boundMRNAs'               ':'  ':'
                'Ribosome'        'mRNAPositions'            ':'  ':'
                'Ribosome'        'tmRNAPositions'           ':'  ':'
                };
            synthesisTimes = zeros(numel(pm.matureIndexs), 6, nSims);
            fracBounds = zeros(numel(pm.matureIndexs), 1, nSims);
            trnaPauses = zeros(36, 1000);
            aaPauses = zeros(21, 1000);
            for i = 1:nSims
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations(i));
                states.Ribosome.states = permute(states.Ribosome.states, [1 3 2]);
                states.Ribosome.boundMRNAs = permute(states.Ribosome.boundMRNAs, [1 3 2]);
                states.Ribosome.mRNAPositions = permute(states.Ribosome.mRNAPositions, [1 3 2]);
                states.Ribosome.tmRNAPositions = permute(states.Ribosome.tmRNAPositions, [1 3 2]);
                
                monCounts = permute(sum(states.ProteinMonomer.counts, 2), [1 3 2]);
                cpxCounts = full(permute(sum(reshape(states.ProteinComplex.counts, ...
                    [numel(pc.matureIndexs)  numel(pc.wholeCellModelIDs)/numel(pc.matureIndexs)  size(states.ProteinComplex.counts, 3)]), 2), [1 3 2]));
                nascs = full(monCounts(pm.nascentIndexs, :));
                procsI = zeros(numel(pm.matureIndexs), size(states.ProteinMonomer.counts, 3));
                procsI(pm.compartments(pm.processedIIIndexs) ~= comp.cytosolIndexs, :) = ...
                    full(permute(sum(states.ProteinMonomer.counts(...
                    pm.processedIIndexs(pm.compartments(pm.processedIIIndexs) ~= comp.cytosolIndexs), ...
                    comp.cytosolIndexs, :), 2), [1 3 2]));
                trans = full(monCounts(pm.processedIIndexs, :)) - procsI;
                procsII = full(monCounts(pm.processedIIIndexs, :));
                folds = full(monCounts(pm.foldedIndexs, :));
                signals = full(monCounts(pm.signalSequenceIndexs, :));
                matrs = full(permute(sum(reshape(permute(monCounts, [1 3 2]), ...
                    [numel(pm.matureIndexs)  numel(pm.wholeCellModelIDs)/numel(pm.matureIndexs)  size(states.ProteinMonomer.counts, 3)]), 2), [1 3 2])) ...
                    - nascs ...
                    - procsI ...
                    - trans ...
                    - procsII ...
                    - folds ...
                    - signals ...
                    + pcComp * cpxCounts;
                bounds = ...
                    + monCounts(pm.boundIndexs, :) ...
                    + pcComp * full(permute(states.ProteinComplex.counts(pc.boundIndexs, :, :), [1 3 2]));
                
                %synthesis times
                [startTimes, startRibs] = find(cat(2, ...
                    states.Ribosome.mRNAPositions(:, 1) > 0 | ...
                    states.Ribosome.tmRNAPositions(:, 1) > 0, ...
                    states.Ribosome.mRNAPositions(:, 1:end-1) == 0 & ...
                    states.Ribosome.tmRNAPositions(:, 1:end-1) == 0 & ...
                    (states.Ribosome.mRNAPositions(:, 2:end) > 0 | ...
                    states.Ribosome.tmRNAPositions(:, 2:end) > 0))');
                [endTimes, endRibs] = find(cat(2, ...
                    (states.Ribosome.mRNAPositions(:, 1:end-1) > 0 | ...
                    states.Ribosome.tmRNAPositions(:, 1:end-1) > 0) & ...
                    states.Ribosome.mRNAPositions(:, 2:end) == 0 & ...
                    states.Ribosome.tmRNAPositions(:, 2:end) == 0, ...
                    states.Ribosome.mRNAPositions(:, end) > 0 | ...
                    states.Ribosome.tmRNAPositions(:, end) > 0)');
                assert(all(startRibs == endRibs));
                assert(all(endTimes >= startTimes));
                boundMRNAs = states.Ribosome.boundMRNAs(sub2ind(...
                    size(states.Ribosome.boundMRNAs), ...
                    startRibs, startTimes));
                synthesisTimes(:, 1, i) = ...
                    accumarray(boundMRNAs, endTimes - startTimes + 1, ...
                    [numel(pm.matureIndexs) 1]) ./ histc(boundMRNAs, (1:numel(pm.matureIndexs))');
                
                %translation pauses
                pause = ...
                    states.Ribosome.states(:, 1:end-1) ~= rib.notExistValue & ...
                    diff(states.Ribosome.mRNAPositions + states.Ribosome.tmRNAPositions, 1, 2) == 0;
                [startTimes, startRibs] = find([pause(:, 1)  diff(pause, 1, 2) > 0]');
                [endTimes, endRibs] = find([diff(pause, 1, 2) < 0  pause(:, end)]');
                assert(all(startRibs == endRibs));
                assert(all(endTimes >= startTimes));
                
                for j = 1:numel(startRibs)
                    if ...
                            states.Ribosome.mRNAPositions(startRibs(j),  startTimes(j)) == pol.monomerLengths(states.Ribosome.boundMRNAs(startRibs(j),  startTimes(j))) || ...
                            states.Ribosome.tmRNAPositions(startRibs(j), startTimes(j)) == pol.proteolysisTagLength
                        continue;
                    end
                    
                    if states.Ribosome.tmRNAPositions(startRibs(j), startTimes(j)) == 0
                        aa = pol.monomerAASequences{states.Ribosome.boundMRNAs(startRibs(j), startTimes(j))}(states.Ribosome.mRNAPositions(startRibs(j),  startTimes(j))+1);
                        tnra = pol.monomerTRNASequences{states.Ribosome.boundMRNAs(startRibs(j), startTimes(j))}(states.Ribosome.mRNAPositions(startRibs(j),  startTimes(j))+1);
                    else
                        aa = pol.proteolysisTagAASequence(states.Ribosome.tmRNAPositions(startRibs(j), startTimes(j))+1);
                        tnra = pol.proteolysisTagTRNASequence(states.Ribosome.tmRNAPositions(startRibs(j), startTimes(j))+1);
                    end
                    if states.Ribosome.mRNAPositions(startRibs(j), startTimes(j)) == 0
                        aaIdx = 21;
                    else
                        aaIdx = find(aa == edu.stanford.covert.cell.kb.ProteinMonomer.bases);
                    end
                    if endTimes(j) - startTimes(j) + 1 > size(aaPauses, 2)
                        aaPauses = [aaPauses zeros(size(aaPauses, 1), (endTimes(j) - startTimes(j) + 1) - size(aaPauses, 2))]; %#ok<AGROW>
                        trnaPauses = [trnaPauses zeros(size(trnaPauses, 1), (endTimes(j) - startTimes(j) + 1) - size(trnaPauses, 2))]; %#ok<AGROW>
                    end
                    aaPauses(aaIdx, endTimes(j) - startTimes(j) + 1) = ...
                        aaPauses(aaIdx, endTimes(j) - startTimes(j) + 1) + 1;
                    trnaPauses(tnra, endTimes(j) - startTimes(j) + 1) = ...
                        trnaPauses(tnra, endTimes(j) - startTimes(j) + 1) + 1;
                end
                
                %processing I time
                synthesisTimes(:, 2, i) = log(2) ./ (sum(max(0, -diff(nascs, 1, 2)), 2) ./ sum(nascs(:, 1:end-1), 2));
                
                %translocation time
                synthesisTimes(:, 3, i) = log(2) ./ (sum(max(0, -diff(procsI, 1, 2)), 2) ./ sum(procsI(:, 1:end-1), 2));
                
                %processing II time
                synthesisTimes(:, 4, i) = log(2) ./ (sum(max(0, -diff(trans, 1, 2)), 2) ./ sum(trans(:, 1:end-1), 2));
                
                %folding time
                synthesisTimes(:, 5, i) = log(2) ./ (sum(max(0, -diff(procsII, 1, 2)), 2) ./ sum(procsII(:, 1:end-1), 2));
                
                %modification time
                synthesisTimes(:, 6, i) = log(2) ./ (sum(max(0, -diff(folds, 1, 2)), 2) ./ sum(folds(:, 1:end-1), 2));
                
                %fraction bound
                fracBounds(:, 1, i) = sum(bounds, 2) ./ sum(matrs, 2);
                
                clear states monCounts cpxCounts nascs procsI trans procsII folds signals matrs bounds;
                clear boundMRNAs startTimes endTimes startRibs endRibs pause;
            end
            
            synthesisTimes(isinf(synthesisTimes)) = NaN;
            fracBounds(isinf(fracBounds)) = NaN;
            
            %% plot
            clf(figHandle);
            axesHandles = zeros(3, 2);
            
            %synthesis times
            axesHandle = subplot(3, 2, 1, 'Parent', figHandle);
            axesHandles(1, 1) = axesHandle;
            tmp = nanmean(synthesisTimes, 3);
            [~, order] = sort(nansum(tmp, 2));
            h = bar(axesHandle, 1:numel(pm.matureIndexs),  tmp(order, :), 'stacked');
            xlim(axesHandle, [0.5 numel(pm.matureIndexs) + 0.5]);
            ylim(axesHandle, [0 max(nansum(tmp, 2))]);
            legend(h, {'Translation', 'Processing-I', 'Translocation', 'Processing-II', 'Folding', 'Modification'}, ...
                'Location', 'NorthWest', 'FontSize', 4);
            xlabel(axesHandle, 'Monomer', 'FontSize', 8);
            ylabel(axesHandle, {'Synthesis' 'Duration (s)'}, 'FontSize', 8);
            
            %transcription time vs length
            axesHandle = subplot(3, 2, 3, 'Parent', figHandle);
            axesHandles(2, 1) = axesHandle;
            meanVal = nanmean(synthesisTimes(:, 1, :), 3);
            stdVal = nanstd(synthesisTimes(:, 1, :), [], 3);
            errorbar(pm.lengths(pm.nascentIndexs), meanVal, stdVal, ...
                'LineStyle', 'none', 'Parent', axesHandle);
            xlim(axesHandle, [min(pm.lengths(pm.nascentIndexs)) max(pm.lengths(pm.nascentIndexs))]);
            ylim(axesHandle, [min(meanVal - stdVal) max(meanVal + stdVal)]);
            xlabel(axesHandle, 'Length', 'FontSize', 8);
            ylabel(axesHandle, {'Translation' 'Duration (s)'}, 'FontSize', 8)
            
            %transcription time vs G/C content
            axesHandle = subplot(3, 2, 5, 'Parent', figHandle);
            axesHandles(3, 1) = axesHandle;
            [~, idx] = min(m.biomassComposition(m.aminoAcidIndexs(1:20)));
            aaContent = pm.baseCounts(pm.nascentIndexs, m.aminoAcidIndexs(idx)) ./ pm.lengths(pm.nascentIndexs);
            meanVal = nanmean(synthesisTimes(:, 1, :), 3);
            stdVal = nanstd(synthesisTimes(:, 1, :), [], 3);
            errorbar(aaContent, meanVal, stdVal, ...
                'LineStyle', 'none', 'Parent', axesHandle);
            xlim(axesHandle, [min(aaContent) max(aaContent)]);
            ylim(axesHandle, [min(meanVal - stdVal) max(meanVal + stdVal)]);
            xlabel(axesHandle, [m.names{m.aminoAcidIndexs(idx)} ' Content'], 'FontSize', 8);
            ylabel(axesHandle, {'Translation' 'Duration'}, 'FontSize', 8)
            
            %pauses -- tRNA
            axesHandle = subplot(3, 2, 2, 'Parent', figHandle);
            axesHandles(1, 2) = axesHandle;
            tmp = trnaPauses ./ repmat(sum(trnaPauses, 2), 1, size(trnaPauses, 2));
            bar(axesHandle, 1:size(trnaPauses, 2), tmp', 'group');
            xlim(axesHandle, [0.5 find(any(tmp > 0.01, 1), 1, 'last')+0.5]);
            ylim(axesHandle, [0 max(tmp(:))]);
            xlabel(axesHandle, 'Pause Duration (s)', 'FontSize', 8);
            ylabel(axesHandle, 'Freq', 'FontSize', 8)
            
            %pauses -- AA
            axesHandle = subplot(3, 2, 4, 'Parent', figHandle);
            axesHandles(2, 2) = axesHandle;
            tmp = aaPauses ./ repmat(sum(aaPauses, 2), 1, size(aaPauses, 2));
            bar(axesHandle, 1:size(aaPauses, 2), tmp', 'group');
            xlim(axesHandle, [0.5 find(any(tmp > 0.01, 1), 1, 'last') + 0.5]);
            ylim(axesHandle, [0 max(tmp(:))]);
            xlabel(axesHandle, 'Pause Duration (s)', 'FontSize', 8);
            ylabel(axesHandle, 'Freq', 'FontSize', 8)
            
            %fraction bound
            axesHandle = subplot(3, 2, 6, 'Parent', figHandle);
            axesHandles(3, 2) = axesHandle;
            [~, order] = sort(nanmean(fracBounds, 3));
            order = order(any(fracBounds(order, :) > 0, 2));
            errorbar((1:numel(order))', ...
                nanmean(fracBounds(order, :, :), 3), ...
                nanstd(fracBounds(order, :, :), [], 3), ...
                'LineStyle', 'none', 'Parent', axesHandle);
            xlim(axesHandle, [0.5 numel(order) + 0.5])
            ylim(axesHandle, [
                min(nanmean(fracBounds, 3) - nanstd(fracBounds, [], 3))  ...
                max(nanmean(fracBounds, 3) + nanstd(fracBounds, [], 3))
                ]);
            xlabel(axesHandle, 'Monomer', 'FontSize', 8)
            ylabel(axesHandle, {'Frac' 'Bound'}, 'FontSize', 8);
            
            %align axes labels
            set(axesHandles, 'FontSize', 6);
            PlotUtil.alignYAxesLabels(axesHandles(:, 1));
            PlotUtil.alignYAxesLabels(axesHandles(:, 2));
            
            %% organize figure data
            figData = struct;
            figData.synthesisTimes = synthesisTimes;
            figData.fracBounds = fracBounds;
            figData.trnaPauses = trnaPauses;
            figData.aaPauses = aaPauses;
        end
        
        function figData = macromolecularComplexes(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            %% constants
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            g = sim.gene;
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            pcComp = sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3);
            monIdxs = find(any(pcComp, 2));
            pcComp = pcComp(monIdxs, :);
            
            %% load data
            stateNames = {
                'Geometry'       'volume'  ':'  ':'
                'ProteinMonomer' 'counts'  ':'  '-sum'
                'ProteinComplex' 'counts'  ':'  '-sum'
                };
            monConc = zeros(numel(monIdxs), numel(selectedSimulations));
            cpxConc = zeros(numel(pc.matureIndexs), numel(selectedSimulations));
            for i = 1:numel(selectedSimulations)
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations(i));
                states.ProteinMonomer.counts = full(sum(reshape(states.ProteinMonomer.counts, ...
                    [numel(pm.matureIndexs)  numel(pm.wholeCellModelIDs)/numel(pm.matureIndexs)  size(states.ProteinMonomer.counts, 3)]), 2));
                states.ProteinMonomer.counts = states.ProteinMonomer.counts(monIdxs, :, :);
                states.ProteinComplex.counts = full(sum(reshape(states.ProteinComplex.counts, ...
                    [numel(pc.matureIndexs)  numel(pc.wholeCellModelIDs)/numel(pc.matureIndexs)  size(states.ProteinComplex.counts, 3)]), 2));
                
                monConc(:, i) = mean((...
                    + permute(states.ProteinMonomer.counts, [1 3 2]) ...
                    + pcComp * permute(states.ProteinComplex.counts, [1 3 2])) ...
                    ./ permute(states.Geometry.volume(ones(numel(monIdxs), 1), :, :) / ConstantUtil.nAvogadro * 1e6, [1 3 2]), ...
                    2);
                cpxConc(:, i) = mean((...
                    permute(states.ProteinComplex.counts, [1 3 2])) ...
                    ./ permute(states.Geometry.volume(ones(numel(pc.matureIndexs), 1), :, :) / ConstantUtil.nAvogadro * 1e6, [1 3 2]), ...
                    2);
                
                clear states;
            end
            
            meanMonConc = mean(monConc, 2);
            stdMonConc = std(monConc, [], 2);
            meanCpxConc = pcComp * mean(cpxConc, 2);
            stdCpxConc = pcComp * std(cpxConc, [], 2);
            
            clear monConc cpxConc;
            
            %% plot
            clf(figHandle);
            [axesHandle, xAxesHandle] = PlotUtil.multiElementPlot(figHandle, 50, [min(meanCpxConc-stdCpxConc) max(meanCpxConc+stdCpxConc)], ...
                struct('xlabelStr', 'Complex Concentration ({\mu}M)'));
            xlabel(xAxesHandle, 'Complex Concentration ({\mu}M)', 'FontSize', 12, 'Units', 'normalized');
            
            PlotUtil.plotBubbles(axesHandle, [meanCpxConc stdCpxConc/3], [meanMonConc stdMonConc/3]);
            
            set(axesHandle, 'YTick', [0 100]);
            ylim(axesHandle, [min(meanMonConc-stdMonConc) max(meanMonConc+stdMonConc)]);
            ylabel(axesHandle, 'Monomer Concentration ({\mu}M)', 'FontSize', 12);
            
            PlotUtil.offsetYAxes(axesHandle, 0.02);
            
            %% organize figure data
            figData = struct;
            figData.meanMonConc = meanMonConc;
            figData.stdMonConc = stdMonConc;
            figData.meanCpxConc = meanCpxConc;
            figData.stdCpxConc = stdCpxConc;
        end
        
        function figData = replication(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            %% constants
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            comp = sim.compartment;
            c = sim.state('Chromosome');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            rep = sim.process('Replication');
            
            cpxIdxs = [
                pc.matureIndexs(rep.enzymeGlobalIndexs([rep.enzymeIndexs_ssb4mer; rep.enzymeIndexs_ssb8mer]))
                pc.boundIndexs(rep.enzymeGlobalIndexs([rep.enzymeIndexs_ssb4mer; rep.enzymeIndexs_ssb8mer]))
                ]';
            
            nSims = numel(selectedSimulations);
            
            %% get data
            stateNames = {
                'mass'         'Mass (fg)'
                'growth_rate'  'Growth (fg h^{-1})'
                'leadingPol1'  'Leading Pol 1'
                'leadingPol2'  'Leading Pol 2'
                'laggingPol1'  'Lagging Pol 1'
                'laggingPol2'  'Lagging Pol 2'
                };
            ensemble = SimulationEnsemble(simBatchDir, stateNames, [], selectedSimulations);
            mass = permute(ensemble.stateData.values(ensemble.getPropertyIndices('mass'), :, :, :) * sim.state('Mass').cellInitialDryWeight / (1-sim.state('Mass').fractionWetWeight) * 1e15, [4 3 1 2]);
            growth = permute(ensemble.stateData.values(ensemble.getPropertyIndices('growth_rate'), :, :, :) * 3600 * sim.state('Mass').cellInitialDryWeight / (1 - sim.state('Mass').fractionWetWeight) * 1e15, [4 3 1 2]);
            ledPol1 = permute(ensemble.stateData.values(ensemble.getPropertyIndices('leadingPol1'), :, :, :), [4 3 1 2]);
            ledPol2 = permute(ensemble.stateData.values(ensemble.getPropertyIndices('leadingPol2'), :, :, :), [4 3 1 2]);
            lagPol1 = permute(ensemble.stateData.values(ensemble.getPropertyIndices('laggingPol1'), :, :, :), [4 3 1 2]);
            lagPol2 = permute(ensemble.stateData.values(ensemble.getPropertyIndices('laggingPol2'), :, :, :), [4 3 1 2]);
            
            stateNames = {
                'ProteinMonomer' 'counts' pm.matureIndexs(rep.enzymeGlobalIndexs(rep.enzymeIndexs_ligase)) comp.cytosolIndexs
                'ProteinComplex' 'counts' cpxIdxs comp.cytosolIndexs
                'Chromosome'     'strandBreaks' ':' ':'
                };
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations);
            ssbs = permute(states.ProteinComplex.counts, [4 3 1 2]);
            [~, ligationTimes] = find(sum(states.Chromosome.strandBreaks, 3));
            
            repStarts = NaN(nSims, 1);
            repEnds = NaN(nSims, 1);
            kineticallyLimitedEnds = NaN(nSims, 1);
            repStartMass = NaN(nSims, 1);
            repStartGrowth = NaN(nSims, 1);
            repStartSSbs = NaN(nSims, 1);
            leadLeadCorr = NaN(nSims, 1);
            okFragSynTimes = [];
            for i = 1:nSims
                if ~all(isnan(ledPol1(i, :)))
                    repStarts(i) = find(~isnan(ledPol1(i, :)), 1, 'first');
                    repEnds(i) = find(~isnan(ledPol1(i, :)), 1, 'last');
                    repStartMass(i) = mass(i, repStarts(i));
                    repStartGrowth(i) = growth(i, repStarts(i));
                    repStartSSbs(i) = ssbs(i, repStarts(i) - 1);
                    
                    tmp = find(conv(diff(ledPol2(i, :), 1, 2), ones(1, 100), 'same') > 0.75 * 100 * rep.dnaPolymeraseElongationRate, 1, 'last') + 1;
                    if ~isempty(tmp)
                        kineticallyLimitedEnds(i) = tmp;
                    else
                        kineticallyLimitedEnds(i) = repEnds(i);
                    end
                    
                    leadLeadCorr(i) = -corr(ledPol1(i, ~isnan(ledPol1(i, :)))', ledPol2(i, ~isnan(ledPol1(i, :)))');
                    
                    okFragSynTimes = [okFragSynTimes diff(find(diff(lagPol1(i, :), 1, 2) < 0)) diff(find(diff(lagPol2(i, :), 1, 2) > 0))]; %#ok<AGROW>
                end
            end
            
            clear ensemble states mass growth ssbs;
            
            %% plot data
            clf(figHandle);
            axesHandles = zeros(5, 2);
            
            axesHandle = subplot(5, 2, 1, 'Parent', figHandle);
            axesHandles(1, 1) = axesHandle;
            h = plot(axesHandle, repStartMass, [(repEnds - repStarts + 1)  (kineticallyLimitedEnds - repStarts + 1)  (repEnds - kineticallyLimitedEnds + 1)] / 3600, '.');
            if any(repStartMass) && range(repStartMass)
                xlim(axesHandle, [min(repStartMass) max(repStartMass)])
            end
            if any(repEnds - repStarts) && range(repEnds - repStarts)
                ylim(axesHandle, [min(repEnds - repStarts + 1) max(repEnds - repStarts + 1)] / 3600)
            end
            legend(h, {'Replication', 'Kin Lim', 'Met Lim'}, 'FontSize', 6)
            xlabel(axesHandle, 'Mass (fg)', 'FontSize', 8);
            ylabel(axesHandle, 'Duration (h)', 'FontSize', 8);
            
            axesHandle = subplot(5, 2, 3, 'Parent', figHandle);
            axesHandles(2, 1) = axesHandle;
            plot(axesHandle, repStartGrowth, [(repEnds - repStarts + 1)  (kineticallyLimitedEnds - repStarts + 1)  (repEnds - kineticallyLimitedEnds + 1)] / 3600, '.');
            if any(repStartGrowth) && range(repStartGrowth)
                xlim(axesHandle, [min(repStartGrowth) max(repStartGrowth)])
            end
            if any(repEnds - repStarts) && range(repEnds - repStarts)
                ylim(axesHandle, [min(repEnds - repStarts + 1) max(repEnds - repStarts + 1)] / 3600)
            end
            xlabel(axesHandle, 'Growth (fg h^{-1})', 'FontSize', 8);
            ylabel(axesHandle, 'Duration (h)', 'FontSize', 8);
            
            axesHandle = subplot(5, 2, 5, 'Parent', figHandle);
            axesHandles(3, 1) = axesHandle;
            plot(axesHandle, repStartSSbs, [(repEnds - repStarts + 1)  (kineticallyLimitedEnds - repStarts + 1)  (repEnds - kineticallyLimitedEnds + 1)] / 3600, '.');
            if any(repStartSSbs) && range(repStartSSbs)
                xlim(axesHandle, [min(repStartSSbs) max(repStartSSbs)])
            end
            if any(repEnds - repStarts) && range(repEnds - repStarts)
                ylim(axesHandle, [min(repEnds - repStarts + 1) max(repEnds - repStarts + 1)] / 3600)
            end
            xlabel(axesHandle, 'SSBs', 'FontSize', 8);
            ylabel(axesHandle, 'Duration (h)', 'FontSize', 8);
            
            axesHandle = subplot(5, 2, 2, 'Parent', figHandle);
            axesHandles(1, 2) = axesHandle;
            tmp = [abs(ledPol1(:) - lagPol1(:)); abs(ledPol2(:) - lagPol2(:))] * 1e-3;
            hist(tmp)
            if any(tmp)
                xlim(axesHandle, [0.5 max(tmp)+0.5])
            end
            ylabel(axesHandle, 'Freq', 'FontSize', 8);
            xlabel(axesHandle, {'Lead-Lag Distance (knt)'}, 'FontSize', 8);
            
            axesHandle = subplot(5, 2, 4, 'Parent', figHandle);
            axesHandles(2, 2) = axesHandle;
            tmp = abs((c.sequenceLen - ledPol1(:)) - ledPol2(:)) * 1e-3;
            hist(tmp)
            if any(tmp)
                xlim(axesHandle, [0.5 max(tmp)+0.5])
            end
            ylabel(axesHandle, 'Freq', 'FontSize', 8);
            xlabel(axesHandle, {'Lead-Lead Distance (knt)'}, 'FontSize', 8);
            
            axesHandle = subplot(5, 2, 6, 'Parent', figHandle);
            axesHandles(3, 2) = axesHandle;
            hist(leadLeadCorr)
            if any(leadLeadCorr) && range(leadLeadCorr)
                xlim(axesHandle, [min(leadLeadCorr) max(leadLeadCorr)])
            end
            ylabel(axesHandle, 'Freq', 'FontSize', 8);
            xlabel(axesHandle, {'Lead-Lead Corr'}, 'FontSize', 8);
            
            axesHandle = subplot(5, 2, 8, 'Parent', figHandle);
            axesHandles(4, 2) = axesHandle;
            hist(okFragSynTimes)
            if ~isempty(okFragSynTimes) && range(okFragSynTimes)
                xlim(axesHandle, [min(okFragSynTimes) max(okFragSynTimes)])
            end
            ylabel(axesHandle, 'Freq', 'FontSize', 8);
            xlabel(axesHandle, {'Okazaki Synthesis Time (s)'}, 'FontSize', 8);
            
            axesHandle = subplot(5, 2, 10, 'Parent', figHandle);
            axesHandles(5, 2) = axesHandle;
            hist(ligationTimes)
            if ~isempty(ligationTimes) && range(ligationTimes)
                xlim(axesHandle, [min(ligationTimes) max(ligationTimes)])
            end
            ylabel(axesHandle, 'Freq', 'FontSize', 8);
            xlabel(axesHandle, {'Okazaki Synthesis Time (s)'}, 'FontSize', 8);
            
            %format axes
            for i = 1:numel(axesHandles)
                if ~axesHandles(i)
                    continue;
                end
                set(axesHandles(i), 'FontSize', 6);
            end
            
            %% organize figure data
            figData = struct;
            figData.ledPol1 = ledPol1;
            figData.ledPol2 = ledPol2;
            figData.lagPol1 = lagPol1;
            figData.lagPol2 = lagPol2;
            figData.ligationTimes = ligationTimes;
            figData.repStarts = repStarts;
            figData.repEnds = repEnds;
            figData.kineticallyLimitedEnds = kineticallyLimitedEnds;
            figData.repStartMass = repStartMass;
            figData.repStartGrowth = repStartGrowth;
            figData.repStartSSbs = repStartSSbs;
            figData.leadLeadCorr = leadLeadCorr;
            figData.okFragSynTimes = okFragSynTimes;
        end
        
        function figData = SSBs(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            %% constants
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            comp = sim.compartment;
            pc = sim.state('ProteinComplex');
            
            cytosolIndexs = comp.cytosolIndexs;
            cpxIdxs = [
                pc.matureIndexs(pc.getIndexs('MG_091_TETRAMER'))
                pc.matureIndexs(pc.getIndexs('MG_091_OCTAMER'))
                pc.boundIndexs(pc.getIndexs('MG_091_OCTAMER'))
                ];
            
            %% get data
            stateNames = {
                'Time'            'values'  ':'      ':'
                'Mass'            'cell'    ':'      '-sum'
                'ProteinComplex'  'counts'  cpxIdxs  cytosolIndexs
                };
            tmpStates = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations);
            
            mass = permute(tmpStates.Mass.cell, [3 4 1 2]);
            simEndTimes = permute(max(tmpStates.Time.values, [], 3), [4 3 1 2]);
            minSimEndTime = min(simEndTimes);
            [~, simOrder] = sort(mass(minSimEndTime, :));
            
            states.ProteinComplex.free4mer = permute(tmpStates.ProteinComplex.counts(1, :, :, :), [3 4 1 2]);
            states.ProteinComplex.free8mer = permute(tmpStates.ProteinComplex.counts(2, :, :, :), [3 4 1 2]);
            states.ProteinComplex.bound8mer = permute(tmpStates.ProteinComplex.counts(3, :, :, :), [3 4 1 2]);
            
            time = (1:size(states.ProteinComplex.bound8mer, 1))' / 3600;
            
            clear tmpStates mass;
            
            %% plot
            clf(figHandle);
            options = struct;
            options.colorOrder = PlotUtil.getRedGreenColorOrder(simOrder);
            options.xdata = time;
            options.ylabelStr = {
                {'Bound' '8mer'}
                {'Free' '8mer'}
                {'Free' '4mer'}
                };
            options.ydata = {
                states.ProteinComplex.bound8mer
                states.ProteinComplex.free8mer
                states.ProteinComplex.free4mer
                };
            PlotUtil.multiElementPlot(figHandle, [10; 10; 10], [0 time(end)], options);
            
            clear options;
            
            %% organize figure data
            figData = struct;
            figData.time = time;
            figData.free4mer = states.ProteinComplex.free4mer;
            figData.free8mer = states.ProteinComplex.free8mer;
            figData.bound8mer = states.ProteinComplex.bound8mer;
        end
        
        function figData = dnaRepair(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.util.CircularSparseMat;
            import edu.stanford.covert.util.ConstantUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            dnaRepair = sim.process('DNARepair');
            
            monIdxs = pm.matureIndexs(dnaRepair.enzymeMonomerGlobalIndexs);
            cpxIdxs = pc.matureIndexs(dnaRepair.enzymeComplexGlobalIndexs);
            
            %% get data
            stateNames = {
                'Time'            'values'                  ':'      ':'
                'Chromosome'      'abasicSites'             '-sum'   '-sum'
                'Chromosome'      'gapSites'                '-sum'   '-sum'
                'Chromosome'      'damagedSugarPhosphates'  '-sum'   '-sum'
                'Chromosome'      'intrastrandCrossLinks'   '-sum'   '-sum'
                'Chromosome'      'strandBreaks'            '-sum'   '-sum'
                'Chromosome'      'hollidayJunctions'       '-sum'   '-sum'
                'Mass'            'cell'                    ':'      '-sum'
                'ProteinMonomer'  'counts'                  monIdxs  '-sum'
                'ProteinComplex'  'counts'                  cpxIdxs  '-sum'
                };
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations);
            
            time = permute(max(states.Time.values, [], 4), [3 4 1 2]) / 3600;
            mass = permute(states.Mass.cell, [3 4 1 2]);
            
            simEndTimes = permute(max(states.Time.values, [], 3), [4 3 1 2]);
            minSimEndTime = min(simEndTimes);
            [~, simOrder] = sort(mass(minSimEndTime, :));
            
            clear mass;
            
            %damage
            gapSites = full(states.Chromosome.gapSites);
            abasicSites = full(states.Chromosome.abasicSites);
            damagedSugarPhosphates = full(states.Chromosome.damagedSugarPhosphates);
            intrastrandCrossLinks = full(states.Chromosome.intrastrandCrossLinks);
            strandBreaks = full(states.Chromosome.strandBreaks);
            hollidayJunctions = full(states.Chromosome.hollidayJunctions);
            
            %enzymes
            enzymes(dnaRepair.enzymeMonomerLocalIndexs, :, :, :) = states.ProteinMonomer.counts;
            enzymes(dnaRepair.enzymeComplexLocalIndexs, :, :, :) = states.ProteinComplex.counts;
            
            %cleanup
            clear states;
            
            %damaged bases -- calculate separately to minimize memory profile
            damagedBases = zeros(size(gapSites));
            for i = 1:numel(selectedSimulations)
                tmpStates = SimulationEnsemble.load(simBatchDir, {'Chromosome' 'damagedBases'}, [], [], 1, 'extract', selectedSimulations(i));
                [posStrnds, vals] = find(tmpStates.Chromosome.damagedBases);
                tfs = ~ismember(vals, dnaRepair.substrateGlobalIndexs(dnaRepair.substrateIndexs_m6AD));
                damagedBases(:, :, 1:size(tmpStates.Chromosome.damagedBases, 3), i) = ...
                    full(sum(sum(CircularSparseMat(posStrnds(tfs, :), ones(sum(tfs), 1), size(tmpStates.Chromosome.damagedBases), 1), 2), 1));
                clear tmpStates posStrnds vals tfs;
            end
            
            %NaN
            for i = 1:numel(simEndTimes)
                gapSites(:, :, simEndTimes(i)+1:end, i) = NaN;
                abasicSites(:, :, simEndTimes(i)+1:end, i) = NaN;
                damagedSugarPhosphates(:, :, simEndTimes(i)+1:end, i) = NaN;
                intrastrandCrossLinks(:, :, simEndTimes(i)+1:end, i) = NaN;
                strandBreaks(:, :, simEndTimes(i)+1:end, i) = NaN;
                hollidayJunctions(:, :, simEndTimes(i)+1:end, i) = NaN;
                enzymes(:, :, simEndTimes(i)+1:end, i) = NaN;
            end
            
            %% plot data
            clf(figHandle);
            
            options = struct();
            options.titleStr = {'Damage'; 'Repair'};
            options.colorOrder = PlotUtil.getRedGreenColorOrder(simOrder);
            options.xdata = time;
            options.ylabelStr = {
                {
                {'Gap' 'Sites'}
                {'Abasic' 'Sites'}
                {'Damaged' 'Sugar-Phosphates'}
                {'Damaged' 'Bases'}
                {'Intrastrand' 'Crosslinks'}
                {'Strand' 'Breaks'}
                {'Holliday' 'Junctions'}
                }
                dnaRepair.enzymeWholeCellModelIDs
                };
            options.ydata = {
                {
                permute(gapSites, [3 4 1 2])
                permute(abasicSites, [3 4 1 2])
                permute(damagedSugarPhosphates, [3 4 1 2])
                permute(damagedBases, [3 4 1 2])
                permute(intrastrandCrossLinks, [3 4 1 2])
                permute(strandBreaks, [3 4 1 2])
                permute(hollidayJunctions, [3 4 1 2])
                }
                permute(mat2cell(permute(enzymes, [3 4 1 2]), size(enzymes, 3), size(enzymes, 4), ones(size(enzymes, 1), 1)), [3 1 2])
                };
            PlotUtil.multiElementPlot(figHandle, cellfun(@(x) 3 * ones(size(x)), options.ylabelStr, 'UniformOutput', false), [0 time(end)], options);
            
            clear options;
            
            %% organize figure data
            figData = struct;
            figData.simOrder = simOrder;
            figData.time = time;
            figData.gapSites = gapSites;
            figData.abasicSites = abasicSites;
            figData.damagedSugarPhosphates = damagedSugarPhosphates;
            figData.damagedBases = damagedBases;
            figData.intrastrandCrossLinks = intrastrandCrossLinks;
            figData.strandBreaks = strandBreaks;
            figData.hollidayJunctions = hollidayJunctions;
            figData.enzymes = enzymes;
        end
        
        function figData = immuneActivation(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            comp = sim.compartment;
            pm = sim.state('ProteinMonomer');
            host = sim.state('Host');
            
            %% get data
            monIdxs = pm.matureIndexs(ismember(pm.compartments(pm.matureIndexs), [comp.terminalOrganelleCytosolIndexs comp.terminalOrganelleMembraneIndexs]));
            stateNames = {
                'Time'            'values'                           ':' ':'
                'Mass'            'cell'                             ':' '-sum'
                'Host'            'isBacteriumAdherent'              ':' ':'
                'Host'            'isTLRActivated'                   ':' ':'
                'Host'            'isNFkBActivated'                  ':' ':'
                'Host'            'isInflammatoryResponseActivated'  ':' ':'
                'ProteinMonomer'  'counts'                           monIdxs [comp.terminalOrganelleCytosolIndexs comp.terminalOrganelleMembraneIndexs]
                };
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations);
            
            time = permute(max(states.Time.values, [], 4), [3 4 1 2]) / 3600;
            mass = permute(states.Mass.cell, [3 4 1 2]);
            adherence = double(permute(states.Host.isBacteriumAdherent, [3 4 1 2]));
            tlr = double(permute(states.Host.isTLRActivated, [3 4 1 2]));
            nfkb = double(permute(states.Host.isNFkBActivated, [3 4 1 2]));
            inflammation = double(permute(states.Host.isInflammatoryResponseActivated, [3 4 1 2]));
            proteins = permute(states.ProteinMonomer.counts, [3 4 1 2]);
            
            simEndTimes = permute(max(states.Time.values, [], 3), [4 3 1 2]);
            minSimEndTime = min(simEndTimes);
            [~, simOrder] = sort(mass(minSimEndTime, :));
            
            for i = 1:numel(simEndTimes)
                adherence(simEndTimes(i)+1:end, i, :) = NaN;
                tlr(simEndTimes(i)+1:end, i, :) = NaN;
                nfkb(simEndTimes(i)+1:end, i, :) = NaN;
                inflammation(simEndTimes(i)+1:end, i, :) = NaN;
                proteins(simEndTimes(i)+1:end, i, :) = NaN;
            end
            
            clear states mass;
            
            %% plot data
            clf(figHandle);
            
            options = struct();
            options.titleStr = {'{\it M. genitalium}' 'Host Response'};
            options.colorOrder = PlotUtil.getRedGreenColorOrder(simOrder);
            options.xdata = time;
            options.ylabelStr = {
                cellfun(@(x) strrep(x, '_', '\_'), pm.wholeCellModelIDs(monIdxs), 'UniformOutput', false)
                {
                'Adherence'
                'TLR 1'
                'TLR 2'
                'TLR 6'
                'NF-{\kappa}B'
                'Inflammation'
                }
                };
            options.ydata = {
                mat2cell(proteins, size(proteins, 1), size(proteins, 2), ones(size(proteins, 3), 1))
                {
                adherence
                tlr(:, :, host.tlrIndexs_1)
                tlr(:, :, host.tlrIndexs_2)
                tlr(:, :, host.tlrIndexs_6)
                nfkb
                inflammation
                }};
            PlotUtil.multiElementPlot(figHandle, {3 * ones(size(monIdxs)) 3 * ones(6, 1)}, [0 time(end)], options);
            
            clear options;
            
            %% organize figure data
            figData = struct;
            figData.simOrder = simOrder;
            figData.time = time;
            figData.proteins = proteins;
            figData.adherence = adherence;
            figData.tlr = tlr;
            figData.nfkb = nfkb;
            figData.inflammation = inflammation;
        end
        
        function figData = smcSpacing(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            if nargin < 2 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 3 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            [simDir, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            simBatchDir = SimulationDiskUtil.getSimulationBatchDir(simDir);
            pc = sim.state('ProteinComplex');
            ch = sim.state('Chromosome');
            
            stateNames = {...
                'Mass'       'cell'                     ':'     '-sum'
                'Time'       'values'                   ':'     ':'
                };
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations);
            simEndTimes = permute(max(states.Time.values, [], 3), [4 1 2 3]);
            timeIdx = min(simEndTimes);
            [~, simOrder] = sort(states.Mass.cell(:, :, timeIdx, :));
            simOrder = simOrder(:);
            clear('states');
            
            maxTime = 0;
            for i = 1:length(selectedSimulations)
                summary = load([simBatchDir filesep num2str(selectedSimulations(i)) filesep 'summary.mat'], 'time');
                maxTime = max(maxTime, sum(summary.time >= 0));
            end
            
            time = 1:maxTime;
            smcSpacing = zeros(numel(selectedSimulations), numel(time));
            
            stateNames = {...
                'Chromosome' 'complexBoundSites'        ':'     ':'
                };
            
            for selSimIdx = 1:numel(selectedSimulations)
                % Run out of memory when using SimulationEnsemble
                % Instead, load and analyze each simulation one at a time
                states = DiskLogger.load([simBatchDir filesep num2str(selectedSimulations(selSimIdx))], stateNames, [], [], 1, 'extract');
                
                cbs = states.Chromosome.complexBoundSites;
                subs = find(cbs == pc.getIndexs('MG_213_214_298_6MER_ADP'));
                subs = subs(subs(:, 2) <= 2, :);
                
                smcSpacing(selSimIdx, :) = ch.sequenceLen ./  histc(subs(:, 3), time);
                
            end
            
            time = time / 3600;
            smcSpacing(smcSpacing == Inf) = NaN;
            
            %% plot data
            clf(figHandle);
            options = struct;
            options.titleStr = {'SMC Spacing'};
            options.colorOrder = PlotUtil.getRedGreenColorOrder(simOrder);
            options.xdata = time;
            options.ydata = {
                smcSpacing
                };
            options.ylabelStr = {
                {'Average Distance (bp)'}
                };
            PlotUtil.multiElementPlot(figHandle, {30}, [0 time(end)], options);
            
            clear options;
            
            %% organize figure data
            figData = struct;
            figData.simOrder = simOrder;
            figData.time = time;
            figData.smcSpacing = smcSpacing;
        end
        
        function figData = geneCopyNumber(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if ~exist('selectedSimulations', 'var') || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            [simDir, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            simBatchDir = SimulationDiskUtil.getSimulationBatchDir(simDir);
            ch = sim.state('Chromosome');
            
            %% get data
            stateNames = {...
                'Mass'       'cell'                     ':'     '-sum'
                'Time'       'values'                   ':'     ':'
                };
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations);
            simEndTimes = permute(max(states.Time.values, [], 3), [4 1 2 3]);
            timeIdx = min(simEndTimes);
            [~, simOrder] = sort(states.Mass.cell(:, :, timeIdx, :));
            simOrder = simOrder(:);
            clear('states');
            
            stateNames = {...
                'Chromosome' 'polymerizedRegions'        ':'     ':'
                };
            
            copyNumber = zeros(numel(selectedSimulations), ch.sequenceLen);
            
            for i = 1:numel(selectedSimulations)
                % Function is ~5x faster if I do this rather than
                % SimulationEnsemble.load() outside of the loop due to
                % avoiding SparseMat concatenation
                states = DiskLogger.load([simBatchDir filesep num2str(selectedSimulations(i))], stateNames, [], [], 1, 'extract');
                
                [subs vals] = find(states.Chromosome.polymerizedRegions);
                
                startpos = subs(:, 1);
                endpos = startpos + vals(:, 1) - 1;
                
                copyNumber(i, :) = cumsum(histc(startpos, 1:ch.sequenceLen) - histc(endpos + 1, 1:ch.sequenceLen)) / simEndTimes(i) / 2; % 2 strands per chromosome
            end
            
            %% plot data
            clf(figHandle);
            options = struct;
            options.titleStr = {'Copy Number Over Cell Cycle'};
            options.colorOrder = PlotUtil.getRedGreenColorOrder(simOrder);
            options.xdata = 1:ch.sequenceLen;
            options.xlabelStr = {'Position on Chromosome'};
            options.ydata = {
                copyNumber
                };
            options.ylabelStr = {
                {'Average Copy Number'}
                };
            [~, xAxesHandle] = PlotUtil.multiElementPlot(figHandle, {30}, [1 ch.sequenceLen], options);
            
            set(xAxesHandle{1}, 'XTick', [1 ch.terCPosition ch.sequenceLen], 'XTickLabel', {'oriC', 'terC', 'oriC'});
            
            clear options;
            
            %% organize figure data
            figData = struct;
            figData.simOrder = simOrder;
            figData.copyNumber = copyNumber;
        end
        
        function figData = unsynchronizedPopulation(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.analysis.Population;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if ~exist('selectedSimulations', 'var') || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            %% get constants
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            
            %% get data
            stateNames = {
                'growth_rate'       'Growth'
                'mass'              'Mass'
                'ploidy'            {'Chromosome' 'Copy Number'}
                'ftsZRing1st'       'ftsZRing'
                'ftsZRing2st'       'ftsZRing'
                'ftsZRing2bt'       'ftsZRing'
                'ftsZRingRbt'       'ftsZRing'
                'pinchedDiameter'   'pinchedDiameter'
                };
            ensemble = SimulationEnsemble(simBatchDir, stateNames, [], selectedSimulations);
            simEndTimes = ensemble.stateData.simulationEndTimes;
            
            [time, values, cellCount, simOrder] = Population.desynchronizePopulation(sim, ensemble.stateData.time, ensemble.stateData.values, simEndTimes);
            time = time / 3600;
            growth = permute(values(ensemble.getPropertyIndices('growth_rate'), :, :, :), [4 3 1 2]) * 3600 * sim.state('Mass').cellInitialDryWeight / (1 - sim.state('Mass').fractionWetWeight);
            mass = permute(values(ensemble.getPropertyIndices('mass'), :, :, :), [4 3 1 2]) * 1e15;
            ploidy = permute(values(ensemble.getPropertyIndices('ploidy'), :, :, :), [4 3 1 2]);
            ftsZRing1st = permute(values(ensemble.getPropertyIndices('ftsZRing1st'), :, :, :), [4 3 1 2]);
            ftsZRing2st = permute(values(ensemble.getPropertyIndices('ftsZRing2st'), :, :, :), [4 3 1 2]);
            ftsZRing2bt = permute(values(ensemble.getPropertyIndices('ftsZRing2bt'), :, :, :), [4 3 1 2]);
            ftsZRingRbt = permute(values(ensemble.getPropertyIndices('ftsZRingRbt'), :, :, :), [4 3 1 2]);
            pinchedDiameter = permute(values(ensemble.getPropertyIndices('pinchedDiameter'), :, :, :), [4 3 1 2]);
            cellCount = permute(cellCount, [4 3 1 2]);
            
            %% plot data
            clf(figHandle);
            options = struct;
            options.colorOrder = PlotUtil.getRedGreenColorOrder(simOrder);
            options.xdata = time;
            options.ydata = {
                growth
                mass
                ploidy
                ftsZRing1st
                ftsZRing2st
                ftsZRing2bt
                ftsZRingRbt
                pinchedDiameter
                cellCount
                };
            options.ylabelStr = {
                {'Growth' '(fg h^{-1})'}
                {'Mass' '(fg)'}
                {'Chr Copy' 'Number'}
                {'FtsZ' 'Ring' '1-Strt'}
                {'FtsZ' 'Ring' '2-Strt'}
                {'FtsZ' 'Ring' '2-Bent'}
                {'FtsZ' 'Ring' 'R-Bent'}
                {'Pinched' 'Diameter'}
                {'Cell' 'Count'}
                };
            PlotUtil.multiElementPlot(figHandle, 10 * ones(size(options.ydata)), [0 time(end)], options);
            
            clear options;
            
            %% organize figure data
            figData = struct;
            figData.simOrder = simOrder;
            figData.time = time;
            figData.growth = growth;
            figData.mass = mass;
            figData.ploidy = ploidy;
            figData.cellCount = cellCount;
        end
        
        function figData = blockedDecayEvents(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if ~exist('selectedSimulations', 'var') || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            %% get data
            nonDecays_DnaA = zeros(0, 1);
            nonDecays_FtsZ = zeros(0, 1);
            nonDecays_Replisome = zeros(0, 1);
            nonDecayCnts_DnaA = zeros(size(selectedSimulations));
            nonDecayCnts_FtsZ = zeros(size(selectedSimulations));
            nonDecayCnts_Replisome = zeros(size(selectedSimulations));
            for i = 1:numel(selectedSimulations)
                warnings = SimulationDiskUtil.getSimulationWarnings(simBatchDir, i);
                
                idxs = find(strcmp({warnings.message}, 'DnaA complex not decayed'));
                for j = 1:numel(idxs)
                    nonDecays_DnaA = [nonDecays_DnaA; warnings(idxs(j)).times]; %#ok<AGROW>
                    nonDecayCnts_DnaA(i) = nonDecayCnts_DnaA(i) + numel(warnings(idxs(j)).times);
                end
                
                idxs = find(strcmp({warnings.message}, 'FtsZ not decayed'));
                for j = 1:numel(idxs)
                    nonDecays_FtsZ = [nonDecays_FtsZ; warnings(idxs(j)).times]; %#ok<AGROW>
                    nonDecayCnts_FtsZ(i) = nonDecayCnts_FtsZ(i) + numel(warnings(idxs(j)).times);
                end
                
                idxs = find(strcmp({warnings.message}, 'DNA polymerase not decayed'));
                for j = 1:numel(idxs)
                    nonDecays_Replisome = [nonDecays_Replisome; warnings(idxs(j)).times]; %#ok<AGROW>
                    nonDecayCnts_Replisome(i) = nonDecayCnts_Replisome(i) + numel(warnings(idxs(j)).times);
                end
            end
            
            nonDecays_DnaA = nonDecays_DnaA / 3600;
            nonDecays_FtsZ = nonDecays_FtsZ / 3600;
            nonDecays_Replisome = nonDecays_Replisome / 3600;
            
            %% plot data
            clf(figHandle);
            options = struct();
            options.titleStr = {
                'Non Decay Event Timing'
                'Non Decay Events Per Simulation'
                };
            [axesHandles, xAxesHandles] = PlotUtil.multiElementPlot(figHandle, repmat({10 * ones(3, 1)}, 2, 1), [0 max([nonDecays_DnaA; nonDecays_FtsZ; nonDecays_Replisome])], options);
            delete(xAxesHandles{2});
            
            %Non Decay Event Timing
            edges = linspace(0, max([nonDecays_DnaA; nonDecays_FtsZ; nonDecays_Replisome]), 50);
            
            axesHandle = axesHandles{1}(1);
            hist(axesHandle, nonDecays_DnaA, edges);
            set(findobj(axesHandle, 'type', 'patch'), 'EdgeColor', 'b', 'FaceColor', 'b');
            ylabel(axesHandle, 'DnaA');
            
            axesHandle = axesHandles{1}(2);
            hist(axesHandle, nonDecays_FtsZ, edges);
            set(findobj(axesHandle, 'type', 'patch'), 'EdgeColor', 'b', 'FaceColor', 'b');
            ylabel(axesHandle, 'FtsZ');
            
            axesHandle = axesHandles{1}(3);
            hist(axesHandle, nonDecays_Replisome, edges);
            set(findobj(axesHandle, 'type', 'patch'), 'EdgeColor', 'b', 'FaceColor', 'b');
            ylabel(axesHandle, 'Replisome');
            
            %Non Decay Events Per Simulation
            axesHandle = axesHandles{2}(1);
            hist(axesHandle, nonDecayCnts_DnaA, 10);
            set(findobj(axesHandle, 'type', 'patch'), 'EdgeColor', 'b', 'FaceColor', 'b');
            xlabel(axesHandle, 'Non-Decays');
            ylabel(axesHandle, 'DnaA');
            
            axesHandle = axesHandles{2}(2);
            hist(axesHandle, nonDecayCnts_DnaA, 10);
            set(findobj(axesHandle, 'type', 'patch'), 'EdgeColor', 'b', 'FaceColor', 'b');
            xlabel(axesHandle, 'Non-Decays');
            ylabel(axesHandle, 'FtsZ');
            
            axesHandle = axesHandles{2}(3);
            hist(axesHandle, nonDecayCnts_DnaA, 10);
            set(findobj(axesHandle, 'type', 'patch'), 'EdgeColor', 'b', 'FaceColor', 'b');
            xlabel(axesHandle, 'Non-Decays');
            ylabel(axesHandle, 'Replisome');
            
            PlotUtil.alignYAxesLabels(axesHandles);
            PlotUtil.offsetYAxes(axesHandles, 0.04);
            
            %% data
            figData = struct;
            figData.nonDecays_DnaA = nonDecays_DnaA;
            figData.nonDecays_FtsZ = nonDecays_FtsZ;
            figData.nonDecays_Replisome = nonDecays_Replisome;
            figData.nonDecayCnts_DnaA = nonDecayCnts_DnaA;
            figData.nonDecayCnts_FtsZ = nonDecayCnts_FtsZ;
            figData.nonDecayCnts_Replisome = nonDecayCnts_Replisome;
        end
        
        function figData = secondaryReplicationInitiation(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SummaryLogger;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if ~exist('selectedSimulations', 'var') || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            %% get data
            secondaryInitiationTimes = [];
            secondaryInitiationInitTimes = NaN(numel(selectedSimulations), 1);
            secondaryInitiationCnts = zeros(numel(selectedSimulations), 1);
            for i = 1:numel(selectedSimulations)
                warnings = SimulationDiskUtil.getSimulationWarnings(simBatchDir, i);
                
                idxs = find(strcmp({warnings.message}, 'Unable to represent multiple replication initiation events per cell cycle'));
                for j = 1:numel(idxs)
                    secondaryInitiationTimes = [secondaryInitiationTimes; warnings(idxs(j)).times]; %#ok<AGROW>
                    secondaryInitiationInitTimes(i) = min(secondaryInitiationInitTimes(i), warnings(idxs).times(1));
                    secondaryInitiationCnts(i) = secondaryInitiationCnts(i) + numel(warnings(idxs(j)).times);
                end
            end
            secondaryInitiationTimes = secondaryInitiationTimes / 3600;
            secondaryInitiationInitTimes = secondaryInitiationInitTimes / 3600;
            
            simStats = SummaryLogger.getSimulationStatistics([SimulationDiskUtil.getBaseDir() filesep simBatchDir], selectedSimulations);
            repInitTimes = simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME) / 3600;
            cellCycleTimes = simStats(:, SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME) / 3600;
            
            clear simStats;
            
            stateNames = {
                'dnaA_box1'
                'dnaA_box2'
                'dnaA_box3'
                'dnaA_box4'
                'dnaA_box5'
                };
            ensemble = SimulationEnsemble(simBatchDir, stateNames, [], selectedSimulations);
            tmp = permute(sum(ensemble.stateData.values, 1), [4 3 1 2]);
            dnaAComplexSizes = tmp(sub2ind(size(tmp), (1:numel(selectedSimulations))', ...
                floor(ensemble.stateData.simulationEndTimes / ensemble.stateData.downsampleStepSec)));
            
            clear ensemble stateNames;
            
            %% plot data
            clf(figHandle);
            options = struct();
            options.xlabelStr = {'Secondary Replication Initiation Time (h)', 'Frequency'};
            [axesHandles, xAxesHandles] = PlotUtil.multiElementPlot(figHandle, {[25 25 25] [1 1]}, [0 max(cellCycleTimes)], options);
            
            axesHandle = axesHandles{1}(1);
            edges = reshape(linspace(0, max(cellCycleTimes), 50), [], 1);
            cla(axesHandle);
            hist(axesHandle, secondaryInitiationTimes, edges);
            set(findobj(axesHandle, 'type', 'patch'), 'EdgeColor', 'b', 'FaceColor', 'b');
            set(axesHandle, 'YTick', [0 max(get(axesHandle, 'YTick'))]);
            ylabel(axesHandle, {'Attempt Frequency'});
            
            axesHandle = axesHandles{1}(2);
            cla(axesHandle)
            edges = reshape(linspace(0, max(cellCycleTimes), 10), [], 1);
            tmp = [histc(repInitTimes, edges)  histc(secondaryInitiationInitTimes, edges)];
            h = bar(axesHandle, edges, tmp);
            legend(h, {'1st', '2nd'}, 'Location', 'NorthWest');
            set(h, 'EdgeColor', 'none');
            set(h(1), 'FaceColor', 'r');
            set(h(2), 'FaceColor', 'b');
            xlim(axesHandle, [0 max(cellCycleTimes)]);
            set(axesHandle, 'YTick', [0 max(get(axesHandle, 'YTick'))]);
            ylabel(axesHandle, {'Frequency'});
            
            axesHandle = axesHandles{1}(3);
            hold(axesHandle, 'on');
            plot(axesHandle, secondaryInitiationInitTimes, repInitTimes, '.');
            line([0 max(cellCycleTimes)], [0 max(cellCycleTimes)], ...
                'Color', 'r', 'LineStyle', ':', 'Parent', axesHandle);
            xlim(axesHandle, [0 max(cellCycleTimes)]);
            ylim(axesHandle, [0 max(cellCycleTimes)]);
            set(axesHandle, 'YTick', get(xAxesHandles{1}, 'XTick'));
            ylabel(axesHandle, {'Replication Initiation Time'});
            
            axesHandle = axesHandles{2}(1);
            cla(axesHandle);
            edges = linspace(min(secondaryInitiationCnts), max(max(secondaryInitiationCnts), 9), 10);
            freq1 = histc(secondaryInitiationCnts, edges);
            barh(axesHandle, edges, freq1, 'EdgeColor', 'b', 'FaceColor', 'b');
            ylabel(axesHandle, 'Secondary Initiation Attempts');
            ylim(axesHandle, edges([1 end]) + [-0.5 0.5]);
            set(axesHandle, 'YTick', [0 max(get(axesHandle, 'YTick'))]);
            
            axesHandle = axesHandles{2}(2);
            cla(axesHandle);
            edges = (0:29)';
            freq2 = histc(dnaAComplexSizes, edges);
            barh(axesHandle, edges, freq2, 'EdgeColor', 'b', 'FaceColor', 'b');
            ylabel(axesHandle, 'DnaA Complex Size');
            ylim(axesHandle, edges([1 end]) + [-0.5; 0.5]);
            set(axesHandle, 'YTick', 0:4:28);
            set(axesHandle, 'YTickLabel', 0:7);
            
            xlim(xAxesHandles{2}, [min([freq1; freq2]) max([freq1; freq2])]);
            xlim(axesHandles{2}(1), [min([freq1; freq2]) max([freq1; freq2])]);
            xlim(axesHandles{2}(2), [min([freq1; freq2]) max([freq1; freq2])]);
            set(xAxesHandles{2}, 'XTickMode', 'auto');
            
            PlotUtil.offsetYAxes(axesHandles, 0.03);
            
            %% data
            figData = struct;
            figData.repInitTimes = repInitTimes;
            figData.secondaryInitiationTimes = secondaryInitiationTimes;
            figData.secondaryInitiationInitTimes = secondaryInitiationInitTimes;
            figData.secondaryInitiationCnts = secondaryInitiationCnts;
            figData.dnaAComplexSizes = dnaAComplexSizes;
        end
        
        function [figData, tabContent, tabColLabels] = warnings(figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SummaryLogger;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            if ~exist('simBatchDir', 'var') || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if ~exist('selectedSimulations', 'var') || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            %% get data
            warnings = [];
            freq = zeros(10, 1);
            for i = 1:numel(selectedSimulations)
                tmpWarnings = SimulationDiskUtil.getSimulationWarnings(simBatchDir, selectedSimulations(i));
                
                repeatedWarning = 0;
                for k = 1:numel(tmpWarnings)
                    if strcmp(tmpWarnings(k).message, 'Unable to represent multiple replication initiation events per cell cycle')
                        freq(1) = freq(1) + numel(tmpWarnings(k).times);
                    elseif strcmp(tmpWarnings(k).message, 'DnaA complex not decayed')
                        freq(2) = freq(2) + numel(tmpWarnings(k).times);
                    elseif strcmp(tmpWarnings(k).message, 'FtsZ not decayed')
                        freq(3) = freq(3) + numel(tmpWarnings(k).times);
                    elseif strcmp(tmpWarnings(k).message, 'DNA polymerase not decayed')
                        freq(4) = freq(4) + numel(tmpWarnings(k).times);
                    elseif strcmp(tmpWarnings(k).message, 'Cell has divided. Cell shape undefined')
                        freq(5) = freq(5) + numel(tmpWarnings(k).times);
                    elseif length(tmpWarnings(k).message) > 24 && strcmp(tmpWarnings(k).message(1:24), 'Linear programming error')
                        freq(6) = freq(6) + numel(tmpWarnings(k).times);
                    elseif length(tmpWarnings(k).message) > 50 && strcmp(tmpWarnings(k).message(1:50), 'Simulation state side effects should preserve mass')
                        freq(7) = freq(7) + numel(tmpWarnings(k).times);
                    elseif length(tmpWarnings(k).message) > 35 && strcmp(tmpWarnings(k).message(1:35), 'Summary figure could not be saved: ')
                        freq(8) = freq(8) + numel(tmpWarnings(k).times);
                    elseif regexp(tmpWarnings(k).message, 'metabolites have negative counts')
                        freq(9) = freq(9) + numel(tmpWarnings(k).times);
                    else
                        freq(10) = freq(10) + numel(tmpWarnings(k).times);
                    end
                    
                    tmpWarnings(k).simulations = repmat(i, size(tmpWarnings(k).times));
                    
                    for j = 1:numel(warnings)
                        if isequal(warnings(j).message, tmpWarnings(k).message) && ...
                                isequal(warnings(j).files, tmpWarnings(k).files) && ...
                                isequal(warnings(j).lineNumbers, tmpWarnings(k).lineNumbers)
                            repeatedWarning = j;
                            break;
                        end
                    end
                    
                    if isempty(warnings)
                        warnings = tmpWarnings(k);
                    elseif repeatedWarning
                        warnings(repeatedWarning).times = ...
                            [warnings(repeatedWarning).times; tmpWarnings(k).times]; %#ok<AGROW>
                        warnings(repeatedWarning).simulations = ...
                            [warnings(repeatedWarning).simulations; tmpWarnings(k).simulations]; %#ok<AGROW>
                    else
                        warnings = [warnings; tmpWarnings(k)]; %#ok<AGROW>
                    end
                end
            end
            
            %% plot data
            clf(figHandle);
            axesHandle = subplot(1, 1, 1, 'Parent', figHandle);
            barh(axesHandle, freq)
            set(axesHandle, 'YTick', 1:numel(freq),  'YTickLabel', {'2nd Rep Init', 'DnaA Decay', 'FtsZ Decay', 'Replisome Decay', 'Shape Undef', 'Lin Prog', 'Neg Mets', 'Sum Fig Save', 'Mass Bal', 'Other'});
            ylabel(axesHandle, 'Warning');
            xlabel(axesHandle, 'Frequency');
            set(axesHandle, 'XScale', 'log');
            ylim(axesHandle, [0.5 numel(freq) + 0.5]);
            
            %% table
            warnings = warnings(~ismember({warnings.message}, {
                'Unable to represent multiple replication initiation events per cell cycle'
                'DnaA complex not decayed'
                'FtsZ not decayed'
                'DNA polymerase not decayed'
                'Cell has divided. Cell shape undefined'
                }));
            tfs = true(size(warnings));
            for i = 1:numel(warnings)
                tfs(i) = length(warnings(i).message) < 35 || ~isequal(warnings(i).message(1:35), 'Summary figure could not be saved: ');
            end
            warnings = warnings(tfs);
            tabColLabels = {'Warning', 'Stack Trace', 'Simulation.Times'};
            tabContent = cell(numel(warnings), 3);
            for i = 1:numel(warnings)
                tabContent{i, 1} = warnings(i).message;
                
                for j = 1:numel(warnings(i).files)
                    tabContent{i, 2} = [tabContent{i, 2} ...
                        sprintf('%s at %d\n', warnings(i).files{j}, warnings(i).lineNumbers(j))];
                end
                tabContent{i, 2} = tabContent{i, 2}(1:end-1);
                
                for j = 1:numel(warnings(i).simulations)
                    tabContent{i, 3} = [tabContent{i, 3} ...
                        sprintf('%d.%d\n', warnings(i).simulations(j), warnings(i).times(j))];
                end
                tabContent{i, 3} = tabContent{i, 3}(1:end-1);
            end
            
            %% data
            figData = struct;
            figData.freq = freq;
        end
    end
    
    %helper methods
    methods (Static = true)
        function [unsyncTime, unsyncData, cellCount, simOrder] = desynchronizePopulation(sim, syncTime, syncData, simEndTimes)
            finTimeIdx = min(simEndTimes);
            finTime = syncTime(finTimeIdx);
            
            cellCycleLength = sim.state('Time').cellCycleLength;
            randStream = edu.stanford.covert.util.RandStream('mcg16807');
            
            unsyncData = NaN(size(syncData, 1), size(syncData, 2), min(simEndTimes), 2 * size(syncData, 4));
            simStartTimes = zeros(size(syncData, 4), 1);
            for i = 1:ceil(size(syncData, 4) * (1 - exp(-log(2) * finTimeIdx / cellCycleLength)))
                simTimeIdx = ceil(randStream.rand * finTimeIdx);
                unsyncData(:, :, 1:simTimeIdx, i) = syncData(:, :, simEndTimes(i)-simTimeIdx+1:simEndTimes(i), i);
                unsyncData(:, :, simTimeIdx+1:end, i) = syncData(:, :, 1:finTimeIdx - simTimeIdx, i);
            end
            for j = i+1:size(syncData, 4)
                simStartTimes(j) = randStream.stochasticRound(log(1 + randStream.rand * (exp(log(2) * finTime / cellCycleLength) - 1)) * cellCycleLength / log(2));
                simTimeIdx = find(syncTime == simStartTimes(j));
                if randStream.rand > 0.5
                    unsyncData(:, :, simTimeIdx+1:end, j) = syncData(:, :, 1:finTimeIdx - simTimeIdx, j);
                else
                    unsyncData(:, :, simTimeIdx+1:end, j) = syncData(:, :, simEndTimes(j) - (finTimeIdx - simTimeIdx)+1:simEndTimes(j), j);
                end
            end
            [~, simOrder] = sort(simStartTimes);
            
            cellCount = sum(~isnan(unsyncData(1, :, :, :)), 4);
            
            unsyncTime = syncTime(:, 1:finTimeIdx);
        end
    end
end
