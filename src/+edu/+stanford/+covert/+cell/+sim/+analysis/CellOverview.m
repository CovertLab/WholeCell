% CellOverview (Executive summary of a single cell)
%
% Author: Derek Macklin, macklin@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 7/6/2011
classdef CellOverview
    methods (Static = true);
        function run(simBatchDir, fileName, simIdxs)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.analysis.CellOverview;
            
            %% options
            if nargin < 1 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            
            if nargin < 3
                simIdxs = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            %% get data
            states = {
                'runTime'           'Run Time'
                'mass'              'Mass'
                'growth_rate'       {'Growth' 'Rate'}
                'ploidy'            'Ploidy'
                'superhelicity'     'Superhelicity'
                'polysaccharides'   'Polysaccharides'
                'immatureMonomers'  {'Immature' 'Monomers'}
                'matureMonomers'    {'Mature' 'Monomers'}
                'immatureComplexs'  {'Immature' 'Complexes'}
                'matureComplexs'    {'Mature' 'Complexes'}
                'proteins'          'Proteins'
                'rnas'              'RNA'
                'dnaA_boxes'        {'Bound' 'DnaA'}
                'dnaA_box1'         {'DnaA' 'Box 1'}
                'dnaA_box2'         {'DnaA' 'Box 2'}
                'dnaA_box3'         {'DnaA' 'Box 3'}
                'dnaA_box4'         {'DnaA' 'Box 4'}
                'dnaA_box5'         {'DnaA' 'Box 5'}
                'helicase1'         'Helicase 1'
                'helicase2'         'Helicase 2'
                'atp'               'ATP'
                'adp'               'ADP'
                'amp'               'AMP'
                'ntps'              'NTPs'
                'amino_acids'       {'Amino' 'Acids'}
                'rnaPolymerases'    {'RNA' 'Pols'}
                'ribosomes'         'Ribosomes'
                'mrnas'             'mRNA'
                'rrnas'             'rRNA'
                'srnas'             'sRNA'
                'trnas'             'tRNA'
                'immatureRnas'      {'Immature' 'RNA'}
                'ftsZ'              'FtsZ'
                'ftsZRing1st'       '1-Strt'
                'ftsZRing2st'       '2-Strt'
                'ftsZRing2bt'       '2-Bent'
                'ftsZRingRbt'       'R-Bent'
                'pinchedDiameter'   {'Pinched' 'Diameter'}
                };
            s = SimulationEnsemble(simBatchDir, states, [], simIdxs);
            
            [~, simOrder] = sort(s.stateData.values(s.getPropertyIndices('growth_rate'), 1, ....
                floor(min(s.stateData.simulationEndTimes) / s.stateData.downsampleStepSec), :));
            
            %% plot, and save
            figHandles = [];
            
            [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            figHandles = [figHandles; figHandle];
            CellOverview.plotMultipleSimulations(s, figHandle, @PlotUtil.plotLine, {'mass'; 'growth_rate'; 'ploidy'; 'polysaccharides'; 'proteins'; 'rnas'}, ...
                struct('simOrder', simOrder));
            if nargin >= 2
                saveas(figHandle, [fileName '-CellOverview-Lines1.pdf']);
                close(figHandle);
            end
            
            [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            figHandles = [figHandles; figHandle];
            CellOverview.plotMultipleSimulations(s, figHandle, @PlotUtil.plotLine, {'growth_rate'; 'dnaA_boxes'; 'dnaA_box1'; 'dnaA_box2'; 'dnaA_box3'; 'dnaA_box4'; 'dnaA_box5'}, ...
                struct('simOrder', simOrder));
            if nargin >= 2
                saveas(figHandle, [fileName '-CellOverview-Lines2.pdf']);
                close(figHandle);
            end
            
            [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            figHandles = [figHandles; figHandle];
            CellOverview.plotMultipleSimulations(s, figHandle, @PlotUtil.plotLine, {'growth_rate'; 'helicase1'; 'helicase2'}, ...
                struct('simOrder', simOrder));
            if nargin >= 2
                saveas(figHandle, [fileName '-CellOverview-Lines3.pdf']);
                close(figHandle);
            end
            
            [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            figHandles = [figHandles; figHandle];
            CellOverview.plotMultipleSimulations(s, figHandle, @PlotUtil.plotLine, {'atp'; 'adp'; 'amp'; 'ntps'; 'amino_acids'; 'rnaPolymerases'; 'ribosomes'}, ...
                struct('simOrder', simOrder));
            if nargin >= 2
                saveas(figHandle, [fileName '-CellOverview-Lines4.pdf']);
                close(figHandle);
            end
            
            [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            figHandles = [figHandles; figHandle];
            CellOverview.plotMultipleSimulations(s, figHandle, @PlotUtil.plotLine, {'rnas'; 'mrnas'; 'rrnas'; 'srnas'; 'trnas'; 'immatureRnas'}, ...
                struct('simOrder', simOrder));
            if nargin >= 2
                saveas(figHandle, [fileName '-CellOverview-Lines5.pdf']);
                close(figHandle);
            end
            
            [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            figHandles = [figHandles; figHandle];
            CellOverview.plotMultipleSimulations(s, figHandle, @PlotUtil.plotLine, {'immatureMonomers'; 'matureMonomers'; 'immatureComplexs'; 'matureComplexs'}, ...
                struct('simOrder', simOrder));
            if nargin >= 2
                saveas(figHandle, [fileName '-CellOverview-Lines6.pdf']);
                close(figHandle);
            end
            
            [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            figHandles = [figHandles; figHandle];
            CellOverview.plotMultipleSimulations(s, figHandle, @PlotUtil.plotLine, {'ftsZ'; 'ftsZRing1st'; 'ftsZRing2st'; 'ftsZRing2bt'; 'ftsZRingRbt'; 'pinchedDiameter'}, ...
                struct('simOrder', simOrder));
            if nargin >= 2
                saveas(figHandle, [fileName '-CellOverview-Lines7.pdf']);
                close(figHandle);
            end
            
            [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            figHandles = [figHandles; figHandle];
            CellOverview.plotMultipleSimulations(s, figHandle, @PlotUtil.plotLine, {'ploidy'; 'superhelicity'}, ...
                struct('simOrder', simOrder));
            if nargin >= 2
                saveas(figHandle, [fileName '-CellOverview-Lines8.pdf']);
                close(figHandle);
            end
            
            [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            figHandles = [figHandles; figHandle];
            CellOverview.plotMultipleSimulations(s, figHandle, @PlotUtil.plotLine, {'runTime'}, ...
                struct('simOrder', simOrder));
            if nargin >= 2
                saveas(figHandle, [fileName '-CellOverview-Lines9.pdf']);
                close(figHandle);
            end
            
            [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            figHandles = [figHandles; figHandle];
            CellOverview.plotInitialGrowthVsEndTime(s, figHandle);
            if nargin >= 2
                saveas(figHandle, [fileName '-CellOverview-InitialGrowthRateVsEndTimes.pdf']);
                close(figHandle);
            end
            
            [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            figHandles = [figHandles; figHandle];
            CellOverview.plotFinalGrowthVsEndTime(s, figHandle);
            if nargin >= 2
                saveas(figHandle, [fileName '-CellOverview-FinalGrowthRateVsEndTimes.pdf']);
                close(figHandle);
            end
            
            if nargin < 2
                figHandles = figHandles(1:min(end, 9));
                for i = 1:numel(figHandles)
                    setappdata(figHandles(i), 'figHandles', figHandles);
                end
                if numel(figHandles) == 9
                    spanFigures(fliplr(reshape(figHandles, [3 3]))');
                end
            end
        end
    end
    
    %multiple simulations
    methods (Static = true)
        function axesHandles = plotMultipleSimulations(simEnsemble, figHandle, doPlot, statesToPlot, options)
            import edu.stanford.covert.cell.sim.analysis.CellOverview;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            
            %% options
            if nargin < 3 || isempty(doPlot)
                doPlot = @PlotUtil.plotLine;
            end
            
            if nargin < 4
                statesToPlot = simEnsemble.stateData.properties;
            end
            
            if nargin < 5
                options = struct();
            end
            if ~isfield(options, 'kaplanMeier')
                options.kaplanMeier = true;
            end
            if ~isfield(options, 'simulationTicks')
                options.simulationTicks = true;
            end
            if ~isfield(options, 'markEarliestSimulationTermination')
                options.markEarliestSimulationTermination = true;
            end
            if ~isfield(options, 'simOrder')
                [~, options.simOrder] = sort(simEnsemble.stateData.values(simEnsemble.getPropertyIndices(statesToPlot(1)), 1, min(simEnsemble.stateData.simulationEndTimes), :));
            end
            
            %% get data
            time = simEnsemble.stateData.time / 3600;
            inds = simEnsemble.getPropertyIndices(statesToPlot);
            values = simEnsemble.stateData.values(inds, :, :, :);
            
            %% layout plots
            nSims = size(values, 4);
            nValuePlots = length(inds);
            nTotalPlots = nValuePlots + 2;
            axesSizes = 4 * ones(nValuePlots, 1);
            dataHandles = zeros(nSims, nTotalPlots);
            
            if options.kaplanMeier
                kaplanMeier = SimulationEnsemble.remainingSimulations(simEnsemble.stateData.time, simEnsemble.stateData.simulationEndTimes);
                axesSizes = [axesSizes; 2];
            end
            
            if options.simulationTicks
                axesSizes = [axesSizes; 0.5];
            end
            
            clf(figHandle);            
            options.titleStr = sprintf('Overview: %s', simEnsemble.getTimeStamp());
            options.colorOrder = PlotUtil.getRedGreenColorOrder(options.simOrder);
            axesHandles = PlotUtil.multiElementPlot(figHandle, axesSizes, time([1 end])', ...
                options);
            
            %% value plots            
            for j = 1:nValuePlots
                % subplot
                axesHandle = axesHandles(j);
                
                % plot data
                tmp = doPlot(axesHandle, time, permute(values(j, :, :, :), [4 3 1 2]));
                if numel(tmp) == nSims
                    dataHandles(:, j) = tmp;
                end
                
                % Y-Axis
                ylabel(axesHandle, simEnsemble.stateData.descriptions{inds(j)});
                
                % Indicate where first simulation terminated
                if options.markEarliestSimulationTermination
                    CellOverview.markEarliestSimulationTermination(simEnsemble, axesHandle);
                end
            end
            
            %% Kaplan-Meier Curve
            if options.kaplanMeier
                axesHandle = axesHandles(nValuePlots + 1);
                
                h = PlotUtil.plotLine(axesHandle, time, 100 * kaplanMeier);
                set(h, 'Color', 'k');
                text(0.1, 0.1, 'Remaining Simulations', ...
                    'HorizontalAlignment', 'Left', ...
                    'VerticalAlignment', 'Bottom', ...
                    'Parent', axesHandle, 'Units', 'data', ...
                    'FontSize', 8);
                set(axesHandle, 'YTick', [0 100], 'YtickLabel', {'0%', '100%'});
                
                simIdxs = find(~isnan(simEnsemble.stateData.simulationEndTimes));
                dataHandles(:, end-1) = scatter(axesHandle, ...
                    simEnsemble.stateData.simulationEndTimes(simIdxs) / 3600, ...
                    100 * kaplanMeier(floor(simEnsemble.stateData.simulationEndTimes(simIdxs) / simEnsemble.stateData.downsampleStepSec)), ...
                    'LineWidth', 2, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'Visible', 'off');
                
                % Indicate where first simulation terminated
                if options.markEarliestSimulationTermination
                    CellOverview.markEarliestSimulationTermination(simEnsemble, axesHandle);
                end
            end
            
            %% Simulation completion tick marks
            if options.simulationTicks
                axesHandle = axesHandles(nValuePlots + 2);
                
                dataHandles(:, end) = PlotUtil.plotColoredTicks(axesHandle, ...
                    simEnsemble.stateData.simulationEndTimes / 3600, options.colorOrder);
                
                set(axesHandle, 'YColor', get(figHandle, 'Color'));
                text(0.1, 0.1, 'Simulation Completions', ...
                    'HorizontalAlignment', 'Left', ...
                    'VerticalAlignment', 'Bottom', ...
                    'Parent', axesHandle, 'Units', 'data', ...
                    'FontSize', 8);
                
                % Indicate where first simulation terminated
                if options.markEarliestSimulationTermination
                    CellOverview.markEarliestSimulationTermination(simEnsemble, axesHandle);
                end
            end
            
            %% Format Y-axis
            %-Align Y-axis labels
            %-offset
            PlotUtil.alignYAxesLabels(axesHandles);
            axesHandles2 = PlotUtil.offsetYAxes(axesHandles, 0.015);
            
            %% callbacks
            setappdata(figHandle, 'nSimulations', nSims);
            setappdata(figHandle, 'nTotalPlots', nTotalPlots);
            setappdata(figHandle, 'axesHandles', axesHandles);
            setappdata(figHandle, 'axesHandles2', axesHandles2);
            setappdata(figHandle, 'dataHandles', dataHandles);
            setappdata(figHandle, 'colorOrder', options.colorOrder);
            for i = 1:nSims
                for j = 1:nTotalPlots
                    if dataHandles(i, j) ~= 0
                        setappdata(dataHandles(i, j), 'simulationIdx', i);
                        setappdata(dataHandles(i, j)', 'propertyIdx', j);
                        set(dataHandles(i, j), 'ButtonDownFcn', @CellOverview.plotMultipleSimulations_Callback);
                    end
                end
            end
        end
        
        function plotMultipleSimulations_Callback(hObject, ~)
            import edu.stanford.covert.cell.sim.analysis.CellOverview;
            
            %get selected figure, simulation, property
            figHandle = get(get(hObject, 'parent'), 'parent');
            simulationIdx = getappdata(hObject, 'simulationIdx');
            propertyIdx = getappdata(hObject, 'propertyIdx');
            
            %print selected simulation
            fprintf('Simulation #%d selected\n', simulationIdx);
            
            %highlight selected simulation
            CellOverview.unselectSimulations(figHandle);
            CellOverview.selectSimulation(simulationIdx, propertyIdx, figHandle)
        end
        
        function unselectSimulations(figHandle)
            %get property of simulations
            nSims = getappdata(figHandle, 'nSimulations');
            colorOrder = getappdata(figHandle, 'colorOrder') * 0.20 + 0.80 * 1;
            lineWidth = ones(nSims, 1);
            
            %get handles of related figures
            figHandles = getappdata(figHandle, 'figHandles');
            if isempty(figHandles)
                figHandles = figHandle;
            end
            
            %highlight selected simulation, dull others
            for k = 1:numel(figHandles)
                try
                    dataHandles = getappdata(figHandles(k), 'dataHandles');
                catch %#ok<CTCH>
                    continue;
                end
                dataHandles = dataHandles(:, any(dataHandles, 1));
                nTotalPlots = size(dataHandles, 2);
                
                set(dataHandles, 'Selected', 'off');
                set(dataHandles(:, end-1), 'visible', 'off');
                
                for j = [1:nTotalPlots-2 nTotalPlots]
                    for i = 1:nSims
                        set(dataHandles(i, j), 'Color', colorOrder(i, :), 'LineWidth', lineWidth(i));
                    end
                end
            end
        end
        
        function selectSimulation(simulationIdx, propertyIdx, figHandle)
            %get property of simulations
            colorOrder = getappdata(figHandle, 'colorOrder');
            
            %get handles of related figures
            figHandles = getappdata(figHandle, 'figHandles');
            if isempty(figHandles)
                figHandles = figHandle;
            end
            
            %highlight selected simulation, dull others
            for k = 1:numel(figHandles)
                try
                    dataHandles = getappdata(figHandles(k), 'dataHandles');
                catch %#ok<CTCH>
                    continue;
                end
                dataHandles = dataHandles(:, any(dataHandles, 1));
                nTotalPlots = size(dataHandles, 2);
                
                set(dataHandles(simulationIdx, end-1), 'visible', 'on');
                
                for j = [1:nTotalPlots-2 nTotalPlots]
                    set(dataHandles(simulationIdx, j), 'Color', colorOrder(simulationIdx, :), 'LineWidth', 2);
                    uistack(dataHandles(simulationIdx, j), 'top');
                end
            end
            
            %briefly set selected plot selected='on'
            dataHandles = getappdata(figHandle, 'dataHandles');
            set(dataHandles(simulationIdx, propertyIdx), 'Selected', 'on');
            drawnow;
            pause(0.1);
            set(dataHandles(simulationIdx, propertyIdx), 'Selected', 'off');
            drawnow;
        end
        
        function restoreSimulations()
            %get property of simulations
            nSims = getappdata(figHandle, 'nSimulations');
            colorOrder = getappdata(figHandle, 'colorOrder');
            lineWidth = ones(nSims, 1);
            
            %get handles of related figures
            figHandles = getappdata(figHandle, 'figHandles');
            if isempty(figHandles)
                figHandles = figHandle;
            end
            
            %highlight selected simulation, dull others
            for k = 1:numel(figHandles)
                try
                    dataHandles = getappdata(figHandles(k), 'dataHandles');
                catch %#ok<CTCH>
                    continue;
                end
                dataHandles = dataHandles(:, any(dataHandles, 1));
                nTotalPlots = size(dataHandles, 2);
                
                set(dataHandles, 'Selected', 'off');
                set(dataHandles(:, end-1), 'visible', 'off');
                
                for j = [1:nTotalPlots-2 nTotalPlots]
                    for i = 1:nSims
                        set(dataHandles(i, j), 'Color', colorOrder(i, :), 'LineWidth', lineWidth(i));
                    end
                end
            end
        end
        
        % Indicate where first simulation terminated
        function markEarliestSimulationTermination(simEnsemble, axesHandle)            
            valid = simEnsemble.getEarliestTermination() / 3600;
            yr = ylim(axesHandle);
            line([valid valid], yr, ...
                'Parent', axesHandle, ...
                'LineStyle', '--',...
                'Color', 'c', ...
                'LineWidth', 1);
        end
        
        function plotInitialGrowthVsEndTime(simEnsemble, figHandle)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            ind = simEnsemble.getPropertyIndices('growth_rate');
            if isnan(ind)
                return
            end
            clf(figHandle)
            
            axesHandle = subplot(1, 1, 1);
            
            growth_rates = permute(simEnsemble.stateData.values(ind, :, :, :), [4 3 1 2]);
            initial_growth_rates = growth_rates(:, 1);
            PlotUtil.plotScatter(axesHandle, initial_growth_rates * 3600, simEnsemble.stateData.simulationEndTimes, 100, 1:length(initial_growth_rates), 'x');
            xlabel(axesHandle, 'Initial Growth Rate (cell/h)', 'FontSize', 14)
            ylabel(axesHandle, 'Simulation End Time (s)', 'FontSize', 14)
            title(axesHandle, sprintf('Overview: %s',...
                simEnsemble.getTimeStamp()), ...
                'Interpreter', 'None', 'FontSize', 14);
        end
        
        function plotFinalGrowthVsEndTime(simEnsemble, figHandle)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            ind = simEnsemble.getPropertyIndices('growth_rate');
            if isnan(ind)
                return
            end
            clf(figHandle)
            
            axesHandle = subplot(1, 1, 1);
            
            growth_rates = permute(simEnsemble.stateData.values(ind, :, :, :), [4 3 1 2]);
            final_growth_rates = zeros(size(simEnsemble.stateData.simulationEndTimes));
            for i = 1:length(final_growth_rates)
                endLoc = find(simEnsemble.stateData.time == simEnsemble.stateData.simulationEndTimes(i), 1);
                if ~isempty(endLoc)
                    final_growth_rates(i) = growth_rates(i, endLoc);
                end
            end
            PlotUtil.plotScatter(axesHandle, final_growth_rates * 3600, simEnsemble.stateData.simulationEndTimes, 100, 1:length(final_growth_rates), 'x');
            xlabel(axesHandle, 'Final Growth Rate (cell/h)', 'FontSize', 14)
            ylabel(axesHandle, 'Simulation End Time (s)', 'FontSize', 14)
            title(axesHandle, sprintf('Overview: %s',...
                simEnsemble.getTimeStamp()), ...
                'Interpreter', 'None', 'FontSize', 14);
        end
    end
end
