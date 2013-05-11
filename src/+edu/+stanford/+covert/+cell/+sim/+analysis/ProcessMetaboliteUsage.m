% ProcessMetaboliteUsage (Examine metabolite usage by each process)
%
% Author: Derek Macklin, macklin@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Created: 8/11/2011
classdef ProcessMetaboliteUsage
    methods (Static = true)
        function run(simBatchDir, selectedSimulations, fileName)
            import edu.stanford.covert.cell.sim.analysis.ProcessMetaboliteUsage;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            if nargin < 1 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            if nargin < 2 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            metIDs = {'ATP'; 'GTP'};
            data = {'Requirements', 'Allocations', 'Usages'};
            
            for i = 1:numel(metIDs)
                for j = 1: numel(data)
                    [~, figHandle] = PlotUtil.newAxesHandle();
                    ProcessMetaboliteUsage.processPopulationMetaboliteHandler(metIDs{i}, data{j}, figHandle, simBatchDir, selectedSimulations);
                    if nargin >= 3
                        saveas(figHandle, [fileName '-' metIDs{i} '-' data{j} '.pdf']);
                        close(figHandle);
                    end
                end
                
                [~, figHandle] = PlotUtil.newAxesHandle();
                ProcessMetaboliteUsage.processPopulationExpecationMetaboliteHandler(metIDs{i}, 'Usages', figHandle, simBatchDir, selectedSimulations);
                if nargin >= 3
                    saveas(figHandle, [fileName '-Expected-' metIDs{i} '-' data{j} '.pdf']);
                    close(figHandle);
                end
            end
        end
        
        function processPopulationMetaboliteHandler(wIDToShow, desiredInfo, figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            
            if ~exist('desiredInfo', 'var') || isempty(desiredInfo)
                desiredInfo = 'Usages';
            end
            if ~exist('selectedSimulations', 'var') || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            [simIndivDir, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            comp = sim.compartment;
            met = sim.state('Metabolite');
            simTimeStamp = SimulationDiskUtil.getSimulationTimeStamp(simIndivDir);
            
            if ~all(ismember(wIDToShow, met.wholeCellModelIDs))
                throw(MException('Population:error', 'At least one element of ''wIDsToShow'' is invalid'));
            end
            
            desiredInfo(1) = upper(desiredInfo(1));
            
            matfile = sprintf('%s%sProcessMetabolite%s-%s.mat', simBatchDir, filesep, desiredInfo, wIDToShow);
            if exist(matfile, 'file')
                try %#ok<TRYNC>
                    load(matfile);
                end
            end
            
            if ~exist('desInfo', 'var')
                [~, idxsToLoad] = ismember(wIDToShow, met.wholeCellModelIDs);
                if ~any(idxsToLoad == met.hydrophobicIndexs)
                    idxsToLoad = sub2ind([numel(met.wholeCellModelIDs) comp.count], ...
                        idxsToLoad, comp.cytosolIndexs);
                else
                    idxsToLoad = sub2ind([numel(met.wholeCellModelIDs) comp.count], ...
                        idxsToLoad, comp.membraneIndexs);
                end
                
                stateNames = {...
                    'Metabolite' ['process' desiredInfo] idxsToLoad     ':'
                    'Time'       'values'                    ':'        ':'
                    'Mass'       'cell'                      ':'        '-sum'
                    };
                
                states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations);
                
                mass = permute(states.Mass.cell, [3 4 1 2]);
                simEndTimes = permute(max(states.Time.values, [], 3), [4 3 1 2]);
                minSimEndTime = min(simEndTimes);
                [~, simOrder] = sort(mass(minSimEndTime, :));
                
                desInfo = states.Metabolite.(['process' desiredInfo]);  % sum along first dimension of: nCompartments x nProcesses x nTimeSteps x nSimulations
                for i = 1:numel(simEndTimes)
                    desInfo(:, :, simEndTimes(i)+1:end, i) = NaN;
                end
                
                time = squeeze(max(states.Time.values, [], 4)) / 3600;
                try %#ok<TRYNC>
                    save(matfile, 'desInfo', 'time', 'simOrder');
                end
            end
            
            clf(figHandle);
            nProcesses = size(desInfo, 2);
            nCols = 4;
            nRows = nProcesses / nCols;
            axesHandles = zeros(nProcesses, 1);
            
            posL = 0.05;
            posB = 0.05;
            posW = 0.19;
            posH = 0.80;
            
            figTitle = sprintf('%s [%s]\nSimulation Set: %s', desiredInfo, wIDToShow, simTimeStamp);
            annotation(figHandle, 'TextBox', [0.25 0.95 0.5 0.025], ...
                'String', figTitle, 'HorizontalAlignment', 'Center', ...
                'VerticalAlignment', 'Middle', 'EdgeColor', 'None', ...
                'FontSize', 14, 'FontWeight', 'Normal');
            
            for j = 1:nCols
                axesHandles((j - 1) * nRows + 1 : j * nRows) = PlotUtil.multiElementPlot(...
                    figHandle, 3.5 * ones(nRows, 1), time([1 end]), struct(...
                    'position', [(j * posL + (j - 1) * posW) posB posW posH], ...
                    'colorOrder', PlotUtil.getRedGreenColorOrder(simOrder)));
            end
            
            th = zeros(nProcesses, 1);
            for j = 1:nProcesses
                procMet = full(permute(desInfo(:, j, :, :), [3 4 2 1]));
                th(j) = title(axesHandles(j), sim.processMetadata.names{j});
                PlotUtil.plotLine(axesHandles(j), time, procMet, true, true, false);
            end
        end
        
        function processPopulationExpecationMetaboliteHandler(wIDToShow, desiredInfo, figHandle, simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            
            if ~isequal(desiredInfo, 'Usages')
                throw(MException('ProcessMetaboliteUsage:error', 'Usages is the only valid desiredInfo'));
            end
            if nargin < 5 || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end
            
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            comp = sim.compartment;
            met = sim.state('Metabolite');
            metIdx = met.getIndexs(wIDToShow);
            if any(metIdx == met.hydrophobicIndexs)
                compIdx = comp.membraneIndexs;
            else
                compIdx = comp.cytosolIndexs;
            end
            metCompIdx = sub2ind([numel(met.wholeCellModelIDs) comp.count], metIdx, compIdx);
            
            expected = met.processBiomassProduction(metIdx, :) - met.processByproduct(metIdx, :);
            
            stateNames = {...
                'Metabolite' 'processUsages'  metCompIdx  ':'
                };
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulations);
            states.Metabolite.processUsages = full(sum(states.Metabolite.processUsages, 3));
            
            clf(figHandle);
            nProcesses = numel(sim.processes);
            nCols = 4;
            nRows = nProcesses / nCols;
            
            figTitle = sprintf('%s [%s]\nSimulation Set: %s', desiredInfo, wIDToShow, 'Population Expectation');
            annotation(figHandle, 'TextBox', [0.1 0.95 0.8 0.025], ...
                'String', figTitle, 'HorizontalAlignment', 'Center', ...
                'VerticalAlignment', 'Middle', 'EdgeColor', 'None', ...
                'FontSize', 14, 'FontWeight', 'Normal');
            
            [axesHandles, xAxisHandles] = PlotUtil.multiElementPlot(...
                figHandle, repmat({1 * ones(nRows, 1)}, nCols, 1), [0 1]);
            
            for i = 1:nProcesses
                axesHandle = axesHandles{ceil(i/nRows)}(mod(i - 1, nRows) + 1);
                set(axesHandle, 'FontSize', 6, 'XColor', 'k', 'XLimMode', 'auto', 'XTickMode', 'auto', 'YLimMode', 'auto', 'YTickMode', 'auto');
                hist(axesHandle, squeeze(states.Metabolite.processUsages(:, i, :, :)));
                set(axesHandle, 'FontSize', 6, 'XColor', 'k', 'XLimMode', 'auto', 'XTickMode', 'auto', 'YLimMode', 'auto', 'YTickMode', 'auto');
                line(expected(:, [i i]), ylim(axesHandle), 'Color', 'r', 'Parent', axesHandle);
                title(axesHandle, sim.processMetadata.names{i});
                ylabel(axesHandle, 'Freq', 'FontSize', 6);
            end
            
            set(cell2mat(xAxisHandles), 'Visible', 'off')
        end
        
        function processSingleCellMetaboliteHandler(wIDsToShow, desiredInfo, figHandle, simBatchDir, selectedSimulation)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            
            if ~exist('desiredInfo', 'var') || isempty(desiredInfo)
                desiredInfo = 'Usages';
            end
            if ~exist('selectedSimulation', 'var') || isempty(selectedSimulation)
                selectedSimulation = SimulationDiskUtil.getCompleteSimulations(simBatchDir);
            end            
            
            [simIndivDir, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            comp = sim.compartment;
            met = sim.state('Metabolite');
            simTimeStamp = SimulationDiskUtil.getSimulationTimeStamp(simIndivDir);
            
            if ~all(ismember(wIDsToShow, met.wholeCellModelIDs))
                throw(MException('Population:error', 'At least one element of ''wIDsToShow'' is invalid'));
            end
            
            desiredInfo(1) = upper(desiredInfo(1));
            wIDsStr = wIDsToShow{1};
            
            if numel(wIDsToShow) > 1
                for i = 2:numel(wIDsToShow)
                    wIDsStr = [wIDsStr ', ' wIDsToShow{i}]; %#ok<AGROW>
                end
            end
            
            
            [~, idxsToLoad] = ismember(wIDsToShow, met.wholeCellModelIDs);
            for i = 1:numel(idxsToLoad)
                if ~any(idxsToLoad(i) == met.hydrophobicIndexs)
                    idxsToLoad(i) = sub2ind([numel(met.wholeCellModelIDs) comp.count], ...
                        idxsToLoad(i), comp.cytosolIndexs);
                else
                    idxsToLoad(i) = sub2ind([numel(met.wholeCellModelIDs) comp.count], ...
                        idxsToLoad(i), comp.membraneIndexs);
                end
            end
            
            stateNames = {...
                'Metabolite' ['process' desiredInfo] idxsToLoad'     ':'
                'Time'       'values'                    ':'        ':'
                };
            
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', selectedSimulation);
            
            
            desInfo = states.Metabolite.(['process' desiredInfo]);  % sum along first dimension of: nCompartments x nProcesses x nTimeSteps x nSimulations
            
            time = squeeze(max(states.Time.values, [], 4)) / 3600;

           
            
            clf(figHandle);
            nProcesses = size(desInfo, 2);
            nCols = 4;
            nRows = nProcesses / nCols;
            axesHandles = zeros(nProcesses, 1);
            
            posL = 0.05;
            posB = 0.05;
            posW = 0.19;
            posH = 0.80;
            
            figTitle = sprintf('%s [%s]\nSimulation Set: %s', desiredInfo, wIDsStr, simTimeStamp);
            annotation(figHandle, 'TextBox', [0.25 0.95 0.5 0.025], ...
                'String', figTitle, 'HorizontalAlignment', 'Center', ...
                'VerticalAlignment', 'Middle', 'EdgeColor', 'None', ...
                'FontSize', 14, 'FontWeight', 'Normal');
            
            for j = 1:nCols
                axesHandles((j - 1) * nRows + 1 : j * nRows) = PlotUtil.multiElementPlot(...
                    figHandle, 3.5 * ones(nRows, 1), time([1 end]), struct(...
                    'position', [(j * posL + (j - 1) * posW) posB posW posH]));
            end
            
            th = zeros(nProcesses, 1);
            for j = 1:nProcesses
                procMet = full(permute(desInfo(:, j, :), [1 3 2]));
                th(j) = title(axesHandles(j), sim.processMetadata.names{j});
                PlotUtil.plotLine(axesHandles(j), time, procMet, true, true, false);
            end
        end
    end
end