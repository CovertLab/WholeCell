% ChromosomeSpaceTimePlot
% Plot space-time densities
%
% Author: Derek Macklin, macklin@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 8/14/2011
classdef ChromosomeSpaceTimePlot
    properties(Constant = true)
        cmap = [...
            1.0000    1.0000    1.0000
            0.3059    0.3961    0.5804
            0.2994    0.4013    0.5688
            0.2928    0.4066    0.5572
            0.2863    0.4118    0.5457
            0.2797    0.4170    0.5341
            0.2732    0.4222    0.5225
            0.2667    0.4275    0.5109
            0.2601    0.4327    0.4994
            0.2536    0.4379    0.4878
            0.2471    0.4432    0.4762
            0.2405    0.4484    0.4646
            0.2340    0.4536    0.4531
            0.2274    0.4588    0.4415
            0.2209    0.4641    0.4299
            0.2144    0.4693    0.4183
            0.2078    0.4745    0.4068
            0.2013    0.4798    0.3952
            0.1948    0.4850    0.3836
            0.1882    0.4902    0.3720
            0.1817    0.4954    0.3605
            0.1751    0.5007    0.3489
            0.1686    0.5059    0.3373
            0.2043    0.4877    0.3195
            0.2400    0.4696    0.3018
            0.2757    0.4514    0.2840
            0.3114    0.4332    0.2663
            0.3472    0.4151    0.2485
            0.3829    0.3969    0.2308
            0.4186    0.3788    0.2130
            0.4543    0.3606    0.1953
            0.4900    0.3424    0.1775
            0.5257    0.3243    0.1598
            0.5614    0.3061    0.1420
            0.5971    0.2879    0.1243
            0.6328    0.2698    0.1065
            0.6685    0.2516    0.0888
            0.7043    0.2335    0.0710
            0.7400    0.2153    0.0533
            0.7757    0.1971    0.0355
            0.8114    0.1790    0.0178
            0.8471    0.1608         0
            0.8486    0.1814         0
            0.8500    0.2020         0
            0.8515    0.2226         0
            0.8530    0.2431         0
            0.8544    0.2637         0
            0.8559    0.2843         0
            0.8574    0.3049         0
            0.8589    0.3255         0
            0.8603    0.3461         0
            0.8618    0.3667         0
            0.8633    0.3873         0
            0.8647    0.4078         0
            0.8662    0.4284         0
            0.8677    0.4490         0
            0.8691    0.4696         0
            0.8706    0.4902         0
            0.7255    0.4085         0
            0.5804    0.3268         0
            0.4353    0.2451         0
            0.2902    0.1634         0
            0.1451    0.0817         0
            0         0         0
            ];
        
        rings = struct(...
            'RNAPolymerase', 7, ...
            'SMC', 2, ...
            'Helicase', 1, ...
            'DnaA', 4, ...
            'Gyrase', 5, ...
            'TopoisomeraseIV', 6, ...
            'Regulation', 3);
        
        chromProteins = {...
            % chromosome property    whole cell ids               strands    color (rgb)               name
            'complexBoundSites',    {'MG_213_214_298_6MER_ADP'},   1:2,     [0 1 0],                'SMC';
            'complexBoundSites',    {'MG_094_HEXAMER'},            1:4,     [0.6980 0.1333 0.1333], 'Helicase';
            'complexBoundSites',    {'MG_469_1MER_ADP', ...
                                     'MG_469_1MER_ATP', ...
                                     'MG_469_2MER_ATP', ...
                                     'MG_469_3MER_ATP', ...
                                     'MG_469_4MER_ATP', ...
                                     'MG_469_5MER_ATP', ...
                                     'MG_469_6MER_ATP', ...
                                     'MG_469_7MER_ATP'},           1:2,     [1.0000 0.6471 0],      'DnaA';
            'complexBoundSites',     {'MG_469_2MER_ATP', ...
                                     'MG_469_3MER_ATP', ...
                                     'MG_469_4MER_ATP', ...
                                     'MG_469_5MER_ATP', ...
                                     'MG_469_6MER_ATP', ...
                                     'MG_469_7MER_ATP'},           1:2,     [1.0000 0.6471 0],      'DnaA Complex';
            'complexBoundSites',    {'DNA_GYRASE'},                1:2,     [0.9 0.9 0.9],          'Gyrase';
            'complexBoundSites',    {'MG_203_204_TETRAMER'},       1:2,     [0.6902 0.8784 0.9020], 'Topoisomerase IV';
            'complexBoundSites',    {'MG_428_DIMER'},              1:2,     [1 1 0],                'LuxR';
            'monomerBoundSites',    {'MG_101_MONOMER'},            1:2,     [1 0.5 0],              'HTH Regulator';
            'monomerBoundSites',    {'MG_236_MONOMER'},            1:2,     [0.2 1 0.8],            'Ferric uptake repressor';
            };
        
        regulationWIDs = {'MG_428_DIMER' 'MG_101_MONOMER' 'MG_236_MONOMER'};        
    end
    
    methods (Static = true);
        function run(relSimPath, selectedTime, fileName)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.analysis.ChromosomeSpaceTimePlot;
            
            if ~exist('relSimPath', 'var') || isempty(relSimPath)
                relSimPath = [SimulationDiskUtil.getLatestSimulationGroup() filesep '1'];
            end
            if ~exist('selectedTime', 'var') || isempty(selectedTime)
                selectedTime = 4 * 3600;
            end
            
            [simDir, ~, sim] = SimulationDiskUtil.getSimulation(relSimPath);
            
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
            
            %% space-time density plots
            
            % Active RNA Polymerase
            RNAPSparseMat = ChromosomeSpaceTimePlot.makeRNAPSparseMat(states, sim);
            
            % Space-Time Overlay Plot
            [ax1, figHandle1] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            whatToShow = {'Gyrase', 'Topoisomerase IV', 'DnaA Complex', 'LuxR', 'HTH Regulator', 'Ferric uptake repressor', 'Helicase'};
            ChromosomeSpaceTimePlot.plotSpaceTimeOverlay(...
                ax1, sim, states, RNAPSparseMat, ...
                ChromosomeSpaceTimePlot.chromProteins, whatToShow, ...
                sprintf('Chromosome 1\nSimulation: %s #%s', simTimeStamp, simIdx));
            if exist('fileName', 'var')
                saveas(figHandle1, [fileName '-SpaceTimeOverlay.pdf']);
                close(figHandle1);
            end
            
            % Circular Density Plot
            [ax3, figHandle3] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            set(figHandle3, 'Color', [1 1 1]);
            [densityBEMatrix, maxTime, maxSpace] = ChromosomeSpaceTimePlot.getProteinBindingDensityMatrix(...
                sim, states, RNAPSparseMat, ...
                ChromosomeSpaceTimePlot.chromProteins, [], 1000, 1000);
            ChromosomeSpaceTimePlot.plotCircularDensity(...
                ax3, densityBEMatrix, maxTime, sim, false, ...
                sprintf('Circular Chromosome Density Plot\nSimulation: %s #%s', simTimeStamp, simIdx));
            if exist('fileName', 'var')
                saveas(figHandle3, [fileName '-CircularDensity.pdf']);
                close(figHandle3);
            end
            
            % Ring Plot with all proteins
            [ax4, figHandle4] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            ChromosomeSpaceTimePlot.makeRingPlot(...
                ax4, sim, states, RNAPSparseMat, ...
                ChromosomeSpaceTimePlot.chromProteins, [], selectedTime, ...
                maxSpace);
            if exist('fileName', 'var')
                saveas(figHandle4, [fileName '-Ring.pdf']);
                close(figHandle4);
            end
        end
        
        function [densityMatrix, pxWidth, pxHeight] = calcRnaPDensityMatrix(activePos, width, height, maxTime, maxSpace)
            pxWidth = maxTime / width;
            pxHeight = maxSpace / height;
            densityMatrix = zeros(height, width);
            for i = 1:height
                for j = 1:width
                    densityMatrix(i, j) = size(...
                        activePos(...
                        activePos(:, round(pxWidth * (j-1) + 1) : round(pxWidth * j)) > round(pxHeight * (i-1)) & ...
                        activePos(:, round(pxWidth * (j-1) + 1) : round(pxWidth * j)) <= round(pxHeight * i)), ...
                        1) / (pxWidth * pxHeight) * 3600;
                end
            end
        end
        
        function [densityMatrix, pxWidth, pxHeight] = calcChromBoundDensityMatrix(subs, width, height, maxTime, maxSpace)
            pxWidth = maxTime / width;
            pxHeight = maxSpace / height;
            
            w = round(pxWidth * (1:width));
            h = round(pxHeight * (1:height));
            
            densityMatrix =  hist3(subs(:, [1 3]), 'Edges', {h w}) / (pxWidth * pxHeight) * 3600;
        end
        
        function plotSpaceTimeOverlay(ax, sim, states, RNAPSparseMat, chProteins, whatToShow, figTitle)
            %% Initialization %%
            if exist('whatToShow', 'var') && ~isempty(whatToShow)
                [~, inds] = ismember(whatToShow, chProteins(:, 5));
                chProteins = chProteins(inds(inds ~= 0), :);
            end
            
            c = sim.state('Chromosome');
            hold(ax, 'on');
            
            % Handle presence/absence of RNAPSparseMat (plot it or don't)
            if exist('RNAPSparseMat', 'var') && ~isempty(RNAPSparseMat)
                pol = sim.state('RNAPolymerase');
                subs = find(RNAPSparseMat >= pol.activelyTranscribingValue);
                subs(subs(:, 2) >= 3, :) = [];
                
                % Flip positions above/below the terC
                subs_old = subs;
                subs(subs_old(:, 1) > c.terCPosition, 1) = subs_old(subs_old(:, 1) > c.terCPosition, 1) - c.terCPosition;
                subs(subs_old(:, 1) <= c.terCPosition, 1) = subs_old(subs_old(:, 1) <= c.terCPosition, 1) + c.terCPosition;
                clear subs_old;
               
                handles = zeros(size(chProteins, 1) + 1, 1);
                names = ['RNA Polymerase'; chProteins(:, 5)];
                handles(1) = patch([...
                    1 / 3600 * subs(:, 3)'
                    1 / 3600 * (subs(:, 3) + 1)'
                    1 / 3600 * (subs(:, 3) + 1)'
                    1 / 3600 * subs(:, 3)'
                    ], [
                    subs(:, 1)'
                    subs(:, 1)'
                    subs(:, 1)' + 1
                    subs(:, 1)' + 1
                    ], ones(1, size(subs, 1)), 'Parent', ax, 'EdgeColor', 'b', 'FaceColor', 'b');
                handle_offset = 1;
            else
                handles = zeros(size(chProteins, 1), 1);
                names = chProteins(:, 5);
                handle_offset = 0;
            end
            
            %% Loop over chromosomal proteins %%
            for i = 1:size(chProteins, 1)
                chromProp = chProteins{i, 1};
                wIDs = chProteins{i, 2};
                strands = chProteins{i, 3};
                color = chProteins{i, 4};
                
                % Determine if it's a complex or monomer
                if strcmp(chromProp, 'complexBoundSites')
                    mod = sim.state('ProteinComplex');
                elseif strcmp(chromProp, 'monomerBoundSites')
                    mod = sim.state('ProteinMonomer');
                end
                
                idxs = find(ismember(mod.wholeCellModelIDs(mod.matureIndexs), wIDs));
                subs = [];
                
                % If we lump multiple proteins together to consider them as
                % one class of proteins (e.g., DnaA), then get all subs
                for j = 1:numel(wIDs)
                    subs = [subs; find(states.Chromosome.(chromProp) == idxs(j))]; %#ok<AGROW>
                end
                
                % Get rid of subs that aren't on the proper strand
                subs = subs(ismember(subs(:, 2), strands), :);
                
                % Flip positions above/below the terC
                subs_old = subs;
                subs(subs_old(:, 1) > c.terCPosition, 1) = subs_old(subs_old(:, 1) > c.terCPosition, 1) - c.terCPosition;
                subs(subs_old(:, 1) <= c.terCPosition, 1) = subs_old(subs_old(:, 1) <= c.terCPosition, 1) + c.terCPosition;
                clear subs_old;
                
                % Plot
                if isempty(subs)
                    continue;
                end
                handles(i + handle_offset) = patch([...
                    1 / 3600 * subs(:, 3)'
                    1 / 3600 * (subs(:, 3) + 1)'
                    1 / 3600 * (subs(:, 3) + 1)'
                    1 / 3600 * subs(:, 3)'
                    ], [
                    subs(:, 1)'
                    subs(:, 1)'
                    subs(:, 1)' + 1
                    subs(:, 1)' + 1
                    ], ones(1, size(subs, 1)), 'Parent', ax, 'EdgeColor', color, 'FaceColor', color);                
            end
            hold(ax, 'off');
            
            % Legend formatting
            lh = legend(handles(handles ~= 0), names(handles ~= 0), 'Location', 'NorthEastOutside', 'FontSize', 8);
            c_lh = get(lh, 'Children');
            for i = 1:numel(c_lh)
                set(c_lh(i), 'LineStyle', '-');
            end
            
            % Labeling
            xlabel(ax, 'Time (hr)', 'FontSize', 16);
            set(ax, 'YTick', [1 c.sequenceLen/2 c.sequenceLen]);
            set(ax, 'YTickLabel', {'TerC', 'OriC', 'TerC'});
            ylim(ax, [0.5 c.sequenceLen + 0.5]);
            ylabel(ax, 'Position', 'FontSize', 16)
            set(ax, 'FontSize', 8);
            xlim(ax, states.Time.values([1 end]) / 3600);
        end
        
        function makeRingPlot(ax, sim, states, RNAPSparseMat, chProteins, whatToMake, time, maxSpace)
            import edu.stanford.covert.cell.sim.analysis.ChromosomeSpaceTimePlot;
            
            if exist('whatToMake', 'var') && ~isempty(whatToMake)
                [~, inds] = ismember(whatToMake, chProteins(:, 5));
                chProteins = chProteins(inds(inds ~= 0), :);
            end
            
            if exist('RNAPSparseMat', 'var') && ~isempty(RNAPSparseMat)
                pol = sim.state('RNAPolymerase');
                subs = find(RNAPSparseMat >= pol.activelyTranscribingValue);
                subs(subs(:, 2) >= 3, :) = [];
                pos = ChromosomeSpaceTimePlot.getPos(subs, time);
                lineHandles = zeros(size(chProteins, 1) + 1, 1);
                names = ['RNA Polymerase'; chProteins(:, 5)];
                h = ChromosomeSpaceTimePlot.plotRing(ax, pos, maxSpace, [1 2], [0 0 1]);
                if ~isempty(h)
                    lineHandles(1) = h(1);
                end
                lineHandle_offset = 1;
            else
                lineHandles = zeros(size(chProteins, 1), 1);
                names = chProteins(:, 5);
                lineHandle_offset = 0;
            end
            
            for i = 1:size(chProteins, 1)
                chromProp = chProteins{i, 1};
                wIDs = chProteins{i, 2};
                strands = chProteins{i, 3};
                color = chProteins{i, 4};
                if strcmp(chromProp, 'complexBoundSites')
                    mod = sim.state('ProteinComplex');
                elseif strcmp(chromProp, 'monomerBoundSites')
                    mod = sim.state('ProteinMonomer');
                end
                idxs = find(ismember(mod.wholeCellModelIDs(mod.matureIndexs), wIDs));
                subs = [];
                for j = 1:numel(wIDs)
                    subs = [subs; find(states.Chromosome.(chromProp) == idxs(j))]; %#ok<AGROW>
                end
                subs = subs(ismember(subs(:, 2), strands), :);
                
                pos = ChromosomeSpaceTimePlot.getPos(subs, time);
                h = ChromosomeSpaceTimePlot.plotRing(ax, pos, maxSpace, [1 2], color);
                if ~isempty(h)
                    lineHandles(i + lineHandle_offset) = h(1);
                end
            end
            
            if any(lineHandles)
                legend(lineHandles, names, 'Location', 'NorthEastOutside');
            end
            
            axis(ax, 'square');
        end
        
        function plotSpaceTimeDensity(axesHandle, densityMatrix, maxTime, pxWidth, figTitle)
            imagesc(densityMatrix, 'Parent', axesHandle);
            hrPos = (0:ceil(maxTime / 3600)) * (3600 / pxWidth);
            set(axesHandle, 'XTick', hrPos);
            set(axesHandle, 'XTickLabel', hrPos * pxWidth/3600);
            set(axesHandle, 'YTick', [0 size(densityMatrix, 1) / 2 size(densityMatrix, 1)]);
            set(axesHandle, 'YTickLabel', {'OriC', 'TerC', 'OriC'});
            set(axesHandle, 'YDir', 'normal');
            axis(axesHandle, [0 size(densityMatrix, 2) 0 size(densityMatrix, 1)]);
            xlabel(axesHandle, 'Time (hr)');
            ylabel(axesHandle, 'Position on Chromosome 1');
            colormap(axesHandle, edu.stanford.covert.cell.sim.analysis.ChromosomeSpaceTimePlot.cmap);
            c_handle = colorbar('Peer', axesHandle);
            ylbl = get(c_handle, 'YLabel');
            set(ylbl, 'String', 'Density (Proteins/bp/hr)');
            title(figTitle, 'Parent', axesHandle);
            xlim(axesHandle, [0 maxTime] / 3600);
            box('off');
        end
        
        function [lineHandles, circleHandles] = plotRing(axesHandle, pos, maxPos, radii, color)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            lineHandles = zeros(size(pos));
            circleHandles = zeros(2, 1);
            
            theta = ((pos / maxPos) * (2 * pi)) + (pi / 2);
            
            for i = 1:size(pos, 1)
                lineXData = radii .* cos(theta(i));
                lineYData = radii .* sin(theta(i));
                lineHandles(i) = line(lineXData, lineYData, 'Color', color, 'Parent', axesHandle);
            end
            
            circleHandles(1) = PlotUtil.drawCircle(0, 0, radii(1), axesHandle);
            circleHandles(2) = PlotUtil.drawCircle(0, 0, radii(2), axesHandle);
            set(axesHandle, 'Visible', 'off');
        end
        
        % Plot circular density
        % TODO: Use loops rather than spaghetti code!
        function plotCircularDensity(axesHandle, densityBEMatrix, maxTime, simulation, indepRingNormalize, figTitle)
            import edu.stanford.covert.cell.sim.analysis.ChromosomeSpaceTimePlot;
            
            ids = simulation.gene.wholeCellModelIDs;
            colormap(flipud(colormap('bone')));
            if indepRingNormalize
                h = bullseye(log10((densityBEMatrix - repmat(min(densityBEMatrix, [], 1), [size(densityBEMatrix, 1) 1])) ./ repmat(range(densityBEMatrix, 1), [size(densityBEMatrix, 1) 1])), 'N', 10, 'tht0', 180, 'axesHandle', axesHandle);
                colorbar_label = 'Relative Protein Count On Each Ring';
                c_handle = colorbar();
            else
                scaleFactor = 1000;
                h = bullseye(log10(scaleFactor * densityBEMatrix * (maxTime/3600) + 1), 'N', 10, 'tht0', 180, 'axesHandle', axesHandle);
                colorbar_label = 'Protein / nt / cell cycle';
                c_handle = colorbar();
                yticks = get(c_handle, 'YTick');
                yticklabels = ((10 .^ yticks) - 1) / (scaleFactor);
                yticklabels_str = cell(size(yticklabels));
                for i = 1:numel(yticklabels_str)
                    yticklabels_str{i} = sprintf('%0.02g', yticklabels(i));
                end
                set(c_handle, 'YTickLabel', yticklabels_str);
            end
            set(h, 'ButtonDownFcn', 'edu.stanford.covert.cell.sim.analysis.ChromosomeSpaceTimePlot.getClickPos(maxSpace)');
            ylbl = get(c_handle, 'YLabel');
            set(ylbl, 'String', colorbar_label);
            line([0 0], [5 5.1], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'Parent', axesHandle);
            text(0, 5.15, 'OriC', 'EdgeColor', 'none', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Rotation', 270, 'Parent', axesHandle);
            line([0 0], [-5 -5.1], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'Parent', axesHandle);
            text(0, -5.15, 'TerC', 'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Rotation', 270, 'Parent', axesHandle);
            
            thisIdx = 297;
            theta = (pi / 2) + simulation.gene.startCoordinates(thisIdx) / 580076 * 2 * pi;
            [x, y] = pol2cart(theta, [5 5.1]);
            line(x, y, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'Parent', axesHandle);
            [x, y] = pol2cart(theta, 5.15);
            halign = 'left';
            if theta >= (pi / 2) && theta < 3 * (pi / 2)
                theta = theta + pi;
                halign = 'right';
            end
            text(x, y, ids{thisIdx}, 'EdgeColor', 'none', 'HorizontalAlignment', halign, 'VerticalAlignment', 'middle', 'Interpreter', 'None', 'Rotation', theta * 180 / pi, 'Parent', axesHandle);
            
            thisIdx = 302;
            theta = (pi / 2) + simulation.gene.startCoordinates(thisIdx) / 580076 * 2 * pi;
            [x, y] = pol2cart(theta, [5 5.1]);
            line(x, y, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'Parent', axesHandle);
            [x, y] = pol2cart(theta, 5.15);
            halign = 'left';
            if theta >= (pi / 2) && theta < 3 * (pi / 2)
                theta = theta + pi;
                halign = 'right';
            end
            text(x, y, ids{thisIdx}, 'EdgeColor', 'none', 'HorizontalAlignment', halign, 'VerticalAlignment', 'middle', 'Interpreter', 'None', 'Rotation', theta * 180 / pi, 'Parent', axesHandle);
            
            thisIdx = [229 236];
            theta = (pi / 2) + simulation.gene.startCoordinates(thisIdx(1)) / 580076 * 2 * pi;
            [x, y] = pol2cart(theta, [5 5.1]);
            line(x, y, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'Parent', axesHandle);
            [x, y] = pol2cart(theta, 5.15);
            halign = 'left';
            if theta >= (pi / 2) && theta < 3 * (pi / 2)
                theta = theta + pi;
                halign = 'right';
            end
            text(x, y, [ids{thisIdx(1)} '-' ids{thisIdx(2)}], 'EdgeColor', 'none', 'HorizontalAlignment', halign, 'VerticalAlignment', 'middle', 'Interpreter', 'None', 'Rotation', theta * 180 / pi, 'Parent', axesHandle);
            
            thisIdx = [144 148];
            theta = (pi / 2) + simulation.gene.startCoordinates(thisIdx(1)) / 580076 * 2 * pi;
            [x, y] = pol2cart(theta, [5 5.1]);
            line(x, y, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'Parent', axesHandle);
            [x, y] = pol2cart(theta, 5.15);
            halign = 'left';
            if theta >= (pi / 2) && theta < 3 * (pi / 2)
                theta = theta + pi;
                halign = 'right';
            end
            text(x, y, [ids{thisIdx(1)} '-' ids{thisIdx(2)}], 'EdgeColor', 'none', 'HorizontalAlignment', halign, 'VerticalAlignment', 'middle', 'Interpreter', 'None', 'Rotation', theta * 180 / pi, 'Parent', axesHandle);
            
            thisIdx = 400;
            theta = (pi / 2) + simulation.gene.startCoordinates(thisIdx) / 580076 * 2 * pi;
            [x, y] = pol2cart(theta, [5 5.1]);
            line(x, y, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'Parent', axesHandle);
            [x, y] = pol2cart(theta, 5.15);
            halign = 'left';
            if theta >= (pi / 2) && theta < 3 * (pi / 2)
                theta = theta + pi;
                halign = 'right';
            end
            text(x, y, ids{thisIdx}, 'EdgeColor', 'none', 'HorizontalAlignment', halign, 'VerticalAlignment', 'middle', 'Interpreter', 'None', 'Rotation', theta * 180 / pi, 'Parent', axesHandle);
            
            thisIdx = 501;
            theta = (pi / 2) + simulation.gene.startCoordinates(thisIdx) / 580076 * 2 * pi;
            [x, y] = pol2cart(theta, [5 5.1]);
            line(x, y, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'Parent', axesHandle);
            [x, y] = pol2cart(theta, 5.15);
            halign = 'left';
            if theta >= (pi / 2) && theta < 3 * (pi / 2)
                theta = theta + pi;
                halign = 'right';
            end
            text(x, y, ids{thisIdx}, 'EdgeColor', 'none', 'HorizontalAlignment', halign, 'VerticalAlignment', 'middle', 'Interpreter', 'None', 'Rotation', theta * 180 / pi, 'Parent', axesHandle);
            
            thisIdx = [16 17];
            theta = (pi / 2) + simulation.gene.startCoordinates(thisIdx(1)) / 580076 * 2 * pi;
            [x, y] = pol2cart(theta, [5 5.1]);
            line(x, y, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'Parent', axesHandle);
            [x, y] = pol2cart(theta, 5.15);
            halign = 'left';
            if theta >= (pi / 2) && theta < 3 * (pi / 2)
                theta = theta + pi;
                halign = 'right';
            end
            text(x, y, [ids{thisIdx(1)} '-' ids{thisIdx(2)}], 'EdgeColor', 'none', 'HorizontalAlignment', halign, 'VerticalAlignment', 'middle', 'Interpreter', 'None', 'Rotation', theta * 180 / pi, 'Parent', axesHandle);
            
            thisIdx = [104 106];
            theta = (pi / 2) + simulation.gene.startCoordinates(thisIdx(1)) / 580076 * 2 * pi;
            [x, y] = pol2cart(theta, [5 5.1]);
            line(x, y, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'Parent', axesHandle);
            [x, y] = pol2cart(theta, 5.15);
            halign = 'left';
            if theta >= (pi / 2) && theta < 3 * (pi / 2)
                theta = theta + pi;
                halign = 'right';
            end
            text(x, y, [ids{thisIdx(1)} '-' ids{thisIdx(2)}], 'EdgeColor', 'none', 'HorizontalAlignment', halign, 'VerticalAlignment', 'middle', 'Interpreter', 'None', 'Rotation', theta * 180 / pi, 'Parent', axesHandle);
            
            thisIdx = 139;
            theta = (pi / 2) + simulation.gene.startCoordinates(thisIdx) / 580076 * 2 * pi;
            [x, y] = pol2cart(theta, [5 5.1]);
            line(x, y, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'Parent', axesHandle);
            [x, y] = pol2cart(theta, 5.15);
            halign = 'left';
            if theta >= (pi / 2) && theta < 3 * (pi / 2)
                theta = theta + pi;
                halign = 'right';
            end
            text(x, y, ids{thisIdx}, 'EdgeColor', 'none', 'HorizontalAlignment', halign, 'VerticalAlignment', 'middle', 'Interpreter', 'None', 'Rotation', theta * 180 / pi, 'Parent', axesHandle);
            
            r = ChromosomeSpaceTimePlot.rings;
            text(0, 1 + (r.RNAPolymerase - 1) * (4 / 7), 'RNA Polymerase', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', 'White', 'Parent', axesHandle);
            text(0, 1 + (r.SMC - 1) * (4 / 7), 'SMC', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', 'White', 'Parent', axesHandle);
            text(0, 1 + (r.Helicase - 1) * (4 / 7), 'Helicase', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', 'Black', 'Parent', axesHandle);
            text(0, 1 + (r.DnaA - 1) * (4 / 7), 'DnaA', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', 'White', 'Parent', axesHandle);
            text(0, 1 + (r.Gyrase - 1) * (4 / 7), 'Gyrase', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', 'White', 'Parent', axesHandle);
            text(0, 1 + (r.TopoisomeraseIV - 1) * (4 / 7), 'Topoisomerase IV', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', 'Black', 'Parent', axesHandle);
            text(0, 1 + (r.Regulation - 1) * (4 / 7), 'Regulation', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', 'Black', 'Parent', axesHandle);
            
            title(axesHandle, figTitle);
        end
        
        function pos = getRnaPPos(activePos, time)
            pos = activePos(:, time);
            pos(isnan(pos)) = [];
        end
        
        function pos = getPos(subs, time)
            pos = subs(subs(:, 3) == time, 1);
        end
        
        function getClickPos(maxVal)
            pt = get(gca, 'CurrentPoint');
            theta = cart2pol(pt(1, 1), pt(1, 2));
            theta = theta - (pi / 2);
            if theta < 0
                theta = theta + (2 * pi);
            elseif theta > (2 * pi)
                theta = theta - (2 * pi);
            end
            fprintf('Clicked on nucleotide: %05g\n', (theta * maxVal / (2 * pi)));
        end
        
        function RNAPSparseMat = makeRNAPSparseMat(states, sim)
            c = sim.state('Chromosome');
            
            nTimePoints = size(states.RNAPolymerase.positionStrands, 3);
            nRnaPol = size(states.RNAPolymerase.positionStrands, 1);
            polStates = reshape(states.RNAPolymerase.states, [nRnaPol * nTimePoints 1]);
            posStrndTimes = reshape(permute([states.RNAPolymerase.positionStrands repmat(permute(1:nTimePoints, [1 3 2]), [nRnaPol 1])], [1 3 2]), [nRnaPol * nTimePoints 3]);
            
            polStates = polStates(posStrndTimes(:, 1) ~= 0, :);
            posStrndTimes = posStrndTimes(posStrndTimes(:, 1) ~= 0, :);
            
            RNAPSparseMat = edu.stanford.covert.util.CircularSparseMat(posStrndTimes, polStates, [c.sequenceLen 4 nTimePoints], 1);
        end
        
        function [ensemble, states, dnaABoxBound, legendHandle] = plotSpaceTime(...
                axesHandle, simBatchDir, simNum, nRNAPols, ...
                nonSpecBoundPols, posSpan, showSingleDnaALine, showLegend, ...
                ensemble, states, dnaABoxBound, ...
                dnaPolLineWidth, rnaPolLineWidth, showRNAPolDetails)
            import edu.stanford.covert.cell.sim.analysis.ChromosomeSpaceTimePlot;
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            
            [simDir, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep num2str(simNum)]);
            c = sim.state('Chromosome');
            rnaPol = sim.state('RNAPolymerase');
            repInit = sim.process('ReplicationInitiation');
            
            if nargin < 11
                stateNames = {
                    'helicase1'         'Helicase 1'
                    'helicase2'         'Helicase 2'
                    'leadingPol1'       'Leading Polymerase 1'
                    'leadingPol2'       'Leading Polymerase 2'
                    'laggingPol1'       'Lagging Polymerase 1'
                    'laggingPol2'       'Lagging Polymerase 2'
                    'dnaA_box1'         'DnaA Box R1'
                    'dnaA_box2'         'DnaA Box R2'
                    'dnaA_box3'         'DnaA Box R3'
                    'dnaA_box4'         'DnaA Box R4'
                    'dnaA_box5'         'DnaA Box R5'
                    };
                ensemble = SimulationEnsemble(simBatchDir, stateNames, [], simNum);
                
                stateNames = {
                    'RNAPolymerase'      'states'
                    'RNAPolymerase'      'positionStrands'
                    'Transcript'         'boundTranscriptionUnits'
                    'Transcript'         'boundTranscriptProgress'
                    'Transcript'         'boundTranscriptChromosome'
                    };
                states = DiskLogger.load(simDir, stateNames, [], [], [], 'extract');
                
                stateNames = {
                    'Chromosome'         'complexBoundSites'
                    };
                tmp = DiskLogger.load(simDir, stateNames, [], [], [], 'extract');
                [subs, vals] = find(tmp.Chromosome.complexBoundSites(repInit.dnaABoxStartPositions, :, :));
                tfs = subs(:, 2) <= 2 & ismember(vals, repInit.enzymeComplexGlobalIndexs);
                isDnaABoxBound = zeros(size(repInit.dnaABoxStartPositions, 1), size(tmp.Chromosome.complexBoundSites, 3));
                isDnaABoxBound(sub2ind(size(isDnaABoxBound), ...
                    subs(tfs, 1), ...
                    subs(tfs, 3))) = 1;
                startIdxs = find([isDnaABoxBound(:, 1)  diff(isDnaABoxBound, [], 2) > 0]);
                endIdxs = find([diff(isDnaABoxBound, [], 2) < 0  isDnaABoxBound(:, end)]);
                [tmpx1, tmpy1] = ind2sub(size(isDnaABoxBound), startIdxs);
                [tmpx2, tmpy2] = ind2sub(size(isDnaABoxBound), endIdxs);
                tmp1 = sortrows([tmpx1, tmpy1], [1 2]);
                tmp2 = sortrows([tmpx2, tmpy2], [1 2]);
                dnaABoxBound = [tmp1(:, 1) tmp1(:, 2) tmp2(:, 2)];
                
                clear tmp subs tfs vals tmpx1 tmpx2 tmpy1 tmpy2 tmp1 tmp2;
            end
            
            time = ensemble.stateData.time / 3600;
            replisome = [
                ensemble.stateData.values(ensemble.getPropertyIndices('helicase1'), :, :, :)
                ensemble.stateData.values(ensemble.getPropertyIndices('helicase2'), :, :, :)
                ensemble.stateData.values(ensemble.getPropertyIndices('leadingPol1'), :, :, :)
                ensemble.stateData.values(ensemble.getPropertyIndices('leadingPol2'), :, :, :)
                ensemble.stateData.values(ensemble.getPropertyIndices('laggingPol1'), :, :, :)
                ensemble.stateData.values(ensemble.getPropertyIndices('laggingPol2'), :, :, :)
                ];
            
            if nargin < 6 || isempty(posSpan)
                posSpan = [1 c.sequenceLen];
                set(axesHandle, 'YTick', [1 c.terCPosition c.sequenceLen]);
                set(axesHandle, 'YTickLabel', {'terC', 'oriC', 'terC'});
            else
                yTick = [1; c.terCPosition; c.sequenceLen];
                yTickLabel = {'oriC'; 'terC'; 'oriC'};
                if posSpan(1) > c.terCPosition
                    yTick = [yTick; posSpan(1)];
                    yTickLabel = [yTickLabel; sprintf('terC+%d', posSpan(1)-c.terCPosition)];
                elseif posSpan(1) < c.terCPosition
                    yTick = [yTick; posSpan(1)];
                    yTickLabel = [yTickLabel; sprintf('terC-%d', c.terCPosition-posSpan(1))];
                end
                if posSpan(2) > c.terCPosition
                    yTick = [yTick; posSpan(2)];
                    yTickLabel = [yTickLabel; sprintf('terC+%d', posSpan(2)-c.terCPosition)];
                elseif posSpan(2) < c.terCPosition
                    yTick = [yTick; posSpan(2)];
                    yTickLabel = [yTickLabel; sprintf('terC-%d', c.terCPosition-posSpan(2))];
                end
                
                [~, order] = unique(yTick);
                set(axesHandle, 'YTick', yTick(order));
                set(axesHandle, 'YTickLabel', yTickLabel(order));
            end
            
            xlim(axesHandle, time([1 end]));
            ylim(axesHandle, posSpan);
            ylabel(axesHandle, 'Position');
            
            %rna polymeraase
            if nargin < 4 || isempty(nRNAPols)
                nRNAPols = size(states.RNAPolymerase.states, 1);
            end
            if nargin < 13
                rnaPolLineWidth = 0.25;
            end
            for i = 1:nRNAPols
                %actively transcribing / specifically bound
                startIdxs = find(cat(3, ...
                    (states.RNAPolymerase.states(i, :, 1) >= rnaPol.activelyTranscribingValue | ...
                    states.RNAPolymerase.states(i, :, 1) == rnaPol.specificallyBoundValue) & ...
                    states.Transcript.boundTranscriptChromosome(i, :, 1) == 1, ...
                    diff((states.RNAPolymerase.states(i, :, :) >= rnaPol.activelyTranscribingValue | ...
                    states.RNAPolymerase.states(i, :, :) == rnaPol.specificallyBoundValue) & ...
                    states.Transcript.boundTranscriptChromosome(i, :, :) == 1, 1, 3) > 0));
                endIdxs = find(cat(3, ...
                    diff((states.RNAPolymerase.states(i, :, :) >= rnaPol.activelyTranscribingValue | ...
                    states.RNAPolymerase.states(i, :, :) == rnaPol.specificallyBoundValue) & ...
                    states.Transcript.boundTranscriptChromosome(i, :, :) == 1, 1, 3) < 0, ...
                    (states.RNAPolymerase.states(i, :, end) >= rnaPol.activelyTranscribingValue | ...
                    states.RNAPolymerase.states(i, :, end) == rnaPol.specificallyBoundValue) & ...
                    states.Transcript.boundTranscriptChromosome(i, :, end) == 1));
                tfs = ...
                    (states.RNAPolymerase.positionStrands(i, 1, startIdxs) >= posSpan(1) & ...
                    states.RNAPolymerase.positionStrands(i, 1, startIdxs) <= posSpan(2)) |  ...
                    (states.RNAPolymerase.positionStrands(i, 1, endIdxs) >= posSpan(1) &  ...
                    states.RNAPolymerase.positionStrands(i, 1, endIdxs) <= posSpan(2)) | ...
                    (states.RNAPolymerase.positionStrands(i, 1, startIdxs) <= posSpan(1) &  ...
                    states.RNAPolymerase.positionStrands(i, 1, endIdxs) >= posSpan(2)) | ...
                    (states.RNAPolymerase.positionStrands(i, 1, endIdxs) <= posSpan(1) &  ...
                    states.RNAPolymerase.positionStrands(i, 1, startIdxs) >= posSpan(2));
                if any(tfs)
                    if nargin >= 14 && showRNAPolDetails
                        tmpIdxs = find(tfs);
                        for j = 1:numel(tmpIdxs)
                            rnaATPolHandle = line((startIdxs(tmpIdxs(j)):endIdxs(tmpIdxs(j))) / 3600, mod(...
                                permute(states.RNAPolymerase.positionStrands(i, 1, startIdxs(tmpIdxs(j)):endIdxs(tmpIdxs(j))), [1 3 2]) ...
                                - c.sequenceLen/2 - 1, c.sequenceLen) + 1, ...
                                'Color', 'b', 'Parent', axesHandle, 'LineWidth', rnaPolLineWidth);
                        end
                    else
                        rnaATPolHandle = line([startIdxs(tfs) endIdxs(tfs)]' / 3600, mod([...
                            permute(states.RNAPolymerase.positionStrands(i, 1, startIdxs(tfs)), [1 3 2])
                            permute(states.RNAPolymerase.positionStrands(i, 1, endIdxs(tfs)), [1 3 2])
                            ] - c.sequenceLen/2 - 1, c.sequenceLen) + 1, ...
                            'Color', 'b', 'Parent', axesHandle, 'LineWidth', rnaPolLineWidth);
                    end
                end
                    
                %non-specifically bound
                if nargin < 5 || nonSpecBoundPols
                    startIdxs = find(states.RNAPolymerase.states(i, :, :) == rnaPol.nonSpecificallyBoundValue & ...
                        states.RNAPolymerase.positionStrands(i, 2, :) <= 2 & ...
                        cat(3, true, diff(states.RNAPolymerase.positionStrands(i, 1, :), 1, 3)));
                    endIdxs = find(states.RNAPolymerase.states(i, :, :) == rnaPol.nonSpecificallyBoundValue & ...
                        states.RNAPolymerase.positionStrands(i, 2, :) <= 2 & ...
                        cat(3, diff(states.RNAPolymerase.positionStrands(i, 1, :), 1, 3), true));
                    tfs = ...
                        (states.RNAPolymerase.positionStrands(i, 1, startIdxs) >= posSpan(1) & ...
                        states.RNAPolymerase.positionStrands(i, 1, startIdxs) <= posSpan(2)) |  ...
                        (states.RNAPolymerase.positionStrands(i, 1, endIdxs) >= posSpan(1) &  ...
                        states.RNAPolymerase.positionStrands(i, 1, endIdxs) <= posSpan(2));
                    if any(tfs)
                        if nargin >= 14 && showRNAPolDetails
                            tmpIdxs = find(tfs);
                            for j = 1:numel(tmpIdxs)
                                rnaNSBPolHandle = line((startIdxs(tmpIdxs(j)):endIdxs(tmpIdxs(j))) / 3600, mod(...
                                    permute(states.RNAPolymerase.positionStrands(i, 1, startIdxs(tmpIdxs(j)):endIdxs(tmpIdxs(j))), [1 3 2]) ...
                                    - c.sequenceLen/2 - 1, c.sequenceLen) + 1, ...
                                    'Color', 'c', 'Parent', axesHandle, 'LineWidth', rnaPolLineWidth);
                            end
                        else
                            rnaNSBPolHandle = line([startIdxs(tfs) endIdxs(tfs)]' / 3600, mod([...
                                permute(states.RNAPolymerase.positionStrands(i, 1, startIdxs(tfs)), [1 3 2])
                                permute(states.RNAPolymerase.positionStrands(i, 1, endIdxs(tfs)), [1 3 2])
                                ] - c.sequenceLen/2 - 1, c.sequenceLen) + 1, ...
                                'Color', 'c', 'Parent', axesHandle, 'LineWidth', rnaPolLineWidth);
                        end
                    end
                end
            end
            
            %DnaA
            dnaABoxPos = mod(repInit.dnaABoxStartPositions(repInit.dnaABoxIndexs_R12345) - c.sequenceLen/2 - 1, c.sequenceLen) + 1;
            dnaABoxMaxOccupancy = [7; 7; 7; 7; 1];
            dnaA = [
                ensemble.stateData.values(ensemble.getPropertyIndices('dnaA_box1'), :, :, :)
                ensemble.stateData.values(ensemble.getPropertyIndices('dnaA_box2'), :, :, :)
                ensemble.stateData.values(ensemble.getPropertyIndices('dnaA_box3'), :, :, :)
                ensemble.stateData.values(ensemble.getPropertyIndices('dnaA_box4'), :, :, :)
                ensemble.stateData.values(ensemble.getPropertyIndices('dnaA_box5'), :, :, :)
                ];
            
            if nargin >= 7 && showSingleDnaALine
                dnaABoxPos = mean(dnaABoxPos);
                dnaABoxMaxOccupancy = sum(dnaABoxMaxOccupancy);
                dnaA = sum(dnaA, 1);
            end
            
            % TODO: Add an argument to the function definition to enable
            % specifying whether or not you want to see ONLY the oriC DnaA
            % complex or all binding events
            for i = 1:numel(dnaABoxPos)
                tmp = find(diff(dnaA(i, :, :)));
                for j = 1:numel(tmp)-1
                    x = [tmp(j)+1; tmp(j+1)+1] / 3600;
                    y = dnaABoxPos(i);
                    dnaAHandle = patch([x(1) x(2) x(2) x(1)]', [y+eps y+eps y-eps y-eps]', 1, ...
                        'Parent', axesHandle, 'LineWidth', 2, 'FaceColor', 'none', 'EdgeColor', [
                        1
                        1 - log(1 + dnaA(i, :, tmp(j)+1) + eps) / log(1+7)
                        1 - log(1 + dnaA(i, :, tmp(j)+1) + eps) / log(1+7)]);
                end
            end
            
%             dnaABoxBound = dnaABoxBound(~ismember(dnaABoxBound(:, 1), repInit.dnaABoxIndexs_R12345), :);
%             dnaAHandle = line(reshape([dnaABoxBound(:, 2:3)/3600 nan(size(dnaABoxBound, 1), 1)]', [], 1), ...
%                 reshape([repInit.dnaABoxStartPositions(dnaABoxBound(:, [1 1]))  nan(size(dnaABoxBound, 1), 1)]', [], 1), ...
%                 'Parent', axesHandle, 'LineWidth', rnaPolLineWidth, 'Color', [1 log(1+1)/log(1+7) log(1+1)/log(1+7)]);
            
            %DNA polymerase
            replisome(1:2:end, :) = replisome(1:2:end, :) - c.sequenceLen/2;
            replisome(2:2:end, :) = replisome(2:2:end, :) + c.sequenceLen/2;
            
            idxs = [
                find(any(~isnan(replisome), 1), 1, 'first')
                find(any(~isnan(replisome), 1), 1, 'last')
                ];
            green = [0 1 0];
            if nargin < 12
                dnaPolLineWidth = 0.5;
            end
            if ~isempty(idxs)
                line(time(:, idxs(1):idxs(2)), ...
                    permute(replisome(3:4, :, idxs(1):idxs(2), :), [1 3 2 4]), ...
                    'Parent', axesHandle, 'Color', green, 'LineWidth', dnaPolLineWidth);
                
                tmpTime = time(:, idxs(1):idxs(2));
                tmp = permute(replisome(5, :, idxs(1):idxs(2), :), [1 3 2 4]);
                idxs = find(diff(tmp) < 0);
                for i = numel(idxs):-1:1
                    tmp = [tmp(1:idxs(i)) NaN tmp(idxs(i)+1:end)];
                    tmpTime = [tmpTime(1:idxs(i)) NaN tmpTime(idxs(i)+1:end)];
                end
                line(tmpTime, tmp, 'Parent', axesHandle, 'Color', green, 'LineWidth', dnaPolLineWidth);
                
                tmpTime = time(:, idxs(1):idxs(2));
                tmp = permute(replisome(6, :, idxs(1):idxs(2), :), [1 3 2 4]);
                idxs = find(diff(tmp) > 0);
                for i = numel(idxs):-1:1
                    tmp = [tmp(1:idxs(i)) NaN tmp(idxs(i)+1:end)];
                    tmpTime = [tmpTime(1:idxs(i)) NaN tmpTime(idxs(i)+1:end)];
                end
                dnaPolHandle = line(tmpTime, tmp, 'Parent', axesHandle, 'Color', green, 'LineWidth', dnaPolLineWidth);
            end
            
            %legend
            if nargin >= 8 && showLegend
                if exist('rnaNSBPolHandle', 'var')
                    handles = [dnaAHandle(1) dnaPolHandle(1) rnaATPolHandle(1) rnaNSBPolHandle(1)];
                    labels = {'DnaA', 'DNA pol', 'RNA pol - SB/AT', 'RNA pol - NSB'};
                else
                    handles = [dnaAHandle(1) dnaPolHandle(1) rnaATPolHandle(1)];
                    labels = {'DnaA Complex', 'DNA pol', 'RNA pol'};
                end
                
                legendHandle = legend(handles, labels, 'FontSize', 7, 'Location', 'NorthWest');
%                 ChromosomeSpaceTimePlot.formatSpaceTimeLegend(axesHandle);
            else
                legendHandle = 0;
            end
        end
        
        function formatSpaceTimeLegend(axesHandle)
            figChildren = get(get(axesHandle, 'parent'), 'children');
            legendHandle = figChildren(1);
            set(legendHandle, 'Location', 'NorthWest');
            legendChildren = get(legendHandle, 'children');
            set(legendChildren(1:3:end), 'visible', 'off')
            set(legendChildren(2:3:end), 'visible', 'on', 'LineWidth', 4)
            set(legendChildren(2:3:end), 'XData', [0.085 0.2])
            for i = 3:3:numel(legendChildren)
                pos = get(legendChildren(i), 'Position');
                set(legendChildren(i), 'Position', [0.25 pos(2) 0])
            end
            set(legendHandle, 'units', 'pixels');
            legendPos = get(legendHandle, 'position');
            set(legendHandle, 'position', [[1.4 1.42] .* legendPos(1:2) 0.57*legendPos(3) legendPos(4)]);
        end
        
        function [densityBEMatrix, maxTime, maxSpace] = getProteinBindingDensityMatrix(...
                sim, states, RNAPSparseMat, chProteins, whatToMake, width, height)
            import edu.stanford.covert.cell.sim.analysis.ChromosomeSpaceTimePlot;
            
            if exist('whatToMake', 'var') && ~isempty(whatToMake)
                [~, inds] = ismember(whatToMake, chProteins(:, 5));
                chProteins = chProteins(inds(inds ~= 0), :);
            end
            
            densityBEMatrix = zeros(height, 7);
            
            if exist('RNAPSparseMat', 'var') && ~isempty(RNAPSparseMat)
                pol = sim.state('RNAPolymerase');
                subs = find(RNAPSparseMat >= pol.activelyTranscribingValue);
                subs(subs(:, 2) >= 3, :) = [];
                
                maxTime = size(RNAPSparseMat, 3);
                maxSpace = size(RNAPSparseMat, 1);
                ringName = 'RNAPolymerase';
                densityMatrix = ChromosomeSpaceTimePlot.calcChromBoundDensityMatrix(subs, width, height, maxTime, maxSpace);
                densityBEMatrix(:, ChromosomeSpaceTimePlot.rings.(ringName)) = ...
                    densityBEMatrix(:, ChromosomeSpaceTimePlot.rings.(ringName)) + (sum(densityMatrix, 2) / size(densityMatrix, 2));
            end
            
            for i = 1:size(chProteins, 1)
                chromProp = chProteins{i, 1};
                wIDs = chProteins{i, 2};
                strands = chProteins{i, 3};
                name = chProteins{i, 5};
                
                if strcmp(chromProp, 'complexBoundSites')
                    mod = sim.state('ProteinComplex');
                elseif strcmp(chromProp, 'monomerBoundSites')
                    mod = sim.state('ProteinMonomer');
                end
                [~, idxs] = ismember(wIDs, mod.wholeCellModelIDs(mod.matureIndexs));
                
                [subs, vals] = find(states.Chromosome.(chromProp));
                subs = subs(ismember(vals, idxs), :);
                subs = subs(ismember(subs(:, 2), strands), :);
                
                maxTime = size(states.Chromosome.complexBoundSites, 3);
                maxSpace = size(states.Chromosome.complexBoundSites, 1);
                densityMatrix = ChromosomeSpaceTimePlot.calcChromBoundDensityMatrix(subs, width, height, maxTime, maxSpace);
                
                if ismember(wIDs{1}, ChromosomeSpaceTimePlot.regulationWIDs)
                    ringName = 'Regulation';
                else
                    ringName = strrep(name, ' ', '');
                end
                
                if isfield(ChromosomeSpaceTimePlot.rings, ringName)
                    densityBEMatrix(:, ChromosomeSpaceTimePlot.rings.(ringName)) = ...
                        densityBEMatrix(:, ChromosomeSpaceTimePlot.rings.(ringName)) + (sum(densityMatrix, 2) / size(densityMatrix, 2));
                end
            end
        end
    end
end
