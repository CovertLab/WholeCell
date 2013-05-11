% Analyzes families of whole-cell simulations
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Created: 11/1/2012
classdef MultiGenerations
    properties (Constant = true)
        ANCESTRY_COLIDX_CELL = 1
        ANCESTRY_COLIDX_GEN = 2
        ANCESTRY_COLIDX_FAMILY = 3
        ANCESTRY_COLIDX_PARENT = 4
        ANCESTRY_COLIDX_CHILD1 = 5
        ANCESTRY_COLIDX_CHILD2 = 6
    end
    
    methods (Static = true)
        function run(simBatchDir, nGen, nCellFirstGen)
            %% import
            import edu.stanford.covert.cell.sim.analysis.MultiGenerations;
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            
            %%
            if nargin < 3
                simBatchDir = '2012_11_15_18_48_23';
                nGen = 3;
                nCellFirstGen = 8;
            end
            
            outDir = [SimulationDiskUtil.getBaseDir() filesep simBatchDir];
            
            if ~exist(outDir, 'dir')
                mkdir(outDir)
            end
            
            %% load constants
            sim = CachedSimulationObjectUtil.load();
            massState = sim.state('Mass');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            %% load data
            ancestry = MultiGenerations.calcAncestry(nGen, nCellFirstGen);
            relations = MultiGenerations.calcAncestryRelations(nGen, nCellFirstGen);
            
            [~, simData, simStartTimes, divSimTfs, finSimTfs, finSimIdxs, propNames, ...
                dnaBndMons, dnaBndCpxs, otherBndCpxs, dnaBndMonTfs, dnaBndCpxTfs, otherBndCpxTfs] = ...
                MultiGenerations.cacheSimData(simBatchDir);
            simData = simData{1};
            simStartTimes = simStartTimes{1};
            divSimTfs = divSimTfs{1};
            finSimTfs = finSimTfs{1};
            finSimIdxs = finSimIdxs{1};
            dnaBndMons = dnaBndMons{1};
            dnaBndCpxs = dnaBndCpxs{1};
            otherBndCpxs = otherBndCpxs{1};
            dnaBndMonTfs = dnaBndMonTfs{1};
            dnaBndCpxTfs = dnaBndCpxTfs{1};
            otherBndCpxTfs = otherBndCpxTfs{1};
            
            %% population growth
            %growth traces
            for iFamily = 1:nCellFirstGen
                famIdxs = find(finSimTfs & ancestry(:, MultiGenerations.ANCESTRY_COLIDX_FAMILY) == iFamily);
                colors = MultiGenerations.calcRedGreenColors(numel(famIdxs));
                [~, order] = sort(simStartTimes(famIdxs));
                colors(order, :) = colors;
                tStep = 10;
                
                [axesHandle, figHandle] = PlotUtil.newAxesHandle();
                hold(axesHandle, 'on');
                for i = 1:numel(famIdxs)
                    tmp = SimulationEnsemble.load(simBatchDir, {'Time' 'values'; 'MetabolicReaction', 'growth'}, [], [], tStep, 'extract', famIdxs(i));
                    tmpTime = permute(tmp.Time.values, [1 3 2]) / 3600;
                    tmpGrowth = permute(tmp.MetabolicReaction.growth, [1 3 2]) * ...
                        massState.cellInitialDryWeight / (1 - massState.fractionWetWeight) * 3600 * 1e15;
                    plot(axesHandle, simStartTimes(famIdxs(i)) + tmpTime, tmpGrowth, ...
                        'Color', colors(i, :));
                end
                xlabel(axesHandle, 'Time (h)');
                ylabel(axesHandle, 'Growth (fg h^{-1})');
                saveas(figHandle, [outDir filesep 'Growth-Family-' num2str(iFamily) '.pdf']);
                close(figHandle);
            end
            
            %number cells
            nCells = NaN(nGen, 1);
            nSurvive = NaN(nGen, 1);
            for iGen = 0:nGen-1
                nCells(iGen + 1) = sum(ancestry(:, MultiGenerations.ANCESTRY_COLIDX_GEN) == iGen & finSimTfs);
                nSurvive(iGen + 1) = sum(ancestry(:, MultiGenerations.ANCESTRY_COLIDX_GEN) == iGen & divSimTfs);
            end
            fracSurvive = nSurvive ./ nCells;
            
            [axesHandle, figHandle] = PlotUtil.newAxesHandle();
            hold(axesHandle, 'on');
            
            plot(axesHandle, 0:nGen-1, nCells, 'b.')
            xlabel(axesHandle, 'Generation');
            ylabel(axesHandle, 'Population');
            xlim(axesHandle, [-0.5 nGen-0.5]);
            ylim(axesHandle, [0 75]);
            
            expFunc = @(a, b, x) a * exp(x / b);
            [f, gof] = fit((0:nGen-1)', nCells, ...
                fittype(expFunc),  ...
                fitoptions('Method', 'NonlinearLeastSquares', 'lower', [0.5 * nCellFirstGen  0.5], 'upper', [2.0 * nCellFirstGen  3], 'StartPoint', [nCellFirstGen 1/log(2)]));
            plot(axesHandle, (0:0.1:nGen-1)', expFunc(f.a, f.b, (0:0.1:nGen-1)'), 'Color', 'r');
            
            saveas(figHandle, [outDir filesep 'PopulationSize.pdf']);
            close(figHandle);
            
            %% survival
            %all families
            [axesHandle, figHandle] = PlotUtil.newAxesHandle();
            hold(axesHandle, 'on');
            
            errorbar(0:nGen-1, fracSurvive * 100, sqrt(fracSurvive .* (1 - fracSurvive)) ./ sqrt(nCells) * 100, 'b.', 'Parent', axesHandle);
            xlabel(axesHandle, 'Generation');
            ylabel(axesHandle, 'Survival (%)');
            xlim(axesHandle, [-0.5 nGen-0.5]);
            ylim(axesHandle, [0 100])
            
            [f, gof] = fit((0:nGen-1)', fracSurvive * 100, ...
                fittype('poly1'), ...
                fitoptions('Method', 'LinearLeastSquares', 'Weights', sqrt(fracSurvive .* (1 - fracSurvive)) ./ sqrt(nCells) * 100));
            plot(axesHandle, (0:0.1:nGen-1)', (0:0.1:nGen-1)' * f.p1 + f.p2, 'Color', 'r');
            
            saveas(figHandle, [outDir filesep 'Survival.pdf']);
            close(figHandle);
            
            %left/right
            nCellsLR = NaN(nGen, 2);
            nSurviveLR = NaN(nGen, 2);
            for iGen = 0:nGen-1
                nCellsLR(iGen + 1, 1) = sum(ancestry(1:2:end, MultiGenerations.ANCESTRY_COLIDX_GEN) == iGen & finSimTfs(1:2:end, 1));
                nCellsLR(iGen + 1, 2) = sum(ancestry(2:2:end, MultiGenerations.ANCESTRY_COLIDX_GEN) == iGen & finSimTfs(2:2:end, 1));
                nSurviveLR(iGen + 1, 1) = sum(ancestry(1:2:end, MultiGenerations.ANCESTRY_COLIDX_GEN) == iGen & divSimTfs(1:2:end, 1));
                nSurviveLR(iGen + 1, 2) = sum(ancestry(2:2:end, MultiGenerations.ANCESTRY_COLIDX_GEN) == iGen & divSimTfs(2:2:end, 1));
            end
            fracSurviveLR = nSurviveLR ./ nCellsLR;
            
            [axesHandle, figHandle] = PlotUtil.newAxesHandle();
            
            h = plot(axesHandle, 1:nGen-1, fracSurviveLR(2:end, :) * 100);
            legend(h, {'Left', 'Right'}, 'Location', 'NorthEastOutside');
            xlabel(axesHandle, 'Generation');
            ylabel(axesHandle, 'Survival (%)');
            xlim(axesHandle, [0.5 nGen-0.5]);
            ylim(axesHandle, [0 100])
            
            saveas(figHandle, [outDir filesep 'Survival-LeftRight.pdf']);
            close(figHandle);
            
            %by family
            [axesHandle, figHandle] = PlotUtil.newAxesHandle();
            hold(axesHandle, 'on');
            
            colors = MultiGenerations.calcRedGreenColors(nCellFirstGen);
            labels = cell(nCellFirstGen, 1);
            h = zeros(nCellFirstGen, 1);
            for i = 1:nCellFirstGen
                tfs = finSimTfs & ancestry(:, MultiGenerations.ANCESTRY_COLIDX_FAMILY) == i;
                fracSurvive = zeros(nGen, 1);
                for iGen = 0:nGen-1
                    fracSurvive(iGen + 1) = ...
                        sum(ancestry(divSimTfs & tfs, MultiGenerations.ANCESTRY_COLIDX_GEN) == iGen) / ...
                        sum(ancestry(tfs, MultiGenerations.ANCESTRY_COLIDX_GEN) == iGen);
                end
                h(i) = plot(axesHandle, 0:nGen-1, fracSurvive * 100, 'Color', colors(i, :));
                labels{i} = num2str(i);
            end
            
            legend(h, labels, 'Location', 'NorthEastOutside');
            xlabel(axesHandle, 'Generation');
            ylabel(axesHandle, 'Survival (%)');
            xlim(axesHandle, [-0.5 nGen-0.5]);
            ylim(axesHandle, [0 100])
            
            saveas(figHandle, [outDir filesep 'Survival-ByFamily.pdf']);
            close(figHandle);
            
            %% left/right chromosome bias
            %- bound protein
            %- methylation
            %- superhelical density
            
            gen_dnaBndMons = NaN(nGen, 2, sum(dnaBndMonTfs));
            gen_dnaBndCpxs = NaN(nGen, 2, sum(dnaBndCpxTfs));
            gen_otherBndCpxs = NaN(nGen, 2, sum(otherBndCpxTfs));
            
            for iGen = 0:nGen - 1
                lefts  = find( isodd(finSimIdxs) & ancestry(finSimIdxs, MultiGenerations.ANCESTRY_COLIDX_GEN) == iGen);
                rights = find(~isodd(finSimIdxs) & ancestry(finSimIdxs, MultiGenerations.ANCESTRY_COLIDX_GEN) == iGen);
                
                gen_dnaBndMons(iGen+1, 1, :) = mean(dnaBndMons(lefts, :), 1);
                gen_dnaBndMons(iGen+1, 2, :) = mean(dnaBndMons(rights, :), 1);
                
                gen_dnaBndCpxs(iGen+1, 1, :) = mean(dnaBndCpxs(lefts, :), 1);
                gen_dnaBndCpxs(iGen+1, 2, :) = mean(dnaBndCpxs(rights, :), 1);
                
                gen_otherBndCpxs(iGen+1, 1, :) = mean(otherBndCpxs(lefts, :), 1);
                gen_otherBndCpxs(iGen+1, 2, :) = mean(otherBndCpxs(rights, :), 1);
            end
            
            [axesHandle, figHandle] = PlotUtil.newAxesHandle();
            hold(axesHandle, 'on');
            h = plot(axesHandle, (1:nGen-1)', permute(gen_dnaBndMons(2:end, 1, :), [1 3 2]));
            g = plot(axesHandle, (1:nGen-1)', permute(gen_dnaBndMons(2:end, 2, :), [1 3 2]));
            colors = MultiGenerations.calcRedGreenColors(size(gen_dnaBndMons, 3));
            for i = 1:size(gen_dnaBndMons, 3)
                set(h(i), 'LineStyle', '-', 'Color', colors(i, :))
                set(g(i), 'LineStyle', ':', 'Color', colors(i, :))
            end
            legend(h, pm.wholeCellModelIDs(pm.boundIndexs(dnaBndMonTfs)), 'Location', 'NorthEastOutside', 'Interpreter', 'none');
            xlim(axesHandle, [0.5 nGen-0.5])
            xlabel(axesHandle, 'Generation');
            ylabel(axesHandle, 'Protein monomers')
            saveas(figHandle, [outDir filesep 'DNABoundMonomerLeftRightBias.pdf']);
            close(figHandle);
            
            [axesHandle, figHandle] = PlotUtil.newAxesHandle();
            hold(axesHandle, 'on');
            h = plot(axesHandle, (1:nGen-1)', permute(gen_dnaBndCpxs(2:end, 1, :), [1 3 2]));
            g = plot(axesHandle, (1:nGen-1)', permute(gen_dnaBndCpxs(2:end, 2, :), [1 3 2]));
            colors = MultiGenerations.calcRedGreenColors(size(gen_dnaBndCpxs, 3));
            for i = 1:size(gen_dnaBndCpxs, 3)
                set(h(i), 'LineStyle', '-', 'Color', colors(i, :))
                set(g(i), 'LineStyle', ':', 'Color', colors(i, :))
            end
            legend(h, pc.wholeCellModelIDs(pc.boundIndexs(dnaBndCpxTfs)), 'Location', 'NorthEastOutside', 'Interpreter', 'none');
            xlim(axesHandle, [0.5 nGen-0.5])
            xlabel(axesHandle, 'Generation');
            ylabel(axesHandle, 'Protein complexes')
            saveas(figHandle, [outDir filesep 'DNABoundComplexLeftRightBias.pdf']);
            close(figHandle);
            
            [axesHandle, figHandle] = PlotUtil.newAxesHandle();
            hold(axesHandle, 'on');
            h = plot(axesHandle, (1:nGen-1)', permute(gen_otherBndCpxs(2:end, 1, :), [1 3 2]));
            g = plot(axesHandle, (1:nGen-1)', permute(gen_otherBndCpxs(2:end, 2, :), [1 3 2]));
            colors = MultiGenerations.calcRedGreenColors(size(gen_otherBndCpxs, 3));
            for i = 1:size(gen_otherBndCpxs, 3)
                set(h(i), 'LineStyle', '-', 'Color', colors(i, :))
                set(g(i), 'LineStyle', ':', 'Color', colors(i, :))
            end
            legend(h, pc.wholeCellModelIDs(pc.boundIndexs(otherBndCpxTfs)), 'Location', 'NorthEastOutside', 'Interpreter', 'none');
            xlim(axesHandle, [0.5 nGen-0.5])
            xlabel(axesHandle, 'Generation');
            ylabel(axesHandle, 'Protein complexes')
            saveas(figHandle, [outDir filesep 'OtherBoundComplexLeftRightBias.pdf']);
            close(figHandle);
            
            %% growth, mass, cell cycle phase durations
            for iProp = 1:size(simData, 2)
                %all families
                [axesHandle, figHandle] = PlotUtil.newAxesHandle();
                boxplot(axesHandle, simData(finSimTfs, iProp), ancestry(finSimTfs, MultiGenerations.ANCESTRY_COLIDX_GEN));
                title(axesHandle, [propNames{iProp, 1} ' - All Cells'])
                xlabel(axesHandle, 'Generation');
                ylabel(axesHandle, propNames{iProp, 2});
                saveas(figHandle, [outDir filesep propNames{iProp, 1} 'BoxPlot-AllCells.pdf']);
                close(figHandle);
                
                [axesHandle, figHandle] = PlotUtil.newAxesHandle();
                boxplot(axesHandle, simData(divSimTfs, iProp), ancestry(divSimTfs, MultiGenerations.ANCESTRY_COLIDX_GEN));
                title(axesHandle, [propNames{iProp, 1} ' - Successful Divisions'])
                xlabel(axesHandle, 'Generation');
                ylabel(axesHandle, propNames{iProp, 2});
                saveas(figHandle, [outDir filesep propNames{iProp, 1} 'BoxPlot-SuccessfulCells.pdf']);
                close(figHandle);
                
                [axesHandle, figHandle] = PlotUtil.newAxesHandle();
                hold(axesHandle, 'on');
                h = zeros(nGen, 1);
                x = simData(finSimTfs, iProp);
                edges = linspace(min(x), max(x), 20);
                labels = cell(nGen, 1);
                colors = MultiGenerations.calcRedGreenColors(nGen);
                for iGen = 0:nGen-1
                    tfs = finSimTfs & ancestry(:, MultiGenerations.ANCESTRY_COLIDX_GEN) == iGen;
                    cnts = histc(simData(tfs, iProp), edges);
                    cnts = cnts / sum(tfs);
                    h(iGen + 1) = plot(axesHandle, edges, cnts, 'Color', colors(iGen + 1, :));
                    labels{iGen + 1} = num2str(iGen);
                end
                title(axesHandle, [propNames{iProp, 1} ' - All Cells'])
                xlabel(axesHandle, propNames{iProp, 2});
                ylabel(axesHandle, 'Frequency');
                legend(h, labels);
                saveas(figHandle, [outDir filesep propNames{iProp, 1} 'Histogram-AllCells.pdf']);
                close(figHandle);
                
                %left/right bias
                [axesHandle, figHandle] = PlotUtil.newAxesHandle();
                lefts = false(size(finSimTfs));
                lefts(1:2:end, 1) = true;
                yData = NaN(nGen, 1);
                for iGen = 0:nGen-1
                    yData(iGen + 1, 1) = nanmean(simData( lefts & finSimTfs & ancestry(:, MultiGenerations.ANCESTRY_COLIDX_GEN) == iGen, iProp));
                    yData(iGen + 1, 2) = nanmean(simData(~lefts & finSimTfs & ancestry(:, MultiGenerations.ANCESTRY_COLIDX_GEN) == iGen, iProp));
                end
                h = plot(axesHandle, 1:nGen-1, yData(2:end, :));
                legend(h, {'Left', 'Right'}, 'Location', 'NorthEastOutside');
                title(axesHandle, [propNames{iProp, 1}])
                xlabel(axesHandle, 'Generation');
                ylabel(axesHandle, propNames{iProp, 2});
                xlim(axesHandle, [0.5 nGen-0.5]);
                saveas(figHandle, [outDir filesep propNames{iProp, 1} 'MeanLeftRight.pdf']);
                close(figHandle);
                
                %by family
                [axesHandle, figHandle] = PlotUtil.newAxesHandle();
                hold(axesHandle, 'on');
                colors = MultiGenerations.calcRedGreenColors(nCellFirstGen);
                labels = cell(nCellFirstGen, 1);
                h = zeros(nCellFirstGen, 1);
                for i = 1:nCellFirstGen
                    tfs = finSimTfs & ancestry(:, MultiGenerations.ANCESTRY_COLIDX_FAMILY) == i;
                    h(i) = plot(axesHandle, ancestry(tfs, MultiGenerations.ANCESTRY_COLIDX_GEN), ...
                        simData(tfs, iProp), '.', 'Color', colors(i, :));
                    
                    labels{i} = num2str(i);
                end
                legend(h, labels, 'Location', 'NorthEastOutside');
                title(axesHandle, [propNames{iProp, 1} ' - All Cells'])
                xlabel(axesHandle, 'Generation');
                ylabel(axesHandle, propNames{iProp, 2});
                xlim(axesHandle, [-0.5 nGen-0.5]);
                saveas(figHandle, [outDir filesep propNames{iProp, 1} 'Distribution.pdf']);
                close(figHandle);
                
                [axesHandle, figHandle] = PlotUtil.newAxesHandle();
                hold(axesHandle, 'on');
                colors = MultiGenerations.calcRedGreenColors(nCellFirstGen);
                labels = cell(nCellFirstGen, 1);
                h = zeros(nCellFirstGen, 1);
                for i = 1:nCellFirstGen
                    tfs = finSimTfs & ancestry(:, MultiGenerations.ANCESTRY_COLIDX_FAMILY) == i;
                    gens = ancestry(tfs, MultiGenerations.ANCESTRY_COLIDX_GEN);
                    x = simData(tfs, iProp);
                    avgs = zeros(nGen, 1);
                    for iGen = 0:nGen-1
                        avgs(iGen + 1) = nanmean(x(gens == iGen));
                    end
                    h(i) = plot(axesHandle, 0:nGen-1, avgs, 'Color', colors(i, :));
                    
                    labels{i} = num2str(i);
                end
                legend(h, labels, 'Location', 'NorthEastOutside');
                title(axesHandle, [propNames{iProp, 1} ' - All Cells'])
                xlabel(axesHandle, 'Generation');
                ylabel(axesHandle, propNames{iProp, 2});
                xlim(axesHandle, [-0.25 nGen-0.75]);
                saveas(figHandle, [outDir filesep propNames{iProp, 1} 'MeanByLineage-AllCells.pdf']);
                close(figHandle);
                
                [axesHandle, figHandle] = PlotUtil.newAxesHandle();
                hold(axesHandle, 'on');
                colors = MultiGenerations.calcRedGreenColors(nCellFirstGen);
                labels = cell(nCellFirstGen, 1);
                h = zeros(nCellFirstGen, 1);
                for i = 1:nCellFirstGen
                    tfs = divSimTfs & ancestry(:, MultiGenerations.ANCESTRY_COLIDX_FAMILY) == i;
                    gens = ancestry(tfs, MultiGenerations.ANCESTRY_COLIDX_GEN);
                    x = simData(tfs, iProp);
                    avgs = zeros(nGen, 1);
                    for iGen = 0:nGen-1
                        avgs(iGen + 1) = nanmean(x(gens == iGen));
                    end
                    h(i) = plot(axesHandle, 0:nGen-1, avgs, 'Color', colors(i, :));
                    
                    labels{i} = num2str(i);
                end
                legend(h, labels, 'Location', 'NorthEastOutside');
                title(axesHandle, [propNames{iProp, 1} ' - Successful Divisions'])
                xlabel(axesHandle, 'Generation');
                ylabel(axesHandle, propNames{iProp, 2});
                xlim(axesHandle, [-0.25 nGen-0.75]);
                saveas(figHandle, [outDir filesep propNames{iProp, 1} 'MeanByLineage-SuccessfulCells.pdf']);
                close(figHandle);
            end
            
            %% mother-daughter correlations within cell cycle
            deltaGenerations = relations(finSimTfs, finSimTfs, 1);
            removal = relations(finSimTfs, finSimTfs, 2);
            
            [idxs1, idxs2] = find(deltaGenerations == 1 & removal == 0);
            
            ensemble = SimulationEnsemble(simBatchDir, cell(0, 1), [1 2], find(finSimTfs));
            simEndTimes = ensemble.stateData.simulationEndTimes;
            states = SimulationEnsemble.load(simBatchDir, {'MetabolicReaction' 'growth'}, [], [], 1, 'extract', find(finSimTfs));
            growth = permute(states.MetabolicReaction.growth, [4 3 1 2]);
            growth(growth == 0) = NaN;
            growth0 = growth(sub2ind(size(growth), (1:numel(simEndTimes))', simEndTimes));
                       
            r = zeros(1, size(growth, 2));
            n = zeros(1, size(growth, 2));
            for i = 1:size(growth, 2)
                tfs = ~isnan(growth(idxs2, i));
                r(i) = corr(growth0(idxs1(tfs)), growth(idxs2(tfs), i), 'type', 'Pearson');
                n(i) = sum(simEndTimes >= i);
            end
            
            [~, figHandle] = PlotUtil.newAxesHandle();
            clf(figHandle);
            
            axesHandle = subplot(2, 1, 1);
            plot(axesHandle, (1:size(growth, 2))/3600, r);
            xlabel(axesHandle, 'Time (h)');
            ylabel(axesHandle, 'Corr');
            
            axesHandle = subplot(2, 1, 2);
            plot(axesHandle, (1:size(growth, 2))/3600, n);
            xlabel(axesHandle, 'Time (h)');
            ylabel(axesHandle, 'No. cells');

            saveas(figHandle, [outDir filesep 'GrowthCorrelationDecay.pdf']);
            close(figHandle);
            
            %% descendents
            deltaGenerations = relations(:, :, 1);
            removal = relations(:, :, 2);
            
            colors = MultiGenerations.calcRedGreenColors(nGen - 1);
            for iProp = 1:size(propNames, 1)
                %all relationships
                [~, figHandle] = PlotUtil.newAxesHandle();
                clf(figHandle);
                
                axesHandle = subplot(2, 1, 1);
                hold(axesHandle, 'on');
                
                r = zeros(nGen - 1, 1);
                h = zeros(nGen - 1, 1);
                labels = cell(nGen - 1, 1);
                for iGen = 1:nGen-1
                    [idxs1, idxs2] = find(deltaGenerations == iGen & removal == 0);
                    tfs = finSimTfs(idxs1) & finSimTfs(idxs2);
                    idxs1 = idxs1(tfs, :);
                    idxs2 = idxs2(tfs, :);
                    
                    r(iGen) = corr(simData(idxs1, iProp), simData(idxs2, iProp), 'type', 'Pearson');
                    h(iGen) = plot(simData(idxs1, iProp), simData(idxs2, iProp), '.', 'Color', colors(iGen, :));
                    labels{iGen} = sprintf('{\\Delta}g=%d (n=%d)', iGen, numel(idxs1));
                end
                
                legend(h, labels, 'Location', 'NorthEastOutside');
                title(axesHandle, ['All descendents - ' propNames{iProp, 1}])
                xlabel(axesHandle, 'Parent');
                ylabel(axesHandle, 'Child');
                
                axesHandle = subplot(2, 1, 2);
                plot(axesHandle, 1:nGen-1, r);
                line([0.5 nGen-0.5], [0 0], 'Color', 0.25 * [1 1 1], 'Parent', axesHandle);
                xlabel(axesHandle, '{\Delta}Generation');
                ylabel(axesHandle, 'Corr');
                xlim(axesHandle, [0.5 nGen-0.5]);
                ylim(axesHandle, [-1 1]);
                set(axesHandle, 'XTick', 1:nGen-1);
                box(axesHandle, 'off');
                
                saveas(figHandle, [outDir filesep 'Descendents-all-' propNames{iProp, 1} '.pdf']);
                close(figHandle);
                
                %relationships from founders
                [~, figHandle] = PlotUtil.newAxesHandle();
                clf(figHandle);
                
                axesHandle = subplot(2, 1, 1);
                hold(axesHandle, 'on');
                
                r = zeros(nGen - 1, 1);
                h = zeros(nGen - 1, 1);
                labels = cell(nGen - 1, 1);
                for iGen = 1:nGen-1
                    [idxs1, idxs2] = find(deltaGenerations == iGen & removal == 0);
                    tfs = finSimTfs(idxs1) & finSimTfs(idxs2);
                    tfs = tfs & ancestry(idxs1, MultiGenerations.ANCESTRY_COLIDX_GEN) == 0;
                    idxs1 = idxs1(tfs, :);
                    idxs2 = idxs2(tfs, :);
                    
                    r(iGen) = corr(simData(idxs1, iProp), simData(idxs2, iProp), 'type', 'Pearson');
                    h(iGen) = plot(simData(idxs1, iProp), simData(idxs2, iProp), '.', 'Color', colors(iGen, :));
                    labels{iGen} = sprintf('{\\Delta}g=%d (n=%d)', iGen, numel(idxs1));
                end
                
                legend(h, labels, 'Location', 'NorthEastOutside');
                title(axesHandle, ['Generations from founder - ' propNames{iProp, 1}])
                xlabel(axesHandle, 'Parent');
                ylabel(axesHandle, 'Child');
                
                axesHandle = subplot(2, 1, 2);
                plot(axesHandle, 1:nGen-1, r);
                line([0.5 nGen-0.5], [0 0], 'Color', 0.25 * [1 1 1], 'Parent', axesHandle);
                xlabel(axesHandle, '{\Delta}Generation');
                ylabel(axesHandle, 'Corr');
                xlim(axesHandle, [0.5 nGen-0.5]);
                ylim(axesHandle, [-1 1]);
                set(axesHandle, 'XTick', 1:nGen-1);
                box(axesHandle, 'off');
                
                saveas(figHandle, [outDir filesep 'Descendents-FromFounders-' propNames{iProp, 1} '.pdf']);
                close(figHandle);
                
                %relationships from founders
                [~, figHandle] = PlotUtil.newAxesHandle();
                clf(figHandle);
                
                axesHandle = subplot(2, 1, 1);
                hold(axesHandle, 'on');
                
                r = zeros(nGen - 1, 1);
                h = zeros(nGen - 1, 1);
                labels = cell(nGen - 1, 1);
                for iGen = 1:nGen-1
                    [idxs1, idxs2] = find(deltaGenerations == 1 & removal == 0);
                    tfs = finSimTfs(idxs1) & finSimTfs(idxs2);
                    tfs = tfs & ancestry(idxs1, MultiGenerations.ANCESTRY_COLIDX_GEN) == (iGen - 1);
                    idxs1 = idxs1(tfs, :);
                    idxs2 = idxs2(tfs, :);
                    
                    r(iGen) = corr(simData(idxs1, iProp), simData(idxs2, iProp), 'type', 'Pearson');
                    h(iGen) = plot(simData(idxs1, iProp), simData(idxs2, iProp), '.', 'Color', colors(iGen, :));
                    labels{iGen} = sprintf('g_1=%d (n=%d)', iGen, numel(idxs1));
                end
                
                legend(h, labels, 'Location', 'NorthEastOutside');
                title(axesHandle, propNames{iProp, 1})
                xlabel(axesHandle, 'Parent');
                ylabel(axesHandle, 'Child');
                
                axesHandle = subplot(2, 1, 2);
                plot(axesHandle, 1:nGen-1, r);
                line([0.5 nGen-0.5], [0 0], 'Color', 0.25 * [1 1 1], 'Parent', axesHandle);
                xlabel(axesHandle, '{\Delta}Generation');
                ylabel(axesHandle, 'Corr');
                xlim(axesHandle, [0.5 nGen-0.5]);
                ylim(axesHandle, [-1 1]);
                set(axesHandle, 'XTick', 1:nGen-1);
                box(axesHandle, 'off');
                
                saveas(figHandle, [outDir filesep 'Descendents-Children-' propNames{iProp, 1} '.pdf']);
                close(figHandle);
            end
            
            %% siblings, cousins
            deltaGenerations = relations(:, :, 1);
            removal = relations(:, :, 2);
            
            colors = MultiGenerations.calcRedGreenColors(nGen - 1);
            for iProp = 1:size(propNames, 1)
                %all cousins
                [~, figHandle] = PlotUtil.newAxesHandle();
                clf(figHandle);
                
                axesHandle = subplot(2, 1, 1);
                hold(axesHandle, 'on');
                
                r = zeros(nGen - 1, 1);
                h = zeros(nGen - 1, 1);
                labels = cell(nGen - 1, 1);
                for iGen = 1:nGen-1
                    [idxs1, idxs2] = find(deltaGenerations == 0 & removal == iGen);
                    tfs = finSimTfs(idxs1) & finSimTfs(idxs2);
                    idxs1 = idxs1(tfs, :);
                    idxs2 = idxs2(tfs, :);
                    
                    r(iGen) = corr(simData(idxs1, iProp), simData(idxs2, iProp), 'type', 'Pearson');
                    h(iGen) = plot(simData(idxs1, iProp), simData(idxs2, iProp), '.', 'Color', colors(iGen, :));
                    labels{iGen} = sprintf('{\\Delta}r=%d (n=%d)', iGen, numel(idxs1));
                end
                
                legend(h, labels, 'Location', 'NorthEastOutside');
                title(axesHandle, ['All cousins - ' propNames{iProp, 1}])
                xlabel(axesHandle, 'Parent');
                ylabel(axesHandle, 'Child');
                
                axesHandle = subplot(2, 1, 2);
                plot(axesHandle, 1:nGen-1, r);
                line([0.5 nGen-0.5], [0 0], 'Color', 0.25 * [1 1 1], 'Parent', axesHandle);
                xlabel(axesHandle, 'Removal');
                ylabel(axesHandle, 'Corr');
                xlim(axesHandle, [0.5 nGen-0.5]);
                ylim(axesHandle, [-1 1]);
                set(axesHandle, 'XTick', 1:nGen-1);
                box(axesHandle, 'off');
                
                saveas(figHandle, [outDir filesep 'Cousins-all-' propNames{iProp, 1} '.pdf']);
                close(figHandle);
                
                %cousins of last generation
                [~, figHandle] = PlotUtil.newAxesHandle();
                clf(figHandle);
                
                axesHandle = subplot(2, 1, 1);
                hold(axesHandle, 'on');
                
                r = zeros(nGen - 1, 1);
                h = zeros(nGen - 1, 1);
                labels = cell(nGen - 1, 1);
                for iGen = 1:nGen-1
                    [idxs1, idxs2] = find(deltaGenerations == 0 & removal == iGen);
                    tfs = finSimTfs(idxs1) & finSimTfs(idxs2) & ancestry(idxs1, MultiGenerations.ANCESTRY_COLIDX_GEN) == iGen;
                    idxs1 = idxs1(tfs, :);
                    idxs2 = idxs2(tfs, :);
                    
                    r(iGen) = corr(simData(idxs1, iProp), simData(idxs2, iProp), 'type', 'Pearson');
                    h(iGen) = plot(simData(idxs1, iProp), simData(idxs2, iProp), '.', 'Color', colors(iGen, :));
                    labels{iGen} = sprintf('{\\Delta}r=%d (n=%d)', iGen, numel(idxs1));
                end
                
                legend(h, labels, 'Location', 'NorthEastOutside');
                title(axesHandle, ['Cousins of last generation - ' propNames{iProp, 1}])
                xlabel(axesHandle, 'Parent');
                ylabel(axesHandle, 'Child');
                
                axesHandle = subplot(2, 1, 2);
                plot(axesHandle, 1:nGen-1, r);
                line([0.5 nGen-0.5], [0 0], 'Color', 0.25 * [1 1 1], 'Parent', axesHandle);
                xlabel(axesHandle, 'Removal');
                ylabel(axesHandle, 'Corr');
                xlim(axesHandle, [0.5 nGen-0.5]);
                ylim(axesHandle, [-1 1]);
                set(axesHandle, 'XTick', 1:nGen-1);
                box(axesHandle, 'off');
                
                saveas(figHandle, [outDir filesep 'Cousins-LastGeneration-' propNames{iProp, 1} '.pdf']);
                close(figHandle);
                
                %first cousins by generation
                [~, figHandle] = PlotUtil.newAxesHandle();
                clf(figHandle);
                
                axesHandle = subplot(2, 1, 1);
                hold(axesHandle, 'on');
                
                r = NaN(nGen - 1, 1);
                h = NaN(nGen - 1, 1);
                labels = cell(nGen - 1, 1);
                for iGen = 1:nGen-1
                    [idxs1, idxs2] = find(deltaGenerations == 0 & removal == 1);
                    tfs = finSimTfs(idxs1) & finSimTfs(idxs2);
                    tfs = tfs & ancestry(idxs1, MultiGenerations.ANCESTRY_COLIDX_GEN) == iGen;
                    if ~any(tfs)
                        continue;
                    end
                    idxs1 = idxs1(tfs, :);
                    idxs2 = idxs2(tfs, :);
                    
                    r(iGen) = corr(simData(idxs1, iProp), simData(idxs2, iProp), 'type', 'Pearson');
                    h(iGen) = plot(simData(idxs1, iProp), simData(idxs2, iProp), '.', 'Color', colors(iGen, :));
                    labels{iGen} = sprintf('g=%d (n=%d)', iGen, numel(idxs1));
                end
                
                legend(h(~isnan(h)), labels(~isnan(h)), 'Location', 'NorthEastOutside');
                title(axesHandle, ['All First Cousins - ' propNames{iProp, 1}])
                xlabel(axesHandle, 'Parent');
                ylabel(axesHandle, 'Child');
                
                axesHandle = subplot(2, 1, 2);
                plot(axesHandle, 1:nGen-1, r);
                line([0.5 nGen-0.5], [0 0], 'Color', 0.25 * [1 1 1], 'Parent', axesHandle);
                xlabel(axesHandle, '{\Delta}Generation');
                ylabel(axesHandle, 'Corr');
                xlim(axesHandle, [0.5 nGen-0.5]);
                ylim(axesHandle, [-1 1]);
                set(axesHandle, 'XTick', 1:nGen-1);
                box(axesHandle, 'off');
                
                saveas(figHandle, [outDir filesep 'FirstCousins-all-' propNames{iProp, 1} '.pdf']);
                close(figHandle);
            end
            
            %% Varying effective DnaA copy number
            MultiGenerations.runEffectiveDnaACopyNumberAnalysis();
        end
        
        function runEffectiveDnaACopyNumberAnalysis(outDir)
            import edu.stanford.covert.cell.sim.analysis.MultiGenerations;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            [simBatchDirs, simData, ~, ~, finSimTfs, ~, propNames] = MultiGenerations.cacheSimData();
            
            nSimBatch = size(simBatchDirs, 1);
            
            %% create output directory
            if nargin < 1
                outDir = [SimulationDiskUtil.getBaseDir() filesep 'multiGenerations'];
            end
            if ~exist(outDir, 'dir')
                mkdir(outDir)
            end
            
            %% analyze
            colors = MultiGenerations.calcRedGreenColors(nSimBatch);
            for iProp = 1:size(propNames, 1)
                [axesHandle, figHandle] = PlotUtil.newAxesHandle();
                hold(axesHandle, 'on');
                
                h = zeros(nSimBatch, 1);
                for iSimBatch = 1:nSimBatch
                    nGen = simBatchDirs{iSimBatch, 2};
                    nCellFirstGen = simBatchDirs{iSimBatch, 3};
                    ancestry = MultiGenerations.calcAncestry(nGen, nCellFirstGen);
                    
                    avgs = zeros(1, nGen);
                    for iGen = 0:nGen - 1
                        avgs(1, iGen + 1) = nanmean(simData{iSimBatch}(finSimTfs{iSimBatch} & ancestry(:, MultiGenerations.ANCESTRY_COLIDX_GEN) == iGen, iProp));
                    end
                    h(iSimBatch) = plot(axesHandle, 0:nGen-1, avgs, 'Color', colors(iSimBatch, :));
                end
                
                xlim(axesHandle, [-0.5 nGen-0.5]);
                title(axesHandle, ['DnaA Copy Number vs ' propNames{iProp, 1}]);
                xlabel(axesHandle, 'Generation');
                ylabel(axesHandle, propNames{iProp, 2});
                legend(h, cellfun(@num2str, simBatchDirs(:, 4), 'UniformOutput', false), 'Location', 'NorthEastOutside');
                
                saveas(figHandle, [outDir filesep 'DnaACopyNumberVs' propNames{iProp, 1} '.pdf']);
                close(figHandle);
            end
        end
        
        function [simBatchDirs, simData, simStartTimes, divSimTfs, finSimTfs, finSimIdxs, propNames, ...
                dnaBndMons, dnaBndCpxs, otherBndCpxs, dnaBndMonTfs, dnaBndCpxTfs, otherBndCpxTfs] = cacheSimData(simBatchDirs)
            import edu.stanford.covert.cell.sim.analysis.MultiGenerations;
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            if nargin >= 1
                if ~iscell(simBatchDirs)
                    simBatchDirs = {simBatchDirs};
                end
            else
                simBatchDirs = {
                    '2012_11_15_18_48_23'  3 8  3
                    '2012_10_24_00_49_53'  6 8  1
                    };
            end
            
            nSimBatch = size(simBatchDirs, 1);
            
            %% load constants
            sim = CachedSimulationObjectUtil.load();
            
            %% get data
            simData = cell(nSimBatch, 1);
            simStartTimes = cell(nSimBatch, 1);
            divSimTfs = cell(nSimBatch, 1);
            finSimTfs = cell(nSimBatch, 1);
            finSimIdxs = cell(nSimBatch, 1);
            dnaBndMons = cell(nSimBatch, 1);
            dnaBndCpxs = cell(nSimBatch, 1);
            otherBndCpxs = cell(nSimBatch, 1);
            dnaBndMonTfs = cell(nSimBatch, 1);
            dnaBndCpxTfs = cell(nSimBatch, 1);
            otherBndCpxTfs = cell(nSimBatch, 1);
            for iSimBatch = 1:nSimBatch
                simBatchDir = simBatchDirs{iSimBatch, 1};
                if exist([SimulationDiskUtil.getBaseDir() filesep simBatchDir filesep 'MultiGenerationsData.mat'], 'file')
                    tmp = load([SimulationDiskUtil.getBaseDir() filesep simBatchDir filesep 'MultiGenerationsData.mat']);
                else
                    nGen = simBatchDirs{iSimBatch, 2};
                    nCellFirstGen = simBatchDirs{iSimBatch, 3};
                    
                    [tmp_simData, tmp_simStartTimes, tmp_divSimTfs, tmp_finSimTfs, tmp_finSimIdxs, tmp_propNames, ...
                        tmp_dnaBndMons, tmp_dnaBndCpxs, tmp_otherBndCpxs, tmp_dnaBndMonTfs, tmp_dnaBndCpxTfs, tmp_otherBndCpxTfs] = ...
                        MultiGenerations.getSimData(simBatchDir, nGen, nCellFirstGen, sim);
                    
                    tmp = struct;
                    tmp.simData = tmp_simData;
                    tmp.simStartTimes = tmp_simStartTimes;
                    tmp.divSimTfs = tmp_divSimTfs;
                    tmp.finSimTfs = tmp_finSimTfs;
                    tmp.finSimIdxs = tmp_finSimIdxs;
                    tmp.propNames = tmp_propNames;
                    tmp.dnaBndMons = tmp_dnaBndMons;
                    tmp.dnaBndCpxs = tmp_dnaBndCpxs;
                    tmp.otherBndCpxs = tmp_otherBndCpxs;
                    tmp.dnaBndMonTfs = tmp_dnaBndMonTfs;
                    tmp.dnaBndCpxTfs = tmp_dnaBndCpxTfs;
                    tmp.otherBndCpxTfs = tmp_otherBndCpxTfs;
                    save([SimulationDiskUtil.getBaseDir() filesep simBatchDir filesep 'MultiGenerationsData.mat'], '-struct', 'tmp');
                end
                simData{iSimBatch} = tmp.simData;
                simStartTimes{iSimBatch} = tmp.simStartTimes;
                divSimTfs{iSimBatch} = tmp.divSimTfs;
                finSimTfs{iSimBatch} = tmp.finSimTfs;
                finSimIdxs{iSimBatch} = tmp.finSimIdxs;
                                
                dnaBndMons{iSimBatch} = tmp.dnaBndMons;
                dnaBndCpxs{iSimBatch} = tmp.dnaBndCpxs;
                otherBndCpxs{iSimBatch} = tmp.otherBndCpxs;
                dnaBndMonTfs{iSimBatch} = tmp.dnaBndMonTfs;
                dnaBndCpxTfs{iSimBatch} = tmp.dnaBndCpxTfs;
                otherBndCpxTfs{iSimBatch} = tmp.otherBndCpxTfs;
                
                propNames = tmp.propNames;
            end
        end
        
        function [simData, simStartTimes, divSimTfs, finSimTfs, finSimIdxs, propNames, ...
                dnaBndMons, dnaBndCpxs, otherBndCpxs, ...
                dnaBndMonTfs, dnaBndCpxTfs, otherBndCpxTfs] = ...
                getSimData(simBatchDir, nGen, nCellFirstGen, sim)
            import edu.stanford.covert.cell.sim.analysis.MultiGenerations;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SummaryLogger;
            
            ancestry = MultiGenerations.calcAncestry(nGen, nCellFirstGen);
            nCells = size(ancestry, 1);
            
            g = sim.gene;
            c = sim.compartment;
            massState = sim.state('Mass');
            met = sim.state('Metabolite');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            simStats = SummaryLogger.getSimulationStatistics([SimulationDiskUtil.getBaseDir() filesep simBatchDir], 1:nCells);
            simStats(:, SummaryLogger.SIM_STATUS_INDEX_INIT_GROWTH) = ...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_INIT_GROWTH)' * ...
                massState.cellInitialDryWeight / (1 - massState.fractionWetWeight) * 3600 * 1e15;
            simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME) = simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME) / 3600;
            simStats(:, SummaryLogger.SIM_STATUS_INDEX_INIT_REP_TIME) = simStats(:, SummaryLogger.SIM_STATUS_INDEX_INIT_REP_TIME) / 3600;
            simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_TIME) = simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_TIME) / 3600;
            simStats(:, SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME) = simStats(:, SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME) / 3600;
            simStats(:, SummaryLogger.SIM_STATUS_INDEX_MASS_DOUBLING_TIME) = simStats(:, SummaryLogger.SIM_STATUS_INDEX_MASS_DOUBLING_TIME) / 3600;
            simStats(:, SummaryLogger.SIM_STATUS_INDEX_FIN_MASS) = simStats(:, SummaryLogger.SIM_STATUS_INDEX_FIN_MASS) * 1e15;
            simStats(:, SummaryLogger.SIM_STATUS_INDEX_INIT_MASS) = simStats(:, SummaryLogger.SIM_STATUS_INDEX_INIT_MASS) * 1e15;
            
            simStartTimes = zeros(size(ancestry, 1), 1);
            for i = 1:size(ancestry, 1)
                if isnan(ancestry(i, MultiGenerations.ANCESTRY_COLIDX_PARENT))
                    simStartTimes(i) = 0;
                else
                    iParent = ancestry(i, MultiGenerations.ANCESTRY_COLIDX_PARENT);
                    simStartTimes(i) = ...
                        simStartTimes(iParent) +  ...
                        simStats(iParent, SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME) + ...
                        1 / 3600;
                end
                
                if ~isnan(ancestry(i, MultiGenerations.ANCESTRY_COLIDX_CHILD1)) &&  ...
                        simStats(i, SummaryLogger.SIM_STATUS_INDEX_STATUS) ~= SummaryLogger.SIM_STATUS_COMPLETED_WITH_DIVISION
                    simStats(ancestry(i, MultiGenerations.ANCESTRY_COLIDX_CHILD1), SummaryLogger.SIM_STATUS_INDEX_STATUS) = SummaryLogger.SIM_STATUS_DIDNT_START;
                    simStats(ancestry(i, MultiGenerations.ANCESTRY_COLIDX_CHILD2), SummaryLogger.SIM_STATUS_INDEX_STATUS) = SummaryLogger.SIM_STATUS_DIDNT_START;
                end
            end
            
            divSimTfs = simStats(:, SummaryLogger.SIM_STATUS_INDEX_STATUS) == SummaryLogger.SIM_STATUS_COMPLETED_WITH_DIVISION;
            nonDivSimTfs = simStats(:, SummaryLogger.SIM_STATUS_INDEX_STATUS) == SummaryLogger.SIM_STATUS_COMPLETED_WITHOUT_DIVISION;
            finSimTfs = divSimTfs | nonDivSimTfs;
            finSimIdxs = find(finSimTfs);
            
            props = {'dnaA_total'; 'dnaA_boxes'; 'rnaPolymerases'; 'ribosomes';
                'amino_acids'; 'ntps'; 'rnas'; 'immatureRnas'; 'matureMonomers'; 'immatureMonomers';
                'matureComplexs'; 'immatureComplexs'; 'atp'; 'adp'; 'amp'};
            ensemble = SimulationEnsemble(simBatchDir, props, [1 2], finSimIdxs);
            
            rnaPolMonIdxs = pm.matureIndexs(any(pc.proteinComplexComposition(g.mRNAIndexs, pc.rnaPolymeraseIndexs(1), :), 3));
            stateNames = {
                'ProteinMonomer'  'counts'  rnaPolMonIdxs   c.cytosolIndexs
                };
            rnaPolSubunitCnts = SimulationEnsemble.load(simBatchDir, stateNames, 1, 1, 1, 'extract', finSimIdxs);
            
            stateNames = {
                'ProteinMonomer'  'counts'                pm.boundIndexs  c.cytosolIndexs
                'ProteinComplex'  'counts'                pc.boundIndexs  c.cytosolIndexs
                'Chromosome'      'damagedBases'          ':'  ':'
                'Chromosome'      'superhelicalDensity'   ':'  ':'
                };
            chromState = SimulationEnsemble.load(simBatchDir, stateNames, 1, 1, 1, 'extract', finSimIdxs);
            
            dnaBndMonTfs = any(any(any(chromState.ProteinMonomer.counts, 2), 3), 4);
            dnaBndCpxTfs = any(any(any(chromState.ProteinComplex.counts, 2), 3), 4) & ~any(any(pc.proteinComplexComposition(g.getIndexs('MG_469'), :, :), 1), 3)';
            otherBndCpxTfs = any(any(any(chromState.ProteinComplex.counts, 2), 3), 4) &  any(any(pc.proteinComplexComposition(g.getIndexs('MG_469'), :, :), 1), 3)';
            
            dnaBndMons = NaN(numel(finSimIdxs), sum(dnaBndMonTfs));
            dnaBndCpxs = NaN(numel(finSimIdxs), sum(dnaBndCpxTfs));
            otherBndCpxs = NaN(numel(finSimIdxs), sum(otherBndCpxTfs));
            methylations = NaN(numel(finSimIdxs), 1);
            sigmas = NaN(numel(finSimIdxs), 1);
            
            for i = 1:numel(finSimIdxs)
                dnaBndMons(i, :) = chromState.ProteinMonomer.counts(dnaBndMonTfs, 1, 1, i);
                dnaBndCpxs(i, :) = chromState.ProteinComplex.counts(dnaBndCpxTfs, 1, 1, i);
                otherBndCpxs(i, :) = chromState.ProteinComplex.counts(otherBndCpxTfs, 1, 1, i);
                methylations(i, 1) = full(sum(sum(chromState.Chromosome.damagedBases(:, :, 1, i) == met.m6ADIndexs, 1), 2));
                sigmas(i, 1) = full(chromState.Chromosome.superhelicalDensity(1, 1, 1, i));
            end
            
            simData = [
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_INIT_GROWTH) ...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_FIN_MASS) ...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME) ...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_TIME)-simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME) ...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME)-simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_TIME) ...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME) ...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_INIT_MASS) ...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_FIN_MASS) ...
                zeros(size(simStats, 1), size(ensemble.stateData.values, 1) + size(rnaPolSubunitCnts.ProteinMonomer.counts, 1) + 2)
                ];
            simData(finSimTfs, 9:end) = [
                permute(ensemble.stateData.values(:, :, 1, :), [4 1 2 3]) ...
                permute(rnaPolSubunitCnts.ProteinMonomer.counts, [4 1 2 3]) ...
                methylations ...
                sigmas
                ];
            propNames = {
                'InitialGrowth'                   'Growth (fg h^{-1})'
                'FinalGrowth'                     'Growth (fg h^{-1})'
                'ReplicationInitiationDuration'   'Time (h)'
                'ReplicationDuration'             'Time (h)'
                'CytokinesisDuration'             'Time (h)'
                'CellCycleLength'                 'Time (h)'
                'InitialMass'                     'Mass (fg)'
                'FinalMass'                       'Mass (fg)'
                'DnaA'                            'Count'
                'BoundDnaA'                       'Count'
                'RNAPolymerases'                  'Count'
                'Ribosomes'                       'Count'
                'AminoAcids'                      'Count'
                'Ntps'                            'Count'
                'Rnas'                            'Count'
                'ImmatureRnas'                    'Count'
                'MatureMonomers'                  'Count'
                'ImmatureMonomers'                'Count'
                'MatureComplexes'                 'Count'
                'ImmatureComplexes'               'Count'
                'Atp'                             'Count'
                'Adp'                             'Count'
                'Amp'                             'Count'
                'RnaPolDelta'                     'Count'
                'RnaPolAlpha'                     'Count'
                'RnaPolBeta1'                     'Count'
                'RnaPolBeta2'                     'Count'
                'Methylations'                    'Count'
                'SuperhelicalDensity'             'dimensionless'
                };
        end
        
        %Returns matrix with five colums:
        %- cell index
        %- generation
        %- parent index
        %- child-1 index
        %- child-2 index
        function ancestry = calcAncestry(nGen, nCellFirstGen)
            import edu.stanford.covert.cell.sim.analysis.MultiGenerations;
            
            nCells = nCellFirstGen * (2^nGen - 1);
            
            ancestry = NaN(nCells, 5);
            
            genOffset = 0;
            iCell = 0;
            for iGen = 0:nGen-1
                for iGenCell = 1:nCellFirstGen * 2^iGen
                    iCell = iCell + 1;
                    
                    if iGen == 0
                        iParent = NaN;
                        iFamily = iGenCell;
                    else
                        iParent = genOffset - nCellFirstGen * 2^(iGen-1) + ceil(iGenCell / 2);
                        iFamily = ancestry(iParent, MultiGenerations.ANCESTRY_COLIDX_FAMILY);
                    end
                    
                    if iGen < nGen-1
                        iChild1 = genOffset + nCellFirstGen * 2^iGen + 2 * iGenCell - 1;
                        iChild2 = genOffset + nCellFirstGen * 2^iGen + 2 * iGenCell;
                    else
                        iChild1 = NaN;
                        iChild2 = NaN;
                    end
                    
                    ancestry(iCell, MultiGenerations.ANCESTRY_COLIDX_CELL) = genOffset + iGenCell;
                    ancestry(iCell, MultiGenerations.ANCESTRY_COLIDX_GEN) = iGen;
                    ancestry(iCell, MultiGenerations.ANCESTRY_COLIDX_FAMILY) = iFamily;
                    ancestry(iCell, MultiGenerations.ANCESTRY_COLIDX_PARENT) = iParent;
                    ancestry(iCell, MultiGenerations.ANCESTRY_COLIDX_CHILD1) = iChild1;
                    ancestry(iCell, MultiGenerations.ANCESTRY_COLIDX_CHILD2) = iChild2;
                end
                genOffset = genOffset + nCellFirstGen * 2^iGen;
            end
        end
        
        function relations = calcAncestryRelations(nGen, nCellFirstGen)
            import edu.stanford.covert.cell.sim.analysis.MultiGenerations;
            
            nCells = 2^nGen - 1;
            
            %intergenerational -- difference in generations
            generations = NaN(nCells, 1);
            genOffset = 0;
            for iGen = 0:nGen-1
                iCell = genOffset + (1:2^iGen);
                genOffset = genOffset + 2^iGen;
                generations(iCell) = iGen;
            end
            deltaGeneration = -generations(:, ones(nCells, 1)) + generations(:, ones(nCells, 1))';
            
            %intragenerational -- number of generations since common ancestor
            ancestry = MultiGenerations.calcAncestry(nGen, 1);
            removal = NaN(nCells, nCells);
            for iCell1 = 1:nCells
                removal(iCell1, iCell1) = 0;
                for iCell2 = iCell1+1:nCells
                    [~, removal(iCell1, iCell2)] = MultiGenerations.calcCommonAncestor(iCell1, iCell2, ancestry);
                    removal(iCell2, iCell1) = removal(iCell1, iCell2);
                end
            end
            
            tmp = cat(3, deltaGeneration, removal);
            relations = NaN(nCellFirstGen * nCells, nCellFirstGen * nCells, 2); %cell-1 x cell-2 x [generation separation,  removal]
            ancestry = MultiGenerations.calcAncestry(nGen, nCellFirstGen);
            for i = 1:nCellFirstGen
                relations(ancestry(:, MultiGenerations.ANCESTRY_COLIDX_FAMILY) == i, ancestry(:, MultiGenerations.ANCESTRY_COLIDX_FAMILY) == i, :) = tmp;
            end
        end
        
        function [iCommonAncestor, removal] = calcCommonAncestor(iCell1, iCell2, ancestry)
            import edu.stanford.covert.cell.sim.analysis.MultiGenerations;
            
            if iCell1 == iCell2
                iCommonAncestor = iCell1;
                removal = 0;
            elseif ancestry(iCell1, MultiGenerations.ANCESTRY_COLIDX_GEN) > ancestry(iCell2, MultiGenerations.ANCESTRY_COLIDX_GEN)
                [iCommonAncestor, removal] = MultiGenerations.calcCommonAncestor(ancestry(iCell1, MultiGenerations.ANCESTRY_COLIDX_PARENT), iCell2, ancestry);
            elseif ancestry(iCell1, MultiGenerations.ANCESTRY_COLIDX_GEN) < ancestry(iCell2, MultiGenerations.ANCESTRY_COLIDX_GEN)
                [iCommonAncestor, removal] = MultiGenerations.calcCommonAncestor(iCell1, ancestry(iCell2, MultiGenerations.ANCESTRY_COLIDX_PARENT), ancestry);
            else
                [iCommonAncestor, removal] = MultiGenerations.calcCommonAncestor(...
                    ancestry(iCell1, MultiGenerations.ANCESTRY_COLIDX_PARENT), ...
                    ancestry(iCell2, MultiGenerations.ANCESTRY_COLIDX_PARENT), ...
                    ancestry);
                removal = removal + 1;
            end
        end
        
        function colors = calcRedGreenColors(n)
            tmp = (0:n-1)';
            r1 = 1;
            r2 = 0;
            g1 = 0;
            g2 = 1;
            b1 = 0;
            b2 = 0;
            colors = [
                r1*2*max(0, 0.5-tmp/(n-1)) + r2*2*max(0, tmp/(n-1)-0.5) ...
                g1*2*max(0, 0.5-tmp/(n-1)) + g2*2*max(0, tmp/(n-1)-0.5) ...
                b1*2*max(0, 0.5-tmp/(n-1)) + b2*2*max(0, tmp/(n-1)-0.5) ...
                ];
        end
    end
end