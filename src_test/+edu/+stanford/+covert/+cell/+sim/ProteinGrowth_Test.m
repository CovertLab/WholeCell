%Protein/Growth medium test
% Test that metabolic growth rate increases exponentially over the cell
% cycle, even considering variation in mRNA and protein expression.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/28/2011
classdef ProteinGrowth_Test < TestCase
    methods
        function this = ProteinGrowth_Test(name)
            this = this@TestCase(name);
        end
    end
    
    methods
        function testExponentialGrowth(this)
            import edu.stanford.covert.util.ConstantUtil;
            
            %% simulate
            iterMax = 31000;
            [initGrowth, initMetabolicEnzymes, growths, cumGrowth, rnaExp, pmExp, pcExp, ...
                matureRNACounts, matureMonCounts, matureCpxCounts] = this.simulateExponentialGrowth(iterMax, 0, 0, 1, 1000);
            
            %% references
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'Metabolism'
                'MacromolecularComplexation'
                });
            g = sim.gene;
            
            time = sim.state('Time');
            mass = sim.state('Mass');
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            %% assertions -- growth
            assertElementsAlmostEqual(initGrowth * exp(log(2) * iterMax / time.cellCycleLength), ...
                growths(end), 'relative', 0.75, 0);
            assertElementsAlmostEqual(initGrowth * time.cellCycleLength / log(2) * exp(log(2) * iterMax / time.cellCycleLength) - 1, ...
                cumGrowth, 'relative', 0.5, 0);
            
            this.calcGrowthLimitingEnzymes(sim, initGrowth, initMetabolicEnzymes);
            
            %% assertions -- RNA
            rnas = matureRNACounts;
            rnas(setdiff(1:end, rna.matureMRNAIndexs)) = ...
                rnas(setdiff(1:end, rna.matureMRNAIndexs)) + ...
                sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * matureCpxCounts;
            assertElementsAlmostEqual(rnas' * rna.molecularWeights(rna.matureIndexs) / ConstantUtil.nAvogadro, ...
                (1 + cumGrowth) * mass.cellInitialDryWeight * mass.dryWeightFractionRNA, ...
                'relative', 75e-2, 0);
            
            %RNA expression matches expectations
            matureRNAExpression = rnaExp;
            matureRNAExpression(setdiff(1:end, rna.matureMRNAIndexs)) = ...
                + matureRNAExpression(setdiff(1:end, rna.matureMRNAIndexs)) ...
                + sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * pcExp;
            assertIn(corr(rna.expression(rna.matureIndexs), matureRNAExpression), [0.95 1]);
            assertIn(180 / pi * acos(rna.expression(rna.matureIndexs)' * matureRNAExpression / ...
                (sqrt(rna.expression(rna.matureIndexs)' * rna.expression(rna.matureIndexs)) * ...
                sqrt(matureRNAExpression' * matureRNAExpression))), ...
                [0 20]);
            
            %mRNA expression matches expectations
            mRNAExp = matureRNAExpression(rna.matureMRNAIndexs);
            assertElementsAlmostEqual(...
                (rna.expression(rna.matureIndexs(rna.matureMRNAIndexs))' * rna.molecularWeights(rna.matureIndexs(rna.matureMRNAIndexs))) / ...
                (rna.expression(rna.matureIndexs)' * rna.molecularWeights(rna.matureIndexs)), ...
                (matureRNACounts(rna.matureMRNAIndexs)' * rna.molecularWeights(rna.matureIndexs(rna.matureMRNAIndexs))) / ...
                (matureRNACounts' * rna.molecularWeights(rna.matureIndexs)),...
                'relative', 60e-2);
            assertIn(corr(rna.expression(rna.matureIndexs(rna.matureMRNAIndexs)), mRNAExp), [0.85 1]);
            assertIn(180 / pi * acos(rna.expression(rna.matureIndexs(rna.matureMRNAIndexs))' * mRNAExp / ...
                (sqrt(rna.expression(rna.matureIndexs(rna.matureMRNAIndexs))' * rna.expression(rna.matureIndexs(rna.matureMRNAIndexs))) * ...
                sqrt(mRNAExp' * mRNAExp))), ...
                [0 25]);
            
            %tRNAs
            assertElementsAlmostEqual(...
                (rna.expression(rna.matureIndexs(rna.matureTRNAIndexs))' * rna.molecularWeights(rna.matureIndexs(rna.matureTRNAIndexs))) / ...
                (rna.expression(rna.matureIndexs)' * rna.molecularWeights(rna.matureIndexs)), ...
                (matureRNACounts(rna.matureTRNAIndexs)' * rna.molecularWeights(rna.matureIndexs(rna.matureTRNAIndexs))) / ...
                (matureRNACounts' * rna.molecularWeights(rna.matureIndexs)),...
                'relative', 30e-2);
            
            %% assertions -- Protein
            expProd = rna.matureRNAGeneComposition(g.mRNAIndexs, :) * matureRNAExpression;
            
            totMonExp = ...
                + pmExp ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * pcExp;
            
            %mature/bound monomer counts increased
            assertElementsAlmostEqual(...
                (1 + cumGrowth) * mass.cellInitialDryWeight * (...
                + mass.initialFractionNTPsInRNAs * mass.dryWeightFractionRNA ...
                + mass.initialFractionAAsInMonomers * mass.dryWeightFractionProtein ...
                ) * edu.stanford.covert.util.ConstantUtil.nAvogadro, ...
                + matureRNACounts' * rna.molecularWeights(rna.matureIndexs) ...
                + matureMonCounts' * pm.molecularWeights(pm.matureIndexs) ...
                + matureCpxCounts' * pc.molecularWeights(pc.matureIndexs), ...
                'relative', 15e-2);
            
            %RNA production / expression matches expectations
            assertIn(corr(expProd, totMonExp), [0.9 1]);
            assertIn(180 / pi * acos((expProd' * totMonExp) / (sqrt(expProd' * expProd) * sqrt(totMonExp' * totMonExp))), ...
                [0 25]);
        end
        
        function testExponentialGrowthDistribution(this)
            nIter = 100;
            iterMax = 31000;
            stepSizeSec = 10;
            metStepSizeSec = 1000;
            initSeed = 2000;
            seed = 5000;
            
            %% calculate growth distribution
            [initGrowths, growths] = this.calcExponentialGrowthDistribution(nIter, iterMax, stepSizeSec, metStepSizeSec, initSeed, seed);
            
            %% calculate doubling time distribution
            doublingTimes = NaN(nIter, 1);
            for  i = 1:nIter
                tmp = metStepSizeSec * find(cumsum(metStepSizeSec * growths(i, :)) > 1, 1, 'first');
                if ~isempty(tmp)
                    doublingTimes(i) = tmp;
                else
                    doublingTimes(i) = iterMax / sum(metStepSizeSec * growths(i, :));
                end
            end
            
            %% load simulation
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], true);
            time = sim.state('Time');
            met = sim.process('Metabolism');
            
            %% assert exponential growth
            assertElementsAlmostEqual(1, ...
                mean(metStepSizeSec * sum(growths, 2)), ...
                'relative', 20e-2, 0);
            assertElementsAlmostEqual(mean(initGrowths), ...
                log(2) / time.cellCycleLength * met.macromoleculeStateInitializationGrowthFactor, ...
                'relative', 5e-2, 0);
            assertEqual(nIter, nnz(initGrowths));
            assertIn(nnz(growths(:, end)), nIter * [0.95 1]);
            assertElementsAlmostEqual(mean(doublingTimes), ...
                time.cellCycleLength, ...
                'relative', 0.2);
            
            %% plot
            %[axesHandle, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            %
            %cla(axesHandle);
            %h = plot(axesHandle, (1:metStepSizeSec:iterMax)/3600, [growths; mean(growths, 1)] / mean(initGrowths));
            %set(h(end), 'LineWidth', 12);
            %xlabel(axesHandle, 'Time (h)')
            %ylabel(axesHandle, 'Growth / Mean initial growth')
            %xlim(axesHandle, [0 iterMax]/3600);
            %saveas(figHandle, 'output/runMediumTests/ExponentialGrowth-GrowthRateDistribution.pdf');
            %
            %cla(axesHandle);
            %h = boundedline((1:metStepSizeSec:iterMax)/3600, repmat(mean(growths, 1) / mean(initGrowths), 4, 1), ...
            %    permute(quantile(growths / mean(initGrowths), 1 - [0.05 0.32 0.68 0.95], 1), [2 3 1]), 'cmap', winter(4), 'transparency', 0.5);
            %for i = 1:4
            %    set(h(i), 'Color', [0 0 0], 'LineWidth', 2);
            %end
            %box(axesHandle, 'on');
            %xlabel(axesHandle, 'Time (h)')
            %ylabel(axesHandle, 'Growth / Mean initial growth')
            %xlim(axesHandle, [0 iterMax]/3600);
            %saveas(figHandle, 'output/runMediumTests/ExponentialGrowth-GrowthRateConfidenceIntervals.pdf');
            %
            %cla(axesHandle);
            %hist(axesHandle, doublingTimes/3600, 40);
            %xlim(axesHandle, [5 20])
            %xlabel(axesHandle, 'Doubling Time (h)');
            %ylabel(axesHandle, 'Frequency');
            %
            %cla(axesHandle);
            %plot(axesHandle, initGrowths, growths(:, end), '.')
            %line([0 sum(initGrowths)], [0 sum(growths(:, end))], 'Color', 'r', 'Parent', axesHandle);
            %xlim(axesHandle, [min(initGrowths) max(initGrowths)])
            %ylim(axesHandle, [min(growths(:, end)) max(growths(:, end))])
            %xlabel(axesHandle, 'Initial growth (cell/s)')
            %ylabel(axesHandle, 'Final growth (cell/s)')
            %saveas(figHandle, 'output/runMediumTests/ExponentialGrowth-InitialVsFinalGrowthRate.pdf');
            %
            %cla(axesHandle);
            %plot(axesHandle, initGrowths, metStepSizeSec * sum(growths, 2), '.')
            %line([0 sum(initGrowths)], [0 metStepSizeSec * sum(sum(growths, 2))], 'Color', 'r', 'Parent', axesHandle);
            %xlim(axesHandle, [min(initGrowths) max(initGrowths)])
            %ylim(axesHandle, metStepSizeSec * [min(sum(growths, 2)) max(sum(growths, 2))])
            %xlabel(axesHandle, 'Initial growth (cell/s)')
            %ylabel(axesHandle, 'Total growth (cell)')
            %saveas(figHandle, 'output/runMediumTests/ExponentialGrowth-InitialVsTotalGrowthRate.pdf');
        end
        
        function testInitialVariance(this)
            nIter = 100;
            iterMax = 31000;
            stepSizeSec = 10;
            metStepSizeSec = Inf;
            initSeed = 3010;
            seed = 5020;
            
            [initCorr, divCorr] = this.calcProteinCopyNumberCorrelationDistribution(...
                nIter, iterMax, stepSizeSec, metStepSizeSec, initSeed, seed);
            
            %% assert initial and final correlations similar
            assertElementsAlmostEqual(1-mean(initCorr(:)), 1-mean(divCorr(:)), 'relative', 0.75);
            
            %% plot
            %h = plot(initCorr, divCorr, '.');
            %legend(h, {'Daughter Cell 1', 'Daughter Cell 2'}, 'Location', 'SouthWest');
            %xlabel('Initial Correlation');
            %ylabel('Final Correlation');
        end
    end
    
    methods (Static = true)
        function fitGrowthRate()
            nIter = 100;
            iterMax = 31000;
            stepSizeSec = 10;
            metStepSizeSec = 1000;
            initSeed = 3000;
            seed = 5010;
            
            [macromoleculeStateInitializationGrowthFactor, ~, exitFlag, output] = edu.stanford.covert.util.ComputationUtil.fzero(...
                @(x) edu.stanford.covert.cell.sim.ProteinGrowth_Test.fitGrowthRate_diff(...
                nIter, iterMax, stepSizeSec, metStepSizeSec, initSeed, seed, x), ...
                [1.2 1.6], optimset('TolX', 1e-2)) %#ok<ASGLU,NOPRT>
            
            if exitFlag ~= 1
                throw(MException('ProteinGrowth_Test:error', output.message));
            end
        end
        
        function value = fitGrowthRate_diff(nIter, iterMax, stepSizeSec, metStepSizeSec, initSeed, seed, macromoleculeStateInitializationGrowthFactor)
            fixtureFileName = 'Simulation_FitGrowthRate.mat';
            
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], true);
            sim.state('ProteinMonomer').minimumAverageExpression = 15;
            sim.process('Metabolism').macromoleculeStateInitializationGrowthFactor = macromoleculeStateInitializationGrowthFactor;
            fitter = edu.stanford.covert.cell.sim.util.FitConstants(sim, struct('method', 'heuristic', 'verbosity', 0));
            fitter.run();
            edu.stanford.covert.cell.sim.SimulationFixture.store(sim, fixtureFileName);
            clear sim fitter ans;
            
            [~, growths] = edu.stanford.covert.cell.sim.ProteinGrowth_Test.calcExponentialGrowthDistribution(...
                nIter, iterMax, stepSizeSec, metStepSizeSec, initSeed, seed, [], 2.9, fixtureFileName);
            value = 1 - mean(metStepSizeSec * sum(growths, 2));
        end
        
        function fitInitialProteinCopyNumberVariance()
            nIter = 100;
            iterMax = 31000;
            stepSizeSec = 10;
            metStepSizeSec = Inf;
            initSeed = 3010;
            seed = 5020;
            
            [macromoleculeStateInitializationVariation, ~, exitFlag, output] = edu.stanford.covert.util.ComputationUtil.fzero(...
                @(x) edu.stanford.covert.cell.sim.ProteinGrowth_Test.fitInitialProteinCopyNumberVariance_diff(...
                nIter, iterMax, stepSizeSec, metStepSizeSec, initSeed, seed, x), ...
                [1 3], optimset('TolX', 1e-2)) %#ok<ASGLU,NOPRT>
            
            if exitFlag ~= 1
                throw(MException('ProteinGrowth_Test:error', output.message));
            end
        end
        
        function value = fitInitialProteinCopyNumberVariance_diff(nIter, iterMax, stepSizeSec, metStepSizeSec, initSeed, seed, macromoleculeStateInitializationVariation)
            [initCorr, divCorr] = edu.stanford.covert.cell.sim.ProteinGrowth_Test.calcProteinCopyNumberCorrelationDistribution(...
                nIter, iterMax, stepSizeSec, metStepSizeSec, initSeed, seed, ...
                macromoleculeStateInitializationVariation);
            value = mean(initCorr(:)) - mean(divCorr(:));
        end
        
        function [initCorr, divCorr] = calcProteinCopyNumberCorrelationDistribution(...
                nIter, iterMax, stepSizeSec, metStepSizeSec, initSeed, seed, ...
                varargin)
            
            %calculate growth distribution
            [~, ~, ~, matureMonCounts, matureCpxCounts, ~, initMatureMonCounts, initMatureCpxCounts] = ...
                edu.stanford.covert.cell.sim.ProteinGrowth_Test.calcExponentialGrowthDistribution(...
                nIter, iterMax, stepSizeSec, metStepSizeSec, initSeed, ...
                seed, @(iter) log(2)/iterMax * exp(log(2) * iter / iterMax), varargin{:});
            
            %load simulation
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], {});
            g = sim.gene;
            time = sim.state('Time');
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            pcComp = sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3);
            mRNAExp = rna.matureRNAGeneComposition(g.mRNAIndexs, :) * rna.expression(rna.matureIndexs) ./ (...
                log(2)/time.cellCycleLength + pm.decayRates(pm.matureIndexs));
            
            %check average correlations very high
            assertIn(corr(mRNAExp, sum(...
                + initMatureMonCounts ...
                + pcComp * initMatureCpxCounts ...
                , 2)), [1-1e-3 1]);
            assertIn(corr(mRNAExp, sum(...
                + matureMonCounts ...
                + pcComp * matureCpxCounts ...
                , 2)), [1-1e-3 1]);
            
            %simulate division
            initCorr = zeros(nIter, 1);
            divCorr = zeros(nIter, 2);
            divMatureMonCounts = zeros([size(matureMonCounts) 2]);
            divMatureCpxCounts = zeros([size(matureCpxCounts) 2]);
            for i = 1:nIter
                tmp = cumsum([matureMonCounts(:, i); matureCpxCounts(:, i)]);
                order = randperm(tmp(end));
                tmpProts = zeros(size(tmp));
                for j = 1:tmp(end)/2
                    idx = find(order(j) <= tmp, 1, 'first');
                    tmpProts(idx) = tmpProts(idx) + 1;
                end
                divMatureMonCounts(:, i, 1) = tmpProts(1:size(matureMonCounts, 1));
                divMatureCpxCounts(:, i, 1) = tmpProts(size(matureMonCounts, 1)+1:end);
                divMatureMonCounts(:, i, 2) = matureMonCounts(:, i) - tmpProts(1:size(matureMonCounts, 1));
                divMatureCpxCounts(:, i, 2) = matureCpxCounts(:, i) - tmpProts(size(matureMonCounts, 1)+1:end);
                
                initCorr(i) = corr(mRNAExp, ...
                    + initMatureMonCounts(:, i) ...
                    + pcComp * initMatureCpxCounts(:, i));
                divCorr(i, 1) = corr(mRNAExp, ...
                    + tmpProts(1:size(matureMonCounts, 1)) + ...
                    + pcComp * tmpProts(size(matureMonCounts, 1)+1:end));
                divCorr(i, 2) = corr(mRNAExp, ...
                    + matureMonCounts(:, i) - tmpProts(1:size(matureMonCounts, 1)) + ...
                    + pcComp * (matureCpxCounts(:, i) - tmpProts(size(matureMonCounts, 1)+1:end)));
            end
        end
        
        function [initGrowths, growths, ...
                matureRNACounts, matureMonCounts, matureCpxCounts, ...
                initMatureRNACounts, initMatureMonCounts, initMatureCpxCounts] = ...
                calcExponentialGrowthDistribution(nIter, iterMax, stepSizeSec, metStepSizeSec, initSeed, seed, ...
                varargin)
            
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'Metabolism'
                'MacromolecularComplexation'
                });
            nRna = numel(sim.state('Rna').matureIndexs);
            nMon = numel(sim.state('ProteinMonomer').matureIndexs);
            nCpx = numel(sim.state('ProteinComplex').matureIndexs);
            verbosity = sim.verbosity;
            clear sim;
            
            growths = zeros(nIter, ceil(iterMax / metStepSizeSec));
            initGrowths = zeros(nIter, 1);
            matureRNACounts = zeros(nRna, nIter);
            matureMonCounts = zeros(nMon, nIter);
            matureCpxCounts = zeros(nCpx, nIter);
            initMatureRNACounts = zeros(nRna, nIter);
            initMatureMonCounts = zeros(nMon, nIter);
            initMatureCpxCounts = zeros(nCpx, nIter);
            
            for i = 1:nIter
                if verbosity > 0
                    fprintf('Iter %d\n', i);
                end
                if isfinite(metStepSizeSec)
                    [initGrowths(i), ~, growths(i, :), ~, ~, ~, ~, ...
                        matureRNACounts(:, i), matureMonCounts(:, i), matureCpxCounts(:, i), ...
                        initMatureRNACounts(:, i), initMatureMonCounts(:, i), initMatureCpxCounts(:, i)] = edu.stanford.covert.cell.sim.ProteinGrowth_Test.simulateExponentialGrowth(...
                        iterMax, initSeed * (1 + i), seed * (1 + i), stepSizeSec, metStepSizeSec, varargin{:});
                else
                    [initGrowths(i), ~, ~, ~, ~, ~, ~, ...
                        matureRNACounts(:, i), matureMonCounts(:, i), matureCpxCounts(:, i), ...
                        initMatureRNACounts(:, i), initMatureMonCounts(:, i), initMatureCpxCounts(:, i)] = edu.stanford.covert.cell.sim.ProteinGrowth_Test.simulateExponentialGrowth(...
                        iterMax, initSeed * (1 + i), seed * (1 + i), stepSizeSec, metStepSizeSec, varargin{:});
                end
            end
        end
        
        function [initGrowth, initMetabolicEnzymes, growths, cumGrowth, rnaExp, pmExp, pcExp, ...
                matureRNACounts, matureMonCounts, matureCpxCounts, ...
                initMatureRNACounts, initMatureMonCounts, initMatureCpxCounts] = ...
                simulateExponentialGrowth(iterMax, initSeed, seed, stepSizeSec, metStepSizeSec, ...
                growthFunc, macromoleculeStateInitializationVariation, fixtureFileName, rnaExpression)
            import edu.stanford.covert.util.ConstantUtil;
            
            if nargin < 6
                growthFunc = [];
            end
            if nargin < 7
                macromoleculeStateInitializationVariation = [];
            end
            if nargin < 8
                fixtureFileName = [];
            end
            if nargin < 9
                rnaExpression = [];
            end
            
            warningState = warning('query', 'WholeCell:warning');
            warning('off', 'WholeCell:warning');
            
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load(fixtureFileName, {
                'Metabolism'
                'MacromolecularComplexation'
                });
            %sim.applyOptions('verbosity', 1);
            
            comp = sim.compartment;
            g = sim.gene;
            
            time = sim.state('Time');
            mass = sim.state('Mass');
            m = sim.state('Metabolite');
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            mr = sim.state('MetabolicReaction');
            
            met = sim.process('Metabolism');
            
            pc.formationProcesses(~ismember(pc.formationProcesses, [9 10])) = 0;
            pc.formationProcesses(pc.formationProcesses == 9) = 1; %Macromolecular complexation
            pc.formationProcesses(pc.formationProcesses == 10) = 2; %Metabolism
            
            if ~isempty(macromoleculeStateInitializationVariation)
                pm.macromoleculeStateInitializationVariation = macromoleculeStateInitializationVariation;
            end
            
            initialGrowthFilterWidth = mr.initialGrowthFilterWidth;
            mr.initialGrowthFilterWidth = Inf;
            
            if ~isempty(rnaExpression)
                rna.expression = rnaExpression;
            end
            
            %% constants
            rnaProd = edu.stanford.covert.util.ComputationUtil.invertCompositionMatrix(rna.nascentRNAMatureRNAComposition) * ...
                (rna.expression(rna.matureIndexs) .* ...
                (log(2) / time.cellCycleLength + rna.decayRates(rna.matureIndexs)));
            rnaProd = rna.nascentRNAMatureRNAComposition * rnaProd / sum(rnaProd) * ...
                (rna.expression(rna.matureIndexs) * mass.cellInitialDryWeight * mass.dryWeightFractionRNA * ConstantUtil.nAvogadro / ...
                (rna.molecularWeights(rna.matureIndexs)' * rna.expression(rna.matureIndexs)))' * (...
                + log(2) / time.cellCycleLength ...
                + rna.decayRates(rna.matureIndexs));
            rnaProd = 0.3929 * rnaProd / sum(rnaProd);
            
            mRNAExp = rna.matureRNAGeneComposition(g.mRNAIndexs, :) * rna.expression(rna.matureIndexs);
            mRNAExp = mRNAExp / sum(mRNAExp);
            nMonProd = mass.cellInitialDryWeight * mass.initialFractionAAsInMonomers * mass.dryWeightFractionProtein * ConstantUtil.nAvogadro / (pm.molecularWeights(pm.matureIndexs)' * mRNAExp) * mRNAExp' * (...
                + log(2) / time.cellCycleLength ...
                + pm.decayRates(pm.matureIndexs));
            
            metGasIdxs = met.substrateIndexs({'O2'; 'CO2'});
            metNTPIdxs = met.substrateIndexs(m.wholeCellModelIDs(m.ntpIndexs));
            metNMPIdxs = met.substrateIndexs(m.wholeCellModelIDs(m.nmpIndexs));
            metAAIdxs = met.substrateIndexs(m.wholeCellModelIDs(m.aminoAcidIndexs));
            metSubIdxs = [metNTPIdxs; metNMPIdxs; metAAIdxs];
            
            matureRNABaseCounts = rna.baseCounts(rna.matureIndexs, m.nmpIndexs)';
            matureRNADecayRates = min(1, rna.decayRates(rna.matureIndexs));
            matureRNAMRNAComposition = rna.matureRNAGeneComposition(g.mRNAIndexs, :);
            
            matureMonBaseCounts = pm.baseCounts(pm.matureIndexs, m.aminoAcidIndexs)';
            matureCpxBaseCounts = pc.baseCounts(pc.matureIndexs, m.aminoAcidIndexs)';
            matureMonDecayRates = pm.decayRates(pm.matureIndexs);
            matureCpxDecayRates = pc.decayRates(pc.matureIndexs);
            
            nComplexs = numel(pc.matureIndexs);
            pcMonComp = sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3);
            notFormComplexTfs = pc.formationProcesses(pc.matureIndexs) ~= 1 | ...
                any(any(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3), 1)';
            
            metSubstrateMonomerLocalIndexs = met.substrateMonomerLocalIndexs(:, 1);
            metSubstrateComplexLocalIndexs = met.substrateComplexLocalIndexs(:, 1);
            metSubstrateMonomerGlobalIndexs = met.substrateMonomerGlobalIndexs(:, 1);
            metSubstrateComplexGlobalIndexs = met.substrateComplexGlobalIndexs(:, 1);
            
            %% initialize
            while true
                if initSeed
                    edu.stanford.covert.cell.sim.ProteinGrowth_Test.seedRandStream(sim, initSeed);
                end
                
                sim.initializeState();
                pc.counts(sub2ind(size(pc.counts), pc.matureIndexs, pc.compartments(pc.matureIndexs))) = ...
                    + pc.counts(sub2ind(size(pc.counts), pc.matureIndexs, pc.compartments(pc.matureIndexs))) ...
                    + pc.counts(sub2ind(size(pc.counts), pc.nascentIndexs, pc.compartments(pc.nascentIndexs)));
                pc.counts(pc.nascentIndexs, :) = 0;
                sim.evolveState();
                
                if abs(mr.growth - mr.meanInitialGrowthRate) / mr.meanInitialGrowthRate < initialGrowthFilterWidth
                    break;
                end
                
                initSeed = mr.randStream.randi([0 2^32-1], 1);
            end
            
            if seed
                edu.stanford.covert.cell.sim.ProteinGrowth_Test.seedRandStream(sim, seed);
            end
            
            %% keep track of initial state
            initGrowth = mr.growth;
            initMetabolicEnzymes = met.enzymes;
            
            initMatureRNACounts = rna.counts(rna.matureIndexs, comp.cytosolIndexs);
            initMatureMonCounts = sum(...
                + pm.counts(pm.matureIndexs, :) ...
                + pm.counts(pm.boundIndexs, :), 2);
            initMatureCpxCounts = sum(...
                + pc.counts(pc.matureIndexs, :) ...
                + pc.counts(pc.boundIndexs, :), 2);
            tfs = any(any(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3), 1);
            initMatureRNACounts(setdiff(1:end, rna.matureMRNAIndexs)) = ...
                initMatureRNACounts(setdiff(1:end, rna.matureMRNAIndexs)) + ...
                sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), tfs, :), 3) * initMatureCpxCounts(tfs);
            initMatureMonCounts = ...
                initMatureMonCounts + ...
                sum(pc.proteinComplexComposition(g.mRNAIndexs, tfs, :), 3) * initMatureCpxCounts(tfs);
            initMatureCpxCounts(tfs) = 0;
            
            initGas = met.substrates(metGasIdxs, met.compartmentIndexs_extracellular);
            
            %% simulate
            met.stepSizeSec = metStepSizeSec;
            
            matureRNACounts = initMatureRNACounts;
            matureMonCounts = initMatureMonCounts;
            matureCpxCounts = initMatureCpxCounts;
            
            ntpCounts = m.counts(m.ntpIndexs, comp.cytosolIndexs);
            nmpCounts = m.counts(m.nmpIndexs, comp.cytosolIndexs);
            aaCounts = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            
            growths = zeros(ceil(iterMax / met.stepSizeSec), 1);
            cumGrowth = 0;
            rnaExp = zeros(size(rna.matureIndexs));
            pmExp = zeros(size(pm.matureIndexs));
            pcExp = zeros(size(pc.matureIndexs));
            
            assertElementsAlmostEqual(...
                mass.cellInitialDryWeight * (...
                + mass.initialFractionNTPsInRNAs * mass.dryWeightFractionRNA ...
                + mass.initialFractionAAsInMonomers * mass.dryWeightFractionProtein ...
                ) * edu.stanford.covert.util.ConstantUtil.nAvogadro, ...
                + matureRNACounts' * rna.molecularWeights(rna.matureIndexs) ...
                + matureMonCounts' * pm.molecularWeights(pm.matureIndexs) ...
                + matureCpxCounts' * pc.molecularWeights(pc.matureIndexs), ...
                'relative', 40e-2);
            
            for iter = 1:stepSizeSec:iterMax
                if sim.verbosity > 0 && mod(iter, 1000) == 1
                    fprintf('Iter %d\n', iter);
                end
                
                %% time
                time.values = iter;
                
                %% mock transcription, RNA maturation, and decay
                newRNAs = rna.randStream.stochasticRound(rnaProd * mr.growth * time.cellCycleLength / log(2) * stepSizeSec);
                decayedRNAs = rna.randStream.stochasticRound(matureRNACounts .* matureRNADecayRates * stepSizeSec);
                matureRNACounts = ...
                    + matureRNACounts ...
                    + newRNAs ...
                    - decayedRNAs;
                ntpCounts = max(0, ...
                    + ntpCounts ...
                    - matureRNABaseCounts * newRNAs ...
                    );
                nmpCounts = ...
                    + nmpCounts ...
                    + matureRNABaseCounts * decayedRNAs;
                
                %% mock translation, protein maturation, and decay
                mRNAExp = matureRNAMRNAComposition * matureRNACounts;
                newMonomers = pm.randStream.stochasticRound(nMonProd * mr.growth * time.cellCycleLength / log(2) * max(0, mRNAExp / sum(mRNAExp)) * stepSizeSec);
                decayedMonomers = pm.randStream.stochasticRound(matureMonCounts .* matureMonDecayRates * stepSizeSec);
                decayedComplexs = pc.randStream.stochasticRound(matureCpxCounts .* matureCpxDecayRates * stepSizeSec);
                
                matureMonCounts = ...
                    + matureMonCounts ...
                    + newMonomers ...
                    - decayedMonomers;
                matureCpxCounts = ...
                    + matureCpxCounts ...
                    - decayedComplexs;
                aaCounts = max(0, ...
                    + aaCounts ...
                    - matureMonBaseCounts * (newMonomers - decayedMonomers) ...
                    - matureCpxBaseCounts * (            - decayedComplexs) ...
                    );
                
                %mock protein complexation
                newComplexs = max(0, floor(min(matureMonCounts(:, ones(nComplexs, 1)) ./ pcMonComp, [], 1))');
                newComplexs(notFormComplexTfs) = 0;
                matureMonCounts = ...
                    + matureMonCounts ...
                    - pcMonComp * newComplexs;
                matureCpxCounts = ...
                    + matureCpxCounts ...
                    + newComplexs;
                
                %simulate metabolism
                if ~isempty(growthFunc)
                    mr.growth = growthFunc(iter);
                elseif mod(iter, met.stepSizeSec) == 1
                    met.substrates(metGasIdxs, met.compartmentIndexs_extracellular) = initGas;
                    met.substrates(metSubIdxs, met.compartmentIndexs_cytosol) = [ntpCounts; nmpCounts; aaCounts];
                    met.substrates(metSubstrateMonomerLocalIndexs, :) = matureMonCounts(metSubstrateMonomerGlobalIndexs, ones(3, 1));
                    met.substrates(metSubstrateComplexLocalIndexs, :) = matureCpxCounts(metSubstrateComplexGlobalIndexs, ones(3, 1));
                    
                    met.enzymes(met.enzymeMonomerLocalIndexs) = matureMonCounts(met.enzymeMonomerGlobalIndexs);
                    met.enzymes(met.enzymeComplexLocalIndexs) = matureCpxCounts(met.enzymeComplexGlobalIndexs);
                    
                    met.evolveState();
                    
                    ntpCounts = met.substrates(metNTPIdxs, met.compartmentIndexs_cytosol);
                    nmpCounts = met.substrates(metNMPIdxs, met.compartmentIndexs_cytosol);
                    aaCounts = met.substrates(metAAIdxs, met.compartmentIndexs_cytosol);
                    
                    growths((iter-1) / met.stepSizeSec + 1) = mr.growth;
                    
                    if mr.growth == 0
                        break;
                    end
                end
                
                %keep track of progres
                cumGrowth = cumGrowth + stepSizeSec * mr.growth;
                if nargout >= 5
                    rnaExp = rnaExp + matureRNACounts;
                    pmExp = pmExp + matureMonCounts;
                    pcExp = pcExp + matureCpxCounts;
                end
            end
            
            warning(warningState.state, 'WholeCell:warning');
        end
        
        function seedRandStream(sim, seed)
            if sim.verbosity > 0
                fprintf('Seed %d\n', seed);
            end
            
            sim.applyOptions('seed', seed);
            sim.seedRandStream();
            for i = 1:numel(sim.states)
                o = sim.states{i};
                o.seed = seed;
                o.seedRandStream();
            end
            
            for i = 1:numel(sim.processes)
                o = sim.processes{i};
                o.seed = seed;
                o.seedRandStream();
            end
        end
        
        function [exponentialGrowthLimitingEnzymeIdxs, growthLimitingEnzymeIdxs] = ...
                calcGrowthLimitingEnzymes(sim, initGrowth, initMetabolicEnzymes)
            %% references
            time = sim.state('Time');
            met = sim.process('Metabolism');
            
            %% enzymes preventing exponential increase in growth rate
            exponentialGrowthLimitingEnzymeIdxs = [];
            enzymes = met.enzymes;
            while true
                [growth, ~, ~, reducedCosts] = met.calcGrowthRate(met.calcFluxBounds(met.substrates, enzymes, met.fbaReactionBounds, met.fbaEnzymeBounds));
                if growth >= initGrowth * exp(log(2) * time.values / time.cellCycleLength)
                    break;
                end
                
                rxnTfs = abs(reducedCosts) > 1e-12;
                enzIdxs = find(any(met.reactionCatalysisMatrix(rxnTfs, :), 1));
                enzymes(enzIdxs) = max(initMetabolicEnzymes(enzIdxs), enzymes(enzIdxs) + 1);
                exponentialGrowthLimitingEnzymeIdxs = unique([exponentialGrowthLimitingEnzymeIdxs enzIdxs]);
            end
            
            for i = numel(exponentialGrowthLimitingEnzymeIdxs):-1:1
                tmp = enzymes;
                tmp(exponentialGrowthLimitingEnzymeIdxs(i)) = met.enzymes(exponentialGrowthLimitingEnzymeIdxs(i));
                if abs(growth - met.calcGrowthRate(met.calcFluxBounds(met.substrates, tmp, met.fbaReactionBounds, met.fbaEnzymeBounds))) / growth < 1e-6
                    exponentialGrowthLimitingEnzymeIdxs(i) = [];
                    enzymes = tmp;
                end
            end
            
            assertIn(met.calcGrowthRate(met.calcFluxBounds(met.substrates, enzymes, met.fbaReactionBounds, met.fbaEnzymeBounds)), ...
                [(1 - 1e-2) * initGrowth * exp(log(2) * time.values / time.cellCycleLength) Inf]);
            
            %% enzymes limiting growth
            [growth, ~, ~, reducedCosts] = met.calcGrowthRate(met.calcFluxBounds(met.substrates, enzymes, met.fbaReactionBounds, met.fbaEnzymeBounds));
            rxnTfs = abs(reducedCosts) > 1e-12;
            growthLimitingEnzymeIdxs = find(any(met.reactionCatalysisMatrix(rxnTfs, :), 1));
            
            for i = numel(growthLimitingEnzymeIdxs):-1:1
                tmp = enzymes;
                tmp(growthLimitingEnzymeIdxs(i)) = enzymes(growthLimitingEnzymeIdxs(i)) - 1;
                if abs(growth - met.calcGrowthRate(met.calcFluxBounds(met.substrates, tmp, met.fbaReactionBounds, met.fbaEnzymeBounds))) / growth < 1e-6
                    growthLimitingEnzymeIdxs(i) = [];
                    enzymes = tmp;
                end
            end
            
            %% print
            if sim.verbosity > 0
                fprintf('Exponential growth-limiting enzymes\n');
                for i = 1:numel(exponentialGrowthLimitingEnzymeIdxs)
                    fprintf('%s\t%d\t%d\t%d\n', ...
                        met.enzymeWholeCellModelIDs{exponentialGrowthLimitingEnzymeIdxs(i)}, ...
                        met.enzymes(exponentialGrowthLimitingEnzymeIdxs(i)), ...
                        enzymes(exponentialGrowthLimitingEnzymeIdxs(i)), ...
                        initMetabolicEnzymes(exponentialGrowthLimitingEnzymeIdxs(i)));
                end
                fprintf('\n');
                
                fprintf('Growth-limiting enzymes\n');
                for i = 1:numel(growthLimitingEnzymeIdxs)
                    fprintf('%s\t%d\t%d\t%d\n', ...
                        met.enzymeWholeCellModelIDs{growthLimitingEnzymeIdxs(i)}, ...
                        met.enzymes(growthLimitingEnzymeIdxs(i)), ...
                        enzymes(growthLimitingEnzymeIdxs(i)), ...
                        initMetabolicEnzymes(growthLimitingEnzymeIdxs(i)));
                end
                fprintf('\n');
            end
        end
    end
end