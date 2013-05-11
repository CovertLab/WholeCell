%RNA
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 4/22/2011
classdef RNA
    methods (Static = true)
        function run(sim, nSample, nIter, fileName)
            import edu.stanford.covert.cell.sim.analysis.RNA;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.PrintUtil;
            
            %simulate
            [expRnaSyn, expRnaExp] = RNA.calcExpectedSynthesisExpression(sim);
            [simRnaSyn, simRnaExp] = RNA.sampleRNASynthesis(sim, nSample, nIter);
            
            %print
            [content, colLabels, indentation] = RNA.tabulateRNASynthesisExpression(sim, expRnaSyn, expRnaExp, simRnaSyn, simRnaExp);
            if nargin == 3
                PrintUtil.printToStdIO(content, colLabels, struct('indentation', indentation));
            else
                PrintUtil.printToFile(content, colLabels, [fileName '.xls'], 'Synthesis-Expression', struct('indentation', indentation));
            end
            
            %plot
            if nargin == 3
                RNA.plotRnaSynthesis(PlotUtil.newAxesHandle(), sim, expRnaSyn, simRnaSyn);
                RNA.plotRnaExpression(PlotUtil.newAxesHandle(), sim, expRnaExp, simRnaExp);
            else
                [axesHandle, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
                
                RNA.plotRnaSynthesis(axesHandle, sim, expRnaSyn, simRnaSyn);
                saveas(figHandle, [fileName '-Synthesis.pdf']);
                
                RNA.plotRnaExpression(axesHandle, sim, expRnaExp, simRnaExp);
                saveas(figHandle, [fileName '-Expression.pdf']);
                
                close(figHandle);
            end
        end
    end
    
    methods (Static = true)
        function [content, colLabels, indentation] = tabulateRNASynthesisExpression(sim, expRnaSyn, expRnaExp, simRnaSyn, simRnaExp)
            rna = sim.state('Rna');
            matIdxs = rna.matureIndexs;
            
            colLabels = [
                'ID', 'Name', 'Transcription Unit', ...
                'Expected Synthesis (molecules / cell cycle)', 'Mean Simulated Synthesis (molecules / cell cycle)', 'Std Simulated Synthesis (molecules / cell cycle)', ...
                cellfun(@(x) ['Synthesis - Simulation #' num2str(x)], num2cell(1:size(simRnaSyn, 2)), 'UniformOutput', false), ...
                'Expected Expression (molecules * s / cell cycle)', 'Mean Simulated Expression (molecules * s / cell cycle)', 'Std Simulated Expression (molecules * s / cell cycle)', ...
                cellfun(@(x) ['Expression - Simulation #' num2str(x)], num2cell(1:size(simRnaExp, 2)), 'UniformOutput', false), ...
                ];
            
            [~, tuIdxs] = find(rna.nascentRNAMatureRNAComposition);
            
            content = [
                rna.wholeCellModelIDs(matIdxs) rna.names(matIdxs) rna.wholeCellModelIDs(tuIdxs) ...
                num2cell([
                expRnaSyn mean(simRnaSyn, 2) std(simRnaSyn, 0, 2) simRnaSyn ...
                expRnaExp mean(simRnaExp, 2) std(simRnaExp, 0, 2) simRnaExp ...
                ])];
            
            indentation = zeros(size(content, 1), 1);
        end
        
        function plotRnaSynthesis(axesHandle, sim, expRnaSyn, simRnaSyn)
            rna = sim.state('Rna');
            
            hold(axesHandle, 'on');
            
            h = [
                plot(axesHandle, expRnaSyn(rna.matureMRNAIndexs), mean(simRnaSyn(rna.matureMRNAIndexs, :), 2), 'r.')
                plot(axesHandle, expRnaSyn(rna.matureRRNAIndexs), mean(simRnaSyn(rna.matureRRNAIndexs, :), 2), 'g.')
                plot(axesHandle, expRnaSyn(rna.matureSRNAIndexs), mean(simRnaSyn(rna.matureSRNAIndexs, :), 2), 'b.')
                plot(axesHandle, expRnaSyn(rna.matureTRNAIndexs), mean(simRnaSyn(rna.matureTRNAIndexs, :), 2), 'c.')
                ];
            line([0 max(max(expRnaSyn), max(mean(simRnaSyn, 2)))], [0 max(max(expRnaSyn), max(mean(simRnaSyn, 2)))], 'Parent', axesHandle);
            xlim(axesHandle, [0 max(max(expRnaSyn), max(mean(simRnaSyn, 2)))]);
            ylim(axesHandle, [0 max(max(expRnaSyn), max(mean(simRnaSyn, 2)))]);
            axis(axesHandle, 'square');
            box(axesHandle, 'on');
            xlabel(axesHandle, 'Expected Synthesis', 'FontSize', 12);
            ylabel(axesHandle, 'Mean Simulated Synthesis', 'FontSize', 12);
            legend(h, {'mRNA', 'rRNA', 'sRNA', 'tRNA'}, 'Location', 'NorthEastOutside');
        end
        
        function plotRnaExpression(axesHandle, sim, expRnaExp, simRnaExp)
            rna = sim.state('Rna');
            
            hold(axesHandle, 'on');
            
            h = [
                plot(axesHandle, expRnaExp(rna.matureMRNAIndexs), mean(simRnaExp(rna.matureMRNAIndexs, :), 2), 'r.')
                plot(axesHandle, expRnaExp(rna.matureRRNAIndexs), mean(simRnaExp(rna.matureRRNAIndexs, :), 2), 'g.')
                plot(axesHandle, expRnaExp(rna.matureSRNAIndexs), mean(simRnaExp(rna.matureSRNAIndexs, :), 2), 'b.')
                plot(axesHandle, expRnaExp(rna.matureTRNAIndexs), mean(simRnaExp(rna.matureTRNAIndexs, :), 2), 'c.')
                ];
            line([0 max(max(expRnaExp), max(mean(simRnaExp, 2)))], [0 max(max(expRnaExp), max(mean(simRnaExp, 2)))], 'Parent', axesHandle);
            xlim(axesHandle, [0 max(max(expRnaExp), max(mean(simRnaExp, 2)))]);
            ylim(axesHandle, [0 max(max(expRnaExp), max(mean(simRnaExp, 2)))]);
            axis(axesHandle, 'square');
            box(axesHandle, 'on');
            xlabel(axesHandle, 'Expected Expression', 'FontSize', 12);
            ylabel(axesHandle, 'Mean Simulated Expression', 'FontSize', 12);
            legend(h, {'mRNA', 'rRNA', 'sRNA', 'tRNA'}, 'Location', 'NorthEastOutside');
        end
        
        function [rnaSyn, rnaExp] = calcExpectedSynthesisExpression(sim)
            import edu.stanford.covert.util.ConstantUtil;
            
            time = sim.state('Time');
            mass = sim.state('Mass');
            rna = sim.state('Rna');
            
            exp = rna.expression(rna.matureIndexs);
            mws = rna.molecularWeights(rna.matureIndexs) / ConstantUtil.nAvogadro;
            drs = rna.decayRates(rna.matureIndexs);
            
            initCnts = exp * mass.cellInitialDryWeight * mass.dryWeightFractionRNA / (exp' * mws);
            
            rnaExp = initCnts * time.cellCycleLength / log(2);
            rnaSyn = initCnts + initCnts .* drs * time.cellCycleLength / log(2);
        end
        
        function [rnaSyn, rnaExp] = sampleRNASynthesis(sim, nSample, nIter)
            import edu.stanford.covert.cell.sim.analysis.RNA;
            
            rna = sim.state('Rna');
            nRna = size(rna.matureIndexs, 1);
            
            rnaSyn = zeros(nRna, nSample);
            rnaExp = zeros(nRna, nSample);
            
            for i = 1:nSample
                [rnaSyn(:, i), rnaExp(:, i)] = RNA.simulateRNASynthesis(nIter, i);
            end
        end
        
        function [rnaSyn, rnaExp] = simulateRNASynthesis(nIter, seed)
            %% references
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'Metabolism'
                'Transcription'
                'RNAProcessing'
                'RNAModification'
                'tRNAAminoacylation'
                'RNADecay'});
            
            g = sim.gene;
            rna = sim.state('Rna');
            time = sim.state('Time');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            %% seed
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
            
            %% turn off decay
            rna.decayRates(setdiff(1:end, rna.intergenicIndexs)) = 0;
            
            %% keep track of initial state
            initRNAs = rna.counts;
            initMonomers = pm.counts;
            initComplexs = pc.counts;
            
            %% simulate
            allRnaSyn = zeros(size(rna.counts));
            allRnaExp = zeros(size(rna.counts));
            for i = 1:nIter
                %mock protein synthesis
                pm.counts = ceil(initMonomers * exp(i * log(2) / time.cellCycleLength));
                pc.counts = ceil(initComplexs * exp(i * log(2) / time.cellCycleLength));
                
                %remember old RNA counts
                oldRnaCounts = rna.counts;
                
                %simulate
                sim.evolveState();
                
                %accounting of new RNAs
                allRnaSyn = allRnaSyn + max(0, rna.counts - oldRnaCounts);
                
                %accounting of expressed RNAs
                allRnaExp = allRnaExp + rna.counts;
                allRnaExp(rna.matureIndexs(setdiff(1:end, rna.matureMRNAIndexs)), sim.compartment.cytosolIndexs) = ...
                    allRnaExp(rna.matureIndexs(setdiff(1:end, rna.matureMRNAIndexs)), sim.compartment.cytosolIndexs) + ...
                    sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * sum(initComplexs(pc.matureIndexs, :) + initComplexs(pc.boundIndexs, :), 2);
            end
            
            rnaSyn = allRnaSyn(rna.matureIndexs, sim.compartment.cytosolIndexs);
            aminoacylatedTfs = allRnaSyn(rna.aminoacylatedIndexs, sim.compartment.cytosolIndexs) ~= 0;
            rnaSyn(aminoacylatedTfs) = allRnaSyn(rna.aminoacylatedIndexs(aminoacylatedTfs), sim.compartment.cytosolIndexs) - ...
                initRNAs(rna.matureIndexs(aminoacylatedTfs), sim.compartment.cytosolIndexs);
            
            rnaExp = ...
                allRnaExp(rna.processedIndexs, sim.compartment.cytosolIndexs) + ...
                allRnaExp(rna.matureIndexs, sim.compartment.cytosolIndexs) + ...
                allRnaExp(rna.aminoacylatedIndexs, sim.compartment.cytosolIndexs) + ...
                allRnaExp(rna.boundIndexs, sim.compartment.cytosolIndexs) + ...
                allRnaExp(rna.misfoldedIndexs, sim.compartment.cytosolIndexs) + ...
                allRnaExp(rna.damagedIndexs, sim.compartment.cytosolIndexs);
            
            %% scale to cell cycle
            rnaSyn = rnaSyn * 1 / (exp(nIter * log(2) / time.cellCycleLength) - 1);
            rnaExp = rnaExp * 1 / (exp(nIter * log(2) / time.cellCycleLength) - 1);
        end
    end
end