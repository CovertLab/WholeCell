%SingleCellVariation
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 3/23/2011
classdef SingleCellVariation
    methods (Static)
        function run(sim, fileName)
            import edu.stanford.covert.cell.sim.analysis.SingleCellVariation;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.PrintUtil;
            
            %sample growth rate distribution
            initialGrowthFilterWidth = sim.state('MetabolicReaction').initialGrowthFilterWidth;
            sim.state('MetabolicReaction').initialGrowthFilterWidth = Inf;
            growthRates = SingleCellVariation.sampleSingleCellVariationDistribution(sim, 100);
            sim.state('MetabolicReaction').initialGrowthFilterWidth = initialGrowthFilterWidth;
            assertElementsAlmostEqual(sim.state('MetabolicReaction').meanInitialGrowthRate, mean(growthRates), 'relative', 0.50, 0);
            
            %% excel file
            [content, colLabels, indentation] = SingleCellVariation.printSingleCellVariationDistribution(sim, growthRates);
            if nargin == 1
                PrintUtil.printToStdIO(content, colLabels, struct('indentation', indentation));
            else
                PrintUtil.printToFile(content, colLabels, [fileName '.xls'], 'GrowthRate', struct('indentation', indentation));
            end
            
            %% plots
            if nargin == 1
                SingleCellVariation.plotSingleCellVariationGrowthRateDistribution(sim, growthRates, PlotUtil.newAxesHandle());
                SingleCellVariation.plotSingleCellVariationDoublingTimeDistribution(sim, growthRates, PlotUtil.newAxesHandle());
                SingleCellVariation.plotExpected_Vs_Simulated_RNA_Weight_Fractions(sim, PlotUtil.newAxesHandle());
                SingleCellVariation.plotExpected_Vs_Simulated_RNA_Expression(sim, PlotUtil.newAxesHandle());
                SingleCellVariation.plotExpected_Vs_Simulated_rRNAExpression(sim, PlotUtil.newAxesHandle());
                SingleCellVariation.plotExpected_Vs_Simulated_sRNA_Expression(sim, PlotUtil.newAxesHandle());
                SingleCellVariation.plotExpected_Vs_Simulated_tRNA_Expression(sim, PlotUtil.newAxesHandle());
                SingleCellVariation.plotExpected_Vs_Simulated_Monomer_Expression(sim, PlotUtil.newAxesHandle());
                SingleCellVariation.plotExpected_Vs_Simulated_NTP_Incorporation(sim, PlotUtil.newAxesHandle());
                SingleCellVariation.plotExpected_Vs_Simulated_AA_Incorporation(sim, PlotUtil.newAxesHandle());
            else
                [axesHandle, figHandle] = PlotUtil.newAxesHandle();
                
                cla(axesHandle);
                SingleCellVariation.plotSingleCellVariationGrowthRateDistribution(sim, growthRates, axesHandle);
                saveas(figHandle, [fileName '-GrowthRate.pdf']);
                
                cla(axesHandle);
                SingleCellVariation.plotSingleCellVariationDoublingTimeDistribution(sim, growthRates, axesHandle);
                saveas(figHandle, [fileName '-DoublingTime.pdf']);
                
                cla(axesHandle);
                SingleCellVariation.plotExpected_Vs_Simulated_RNA_Weight_Fractions(sim, axesHandle);
                saveas(figHandle, [fileName '-Expected_Vs_Simulated_RNA_Weight_Fractions.pdf']);
                
                cla(axesHandle);
                SingleCellVariation.plotExpected_Vs_Simulated_RNA_Expression(sim, axesHandle);
                saveas(figHandle, [fileName '-Expected_Vs_Simulated_RNA_Expression.pdf']);
                
                cla(axesHandle);
                SingleCellVariation.plotExpected_Vs_Simulated_rRNAExpression(sim, axesHandle);
                saveas(figHandle, [fileName '-Expected_Vs_Simulated_rRNA_Expression.pdf']);
                
                cla(axesHandle);
                SingleCellVariation.plotExpected_Vs_Simulated_sRNA_Expression(sim, axesHandle);
                saveas(figHandle, [fileName '-Expected_Vs_Simulated_sRNA_Expression.pdf']);
                
                cla(axesHandle);
                SingleCellVariation.plotExpected_Vs_Simulated_tRNA_Expression(sim, axesHandle);
                saveas(figHandle, [fileName '-Expected_Vs_Simulated_tRNA_Expression.pdf']);
                
                cla(axesHandle);
                SingleCellVariation.plotExpected_Vs_Simulated_Monomer_Expression(sim, axesHandle);
                saveas(figHandle, [fileName '-Expected_Vs_Simulated_Monomer_Expression.pdf']);
                
                cla(axesHandle);
                SingleCellVariation.plotExpected_Vs_Simulated_NTP_Incorporation(sim, axesHandle);
                saveas(figHandle, [fileName '-Expected_Vs_Simulated_NTP_Incorporation.pdf']);
                
                cla(axesHandle);
                SingleCellVariation.plotExpected_Vs_Simulated_AA_Incorporation(sim, axesHandle);
                saveas(figHandle, [fileName '-Expected_Vs_Simulated_AA_Incorporation.pdf']);
                
                close(figHandle);
            end
        end
    end
    
    methods (Static = true)
        function growthRates = sampleSingleCellVariationDistribution(sim, nTrials)
            warnStatus = warning('query', 'WholeCell:warning');
            warning('off', 'WholeCell:warning');
            
            growthRates = zeros(nTrials, 1);
            r = sim.state('MetabolicReaction');
            for i = 1:nTrials
                %seed rand stream
                sim.applyOptions(struct('seed', i));
                for j = 1:numel(sim.states)
                    o = sim.states{j};
                    o.seed = i;
                end
                for j = 1:numel(sim.processes)
                    o = sim.processes{j};
                    o.seed = i;
                end
                
                %calculate growth rate
                sim.initializeState();
                growthRates(i) = r.growth;
            end
            
            warning(warnStatus.state, 'WholeCell:warning');
        end
    end
    
    %printing
    methods (Static = true)
        function [content, colLabels, indentation] = printSingleCellVariationDistribution(sim, growthRates)
            %import classes
            import edu.stanford.covert.util.ConstantUtil;
            import edu.stanford.covert.cell.sim.Simulation;
            
            content = cell(0, 4);
            colLabels = {'Trial', 'Growth Rate (cell/s)', 'Doubling Time (hr)'};
            
            %time
            content = [content;
                num2cell(zeros(size(growthRates)))  cellfun(@(x) num2str(x), num2cell((1:numel(growthRates))'), 'UniformOutput', false)  num2cell(growthRates)  num2cell(1./(growthRates*3600)/sim.state('Mass').timeAveragedCellWeight)
                ];
            
            content = [content;{
                0 'Mean'          mean(growthRates)  1/(mean(growthRates)*3600)/sim.state('Mass').timeAveragedCellWeight
                0 'Min'           min(growthRates)   1/(max(growthRates)*3600)/sim.state('Mass').timeAveragedCellWeight
                0 'Max'           max(growthRates)   1/(min(growthRates)*3600)/sim.state('Mass').timeAveragedCellWeight
                0 'Experimental'  1/(sim.state('Time').cellCycleLength)/sim.state('Mass').timeAveragedCellWeight sim.state('Time').cellCycleLength/3600
                }];
            
            %format output
            indentation = cell2mat(content(:, 1));
            content = content(:, 2:end);
        end
    end
    
    %plotting
    methods (Static = true)
        function plotSingleCellVariationGrowthRateDistribution(sim, growthRates, axesHandle)
            hist(axesHandle, growthRates, 10);
            xlims = xlim(axesHandle);
            ylims = ylim(axesHandle);
            xlabel('Growth (cell/s)', 'FontSize', 12);
            ylabel('Frequency', 'FontSize', 12);
            
            h2 = line(mean(growthRates) * [1 1], ylims);
            set(h2, 'Color', 'g');
            
            h3 = line(median(growthRates) * [1 1], ylims);
            set(h3, 'Color', 'c');
            
            mu = mean(growthRates);
            sigma = std(growthRates);
            h4 = line([mu-sigma; mu-sigma], ylims);
            h5 = line([mu+sigma; mu+sigma], ylims);
            set(h4, 'Color', 'r');
            set(h5, 'Color', 'r');
            
            mr = sim.state('MetabolicReaction');
            mu = mr.meanInitialGrowthRate;
            wd = mr.initialGrowthFilterWidth;
            h6 = line(mu*(1-wd)*[1 1], ylims);
            h7 = line(mu*(1+wd)*[1 1], ylims);
            set(h6, 'Color', 'y');
            set(h7, 'Color', 'y');
            
            h = legend([h2 h3 h4 h6], 'Mean', 'Median', 'Mean \pm 1 Std', 'Filter');
            set(h, 'Location', 'NorthWest');
            
            xlim(axesHandle, xlims);
            ylim(axesHandle, ylims);
        end
        
        function plotSingleCellVariationDoublingTimeDistribution(sim, growthRates, axesHandle)
            thresh = sim.state('Time').cellCycleLength/3600 * 2;
            
            doublingTimes = 1./(growthRates * 3600)/sim.state('Mass').timeAveragedCellWeight;
            hist(axesHandle, doublingTimes(doublingTimes <= thresh), 10);
            h1 = findobj(gca, 'Type', 'patch');
            xlabel('Doubling Time (h)', 'FontSize', 12);
            ylabel('Frequency', 'FontSize', 12);
            xlim([floor(min(doublingTimes(doublingTimes <= thresh))) ceil(max(doublingTimes(doublingTimes <= thresh)))]);
            
            h2 = line(1/(mean(growthRates) * 3600) / sim.state('Mass').timeAveragedCellWeight * [1 1], [0 max(ylim)]);
            set(h2, 'Color', 'g');
            
            h3 = line(1/(median(growthRates) * 3600) / sim.state('Mass').timeAveragedCellWeight * [1 1], [0 max(ylim)]);
            set(h3, 'Color', 'c');
            
            h4 = line(sim.state('Time').cellCycleLength * [1 1] / 3600, [0 max(ylim)]);
            set(h4, 'Color', 'r');
            
            h = legend([h1(1) h2 h3 h4], sprintf('Simulation, growth (%.1f%%)', 100 * sum(doublingTimes <= thresh) / numel(doublingTimes)), 'Sim-Mean', 'Sim-Median', 'Exp-Mean');
            set(h, 'Location', 'NorthEast');
        end
        
        function plotExpected_Vs_Simulated_RNA_Weight_Fractions(sim, axesHandle)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            import edu.stanford.covert.util.ConstantUtil;
            
            r = sim.state('Rna');
            
            weightFractions = [...
                r.molecularWeights(r.matureIndexs(r.matureMRNAIndexs))' * ...
                sum(sum(r.counts(r.matureIndexs(r.matureMRNAIndexs), :, :), 3), 2);...
                r.molecularWeights(r.matureIndexs(r.matureRibosomalRRNAIndexs)) .* ...
                sum(sum(r.counts(r.matureIndexs(r.matureRibosomalRRNAIndexs), :, :), 3), 2);...
                r.molecularWeights(r.matureIndexs(r.matureSRNAIndexs))' * ...
                sum(sum(r.counts(r.matureIndexs(r.matureSRNAIndexs), :, :), 3), 2)
                r.molecularWeights(r.matureIndexs(r.matureTRNAIndexs))' * ...
                sum(sum(r.counts(r.matureIndexs(r.matureTRNAIndexs), :, :), 3), 2)] / ...
                ConstantUtil.nAvogadro;
            
            plot(axesHandle, r.expectedWeightFractions, weightFractions, '.', 'MarkerSize', 20);
            line([0 1], [0 sum(weightFractions)], 'Parent', axesHandle, 'Color', 'r');
            xlim([min(r.expectedGeneDecayRates) max(r.weightFractions)]);
            ylim([min(weightFractions) max(weightFractions)]);
            xlabel(axesHandle, 'Expected Weight Fraction', 'fontSize', 16);
            ylabel(axesHandle, 'Simulated Weight (g)', 'fontSize', 16);
        end
        
        function plotExpected_Vs_Simulated_RNA_Expression(sim, axesHandle)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            r = sim.state('Rna');
            
            expression = r.expression(r.matureIndexs);
            
            RNAs = sum(sum(...
                r.counts(r.processedIndexs,     :, :) + ...
                r.counts(r.matureIndexs,        :, :) + ...
                r.counts(r.boundIndexs,         :, :) + ...
                r.counts(r.misfoldedIndexs,     :, :) + ...
                r.counts(r.damagedIndexs,       :, :) + ...
                r.counts(r.aminoacylatedIndexs, :, :), ...
                3), 2);
            
            plot(axesHandle, expression, RNAs, '.', 'MarkerSize', 20);
            line([0 1], [0 sum(RNAs)], 'Parent', axesHandle, 'Color', 'r');
            xlim([min(expression) max(expression)]);
            ylim([min(RNAs) max(RNAs)]);
            xlabel(axesHandle, 'Expected Expression', 'fontSize', 16);
            ylabel(axesHandle, 'Simulated Counts', 'fontSize', 16);
        end
        
        function plotExpected_Vs_Simulated_rRNAExpression(sim, axesHandle)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            r = sim.state('Rna');
            c = sim.state('ProteinComplex');
            
            expression = r.geneExpression(sim.gene.rRNAIndexs);
            expression = expression / sum(expression);
            
            RNAs = sum(sum(...
                r.counts(r.processedIndexs(    r.matureRRNAIndexs), :, :) + ...
                r.counts(r.matureIndexs(       r.matureRRNAIndexs), :, :) + ...
                r.counts(r.boundIndexs(        r.matureRRNAIndexs), :, :) + ...
                r.counts(r.misfoldedIndexs(    r.matureRRNAIndexs), :, :) + ...
                r.counts(r.damagedIndexs(      r.matureRRNAIndexs), :, :) + ...
                r.counts(r.aminoacylatedIndexs(r.matureRRNAIndexs), :, :), ...
                3), 2) + ...
                sum(c.proteinComplexComposition(sim.gene.rRNAIndexs, :, :), 3) * ...
                sum(sum(...
                c.counts(c.matureIndexs,      :, :) + ...
                c.counts(c.inactivatedIndexs, :, :) + ...
                c.counts(c.boundIndexs,       :, :) + ...
                c.counts(c.misfoldedIndexs,   :, :) + ...
                c.counts(c.damagedIndexs,     :, :), 3), 2);
            
            plot(axesHandle, expression, RNAs, '.', 'MarkerSize', 20);
            line([0 1], [0 sum(RNAs)], 'Parent', axesHandle, 'Color', 'r');
            xlim([min(0.3, min(expression)) max(0.4, max(expression))]);
            ylim([min(RNAs)-1 max(RNAs)+1]);
            xlabel(axesHandle, 'Expected Expression', 'fontSize', 16);
            ylabel(axesHandle, 'Simulated Counts', 'fontSize', 16);
        end
        
        function plotExpected_Vs_Simulated_sRNA_Expression(sim, axesHandle)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            r = sim.state('Rna');
            
            expression = r.geneExpression(sim.gene.sRNAIndexs);
            expression = expression / sum(expression);
            
            RNAs = sum(sum(...
                r.counts(r.processedIndexs(    r.matureSRNAIndexs), :, :) + ...
                r.counts(r.matureIndexs(       r.matureSRNAIndexs), :, :) + ...
                r.counts(r.boundIndexs(        r.matureSRNAIndexs), :, :) + ...
                r.counts(r.misfoldedIndexs(    r.matureSRNAIndexs), :, :) + ...
                r.counts(r.damagedIndexs(      r.matureSRNAIndexs), :, :) + ...
                r.counts(r.aminoacylatedIndexs(r.matureSRNAIndexs), :, :), ...
                3), 2);
            
            plot(axesHandle, expression, RNAs, '.', 'MarkerSize', 20);
            line([0 1], [0 sum(RNAs)], 'Parent', axesHandle, 'Color', 'r');
            xlim([min(expression) max(expression)]);
            ylim([min(RNAs) max(RNAs)]);
            xlabel(axesHandle, 'Expected Expression', 'fontSize', 16);
            ylabel(axesHandle, 'Simulated Counts', 'fontSize', 16);
        end
        
        function plotExpected_Vs_Simulated_tRNA_Expression(sim, axesHandle)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            r = sim.state('Rna');
            
            expression = r.geneExpression(sim.gene.tRNAIndexs);
            expression = expression / sum(expression);
            
            RNAs = sum(sum(...
                r.counts(r.processedIndexs(    r.matureTRNAIndexs), :, :) + ...
                r.counts(r.matureIndexs(       r.matureTRNAIndexs), :, :) + ...
                r.counts(r.boundIndexs(        r.matureTRNAIndexs), :, :) + ...
                r.counts(r.misfoldedIndexs(    r.matureTRNAIndexs), :, :) + ...
                r.counts(r.damagedIndexs(      r.matureTRNAIndexs), :, :) + ...
                r.counts(r.aminoacylatedIndexs(r.matureTRNAIndexs), :, :), ...
                3), 2);
            
            plot(axesHandle, expression, RNAs, '.', 'MarkerSize', 20);
            line([0 1], [0 sum(RNAs)], 'Parent', axesHandle, 'Color', 'r');
            xlim([min(expression) max(expression)]);
            ylim([min(RNAs) max(RNAs)]);
            xlabel(axesHandle, 'Expected Expression', 'fontSize', 16);
            ylabel(axesHandle, 'Copy Number', 'fontSize', 16);
        end
                
        function plotExpected_Vs_Simulated_Monomer_Expression(sim, axesHandle)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            r = sim.state('Rna');
            m = sim.state('ProteinMonomer');
            c = sim.state('ProteinComplex');
            
            expression = r.geneExpression(sim.gene.mRNAIndexs) ./ m.halfLives(m.matureIndexs);
            expression = expression / sum(expression);
            
            monomers = ...
                sum(sum(...
                m.counts(m.matureIndexs,      :, :) + ...
                m.counts(m.inactivatedIndexs, :, :) + ...
                m.counts(m.boundIndexs,       :, :) + ...
                m.counts(m.misfoldedIndexs,   :, :) + ...
                m.counts(m.damagedIndexs,     :, :), 3), 2) + ...
                sum(c.proteinComplexComposition(sim.gene.mRNAIndexs, :, :), 3) * ...
                sum(sum(...
                c.counts(c.matureIndexs,      :, :) + ...
                c.counts(c.inactivatedIndexs, :, :) + ...
                c.counts(c.boundIndexs,       :, :) + ...
                c.counts(c.misfoldedIndexs,   :, :) + ...
                c.counts(c.damagedIndexs,     :, :), 3), 2);
            
            plot(axesHandle, expression, monomers, '.', 'MarkerSize', 10);
            line([0 1], [0 sum(monomers)], 'Parent', axesHandle, 'Color', 'r');
            xlim([min(expression) max(expression)]);
            ylim([min(monomers) max(monomers)]);
            xlabel(axesHandle, 'Expected Expression', 'fontSize', 16);
            ylabel(axesHandle, 'Simulated Counts', 'fontSize', 16);
        end
        
        function plotExpected_Vs_Simulated_NTP_Incorporation(sim, axesHandle)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            r = sim.state('Rna');
            m = sim.state('Metabolite');
            
            nmpComposition = sum(m.nmpComposition, 2);
            
            ntps = sum(sum(multiprod(r.baseCounts(:, sim.state('Metabolite').nmpIndexs)', r.counts, [1 2], [1 2]), 3), 2);
            
            plot(axesHandle, nmpComposition, ntps, '.', 'MarkerSize', 20);
            line([0 1], [0 sum(ntps)], 'Parent', axesHandle, 'Color', 'r');
            xlim([min(nmpComposition) max(nmpComposition)]);
            ylim([min(ntps) max(ntps)]);
            xlabel(axesHandle, 'Expected Incorporation', 'fontSize', 16);
            ylabel(axesHandle, 'Simulated Counts', 'fontSize', 16);
        end
        
        function plotExpected_Vs_Simulated_AA_Incorporation(sim, axesHandle)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            met = sim.state('Metabolite');
            m = sim.state('ProteinMonomer');
            c = sim.state('ProteinComplex');
            
            aaComposition = sum(met.aaComposition, 2);
            
            aas = multiprod(m.baseCounts(:, sim.state('Metabolite').aminoAcidIndexs)', m.counts, [1 2], [1 2]) + ...
                multiprod(c.baseCounts(:, sim.state('Metabolite').aminoAcidIndexs)', c.counts, [1 2], [1 2]);
            aas = sum(sum(aas, 3), 2);
            
            plot(axesHandle, aaComposition, aas, '.', 'MarkerSize', 20);
            line([0 1], [0 sum(aas)],'Parent', axesHandle, 'Color', 'r');
            xlim([min(aaComposition) max(aaComposition)]);
            ylim([min(aas) max(aas)]);
            xlabel(axesHandle, 'Expected Incorporation', 'fontSize', 16);
            ylabel(axesHandle, 'Simulated Counts', 'fontSize', 16);
        end
    end
end