%RNA Expression Analysis
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 12/22/2011
classdef RNAExpression
    methods (Static = true)
        function run(sim, fileName)
            %import classes
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.PrintUtil;
            import edu.stanford.covert.util.ComputationUtil;
            
            %load simulation
            g = sim.gene;
            r = sim.state('Rna');
            
            %get data
            rnaExp = r.expectedGeneExpression(:, 1);
            
            exp_mRNAWt = 0.041;
            exp_rRNAWt = [0.017; 0.271; 0.525];
            exp_sRNAWt = 0;
            exp_tRNAWt = 0.146;
            assertEqual(1, exp_mRNAWt + sum(exp_rRNAWt) + exp_sRNAWt + exp_tRNAWt);
            
            %calculate gene molecular weights approximately
            geneMWs = (r.matureRNAGeneComposition * (r.molecularWeights(r.aminoacylatedIndexs) ./ (r.matureRNAGeneComposition' * g.lengths))) ...
                .* g.lengths;
            
            %check RNA MWs computed correctly
            monocistronicIdxs = find(sum(r.matureRNAGeneComposition(g.mRNAIndexs, :), 1) == 1);
            monocistronicmRNAGeneIdxs = g.mRNAIndexs(any(r.matureRNAGeneComposition(g.mRNAIndexs, monocistronicIdxs), 2));
            assertElementsAlmostEqual(geneMWs(monocistronicmRNAGeneIdxs), r.molecularWeights(r.aminoacylatedIndexs(monocistronicIdxs)), 'relative', 1e-12);
            assertElementsAlmostEqual(geneMWs(g.rRNAIndexs), r.molecularWeights(r.aminoacylatedIndexs(r.matureRRNAIndexs)), 'relative', 1e-12);
            assertElementsAlmostEqual(geneMWs(g.sRNAIndexs), r.molecularWeights(r.aminoacylatedIndexs(r.matureSRNAIndexs)), 'relative', 1e-12);
            assertElementsAlmostEqual(geneMWs(g.tRNAIndexs), r.molecularWeights(r.aminoacylatedIndexs(r.matureTRNAIndexs)), 'relative', 1e-12);
            assertElementsAlmostEqual(r.molecularWeights(r.aminoacylatedIndexs), r.matureRNAGeneComposition' * geneMWs, 'relative', 1e-12);
            
            %reconstruct gene expression
            mRNAWt = rnaExp(g.mRNAIndexs)' * geneMWs(g.mRNAIndexs);
            rRNAWt = rnaExp(g.ribosomalRRNAIndexs) .* geneMWs(g.ribosomalRRNAIndexs);
            sRNAWt = rnaExp(g.sRNAIndexs)' * geneMWs(g.sRNAIndexs);
            tRNAWt = rnaExp(g.tRNAIndexs)' * geneMWs(g.tRNAIndexs);
            
            rnaExp(g.mRNAIndexs) = rnaExp(g.mRNAIndexs) * exp_mRNAWt / mRNAWt;
            rnaExp(g.ribosomalRRNAIndexs) = rnaExp(g.ribosomalRRNAIndexs) .* exp_rRNAWt ./ rRNAWt;
            rnaExp(g.sRNAIndexs) = rnaExp(g.sRNAIndexs) * exp_sRNAWt / sRNAWt;
            rnaExp(g.tRNAIndexs) = rnaExp(g.tRNAIndexs) * exp_tRNAWt / tRNAWt;
            
            rnaExp = rnaExp / sum(rnaExp);
            
            exp_RNAWt = [exp_mRNAWt; exp_rRNAWt; exp_sRNAWt; exp_tRNAWt];
            RNAwt = [mRNAWt; rRNAWt; sRNAWt; tRNAWt];
            exp_RNAWt = exp_RNAWt / sum(exp_RNAWt);
            RNAwt = RNAwt / sum(RNAwt);
            
            %plot RNA weight fractions
            [axesHandle, figHandle] = PlotUtil.newAxesHandle();
            
            cla(axesHandle);
            hold(axesHandle, 'on');
            h = [
                plot(axesHandle, RNAwt(1), exp_RNAWt(1), 'r.')
                plot(axesHandle, RNAwt(2:4), exp_RNAWt(2:4), 'g.')
                plot(axesHandle, RNAwt(5), exp_RNAWt(5), 'b.')
                plot(axesHandle, RNAwt(6), exp_RNAWt(6), 'c.')
                ];            
            legend(h, {'mRNA', 'rRNA', 'sRNA', 'tRNA'}, 'Location', 'NorthWest');
            xlabel(axesHandle, 'Reconstruction - Earliest Stage', 'FontSize', 12)
            ylabel(axesHandle, 'Reconstruction - Later Stage', 'FontSize', 12)
            set(axesHandle, 'xscale', 'log', 'yscale', 'log', 'xlimmode', 'auto',  'ylimmode', 'auto', 'xminortick', 'off', 'yminortick', 'off')
            axis(axesHandle, 'equal');
            xlim(axesHandle, [min([xlim(axesHandle) ylim(axesHandle)]) max([xlim(axesHandle) ylim(axesHandle)])]);
            ylim(axesHandle, xlim(axesHandle));
            line(xlim(axesHandle), ylim(axesHandle), 'Parent', axesHandle, 'LineStyle', ':', 'Color', [0.5 0.5 0.5])
            
            saveas(figHandle, [fileName '-EarlyVsLaterPhaseReconstruction-WeightFractions.pdf']);
            close(figHandle);
            
            %plot all RNA
            [axesHandle, figHandle] = PlotUtil.newAxesHandle();
            
            cla(axesHandle);
            hold(axesHandle, 'on');
            h = [
                plot(axesHandle, rnaExp(g.mRNAIndexs), r.expectedGeneExpression(g.mRNAIndexs, 1), 'r.')
                plot(axesHandle, rnaExp(g.rRNAIndexs), r.expectedGeneExpression(g.rRNAIndexs, 1), 'g.')
                plot(axesHandle, rnaExp(g.sRNAIndexs), r.expectedGeneExpression(g.sRNAIndexs, 1), 'b.')
                plot(axesHandle, rnaExp(g.tRNAIndexs), r.expectedGeneExpression(g.tRNAIndexs, 1), 'c.')
                ];
            legend(h, {'mRNA', 'rRNA', 'sRNA', 'tRNA'}, 'Location', 'NorthWest');
            xlabel(axesHandle, 'Reconstruction - Earliest Stage', 'FontSize', 12)
            ylabel(axesHandle, 'Reconstruction - Later Stage', 'FontSize', 12)
            set(axesHandle, 'xscale', 'log', 'yscale', 'log', 'xlimmode', 'auto',  'ylimmode', 'auto', 'xminortick', 'off', 'yminortick', 'off')
            axis(axesHandle, 'equal');
            xlim(axesHandle, [min([xlim(axesHandle) ylim(axesHandle)]) max([xlim(axesHandle) ylim(axesHandle)])]);
            ylim(axesHandle, xlim(axesHandle));
            line(xlim(axesHandle), ylim(axesHandle), 'Parent', axesHandle, 'LineStyle', ':', 'Color', [0.5 0.5 0.5])
            
            saveas(figHandle, [fileName '-EarlyVsLaterPhaseReconstruction.pdf']);
            close(figHandle);
            
            %save to excel sheet
            content = [g.wholeCellModelIDs num2cell(rnaExp) num2cell(r.expectedGeneExpression(:, 1))];
            colLabels = {'Gene', 'Expression-Early Reconstruction', 'Expression-Later Reconstruction'};
            PrintUtil.printToFile(content, colLabels, [fileName '.xls'], 'RNA Expression');
        end
    end
end