%FBA
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 3/24/2011
classdef FBA
    %printing
    methods (Static)
        function run(sim, fileName)
            import edu.stanford.covert.cell.sim.analysis.FBA;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.PrintUtil;
            
            %excel file
            [content, colLabels, indentation] = FBA.printNetworkReduction(sim);
            if nargin == 1
                PrintUtil.printToStdIO(content, colLabels, struct('indentation', indentation));
            else
                PrintUtil.printToFile(content, colLabels, [fileName '.xls'], 'NetworkReduction', struct('indentation', indentation));
            end
            
            %plots
            if nargin == 1
                FBA.plotNetworkReduction(sim, PlotUtil.newAxesHandle());
            else
                [axesHandle, figHandle] = PlotUtil.newAxesHandle();
                
                cla(axesHandle);
                FBA.plotNetworkReduction(sim, axesHandle);
                saveas(figHandle, [fileName '-NetworkReduction.pdf']);
                
                close(figHandle);
            end
        end
        
        function [content, colLabels, indentation] = printNetworkReduction(sim)
            p = sim.process('Metabolism');
            
            [subCmpRemoved, rxnRemoved] = p.formulateFBA([], [], true);
            
            content = cell(0, 6);
            colLabels = {'ID', 'Compartment', 'Reason', 'Iteration-1', 'Iteration-2'};
            reasons = {
                'Disconnected'
                'Substrate cannot be balanced: only appears in single reaction'
                'Substrate cannot be balanced: move in same direction in all reactions in which it appears'
                'Reaction has zero flux in all solutions of S*v = 0: reaction has no non-zero component in null space of S';
                };
            
            subIDs = p.substrateWholeCellModelIDs;
            cmpIDs = p.compartment.wholeCellModelIDs(p.substrateMetaboliteCompartmentIndexs(1, :));
            content = [content;{
                0 'Substrates' [] [] [] []}];
            for i = 1:size(subCmpRemoved, 1)
                if subCmpRemoved(i, 1) == 0
                    continue;
                end
                [sIdx, cIdx] = ind2sub([size(p.reactionStoichiometryMatrix, 1) size(p.reactionStoichiometryMatrix, 3)], i);
                content = [content;{
                    1 subIDs{sIdx} cmpIDs{cIdx} reasons{subCmpRemoved(i, 3)} subCmpRemoved(i, 1) subCmpRemoved(i, 2)}]; %#ok<AGROW>
            end
            
            rxnIDs = p.reactionWholeCellModelIDs;
            content = [content;{
                0 'Reactions' [] [] [] []}];
            for i = 1:size(rxnRemoved, 1)
                if rxnRemoved(i, 3) == 0
                    continue;
                end
                content = [content;{
                    1 rxnIDs{i} [] reasons{rxnRemoved(i, 3)} rxnRemoved(i, 1) rxnRemoved(i, 2)}]; %#ok<AGROW>
            end
            
            %format output
            indentation = cell2mat(content(:, 1));
            content = content(:, 2:end);
        end
        
        function plotNetworkReduction(sim, axesHandle)
            p = sim.process('Metabolism');
            
            [subCmpRemoved, rxnRemoved] = p.formulateFBA([], [], true);
            subCmpRemoved(subCmpRemoved(:, 1) == 0, 1) = max(subCmpRemoved(:, 1)) + 1;
            rxnRemoved(rxnRemoved(:, 1) == 0, 1) = max(rxnRemoved(:, 1)) + 1;
            
            levels = sortrows(unique([subCmpRemoved; rxnRemoved], 'rows'), [1 2 3]);
            
            [~, subCmpIdxs] = sortrows(subCmpRemoved, 1:3);
            [~, rxnIdxs] = sortrows(rxnRemoved, 1:3);
            
            [~, subCmpLevels] = ismember(subCmpRemoved, levels, 'rows');
            [~, rxnLevels] = ismember(rxnRemoved, levels, 'rows');
            
            hold on;
            
            x = [
                0.5
                numel(rxnIdxs)+0.5
                numel(rxnIdxs)+0.5
                0.5
                ];
            y = [
                0.5
                0.5
                numel(subCmpIdxs)+0.5
                numel(subCmpIdxs)+0.5
                ];
            patch(x, y, 1, 'FaceColor', [1 1 1], 'Parent', axesHandle, 'FaceAlpha', 1);
            patch(x, y, 1, 'CDataMapping', 'direct', 'Parent', axesHandle, 'FaceAlpha', 0.25, 'EdgeAlpha', 1);
            for i = 2:size(levels, 1)
                j = find(subCmpLevels(subCmpIdxs) < i, 1, 'last');
                k = find(rxnLevels(rxnIdxs) < i, 1, 'last');
                if isempty(j), j = 0; end;
                if isempty(k), k = 0; end;
                
                x = [
                    k+0.5
                    numel(rxnIdxs)+0.5
                    numel(rxnIdxs)+0.5
                    k+0.5
                    ];
                y = [
                    j+0.5
                    j+0.5
                    numel(subCmpIdxs)+0.5
                    numel(subCmpIdxs)+0.5
                    ];
                patch(x, y, 1, 'FaceColor', [1 1 1], 'Parent', axesHandle, 'FaceAlpha', 1);                
                patch(x, y, i, 'CDataMapping', 'direct', 'Parent', axesHandle, 'FaceAlpha', 0.25, 'EdgeAlpha', 1);
            end
            
            rxnSMat = reshape(permute(p.reactionStoichiometryMatrix, [2 1 3]), ...
                size(p.reactionStoichiometryMatrix, 2), [])';
            
            [y, x] = find(...
                ((rxnSMat(subCmpIdxs, rxnIdxs) < 0 & repmat(p.reactionBounds(rxnIdxs, 2)' > 0, numel(subCmpIdxs), 1)) | ...
                (rxnSMat(subCmpIdxs, rxnIdxs) > 0 & repmat(p.reactionBounds(rxnIdxs, 1)' < 0, numel(subCmpIdxs), 1))) & ...
                ~((rxnSMat(subCmpIdxs, rxnIdxs) > 0 & repmat(p.reactionBounds(rxnIdxs, 2)' > 0, numel(subCmpIdxs), 1)) | ...
                (rxnSMat(subCmpIdxs, rxnIdxs) < 0 & repmat(p.reactionBounds(rxnIdxs, 1)' < 0, numel(subCmpIdxs), 1))));
            h1 = plot(axesHandle, x, y, 'r.', 'MarkerSize', 6);
            
            [y, x] = find(...
                ~((rxnSMat(subCmpIdxs, rxnIdxs) < 0 & repmat(p.reactionBounds(rxnIdxs, 2)' > 0, numel(subCmpIdxs), 1)) | ...
                (rxnSMat(subCmpIdxs, rxnIdxs) > 0 & repmat(p.reactionBounds(rxnIdxs, 1)' < 0, numel(subCmpIdxs), 1))) &  ...
                ((rxnSMat(subCmpIdxs, rxnIdxs) > 0 & repmat(p.reactionBounds(rxnIdxs, 2)' > 0, numel(subCmpIdxs), 1)) | ...
                (rxnSMat(subCmpIdxs, rxnIdxs) < 0 & repmat(p.reactionBounds(rxnIdxs, 1)' < 0, numel(subCmpIdxs), 1))));
            h2 = plot(axesHandle, x, y, 'g.', 'MarkerSize', 6);
            
            [y, x] = find(...
                ((rxnSMat(subCmpIdxs, rxnIdxs) < 0 & repmat(p.reactionBounds(rxnIdxs, 2)' > 0, numel(subCmpIdxs), 1)) | ...
                (rxnSMat(subCmpIdxs, rxnIdxs) > 0 & repmat(p.reactionBounds(rxnIdxs, 1)' < 0, numel(subCmpIdxs), 1))) &  ...
                ((rxnSMat(subCmpIdxs, rxnIdxs) > 0 & repmat(p.reactionBounds(rxnIdxs, 2)' > 0, numel(subCmpIdxs), 1)) | ...
                (rxnSMat(subCmpIdxs, rxnIdxs) < 0 & repmat(p.reactionBounds(rxnIdxs, 1)' < 0, numel(subCmpIdxs), 1))));
            h3 = plot(axesHandle, x, y, 'b.', 'MarkerSize', 6);
            
            legend([h1 h2 h3], {'Drained'; 'Produced'; 'Bidirectional'}, 'Location', 'NorthEastOutside');
            
            colormap(jet);
            
            xlim([0.5 numel(rxnIdxs)+0.5])
            ylim([0.5 numel(subCmpIdxs)+0.5]);
            
            box('on');
            xlabel('Reactions', 'FontSize', 12);
            ylabel('Substrates', 'FontSize', 12);
            set(axesHandle, 'YDir', 'reverse');
            set(axesHandle, 'XAxisLocation', 'top');
            axis('square')
        end
    end
end