%SimulationStructure
% Outputs analysis of the structure of the simulation:
% - which processes use which metabolites / gene products / enzymes
% - which metabolites / gene products / enzymes aren't used by any process
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/12/2010
classdef SimulationStructure
    methods (Static = true)        
        function run(sim, fileName)
            import edu.stanford.covert.cell.sim.analysis.SimulationStructure;
            import edu.stanford.covert.cell.sim.util.PrintUtil;
            
            % excel file
            [content, colLabels] = SimulationStructure.processMetabolites(sim);
            if nargin == 1
                PrintUtil.printToStdIO(content, colLabels);
            else
                PrintUtil.printToFile(content, colLabels, [fileName '.xls'], 'Process-Metabolites');
            end
            
            [content, colLabels] = SimulationStructure.processGeneProducts(sim);
            if nargin == 1
                PrintUtil.printToStdIO(content, colLabels);
            else
                PrintUtil.printToFile(content, colLabels, [fileName '.xls'], 'Process-Gene Product');
            end
            
            [content, colLabels] = SimulationStructure.unincludedMetabolites(sim);
            if nargin == 1
                PrintUtil.printToStdIO(content, colLabels);
            else
                PrintUtil.printToFile(content, colLabels, [fileName '.xls'], 'Unincluded Metabolites');
            end
            
            [content, colLabels] = SimulationStructure.unincludedGenes(sim);
            if nargin == 1
                PrintUtil.printToStdIO(content, colLabels);
            else
                PrintUtil.printToFile(content, colLabels, [fileName '.xls'], 'Unincluded Genes');
            end
            
            [content, colLabels] = SimulationStructure.unincludedProteins(sim);
            if nargin == 1
                PrintUtil.printToStdIO(content, colLabels);
            else
                PrintUtil.printToFile(content, colLabels, [fileName '.xls'], 'Unincluded Proteins');
            end
            
            % plots
            if nargin == 1
                SimulationStructure.plotProcessMetaboliteSharing(sim);
                SimulationStructure.plotProcessMetabolites(sim);
                SimulationStructure.plotProcessSharedMetabolites(sim);
                SimulationStructure.plotProcessGeneProducts(sim);
            else
                [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
                
                clf(figHandle);
                SimulationStructure.plotProcessMetaboliteSharing(sim, figHandle);
                saveas(figHandle, [fileName '-ProcessMetaboliteSharing.pdf']);
                
                clf(figHandle);
                SimulationStructure.plotProcessMetabolites(sim, figHandle);
                saveas(figHandle, [fileName '-ProcessMetabolites.pdf']);
                
                clf(figHandle);
                SimulationStructure.plotProcessSharedMetabolites(sim, figHandle);
                saveas(figHandle, [fileName '-ProcessSharedMetabolites.pdf']);
                
                clf(figHandle);
                SimulationStructure.plotProcessGeneProducts(sim, figHandle);
                saveas(figHandle, [fileName '-ProcessGeneProducts.pdf']);
                
                close(figHandle);
            end
        end
    end
    
    %analysis
    methods (Static)
        %which processes use which metabolites
        function [content, colLabels] = processMetabolites(sim)
            colLabels = {'Process ID', 'Process Name', 'Metabolite ID', 'Metabolite Name'};
            content = cell(0, numel(colLabels));
            
            for i = 1:length(sim.processes)
                m = sim.processes{i};
                
                idxs = m.substrateMetaboliteLocalIndexs;
                
                content = [content;
                    repmat({m.wholeCellModelID m.name}, numel(idxs), 1) ...
                    m.substrateWholeCellModelIDs(idxs) ...
                    m.substrateNames(idxs)]; %#ok<AGROW>
            end
        end
        
        %which metabolites aren't used in any processes
        function [content, colLabels] = unincludedMetabolites(sim)
            colLabels = {'Metabolite ID', 'Metabolite Name'};
            
            included = edu.stanford.covert.cell.sim.analysis.SimulationStructure.processMetabolites(sim);
            
            [~, idxs] = setdiff(sim.state('Metabolite').wholeCellModelIDs, unique(included(:,3)));
            content = [
                sim.state('Metabolite').wholeCellModelIDs(idxs) ...
                sim.state('Metabolite').names(idxs)];
        end
        
        %which processes use which gene products, and in which compartments, as enzymes
        function [content, colLabels] = processGeneProducts(sim)
            colLabels = {'Process ID', 'Process Name', 'Gene ID', 'Gene Name', 'Product ID', 'Product Name', 'Product Compartment ID', 'Product Compartment Name', 'Implemented'};
            content = cell(0, numel(colLabels));
            
            for i = 1:length(sim.processes)
                m = sim.processes{i};
                
                %stimuli
                stimuliGeneComposition = m.stimuliGeneComposition();
                for j = 1:length(sim.gene.wholeCellModelIDs)
                    for k = 1:length(m.stimuliWholeCellModelIDs)
                        if ~stimuliGeneComposition(j, k); continue; end;
                        
                        if size(m.stimuliCompartments, 2) == 1
                            compartmentWholeCellModelID = sim.compartment.wholeCellModelIDs{m.stimuliCompartments(k)};
                            compartmentName = sim.compartment.names{m.stimuliCompartments(k)};
                        else
                            compartmentWholeCellModelID = cell(1,1);
                            compartmentName = cell(1,1);
                        end
                        
                        content = [content;
                            m.wholeCellModelID            m.name            ...
                            sim.gene.wholeCellModelIDs(j) sim.gene.names(j)  ...
                            m.stimuliWholeCellModelIDs(k) m.stimuliNames(k) ...
                            compartmentWholeCellModelID   compartmentName   ...
                            'Stimulus']; %#ok<AGROW>
                    end
                end
                
                %substrates
                substrateGeneComposition = m.substrateGeneComposition();
                for j = 1:length(sim.gene.wholeCellModelIDs)
                    for k = 1:length(m.substrateWholeCellModelIDs)
                        if ~substrateGeneComposition(j, k); continue; end;
                        
                        if size(m.substrateCompartments, 2) == 1
                            compartmentWholeCellModelID = sim.compartment.wholeCellModelIDs{m.substrateCompartments(k)};
                            compartmentName = sim.compartment.names{m.substrateCompartments(k)};
                        else
                            compartmentWholeCellModelID = cell(1,1);
                            compartmentName = cell(1,1);
                        end
                        
                        content = [content;
                            m.wholeCellModelID              m.name            ...
                            sim.gene.wholeCellModelIDs(j)   sim.gene.names(j)  ...
                            m.substrateWholeCellModelIDs(k) m.substrateNames(k) ...
                            compartmentWholeCellModelID     compartmentName   ...
                            'Substrate']; %#ok<AGROW>
                    end
                end
                
                %enzymes
                enzymeGeneComposition = m.enzymeGeneComposition();
                for j = 1:length(sim.gene.wholeCellModelIDs)
                    for k = 1:length(m.enzymeWholeCellModelIDs)
                        if ~enzymeGeneComposition(j, k); continue; end;
                        
                        if size(m.enzymeCompartments, 2) == 1
                            compartmentWholeCellModelID = sim.compartment.wholeCellModelIDs{m.enzymeCompartments(k)};
                            compartmentName = sim.compartment.names{m.enzymeCompartments(k)};
                        else
                            compartmentWholeCellModelID = cell(1,1);
                            compartmentName = cell(1,1);
                        end
                        
                        content = [content;
                            m.wholeCellModelID            m.name           ...
                            sim.gene.wholeCellModelIDs(j) sim.gene.names(j) ...
                            m.enzymeWholeCellModelIDs(k)  m.enzymeNames(k) ...
                            compartmentWholeCellModelID   compartmentName  ...
                            'Enzyme']; %#ok<AGROW>
                    end
                end
            end
        end
        
        %gene whose products aren't functionally used in any process
        function [content, colLabels] = unincludedGenes(sim)
            colLabels = {'Gene ID','Gene Name'};
            
            included = edu.stanford.covert.cell.sim.analysis.SimulationStructure.processGeneProducts(sim);
            
            [~, idxs] = setdiff(sim.gene.wholeCellModelIDs, unique(included(:,3)));
            content = [
                sim.gene.wholeCellModelIDs(idxs) ...
                sim.gene.names(idxs)];
        end
        
        %protein monomers and complexes which aren't used in any process
        function [content, colLabels] = unincludedProteins(sim)
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            colLabels = {'Protein ID', 'Protein Name'};
            
            includedMonomers = {};
            includedComplexs = {};
            for i = 1:length(sim.processes)
                m = sim.processes{i};
                
                includedMonomers = [includedMonomers; m.stimuliWholeCellModelIDs(m.stimulusMonomerLocalIndexs)]; %#ok<AGROW>
                includedMonomers = [includedMonomers; m.substrateWholeCellModelIDs(m.substrateMonomerLocalIndexs)]; %#ok<AGROW>
                includedMonomers = [includedMonomers; m.enzymeWholeCellModelIDs(m.enzymeMonomerLocalIndexs)]; %#ok<AGROW>
                
                includedComplexs = [includedComplexs; m.stimuliWholeCellModelIDs(m.stimulusComplexLocalIndexs)]; %#ok<AGROW>
                includedComplexs = [includedComplexs; m.substrateWholeCellModelIDs(m.substrateComplexLocalIndexs)]; %#ok<AGROW>
                includedComplexs = [includedComplexs; m.enzymeWholeCellModelIDs(m.enzymeComplexLocalIndexs)]; %#ok<AGROW>
            end
            
            [~, unincludedMonomerIdxs] = setdiff(pm.wholeCellModelIDs(pm.matureIndexs), unique(includedMonomers));
            [~, unincludedComplexIdxs] = setdiff(pc.wholeCellModelIDs(pc.matureIndexs), unique(includedComplexs));
            
            unincludedMonomerIdxs(any(any(pc.proteinComplexComposition(sim.gene.mRNAIndexs(unincludedMonomerIdxs), setdiff(1:end, unincludedComplexIdxs), :),3),2)) = [];
            
            content = [
                pm.wholeCellModelIDs(pm.matureIndexs(unincludedMonomerIdxs)) ...
                pm.names(pm.matureIndexs(unincludedMonomerIdxs));
                pc.wholeCellModelIDs(pc.matureIndexs(unincludedComplexIdxs)) ...
                pc.names(pc.matureIndexs(unincludedComplexIdxs))];
        end
    end
    
    %plots
    methods (Static)
        function plotProcessMetaboliteSharing(sim, figHandle)
            import edu.stanford.covert.cell.sim.analysis.SimulationStructure;
            
            content = SimulationStructure.processMetabolites(sim);
            [~, modIdxs] = ismember(content(:, 1), sim.processWholeCellModelIDs);
            [~, metIdxs] = ismember(content(:, 3), sim.state('Metabolite').wholeCellModelIDs);
            mat = zeros(numel(sim.processWholeCellModelIDs), numel(sim.state('Metabolite').wholeCellModelIDs));
            mat(sub2ind(size(mat), modIdxs, metIdxs)) = 1;
            
            mat2 = ones(numel(sim.processWholeCellModelIDs));
            for i = 1:size(mat,1)
                for j = 1:size(mat,1)
                    mat2(i,j) = 1 + double(any(mat(i, :) & mat(j, :)));
                end
            end
            mat2(sub2ind(size(mat2), 1:size(mat2,1), 1:size(mat2,1))) = ...
                double(mat2(sub2ind(size(mat2), 1:size(mat2,1), 1:size(mat2,1)))==1) + ...
                double(mat2(sub2ind(size(mat2), 1:size(mat2,1), 1:size(mat2,1)))~=1)*3;
            
            dist = pdist(mat2, 'cityblock');
            tree = linkage(dist, 'average');
            modOrder = optimalleaforder(tree, dist);
            
            %plot
            if ~exist('figHandle', 'var')
                axesHandle = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            else
                axesHandle = subplot(1, 1, 1, 'Parent', figHandle);
            end
            
            image(mat2(modOrder, modOrder), 'Parent', axesHandle);
            colormap(axesHandle, [
                1 1 1;
                0 0 0;
                0.5 0.5 0.5;
                ]);
            
            processNames = SimulationStructure.getShortProcessNames(sim);
            
            for i=1:numel(sim.processWholeCellModelIDs)
                text(0, i, processNames{modOrder(i)}, 'Parent', axesHandle, 'FontSize', 10, 'Rotation', 0,   'HorizontalAlignment', 'right', 'VerticalAlignment','middle', 'units', 'data');
                text(i, 0, processNames{modOrder(i)}, 'Parent', axesHandle, 'FontSize', 10, 'Rotation', 270, 'HorizontalAlignment', 'right', 'VerticalAlignment','middle', 'units', 'data');
            end
            
            set(axesHandle, 'XTick', 1:numel(processNames));
            set(axesHandle, 'YTick', 1:numel(processNames));
            set(axesHandle, 'XTickLabel', cell(size(processNames)));
            set(axesHandle, 'YTickLabel', cell(size(processNames)));
            
            ylabel('Process', 'FontSize', 18);
            xlabel('Process', 'FontSize', 18);
            set(axesHandle, 'XAxisLocation','top')
            axesPos = get(axesHandle', 'Position');
            set(axesHandle, 'Position', axesPos + [.11 -.1 -.1 -.1]);
            xLabelPos = get(get(axesHandle,'XLabel'), 'Position');
            set(get(axesHandle,'XLabel'), 'Position', xLabelPos + [0 -7 0]);
            yLabelPos = get(get(axesHandle,'YLabel'), 'Position');
            set(get(axesHandle,'YLabel'), 'Position', yLabelPos + [-6 0 0]);
        end
        
        function plotProcessMetabolites(sim, figHandle)
            import edu.stanford.covert.cell.sim.analysis.SimulationStructure;
            
            content = SimulationStructure.processMetabolites(sim);
            [~, modIdxs] = ismember(content(:, 1), sim.processWholeCellModelIDs);
            [~, metIdxs] = ismember(content(:, 3), sim.state('Metabolite').wholeCellModelIDs);
            mat = zeros(numel(sim.processWholeCellModelIDs), numel(sim.state('Metabolite').wholeCellModelIDs));
            mat(sub2ind(size(mat), modIdxs, metIdxs)) = 1;
            
            [mDNAIdxs, mRNAIdxs, mProtIdxs, mOtherIdxs] = SimulationStructure.processIndexs(sim);
            order = [mDNAIdxs; mRNAIdxs; mProtIdxs; mOtherIdxs];
            
            processWholeCellModelIDs = {sim.processWholeCellModelIDs{order}}'; %#ok<CCAT1>
            processNames = SimulationStructure.getShortProcessNames(sim);
            processNames = processNames(order(end:-1:1));
            mat = mat(order, :);
            
            %plot
            if ~exist('figHandle', 'var')
                [axesHandle, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            else
                axesHandle = subplot(1, 1, 1, 'Parent', figHandle);
            end
            
            image(mat+1, 'Parent', axesHandle);
            colormap(axesHandle, [
                1 1 1;
                0 0 0]);
            
            set(axesHandle, 'XTick', []);
            set(axesHandle, 'YTick', 1:numel(processWholeCellModelIDs))
            set(axesHandle, 'YTickLabel', cell(size(processNames)));
            
            ylabel('Process', 'FontSize', 28, 'Units','normalized');
            xlabel('Metabolite', 'FontSize', 28);
            set(axesHandle, 'XAxisLocation','top')
            axesPos = get(axesHandle', 'Position');
            set(axesHandle, 'Position', axesPos + [.15 -.10 -.08 -.08]);
            xLabelPos = get(get(axesHandle,'XLabel'), 'Position');
            set(get(axesHandle,'XLabel'), 'Position', xLabelPos + [0 -5.5 0]);
            yLabelPos = get(get(axesHandle,'YLabel'), 'Position');
            set(get(axesHandle,'YLabel'), 'Position', yLabelPos + [-.27 0 0]);
            
            %process labels
            axisPos = get(axesHandle,'position');
            props = {'FontSize', 10, 'HorizontalAlignment', 'right', 'VerticalAlignment','middle', 'units', 'normalized', 'Parent', axesHandle};
            for i = 1:numel(processWholeCellModelIDs)
                text(-axisPos(3)/100, (i+0.1)/(numel(order)+1), processNames{i}, props{:});
            end
            
            %process  category labels
            props = {'Rotation', 90, 'FontSize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'units', 'normalized', 'Parent', axesHandle};
            text(-0.27, 1/(numel(order)+1)*(1/2 + numel(mOtherIdxs)/2), 'Other', props{:});
            text(-0.27, 1/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs)/2), 'Protein', props{:});
            text(-0.27, 1/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + numel(mRNAIdxs)/2), 'RNA', props{:});
            text(-0.27, 1/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + numel(mRNAIdxs) + numel(mDNAIdxs)/2), 'DNA', props{:});
            
            %process category brackets
            props = {'LineWidth', 1};
            margin = 0.2;
            tickW = 0.015;
            yTickW = -0.27+0.03;
            
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; 0]),     axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + [margin; numel(mOtherIdxs)-margin]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + [margin; margin]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + [numel(mOtherIdxs)-margin; numel(mOtherIdxs)-margin]), props{:});
            
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; 0]),     axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + [margin; numel(mProtIdxs)-margin]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + [margin; margin]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + [numel(mProtIdxs)-margin; numel(mProtIdxs)-margin]), props{:});
            
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; 0]),     axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + [margin; numel(mRNAIdxs)-margin]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + [margin; margin]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + [numel(mRNAIdxs)-margin; numel(mRNAIdxs)-margin]), props{:});
            
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; 0]),     axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + numel(mRNAIdxs) + [margin; numel(mDNAIdxs)-margin]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + numel(mRNAIdxs) + [margin; margin]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + numel(mRNAIdxs) + [numel(mDNAIdxs)-margin; numel(mDNAIdxs)-margin]), props{:});
        end
        
        function plotProcessSharedMetabolites(sim, figHandle)
            import edu.stanford.covert.cell.sim.analysis.SimulationStructure;
            
            content = SimulationStructure.processMetabolites(sim);
            [~, modIdxs] = ismember(content(:, 1), sim.processWholeCellModelIDs);
            [~, metIdxs] = ismember(content(:, 3), sim.state('Metabolite').wholeCellModelIDs);
            mat = zeros(numel(sim.processWholeCellModelIDs), numel(sim.state('Metabolite').wholeCellModelIDs));
            mat(sub2ind(size(mat), modIdxs, metIdxs)) = 1;
            
            idxs = find(sum(mat, 1) > 4);
            metaboliteWholeCellModelIDs = sim.state('Metabolite').wholeCellModelIDs(idxs);
            mat = mat(:, idxs);
            
            [~,idxs] = ismember({
                'ATP'; 'ADP'; 'AMP'; 'GTP'; 'GDP';
                'PI'; 'PPI';
                'H'; 'H2O';
                'ALA'; 'ARG'; 'ASN'; 'ASP'; 'CYS'; 'FMET'; 'GLN'; 'GLU'; 'GLY'; 'HIS'; 'ILE'; 'LEU'; 'LYS'; 'MET'; 'PHE'; 'PRO'; 'SER'; 'THR'; 'TRP'; 'TYR'; 'VAL'}, ...
                metaboliteWholeCellModelIDs);
            metaboliteWholeCellModelIDs = metaboliteWholeCellModelIDs(idxs);
            mat = mat(:, idxs);
            
            nucleicAcidIndexs = (1:5)';
            phosphateIndexs = (6:7)';
            waterIndexs = (8:9)';
            aminoAcidIndexs = (10:30)';
            
            [mDNAIdxs, mRNAIdxs, mProtIdxs, mOtherIdxs] = SimulationStructure.processIndexs(sim);
            order = [mDNAIdxs; mRNAIdxs; mProtIdxs; mOtherIdxs];
            
            processWholeCellModelIDs = {sim.processWholeCellModelIDs{order}}'; %#ok<CCAT1>
            processNames = SimulationStructure.getShortProcessNames(sim);
            processNames = processNames(order(end:-1:1));
            mat = mat(order, :);
            
            %plot
            if ~exist('figHandle', 'var')
                [axesHandle, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            else
                axesHandle = subplot(1, 1, 1, 'Parent', figHandle);
            end
            
            image(mat+1, 'Parent', axesHandle);
            colormap(axesHandle, [
                1 1 1;
                0 0 0]);
            
            set(axesHandle, 'XTick', []);
            set(axesHandle, 'YTick', 1:numel(processWholeCellModelIDs))
            set(axesHandle, 'YTickLabel', cell(size(processNames)));
            
            ylabel('Process', 'FontSize', 28, 'Units','normalized');
            xlabel('Metabolite', 'FontSize', 28);
            set(axesHandle, 'XAxisLocation','top')
            axesPos = get(axesHandle', 'Position');
            set(axesHandle, 'Position', axesPos + [.15 -.10 -.08 -.08]);
            xLabelPos = get(get(axesHandle,'XLabel'), 'Position');
            set(get(axesHandle,'XLabel'), 'Position', xLabelPos + [0 -5.5 0]);
            yLabelPos = get(get(axesHandle,'YLabel'), 'Position');
            set(get(axesHandle,'YLabel'), 'Position', yLabelPos + [-.27 0 0]);
            
            %metabolite labels
            metIDs = metaboliteWholeCellModelIDs(sum(mat, 1) > 4);
            for i = 1:numel(metIDs)
                x = 1/2 + find(strcmp(metaboliteWholeCellModelIDs, metIDs{i}));
                text(x, 0, metaboliteWholeCellModelIDs{i}, 'Parent', axesHandle, 'FontSize', 10, 'Rotation', 270, 'HorizontalAlignment', 'right', 'VerticalAlignment','middle', 'units', 'data');
            end
            
            %metabolite  category labels
            props = {'FontSize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'units', 'normalized', 'Parent', axesHandle};
            text(1/(numel(metaboliteWholeCellModelIDs)+1)*(1/2 + numel(nucleicAcidIndexs)/2),1.19, 'NXP', props{:});
            text(1/(numel(metaboliteWholeCellModelIDs)+1)*(1/2 + numel(nucleicAcidIndexs) + numel(phosphateIndexs)/2), 1.19, 'Pi', props{:});
            text(1/(numel(metaboliteWholeCellModelIDs)+1)*(1/2 + numel(nucleicAcidIndexs) + numel(phosphateIndexs) + numel(waterIndexs)/2), 1.19, 'H', props{:});
            text(1/(numel(metaboliteWholeCellModelIDs)+1)*(1/2 + numel(nucleicAcidIndexs) + numel(phosphateIndexs) + numel(waterIndexs) + numel(aminoAcidIndexs)/2), 1.19, 'Amino Acid', props{:});
            
            %metabolite category brackets
            props = {'LineWidth', 1};
            margin = 0.2;
            tickW = -0.015;
            yTickW = 1.15;
            axisPos = get(axesHandle,'position');
            
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)/(numel(metaboliteWholeCellModelIDs)+1)*(1/2 + [margin; numel(nucleicAcidIndexs)-margin]), axisPos(2) + axisPos(4)*(yTickW + [0; 0]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)/(numel(metaboliteWholeCellModelIDs)+1)*(1/2 + [margin; margin]), axisPos(2) + axisPos(4)*(yTickW + [0; tickW]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)/(numel(metaboliteWholeCellModelIDs)+1)*(1/2 + [numel(nucleicAcidIndexs)-margin; numel(nucleicAcidIndexs)-margin]), axisPos(2) + axisPos(4)*(yTickW + [0; tickW]), props{:});
            
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)/(numel(metaboliteWholeCellModelIDs)+1)*(1/2 + numel(nucleicAcidIndexs) + [margin; numel(phosphateIndexs)-margin]), axisPos(2) + axisPos(4)*(yTickW + [0; 0]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)/(numel(metaboliteWholeCellModelIDs)+1)*(1/2 + numel(nucleicAcidIndexs) + [margin; margin]), axisPos(2) + axisPos(4)*(yTickW + [0; tickW]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)/(numel(metaboliteWholeCellModelIDs)+1)*(1/2 + numel(nucleicAcidIndexs) + [numel(phosphateIndexs)-margin; numel(phosphateIndexs)-margin]), axisPos(2) + axisPos(4)*(yTickW + [0; tickW]), props{:});
            
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)/(numel(metaboliteWholeCellModelIDs)+1)*(1/2 + numel(nucleicAcidIndexs) + numel(phosphateIndexs) + [margin; numel(waterIndexs)-margin]), axisPos(2) + axisPos(4)*(yTickW + [0; 0]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)/(numel(metaboliteWholeCellModelIDs)+1)*(1/2 + numel(nucleicAcidIndexs) + numel(phosphateIndexs) + [margin; margin]), axisPos(2) + axisPos(4)*(yTickW + [0; tickW]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)/(numel(metaboliteWholeCellModelIDs)+1)*(1/2 + numel(nucleicAcidIndexs) + numel(phosphateIndexs) + [numel(waterIndexs)-margin; numel(waterIndexs)-margin]), axisPos(2) + axisPos(4)*(yTickW + [0; tickW]), props{:});
            
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)/(numel(metaboliteWholeCellModelIDs)+1)*(1/2 + numel(nucleicAcidIndexs) + numel(phosphateIndexs) + numel(waterIndexs) + [margin; numel(aminoAcidIndexs)-margin]), axisPos(2) + axisPos(4)*(yTickW + [0; 0]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)/(numel(metaboliteWholeCellModelIDs)+1)*(1/2 + numel(nucleicAcidIndexs) + numel(phosphateIndexs) + numel(waterIndexs) + [margin; margin]), axisPos(2) + axisPos(4)*(yTickW + [0; tickW]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)/(numel(metaboliteWholeCellModelIDs)+1)*(1/2 + numel(nucleicAcidIndexs) + numel(phosphateIndexs) + numel(waterIndexs) + [numel(aminoAcidIndexs)-margin; numel(aminoAcidIndexs)-margin]), axisPos(2) + axisPos(4)*(yTickW + [0; tickW]), props{:});
            
            %process labels
            props = {'FontSize', 10, 'HorizontalAlignment', 'right', 'VerticalAlignment','middle', 'units', 'normalized', 'Parent', axesHandle};
            for i = 1:numel(processWholeCellModelIDs)
                text(-axisPos(3)/100, (i+0.1)/(numel(order)+1), processNames{i}, props{:});
            end
            
            %process  category labels
            props = {'Rotation', 90, 'FontSize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'units', 'normalized', 'Parent', axesHandle};
            text(-0.27, 1/(numel(order)+1)*(1/2 + numel(mOtherIdxs)/2), 'Other', props{:});
            text(-0.27, 1/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs)/2), 'Protein', props{:});
            text(-0.27, 1/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + numel(mRNAIdxs)/2), 'RNA', props{:});
            text(-0.27, 1/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + numel(mRNAIdxs) + numel(mDNAIdxs)/2), 'DNA', props{:});
            
            %process category brackets
            props = {'LineWidth', 1};
            margin = 0.2;
            tickW = 0.015;
            yTickW = -0.27+0.03;
            
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; 0]),     axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + [margin; numel(mOtherIdxs)-margin]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + [margin; margin]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + [numel(mOtherIdxs)-margin; numel(mOtherIdxs)-margin]), props{:});
            
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; 0]),     axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + [margin; numel(mProtIdxs)-margin]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + [margin; margin]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + [numel(mProtIdxs)-margin; numel(mProtIdxs)-margin]), props{:});
            
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; 0]),     axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + [margin; numel(mRNAIdxs)-margin]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + [margin; margin]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + [numel(mRNAIdxs)-margin; numel(mRNAIdxs)-margin]), props{:});
            
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; 0]),     axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + numel(mRNAIdxs) + [margin; numel(mDNAIdxs)-margin]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + numel(mRNAIdxs) + [margin; margin]), props{:});
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + numel(mRNAIdxs) + [numel(mDNAIdxs)-margin; numel(mDNAIdxs)-margin]), props{:});
        end
        
        function plotProcessGeneProducts(sim, figHandle)
            import edu.stanford.covert.cell.sim.analysis.SimulationStructure;
            
            %process names
            processNames = SimulationStructure.getShortProcessNames(sim);
            
            %order processes
            [mDNAIdxs, mRNAIdxs, mProtIdxs, mOtherIdxs] = SimulationStructure.processIndexs(sim);
            
            order = [mDNAIdxs; mRNAIdxs; mProtIdxs; mOtherIdxs];
            order = order(end:-1:1);
            if numel(order) < numel(sim.processWholeCellModelIDs)
                throw(MException('SimulationStructure:error', 'Not all processes included'));
            end
            
            %count numbers of genes used in each process
            content = SimulationStructure.processGeneProducts(sim);
            [~, mIdxs] = ismember(content(:,1), sim.processWholeCellModelIDs);
            mCnts = histc(mIdxs, 1:numel(sim.processWholeCellModelIDs));
            
            %bar graph
            if ~exist('figHandle', 'var')
                [axesHandle, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            else
                axesHandle = subplot(1, 1, 1, 'Parent', figHandle);
            end
            
            hold on;
            colors = SimulationStructure.barColors();
            %barh((1:numel(sim.processWholeCellModelIDs))', mCnts(order));
            barh((1:numel(mOtherIdxs))', mCnts(mOtherIdxs(end:-1:1)), 'facecolor', colors(1, :), 'Parent', axesHandle);
            barh((1:numel(mProtIdxs))' + numel(mOtherIdxs), mCnts(mProtIdxs(end:-1:1)), 'facecolor', colors(2, :), 'Parent', axesHandle);
            barh((1:numel(mRNAIdxs))' + numel(mOtherIdxs) + numel(mProtIdxs), mCnts(mRNAIdxs(end:-1:1)), 'facecolor', colors(3, :), 'Parent', axesHandle);
            barh((1:numel(mDNAIdxs))' + numel(mOtherIdxs) + numel(mProtIdxs) + numel(mRNAIdxs), mCnts(mDNAIdxs(end:-1:1)), 'facecolor', colors(4, :), 'Parent', axesHandle);
            
            ylim([0 numel(sim.processWholeCellModelIDs)+1]);
            xlim([1 max(mCnts)*1.1]);
            set(axesHandle, 'YTick', []);
            set(axesHandle, 'YTickLabel', []);
            set(axesHandle, 'XScale', 'log');
            set(axesHandle, 'box', 'on');
            
            %axis labels
            ylabel('Process', 'FontSize', 18);
            xlabel('Genes', 'FontSize', 18);
            yLabelMargin  = 0.10;
            
            axisPos = get(axesHandle, 'Position');
            set(axesHandle, 'Position', axisPos + [1.7*yLabelMargin 0 -yLabelMargin 0]);
            axisPos = get(axesHandle, 'Position');
            
            set(get(axesHandle, 'YLabel'), 'Units', 'normalized');
            yLabelPos = get(get(axesHandle,'YLabel'), 'Position');
            set(get(axesHandle, 'YLabel'), 'Position', yLabelPos + [-2.1*yLabelMargin/axisPos(3) 0 0]);
            
            %tick labels
            props = {'FontSize', 10, 'HorizontalAlignment', 'right', 'VerticalAlignment','middle', 'units', 'normalized', 'Parent', axesHandle};
            for i = 1:numel(sim.processWholeCellModelIDs)
                if ismember(order(i), mOtherIdxs)
                    color = colors(1, :);
                elseif ismember(order(i), mProtIdxs)
                    color = colors(2, :);
                elseif ismember(order(i), mRNAIdxs)
                    color = colors(3, :);
                else
                    color = colors(4, :);
                end
                text(-axisPos(3)/100, (i+0.1)/(numel(order)+1), processNames{order(i)}, props{:}, 'Color', color);
            end
            
            %category labels
            props = {'Rotation', 90, 'FontSize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'units', 'normalized', 'Parent', axesHandle};
            yLabelPos = get(get(axesHandle,'ylabel'), 'position');
            text(yLabelPos(1)+yLabelMargin - 0.06, 1/(numel(order)+1)*(1/2 + numel(mOtherIdxs)/2), 'Other', props{:}, 'Color', colors(1, :));
            text(yLabelPos(1)+yLabelMargin - 0.06, 1/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs)/2), 'Protein', props{:}, 'Color', colors(2, :));
            text(yLabelPos(1)+yLabelMargin - 0.06, 1/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + numel(mRNAIdxs)/2), 'RNA', props{:}, 'Color', colors(3, :));
            text(yLabelPos(1)+yLabelMargin - 0.06, 1/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + numel(mRNAIdxs) + numel(mDNAIdxs)/2), 'DNA', props{:}, 'Color', colors(4, :));
            
            %category brackets
            props = {'LineWidth', 1};
            margin = 0.2;
            tickW = 0.015;
            yTickW = yLabelPos(1) + yLabelMargin + 0.03 - 0.05;
            
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; 0]),     axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + [margin; numel(mOtherIdxs)-margin]), props{:}, 'Color', colors(1, :));
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + [margin; margin]), props{:}, 'Color', colors(1, :));
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + [numel(mOtherIdxs)-margin; numel(mOtherIdxs)-margin]), props{:}, 'Color', colors(1, :));
            
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; 0]),     axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + [margin; numel(mProtIdxs)-margin]), props{:}, 'Color', colors(2, :));
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + [margin; margin]), props{:}, 'Color', colors(2, :));
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + [numel(mProtIdxs)-margin; numel(mProtIdxs)-margin]), props{:}, 'Color', colors(2, :));
            
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; 0]),     axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + [margin; numel(mRNAIdxs)-margin]), props{:}, 'Color', colors(3, :));
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + [margin; margin]), props{:}, 'Color', colors(3, :));
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + [numel(mRNAIdxs)-margin; numel(mRNAIdxs)-margin]), props{:}, 'Color', colors(3, :));
            
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; 0]),     axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + numel(mRNAIdxs) + [margin; numel(mDNAIdxs)-margin]), props{:}, 'Color', colors(4, :));
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + numel(mRNAIdxs) + [margin; margin]), props{:}, 'Color', colors(4, :));
            annotation(figHandle, 'line', axisPos(1) + axisPos(3)*(yTickW + [0; tickW]), axisPos(2) + axisPos(4)/(numel(order)+1)*(1/2 + numel(mOtherIdxs) + numel(mProtIdxs) + numel(mRNAIdxs) + [numel(mDNAIdxs)-margin; numel(mDNAIdxs)-margin]), props{:}, 'Color', colors(4, :));
        end
        
        function processNames = getShortProcessNames(sim)
            processNames = cellfun(@(x) x.name, sim.processes, 'UniformOutput', false);
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_ChromosomeCondensation')}     = 'Condensation';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_ChromosomeSegregation')}      = 'Segregation';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_Cytokinesis')}                = 'Cytokinesis';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_DNADamage')}                  = 'Damage';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_DNARepair')}                  = 'Repair';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_FtsZPolymerization')}         = 'FtsZ';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_HostInteraction')}            = 'Host Interaction';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_MacromolecularComplexation')} = 'Complexation';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_Metabolism')}                 = 'Metabolism';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_ProteinActivation')}          = 'Activation';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_ProteinDecay')}               = 'Degradation';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_ProteinFolding')}             = 'Folding';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_ProteinModification')}        = 'Modification';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_ProteinProcessingI')}         = 'Processing I';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_ProteinProcessingII')}        = 'Processing II';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_ProteinTranslocation')}       = 'Translocation';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_RNADecay')}                   = 'Degradation';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_RNAModification')}            = 'Modification';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_RNAProcessing')}              = 'Processing';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_Replication')}                = 'Replication';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_ReplicationInitiation')}      = 'Rep Init';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_RibosomeAssembly')}           = 'Ribosome';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_TerminalOrganelleAssembly')}  = 'Term Org';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_Transcription')}              = 'Transcription';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_TranscriptionalRegulation')}  = 'Trans Reg';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_Translation')}                = 'Translation';
            processNames{strcmp(sim.processWholeCellModelIDs, 'Process_tRNAAminoacylation')}         = 'Aminoacylation';
        end
        
        function [mDNAIdxs, mRNAIdxs, mProtIdxs, mOtherIdxs] = processIndexs(sim)
            mDNAIdxs = find(ismember(sim.processWholeCellModelIDs, {
                'Process_ChromosomeCondensation';
                'Process_ChromosomeSegregation';
                'Process_DNADamage';
                'Process_DNARepair';
                'Process_Replication';
                'Process_ReplicationInitiation';
                'Process_DNASupercoiling';
                'Process_TranscriptionalRegulation'}));
            mRNAIdxs = find(ismember(sim.processWholeCellModelIDs, {
                'Process_RNADecay';
                'Process_RNAModification';
                'Process_RNAProcessing';
                'Process_Transcription';
                'Process_tRNAAminoacylation'}));
            mProtIdxs = find(ismember(sim.processWholeCellModelIDs, {
                'Process_MacromolecularComplexation';
                'Process_ProteinActivation';
                'Process_ProteinDecay';
                'Process_ProteinFolding';
                'Process_ProteinModification';
                'Process_ProteinProcessingI';
                'Process_ProteinProcessingII';
                'Process_ProteinTranslocation';
                'Process_RibosomeAssembly';
                'Process_TerminalOrganelleAssembly';
                'Process_Translation'}));
            mOtherIdxs = find(ismember(sim.processWholeCellModelIDs, {
                'Process_Cytokinesis';
                'Process_FtsZPolymerization';
                'Process_HostInteraction';
                'Process_Metabolism'}));
        end
        
        function colors = barColors()
            colors1 = [
                179 162 199;
                217 150 148;
                195 214 155;
                149 179 215;
                250 192 144;]/255;
            colors2 = [
                96 74 123;
                149 55 53;
                119 147 60;
                55 96 146;
                228 108 10;]/255;
            colors = (colors1+colors2)/2;
        end
    end
end