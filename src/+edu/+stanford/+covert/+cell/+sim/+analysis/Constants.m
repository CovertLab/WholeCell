%Constants
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/12/2010
classdef Constants
    methods (Static)
        function run(sim, fileName)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.PrintUtil;
            
            metaClass = meta.class.fromName('edu.stanford.covert.cell.sim.analysis.Constants');
            
            %% printing
            for i = 1:numel(metaClass.Methods)
                if numel(metaClass.Methods{i}.Name) <= 5 || ~strcmp('print', metaClass.Methods{i}.Name(1:5))
                    continue;
                end
                
                func = str2func(['edu.stanford.covert.cell.sim.analysis.Constants.' metaClass.Methods{i}.Name]);
                [content, colLabels, indentation] = func(sim);
                if nargin == 1
                    PrintUtil.printToStdIO(content, colLabels, struct('indentation', indentation));
                else
                    PrintUtil.printToFile(content, colLabels, [fileName '.xls'], metaClass.Methods{i}.Name(6:min(36, end)), struct('indentation', indentation));
                end
            end
            
            %% plotting
            [axesHandle, figHandle] = PlotUtil.newAxesHandle();
            
            for i = 1:numel(metaClass.Methods)
                if numel(metaClass.Methods{i}.Name) <= 4 || ~strcmp('plot', metaClass.Methods{i}.Name(1:4))
                    continue;
                end
                
                cla(axesHandle);
                func = str2func(['edu.stanford.covert.cell.sim.analysis.Constants.' metaClass.Methods{i}.Name]);
                if nargin == 1
                    func(sim, PlotUtil.newAxesHandle());
                else
                    func(sim, axesHandle, false);
                    saveas(figHandle, [fileName '-' metaClass.Methods{i}.Name(5:end) '.pdf']);
                end
            end
            
            close(figHandle);
        end
    end
    
    %printing
    methods (Static)
        function [content, colLabels, indentation] = printConstants(sim)
            %import classes
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            %initialize
            colLabels = {'Constant Name', 'State/Process ID', 'State/Process Name', 'Parameter', 'Fixed', 'class', 'size', 'numel', 'nnz'};
            content = cell(0, numel(colLabels));
            
            %states
            for i = 1:numel(sim.states)
                s = sim.states{i};
                content = [content;
                    Constants.objectConstants(s, s.fixedConstantNames,  s.wholeCellModelID, s.name, 'Y');
                    Constants.objectConstants(s, s.fittedConstantNames, s.wholeCellModelID, s.name, 'N')]; %#ok<AGROW>
            end
            
            %processes
            for i = 1:numel(sim.processes)
                m = sim.processes{i};
                content = [content;
                    Constants.objectConstants(m, m.fixedConstantNames,  m.wholeCellModelID, m.name, 'Y');
                    Constants.objectConstants(m, m.fittedConstantNames, m.wholeCellModelID, m.name, 'N')]; %#ok<AGROW>
            end
            
            indentation = zeros(size(content, 1), 1);
        end
        
        function content = objectConstants(obj, propertyNames, processID, processName, fitted)
            content = cell(numel(propertyNames), 9);
            fieldNames = fieldnames(obj);
            for i = 1:numel(propertyNames)
                if ~any(strcmp(fieldNames, propertyNames{i}))
                    warning('WholeCell:warning', '%s has no property %s', class(obj), propertyNames{i});
                    continue;
                end
                propVal = obj.(propertyNames{i});
                
                if ismember(propertyNames{i}, obj.parameterNames)
                    parameter = 'Y';
                else
                    parameter = 'N';
                end
                
                siz = size(propVal);
                siz = [sprintf('%d x ', siz(1:end-1)) sprintf('%d', siz(end))];
                
                switch class(propVal)
                    case 'cell'
                        nnzVal = sum(cellfun(@(x) ~isempty(x), propVal));
                    otherwise
                        try
                            nnzVal = nnz(propVal);
                        catch %#ok<CTCH>
                            nnzVal = [];
                        end
                end
                
                content(i, :) = {propertyNames{i}, processID, processName, parameter, fitted, class(propVal), siz, numel(propVal), nnzVal};
            end
        end
        
        function [content, colLabels, indentation] = printExpectedExpressions(sim)
            r = sim.state('Rna');
            pc = sim.state('ProteinComplex');
            fitter = edu.stanford.covert.cell.sim.util.FitConstants(sim);
            
            %% Fitted expressions
            [freeRnas, freeMons, ~, ~, ~, totRnas, totMons, totCpxs] = ...
                fitter.calcMacromolecularCounts(fitter.constructParameterVectorFromSimulation());
            
            geneRNAs     = r.matureRNAGeneComposition * totRnas;
            geneFreeRNAs = r.matureRNAGeneComposition * freeRnas;
            geneMonomers = zeros(size(sim.gene.wholeCellModelIDs));
            geneMonomers(sim.gene.mRNAIndexs) = totMons;
            geneComplexs = (sum(pc.proteinComplexComposition, 3) > 0) * totCpxs;
            geneFreeMonomers = zeros(size(sim.gene.wholeCellModelIDs));
            geneFreeMonomers(sim.gene.mRNAIndexs) = freeMons;
            
            calc = [geneRNAs geneFreeRNAs geneMonomers geneFreeMonomers geneComplexs];
            
            %% Experimentally observed expressions
            [freeRnas, freeMons, ~, ~, ~, totRnas, totMons, totCpxs] = ...
                fitter.calcMacromolecularCounts(fitter.initializeFittedConstants());
            
            geneRNAs     = r.matureRNAGeneComposition * totRnas;
            geneFreeRNAs = r.matureRNAGeneComposition * freeRnas;
            geneMonomers = zeros(size(sim.gene.wholeCellModelIDs));
            geneMonomers(sim.gene.mRNAIndexs) = totMons;
            geneComplexs = (sum(pc.proteinComplexComposition, 3) > 0) * totCpxs;
            geneFreeMonomers = zeros(size(sim.gene.wholeCellModelIDs));
            geneFreeMonomers(sim.gene.mRNAIndexs) = freeMons;
            
            expt = [geneRNAs geneFreeRNAs geneMonomers geneFreeMonomers geneComplexs];
            
            %% format output
            colLabels = {'ID', 'Name', ...
                'Calc. Total RNA', 'Calc. Free RNA', 'Calc. Total Monomer', 'Calc. Free Monomer', 'Calc. Complex', ...
                'Expt. Total RNA', 'Expt. Free RNA', 'Expt. Total Monomer', 'Expt. Free Monomer', 'Expt. Complex'};
            content = [
                sim.gene.wholeCellModelIDs sim.gene.names num2cell(calc) num2cell(expt)
                ];
            indentation = zeros(numel(sim.gene.names), 1);
        end
        
        function [content, colLabels, indentation] = printCalculatedExperimentalComparison(sim)
            m = sim.state('Metabolite');
            r = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            
            aminoAcidIndexs   = m.aminoAcidIndexs;
            
            rnaBaseCnts = r.baseCounts(r.matureIndexs, m.nmpIndexs);
            rnaLengths  = r.lengths(r.matureIndexs);
            monAACnts   = pm.baseCounts(pm.matureIndexs, aminoAcidIndexs);
            monLens     = pm.lengths(pm.matureIndexs);
            
            rnaExpr  = r.expression(r.matureIndexs) / sum(r.expression(r.matureIndexs));
            geneExpr = r.geneExpression;
            monExpr  = geneExpr(sim.gene.mRNAIndexs) ./ pm.halfLives(pm.matureIndexs);
            monExpr  = monExpr / sum(monExpr);
            tRNAExpr = geneExpr(sim.gene.tRNAIndexs) / sum(geneExpr(sim.gene.tRNAIndexs));
            nmpComp  = sum(m.nmpComposition, 2);
            aaComp   = sum(m.aaComposition, 2);
            
            exptGeneExpr = r.expectedGeneExpression(:,1);
            exptMonExpr  = exptGeneExpr(sim.gene.mRNAIndexs) ./ pm.halfLives(pm.matureIndexs);
            exptMonExpr  = exptMonExpr / sum(exptMonExpr);
            exptTRNAExpr = exptGeneExpr(sim.gene.tRNAIndexs) / sum(exptGeneExpr(sim.gene.tRNAIndexs));
            exptNMPComp  = m.experimentalNMPComposition;
            exptAAComp   = m.experimentalAAComposition;
            
            %initialize
            content = cell(0, 10);
            colLabels = {
                'Quantity';
                'Expt. Produced'; 'Expt. Expressed'; 'Expt.  Difference'; 'Expt.  Ratio';
                'Sim. Prod'; 'Sim. Expressed'; 'Sim. Difference'; 'Sim. Ratio'
                }';
            
            %content
            content = [content;{
                0 'Deviation from experimental values'  [] [] [] [] [] [] [] []
                1 'Gene expression'    norm(geneExpr - exptGeneExpr, 'fro') / norm(exptGeneExpr, 'fro') [] [] [] [] [] [] []
                1 'tRNA expression'    norm(tRNAExpr - exptTRNAExpr, 'fro') / norm(exptTRNAExpr, 'fro') [] [] [] [] [] [] []
                1 'Monomer expression' norm(monExpr  - exptMonExpr,  'fro') / norm(exptMonExpr,  'fro') [] [] [] [] [] [] []
                1 'NMP production'     norm(nmpComp  - exptNMPComp,  'fro') / norm(exptNMPComp,  'fro') [] [] [] [] [] [] []
                1 'AA production'      norm(aaComp   - exptAAComp,   'fro') / norm(exptAAComp,   'fro') [] [] [] [] [] [] []
                }];
            
            content = [content;{
                0 'NMPs' [] [] [] [] [] [] [] []
                }];
            for i = 1:length(m.nmpIndexs)
                content = [content;{
                    1
                    m.wholeCellModelIDs{m.nmpIndexs(i)}
                    exptNMPComp(i)
                    NaN
                    NaN
                    NaN
                    nmpComp(i)
                    rnaBaseCnts(:, i)' * rnaExpr / (rnaLengths' * rnaExpr)
                    nmpComp(i) - rnaBaseCnts(:, i)' * rnaExpr / (rnaLengths' * rnaExpr)
                    nmpComp(i) / (rnaBaseCnts(:, i)' * rnaExpr / (rnaLengths' * rnaExpr))
                    }']; %#ok<AGROW>
            end
            
            content = [content;{
                0 'Amino acids' [] [] [] [] [] [] [] []
                }];
            for i = 1:length(aminoAcidIndexs)
                content = [content;{
                    1
                    m.wholeCellModelIDs{aminoAcidIndexs(i)}
                    exptAAComp(i)
                    monAACnts(:,i)' * exptMonExpr / (monLens' * exptMonExpr)
                    exptAAComp(i) - monAACnts(:,i)' * exptMonExpr / (monLens' * exptMonExpr)
                    exptAAComp(i) / (monAACnts(:,i)' * exptMonExpr / (monLens' * exptMonExpr))
                    aaComp(i)
                    monAACnts(:,i)' * monExpr / (monLens' * monExpr)
                    aaComp(i) - monAACnts(:,i)' * monExpr / (monLens' * monExpr)
                    aaComp(i) / (monAACnts(:,i)' * monExpr / (monLens' * monExpr))
                    }']; %#ok<AGROW>
            end
            
            %format output
            indentation = cell2mat(content(:, 1));
            content = content(:, 2:end);
        end
        
        function [content, colLabels, indentation] = printBiomassNormalization(sim)
            import edu.stanford.covert.util.ConstantUtil;
            
            m = sim.state('Metabolite');
            metMWs = m.molecularWeights / ConstantUtil.nAvogadro;
            cllInitDryWt = sim.state('Mass').cellInitialDryWeight;
            
            %initialize
            content = cell(0, 4);
            colLabels = {'Quantity', 'Value', 'Units'};
            
            %content
            content = [content;{
                0 'Metabolism' [] [];
                1 'Metabolism production' sum(m.metabolismProduction, 2)' * metMWs / cllInitDryWt 'cells'
                1 'Biomass production'    sum(m.biomassProduction, 2)'    * metMWs / cllInitDryWt 'cells'
                1 'Byproducts'            sum(m.byproducts, 2)'           * metMWs / cllInitDryWt 'cells'
                }];
            
            %format output
            indentation = cell2mat(content(:, 1));
            content = content(:, 2:end);
        end
    end
    
    %scatter plots of experimental vs calculated constants
    methods (Static = true)
        function plotExperimental_Vs_Calculated_Gene_Expression(sim, axesHandle, showLabelsAndLegend)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            r = sim.state('Rna');
            
            calculated = r.geneExpression;
            experimental = r.expectedGeneExpression(:, 1);
            hold(axesHandle, 'on');
            h = zeros(3, 1);
            h(1) = plot(axesHandle, calculated(sim.gene.mRNAIndexs), experimental(sim.gene.mRNAIndexs), 'r.', 'MarkerSize', 24);
            h(2) = plot(axesHandle, calculated(sim.gene.rRNAIndexs), experimental(sim.gene.rRNAIndexs), 'g.', 'MarkerSize', 24);
            h(3) = plot(axesHandle, calculated(sim.gene.sRNAIndexs), experimental(sim.gene.sRNAIndexs), 'c.', 'MarkerSize', 24);
            h(4) = plot(axesHandle, calculated(sim.gene.tRNAIndexs), experimental(sim.gene.tRNAIndexs), 'b.', 'MarkerSize', 24);            
            hold(axesHandle, 'off');
            Constants.formatAxesExperimentalCalculatedScatterPlot(axesHandle, calculated, experimental, {'Gene' 'Expression'});
            if nargin >= 3 && showLabelsAndLegend
                legend(h, {'mRNA', 'rRNA', 'sRNA', 'tRNA'}, 'Location', 'NorthWest');
                Constants.labelExperimentalCalculatedScatterPlot(axesHandle,...
                    calculated, experimental, ...
                    sim.gene.wholeCellModelIDs, 10);
            end
        end
        
        function plotExperimental_Vs_Calculated_Monomer_Expression(sim, axesHandle, showLabelsAndLegend)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            r = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            
            calculated = r.geneExpression(sim.gene.mRNAIndexs) ./ pm.halfLives(pm.matureIndexs);
            experimental = r.expectedGeneExpression(sim.gene.mRNAIndexs, sim.compartment.cytosolIndexs) ./ pm.halfLives(pm.matureIndexs);
            
            calculated = calculated / sum(calculated);
            experimental = experimental / sum(experimental);
            
            plot(axesHandle, calculated, experimental, 'r.', 'MarkerSize', 24);
            Constants.formatAxesExperimentalCalculatedScatterPlot(axesHandle, calculated, experimental, {'Protein' 'Expression'});
            if nargin >= 3 && showLabelsAndLegend
                Constants.labelExperimentalCalculatedScatterPlot(axesHandle,...
                    calculated, experimental, ...
                    sim.gene.wholeCellModelIDs(sim.gene.mRNAIndexs), 10);
            end
        end
        
        function plotExperimental_Vs_Calculated_TRNA_Expression(sim, axesHandle, showLabelsAndLegend)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            r = sim.state('Rna');
            
            calculated = r.geneExpression(sim.gene.tRNAIndexs);
            experimental = r.expectedGeneExpression(sim.gene.tRNAIndexs, sim.compartment.cytosolIndexs);
            
            calculated = calculated / sum(calculated);
            experimental = experimental / sum(experimental);
            
            plot(axesHandle, calculated, experimental, 'b.', 'MarkerSize', 24);
            Constants.formatAxesExperimentalCalculatedScatterPlot(axesHandle, calculated, experimental, {'tRNA' 'Expression'});
            if nargin >= 3 && showLabelsAndLegend
                Constants.labelExperimentalCalculatedScatterPlot(axesHandle,...
                    calculated, experimental, ...
                    sim.gene.wholeCellModelIDs(sim.gene.tRNAIndexs), 3);
            end
        end
        
        function plotExperimental_Vs_Calculated_NMP_Composition(sim, axesHandle, showLabelsAndLegend)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            m = sim.state('Metabolite');
            calculated = sum(m.nmpComposition, 2);
            experimental = m.experimentalNMPComposition;
            plot(axesHandle, calculated, experimental, 'b.', 'MarkerSize', 24);
            Constants.formatAxesExperimentalCalculatedScatterPlot(axesHandle, calculated, experimental, {'NMP' 'Composition'});
            if nargin >= 3 && showLabelsAndLegend
                Constants.labelExperimentalCalculatedScatterPlot(axesHandle,...
                    calculated, experimental, ...
                    m.wholeCellModelIDs(m.nmpIndexs));
            end
        end
        
        function plotExperimental_Vs_Calculated_AA_Composition(sim, axesHandle, showLabelsAndLegend)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            m = sim.state('Metabolite');
            calculated = sum(m.aaComposition, 2);
            experimental = m.experimentalAAComposition;
            plot(axesHandle, calculated, experimental, 'b.', 'MarkerSize', 24);
            Constants.formatAxesExperimentalCalculatedScatterPlot(axesHandle, calculated, experimental, {'AA' 'Composition'});
            if nargin >= 3 && showLabelsAndLegend
                Constants.labelExperimentalCalculatedScatterPlot(axesHandle,...
                    calculated, experimental, ...
                    m.wholeCellModelIDs(m.aminoAcidIndexs), 3);
            end
        end
        
        function plotExperimental_Vs_Calculated_RNA_Weight_Fractions(sim, axesHandle, showLabelsAndLegend)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            r = sim.state('Rna');
            
            calculated = r.weightFractions;
            experimental = r.expectedWeightFractions;
            plot(axesHandle, calculated, experimental, 'b.', 'MarkerSize', 24);
            Constants.formatAxesExperimentalCalculatedScatterPlot(axesHandle, calculated, experimental, {'RNA' 'Composition'});
            if nargin >= 3 && showLabelsAndLegend
                Constants.labelExperimentalCalculatedScatterPlot(axesHandle,...
                    calculated, experimental, {'mRNA', 'rRNA 5S', 'rRNA 16S', 'rRNA 23S', 'sRNA', 'tRNA'});
            end
        end
        
        function plotExperimental_Vs_Calculated_Gene_Decay_Rates(sim, axesHandle, showLabelsAndLegend)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            r = sim.state('Rna');
            
            calculated = r.geneDecayRates;
            experimental = r.expectedGeneDecayRates;
            tfs = isfinite(experimental);
            plot(axesHandle, experimental(tfs), calculated(tfs), 'b.', 'MarkerSize', 24);
            Constants.formatAxesExperimentalCalculatedScatterPlot(axesHandle, calculated(tfs), experimental(tfs), {'RNA Decay' 'Rates'});
            minVal = min([calculated(tfs); experimental(tfs)]);
            maxVal = max([calculated(tfs); experimental(tfs)]);
            axis(axesHandle, [minVal maxVal minVal maxVal]);
            axis(axesHandle, 'square');
            
            if nargin >= 3 && showLabelsAndLegend
                Constants.labelExperimentalCalculatedScatterPlot(axesHandle,...
                    calculated(tfs), experimental(tfs), sim.gene.wholeCellModelIDs, 10);
            end
        end
        
        function formatAxesExperimentalCalculatedScatterPlot(axesHandle, calculated, experimental, titleStr)
            minVal = min([calculated; experimental]);
            maxVal = max([calculated; experimental]);
            range = maxVal - minVal;
            minVal = max(0, minVal - range * 0);
            maxVal = min(1, maxVal + range * 0);
            line([0 1], [0 1], 'Color', 'k', 'LineStyle', ':');
            axis(axesHandle, [minVal maxVal minVal maxVal]);
            xlabel(axesHandle, 'Model', 'FontSize', 10);
            ylabel(axesHandle, 'Experiment', 'FontSize', 10);
            set(axesHandle, 'XTick', []);
            set(axesHandle, 'YTick', []);
            titleStr = [titleStr{1} ' ' titleStr{2}];
            title(axesHandle, titleStr, 'FontSize', 12);
            
            axis(axesHandle, 'square');
            figHandle = get(axesHandle, 'Parent');
            pos = get(figHandle, 'Position');
            set(figHandle, 'Position', [pos(1:2) max(pos(3:4)) max(pos(3:4))]);
        end
        
        function labelExperimentalCalculatedScatterPlot(axesHandle, calculated, experimental, labels, count)
            if ~exist('count', 'var')
                for i = 1:length(calculated)
                    if(calculated(i) < experimental(i))
                        text(calculated(i), ...
                            experimental(i), ...
                            labels{i}, ...
                            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom',...
                            'FontSize', 10, 'Interpreter', 'none', ...
                            'Parent', axesHandle);
                    else
                        text(calculated(i), ...
                            experimental(i), ...
                            labels{i}, ...
                            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top',...
                            'FontSize', 10, 'Interpreter', 'none', ...
                            'Parent', axesHandle);
                    end
                end
            else
                [~, idxs] = sort(calculated ./ experimental);
                for i = 1:count
                    text(calculated(idxs(i)), ...
                        experimental(idxs(i)), ...
                        labels{idxs(i)},...
                        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
                        'FontSize', 10, 'Interpreter', 'none', ...
                        'Parent', axesHandle);
                    text(calculated(idxs(end - i + 1)), ...
                        experimental(idxs(end - i + 1)), ...
                        labels{idxs(end - i + 1)},...
                        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top',...
                        'FontSize', 10, 'Interpreter', 'none', ...
                        'Parent', axesHandle);
                end
            end
        end
    end
    
    %plot experimental vs calculated constant ratios
    methods (Static = true)
        function plotExperimental_Calculated_Gene_Expression_Ratios(sim, axesHandle, showLabelsAndLegend)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            r = sim.state('Rna');
            
            x = 1:length(sim.gene.wholeCellModelIDs);
            expectedGeneExpression = r.expectedGeneExpression(:, 1);
            ratios = log2(r.geneExpression ./ expectedGeneExpression);
            
            hold(axesHandle,'on');
            h=zeros(3,1);
            h(1) = bar(x(sim.gene.mRNAIndexs), ratios(sim.gene.mRNAIndexs),...
                'Parent', axesHandle, 'FaceColor', 'r', 'EdgeColor', 'r');
            h(2) = bar(x(sim.gene.rRNAIndexs), ratios(sim.gene.rRNAIndexs),...
                'Parent', axesHandle, 'FaceColor', 'g', 'EdgeColor', 'g');
            h(3) = bar(x(sim.gene.sRNAIndexs), ratios(sim.gene.sRNAIndexs),...
                'Parent', axesHandle, 'FaceColor', 'c', 'EdgeColor', 'c');
            h(4) = bar(x(sim.gene.tRNAIndexs), ratios(sim.gene.tRNAIndexs),...
                'Parent', axesHandle, 'FaceColor', 'b', 'EdgeColor', 'b');            
            hold(axesHandle, 'off');
            
            Constants.formatAxesExperimentalCalculatedBarPlot(axesHandle, x, ratios, 'Gene expression');
            if nargin >= 3 && showLabelsAndLegend
                legend(h, {'mRNA', 'rRNA', 'sRNA', 'tRNA'});
                Constants.labelExperimentalCalculatedBarPlot(axesHandle, x, ratios, sim.gene.wholeCellModelIDs, 10);
            end
        end
        
        function plotExperimental_Calculated_Monomer_Expression_Ratios(sim, axesHandle, showLabelsAndLegend)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            r = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            
            calculated = r.geneExpression(sim.gene.mRNAIndexs) ./ pm.halfLives(pm.matureIndexs);
            experimental = r.expectedGeneExpression(sim.gene.mRNAIndexs, sim.compartment.cytosolIndexs) ./ pm.halfLives(pm.matureIndexs);
            
            calculated = calculated / sum(calculated);
            experimental = experimental / sum(experimental);
            
            x = 1:length(sim.gene.mRNAIndexs);
            ratios = log2(calculated ./ experimental);
            bar(x, ratios, 'FaceColor', 'r', 'Parent', axesHandle);
            Constants.formatAxesExperimentalCalculatedBarPlot(axesHandle, x, ratios, 'Monomer expression');
            if nargin >= 3 && showLabelsAndLegend
                Constants.labelExperimentalCalculatedBarPlot(axesHandle, x, ratios, ...
                    sim.gene.wholeCellModelIDs(sim.gene.mRNAIndexs), 10);
            end
        end
        
        function plotExperimental_Calculated_TRNA_Expression_Ratios(sim, axesHandle, showLabelsAndLegend)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            r = sim.state('Rna');
            
            calculatedExpression = r.geneExpression(sim.gene.tRNAIndexs);
            experimentalExpression = r.expectedGeneExpression(sim.gene.tRNAIndexs, sim.compartment.cytosolIndexs);
            
            calculatedExpression = calculatedExpression / sum(calculatedExpression);
            experimentalExpression = experimentalExpression / sum(experimentalExpression);
            
            x = 1:length(sim.gene.tRNAIndexs);
            ratios = log2(calculatedExpression ./ experimentalExpression);
            bar(x,ratios, 'FaceColor', 'b', 'Parent', axesHandle);
            Constants.formatAxesExperimentalCalculatedBarPlot(axesHandle, x, ratios, 'tRNA expression');
            if nargin >= 3 && showLabelsAndLegend
                Constants.labelExperimentalCalculatedBarPlot(axesHandle, x, ratios, ...
                    sim.gene.wholeCellModelIDs(sim.gene.tRNAIndexs), 3);
            end
        end
        
        function plotExperimental_Calculated_NMP_Composition_Ratios(sim, axesHandle, showLabelsAndLegend)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            m = sim.state('Metabolite');
            x = 1:length(m.nmpIndexs);
            ratios = log2(sum(m.nmpComposition, 2) ./ m.experimentalNMPComposition);
            bar(x, ratios, 'FaceColor', 'b', 'Parent', axesHandle);
            Constants.formatAxesExperimentalCalculatedBarPlot(axesHandle, x, ratios, 'NMP composition');
            if nargin >= 3 && showLabelsAndLegend
                Constants.labelExperimentalCalculatedBarPlot(axesHandle, x, ratios, ...
                    m.wholeCellModelIDs(m.nmpIndexs), 2);
            end
        end
        
        function plotExperimental_Calculated_AA_Composition_Ratios(sim, axesHandle, showLabelsAndLegend)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            m = sim.state('Metabolite');
            x = 1:length(m.aminoAcidIndexs);
            ratios = log2(sum(m.aaComposition, 2) ./ m.experimentalAAComposition);
            bar(x, ratios, 'FaceColor', 'b', 'Parent', axesHandle);
            Constants.formatAxesExperimentalCalculatedBarPlot(axesHandle, x, ratios, 'Amino acid composition');
            if nargin >= 3 && showLabelsAndLegend
                Constants.labelExperimentalCalculatedBarPlot(axesHandle, x, ratios, ...
                    m.wholeCellModelIDs(m.aminoAcidIndexs), 3);
            end
        end
        
        function plotExperimental_Calculated_RNA_Weight_Fraction_Ratios(sim, axesHandle, showLabelsAndLegend)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            r = sim.state('Rna');
            
            ratios = log2(r.weightFractions ./ r.expectedWeightFractions);
            x = 1:numel(ratios);
            bar(x, ratios, 'FaceColor', 'b', 'Parent', axesHandle);
            Constants.formatAxesExperimentalCalculatedBarPlot(axesHandle, x, ratios, 'RNA weight fraction');
            if nargin >= 3 && showLabelsAndLegend
                Constants.labelExperimentalCalculatedBarPlot(axesHandle, x, ratios, ...
                    {'mRNA', 'rRNA 5S', 'rRNA 16S', 'rRNA 23S', 'sRNA', 'tRNA'});
            end
        end
        
        function plotExperimental_Calculated_Gene_Decay_Rate_Ratios(sim, axesHandle, showLabelsAndLegend)
            import edu.stanford.covert.cell.sim.analysis.Constants;
            
            r = sim.state('Rna');
            
            x = 1:length(sim.gene.wholeCellModelIDs);
            ratios = log2(r.geneDecayRates ./ r.expectedGeneDecayRates);
            bar(x, ratios, 'FaceColor', 'b', 'Parent', axesHandle);
            Constants.formatAxesExperimentalCalculatedBarPlot(axesHandle, x, ratios, 'Gene decay rate');
            if nargin >= 3 && showLabelsAndLegend
                Constants.labelExperimentalCalculatedBarPlot(axesHandle, x, ratios, ...
                    sim.gene.wholeCellModelIDs, 10);
            end
        end
        
        function formatAxesExperimentalCalculatedBarPlot(axesHandle, x, ratios, xAxisLabel)
            if ~exist('xAxisLabel', 'var')
                xAxisLabel = 'Index';
            end
            
            minVal = min(ratios);
            maxVal = max(ratios);
            range = maxVal - minVal;
            minVal = minVal - range * 0.05;
            maxVal = maxVal + range * 0.05;
            
            if ~isnan(min(x)) && ~isnan(max(x))
                xlim(axesHandle, [min(x)-0.5 max(x)+0.5]);
            end
            if ~isnan(minVal) && ~isnan(maxVal)
                ylim(axesHandle, [minVal maxVal]);
            end
            
            xlabel(axesHandle, xAxisLabel, 'FontSize', 16);
            ylabel(axesHandle, 'log_2(Calculated / Experimental)', 'FontSize', 16);
        end
        
        function labelExperimentalCalculatedBarPlot(axesHandle, x, ratios, labels, count)
            if ~exist('count', 'var')
                for i = 1:length(x)
                    if(ratios(i)>0)
                        text(x(i), ratios(i), labels{i},...
                            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                            'FontSize', 10, 'Interpreter', 'none', ...
                            'Parent', axesHandle);
                    else
                        text(x(i), ratios(i), labels{i},...
                            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                            'FontSize', 10, 'Interpreter', 'none', ...
                            'Parent', axesHandle);
                    end
                end
            else
                [~, idxs] = sort(ratios);
                for i = 1:count
                    text(x(idxs(i)), ratios(idxs(i)), labels{idxs(i)},...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top',...
                        'FontSize', 10, 'Interpreter', 'none', ...
                        'Parent', axesHandle);
                    text(x(idxs(end-i+1)), ratios(idxs(end-i+1)), labels{idxs(end-i+1)},...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom',...
                        'FontSize', 10, 'Interpreter', 'none', ...
                        'Parent', axesHandle);
                end
            end
        end
    end
end