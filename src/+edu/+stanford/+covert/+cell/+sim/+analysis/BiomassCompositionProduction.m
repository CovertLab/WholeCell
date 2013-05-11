%BiomassCompositionProduction
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/12/2010
classdef BiomassCompositionProduction
    properties (Constant)
        nDryWeightFractions = 8;
        dryWeightFractionIndex_DNA          = 1;
        dryWeightFractionIndex_RNA          = 2;
        dryWeightFractionIndex_Protein      = 3;
        dryWeightFractionIndex_Lipid        = 4;
        dryWeightFractionIndex_Polyamine    = 5;
        dryWeightFractionIndex_Carbohydrate = 6;
        dryWeightFractionIndex_Vitamin      = 7;
        dryWeightFractionIndex_Ion          = 8;
    end
    
    methods (Static)
        function run(sim, fileName)
            import edu.stanford.covert.cell.sim.analysis.BiomassCompositionProduction;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.PrintUtil;
            
            [bmComp, bmProd, byProd] = edu.stanford.covert.cell.sim.analysis.BiomassCompositionProduction.getResourceRequirements(sim);
            
            %excel file
            [content, colLabels, indentation] = BiomassCompositionProduction.allMetabolites(sim, bmComp, bmProd, byProd);
            if nargin == 1
                PrintUtil.printToStdIO(content, colLabels, struct('indentation', indentation));
            else
                PrintUtil.printToFile(content, colLabels, [fileName '.xls'], 'All Metabolites', struct('indentation', indentation));
            end
            
            [content, colLabels, indentation] = BiomassCompositionProduction.metabolitesByWeightFraction(sim);
            if nargin == 1
                PrintUtil.printToStdIO(content, colLabels, struct('indentation', indentation));
            else
                PrintUtil.printToFile(content, colLabels, [fileName '.xls'], 'By Weight Fraction', struct('indentation', indentation));
            end
                
            % plots
            if nargin == 1
                BiomassCompositionProduction.plotWeightFractions(sim, PlotUtil.newAxesHandle());
                BiomassCompositionProduction.plotDNMPComposition(sim, PlotUtil.newAxesHandle(), bmComp);
                BiomassCompositionProduction.plotNMPComposition(sim, PlotUtil.newAxesHandle(), bmComp);
                BiomassCompositionProduction.plotAAComposition(sim, PlotUtil.newAxesHandle(), bmComp);
            else
                [axesHandle, figHandle] = PlotUtil.newAxesHandle();
                BiomassCompositionProduction.plotBiomassComposition(sim, axesHandle, bmComp);
                saveas(figHandle, [fileName '.pdf']);
                close(figHandle);
                
                [axesHandle, figHandle] = PlotUtil.newAxesHandle();
                BiomassCompositionProduction.plotWeightFractions(sim, axesHandle);
                saveas(figHandle, [fileName '-WeightFractions.pdf']);
                close(figHandle);
                
                [axesHandle, figHandle] = PlotUtil.newAxesHandle();
                BiomassCompositionProduction.plotDNMPComposition(sim, axesHandle, bmComp);
                saveas(figHandle, [fileName '-dNMPComposition.pdf']);
                close(figHandle);
                
                [axesHandle, figHandle] = PlotUtil.newAxesHandle();
                BiomassCompositionProduction.plotNMPComposition(sim, axesHandle, bmComp);
                saveas(figHandle, [fileName '-NMPComposition.pdf']);
                close(figHandle);
                
                [axesHandle, figHandle] = PlotUtil.newAxesHandle();
                BiomassCompositionProduction.plotAAComposition(sim, axesHandle, bmComp);
                saveas(figHandle, [fileName '-AAComposition.pdf']);                
                close(figHandle);
            end
        end
    end
    
    methods (Static)
        function [content, colLabels, indentation] = allMetabolites(sim, bmComp, bmProd, byProd)
            m = sim.state('Metabolite');
            
            colLabels = {'ID', 'Name', 'Composition', 'Production', 'Byproduction'};
            content = cell(0, 6);
            
            %dNMPs, dNTPs
            content = [content;
                {0 'dNMPs' '' sum(sum(bmComp(m.dnmpIndexs, :))) sum(sum(bmProd(m.dnmpIndexs, :))) sum(sum(byProd(m.dnmpIndexs, :)))}
                num2cell(ones(size(m.dnmpIndexs))) m.wholeCellModelIDs(m.dnmpIndexs) m.names(m.dnmpIndexs) num2cell(sum(bmComp(m.dnmpIndexs, :), 2)) num2cell(sum(bmProd(m.dnmpIndexs, :), 2)) num2cell(sum(byProd(m.dnmpIndexs, :), 2))
                ];
            content = [content;
                {0 'dNTPs' '' sum(sum(bmComp(m.dntpIndexs, :))) sum(sum(bmProd(m.dntpIndexs, :))) sum(sum(byProd(m.dntpIndexs, :)))}
                num2cell(ones(size(m.dntpIndexs))) m.wholeCellModelIDs(m.dntpIndexs) m.names(m.dntpIndexs) num2cell(sum(bmComp(m.dntpIndexs, :), 2)) num2cell(sum(bmProd(m.dntpIndexs, :), 2)) num2cell(sum(byProd(m.dntpIndexs, :), 2))
                ];
            
            %NMPs, NTPs
            content = [content;
                {0 'NMPs' '' sum(sum(bmComp(m.nmpIndexs, :))) sum(sum(bmProd(m.nmpIndexs, :))) sum(sum(byProd(m.nmpIndexs, :)))}
                num2cell(ones(size(m.nmpIndexs))) m.wholeCellModelIDs(m.nmpIndexs) m.names(m.nmpIndexs) num2cell(sum(bmComp(m.nmpIndexs, :), 2)) num2cell(sum(bmProd(m.nmpIndexs, :), 2)) num2cell(sum(byProd(m.nmpIndexs, :), 2))
                ];
            content = [content;
                {0 'NTPs' '' sum(sum(bmComp(m.ntpIndexs, :))) sum(sum(bmProd(m.ntpIndexs, :))) sum(sum(byProd(m.ntpIndexs, :)))}
                num2cell(ones(size(m.ntpIndexs))) m.wholeCellModelIDs(m.ntpIndexs) m.names(m.ntpIndexs) num2cell(sum(bmComp(m.ntpIndexs, :), 2)) num2cell(sum(bmProd(m.ntpIndexs, :), 2)) num2cell(sum(byProd(m.ntpIndexs, :), 2))
                ];
            
            %amino acids
            content = [content;
                {0 'Amino Acids' '' sum(sum(bmComp(m.aminoAcidIndexs, :))) sum(sum(bmProd(m.aminoAcidIndexs, :))) sum(sum(byProd(m.aminoAcidIndexs, :)))}
                num2cell(ones(size(m.aminoAcidIndexs))) m.wholeCellModelIDs(m.aminoAcidIndexs) m.names(m.aminoAcidIndexs) num2cell(sum(bmComp(m.aminoAcidIndexs, :), 2)) num2cell(sum(bmProd(m.aminoAcidIndexs, :), 2)) num2cell(sum(byProd(m.aminoAcidIndexs, :), 2))
                ];
            
            %lipids
            content = [content;
                {0 'Lipids' '' sum(sum(bmComp(m.lipidIndexs, :))) sum(sum(bmProd(m.lipidIndexs, :))) sum(sum(byProd(m.lipidIndexs, :)))}
                num2cell(ones(size(m.lipidIndexs))) m.wholeCellModelIDs(m.lipidIndexs) m.names(m.lipidIndexs) num2cell(sum(bmComp(m.lipidIndexs, :), 2)) num2cell(sum(bmProd(m.lipidIndexs, :), 2)) num2cell(sum(byProd(m.lipidIndexs, :), 2))
                ];
            
            %polyamines
            content = [content;
                {0 'Polyamines' '' sum(sum(bmComp(m.polyamineIndexs, :))) sum(sum(bmProd(m.polyamineIndexs, :))) sum(sum(byProd(m.polyamineIndexs, :)))}
                num2cell(ones(size(m.polyamineIndexs))) m.wholeCellModelIDs(m.polyamineIndexs) m.names(m.polyamineIndexs) num2cell(sum(bmComp(m.polyamineIndexs, :), 2)) num2cell(sum(bmProd(m.polyamineIndexs, :), 2)) num2cell(sum(byProd(m.polyamineIndexs, :), 2))
                ];
            
            %carbohydrates
            content = [content;
                {0 'Carbohydrates' '' sum(sum(bmComp(m.carbohydrateIndexs, :))) sum(sum(bmProd(m.carbohydrateIndexs, :))) sum(sum(byProd(m.carbohydrateIndexs, :)))}
                num2cell(ones(size(m.carbohydrateIndexs))) m.wholeCellModelIDs(m.carbohydrateIndexs) m.names(m.carbohydrateIndexs) num2cell(sum(bmComp(m.carbohydrateIndexs, :), 2)) num2cell(sum(bmProd(m.carbohydrateIndexs, :), 2)) num2cell(sum(byProd(m.carbohydrateIndexs, :), 2))
                ];
            
            %vitamins
            content = [content;
                {0 'Vitamins' '' sum(sum(bmComp(m.vitaminIndexs, :))) sum(sum(bmProd(m.vitaminIndexs, :))) sum(sum(byProd(m.vitaminIndexs, :)))}
                num2cell(ones(size(m.vitaminIndexs))) m.wholeCellModelIDs(m.vitaminIndexs) m.names(m.vitaminIndexs) num2cell(sum(bmComp(m.vitaminIndexs, :), 2)) num2cell(sum(bmProd(m.vitaminIndexs, :), 2)) num2cell(sum(byProd(m.vitaminIndexs, :), 2))
                ];
            
            %ions
            content = [content;
                {0 'Ions' '' sum(sum(bmComp(m.ionIndexs, :))) sum(sum(bmProd(m.ionIndexs, :))) sum(sum(byProd(m.ionIndexs, :)))}
                num2cell(ones(size(m.ionIndexs))) m.wholeCellModelIDs(m.ionIndexs) m.names(m.ionIndexs) num2cell(sum(bmComp(m.ionIndexs, :), 2)) num2cell(sum(bmProd(m.ionIndexs, :), 2)) num2cell(sum(byProd(m.ionIndexs, :), 2))
                ];
            
            %other
            otherIdxs = setdiff((1:numel(m.wholeCellModelIDs))', [
                m.dnmpIndexs
                m.dntpIndexs
                m.nmpIndexs
                m.ntpIndexs
                m.aminoAcidIndexs
                m.lipidIndexs
                m.polyamineIndexs
                m.carbohydrateIndexs
                m.vitaminIndexs
                m.ionIndexs]);
            
            otherIdxs = otherIdxs(any(bmComp(otherIdxs, :), 2) | any(bmProd(otherIdxs, :), 2) | any(byProd(otherIdxs, :), 2));
            
            content = [content;
                {0 'Other' '' sum(sum(bmComp(otherIdxs, :))) sum(sum(bmProd(otherIdxs, :))) sum(sum(byProd(otherIdxs, :)))}
                num2cell(ones(size(otherIdxs))) m.wholeCellModelIDs(otherIdxs) m.names(otherIdxs) num2cell(sum(bmComp(otherIdxs, :), 2)) num2cell(sum(bmProd(otherIdxs, :), 2)) num2cell(sum(byProd(otherIdxs, :), 2))
                ];
            
            %format output
            indentation = cell2mat(content(:, 1));
            content = content(:, 2:end);
        end
        
        function [content, colLabels, indentation] = metabolitesByWeightFraction(sim)
            import edu.stanford.covert.cell.sim.analysis.BiomassCompositionProduction;
            mass = sim.state('Mass');
            m = sim.state('Metabolite');
            
            colLabels = {'ID', 'Name', 'Composition', 'Production', 'Byproduction'};
            content = cell(0, 6);
            
            %DNA
            [bmComp, bmProd, byProd] = edu.stanford.covert.cell.sim.analysis.BiomassCompositionProduction.getResourceRequirements(sim);
            bmComp = bmComp * mass.dryWeightFractionDNA;
            bmProd = bmProd * mass.dryWeightFractionDNA;
            byProd = byProd * mass.dryWeightFractionDNA;
            idxs = find(any(bmComp, 2) | any(bmProd, 2) | any(byProd, 2));
            content = [content;
                {0 'DNA' '' sum(sum(bmComp(idxs, :))) sum(sum(bmProd(idxs, :))) sum(sum(byProd(idxs, :)))}
                num2cell(ones(size(idxs))) m.wholeCellModelIDs(idxs) m.names(idxs) num2cell(sum(bmComp(idxs, :), 2)) num2cell(sum(bmProd(idxs, :), 2)) num2cell(sum(byProd(idxs, :), 2))
                ];
            
            %RNA
            [bmComp, bmProd, byProd] = edu.stanford.covert.cell.sim.analysis.BiomassCompositionProduction.getResourceRequirements(sim);
            bmComp = bmComp * mass.dryWeightFractionRNA;
            bmProd = bmProd * mass.dryWeightFractionRNA;
            byProd = byProd * mass.dryWeightFractionRNA;
            idxs = find(any(bmComp, 2) | any(bmProd, 2) | any(byProd, 2));
            content = [content;
                {0 'RNA' '' sum(sum(bmComp(idxs, :))) sum(sum(bmProd(idxs, :))) sum(sum(byProd(idxs, :)))}
                num2cell(ones(size(idxs))) m.wholeCellModelIDs(idxs) m.names(idxs) num2cell(sum(bmComp(idxs, :), 2)) num2cell(sum(bmProd(idxs, :), 2)) num2cell(sum(byProd(idxs, :), 2))
                ];
            
            %Protein
            [bmComp, bmProd, byProd] = edu.stanford.covert.cell.sim.analysis.BiomassCompositionProduction.getResourceRequirements(sim);
            bmComp = bmComp * mass.dryWeightFractionProtein;
            bmProd = bmProd * mass.dryWeightFractionProtein;
            byProd = byProd * mass.dryWeightFractionProtein;
            idxs = find(any(bmComp, 2) | any(bmProd, 2) | any(byProd, 2));
            content = [content;
                {0 'Protein' '' sum(sum(bmComp(idxs, :))) sum(sum(bmProd(idxs, :))) sum(sum(byProd(idxs, :)))}
                num2cell(ones(size(idxs))) m.wholeCellModelIDs(idxs) m.names(idxs) num2cell(sum(bmComp(idxs, :), 2)) num2cell(sum(bmProd(idxs, :), 2)) num2cell(sum(byProd(idxs, :), 2))
                ];
            
            %Lipid
            [bmComp, bmProd, byProd] = edu.stanford.covert.cell.sim.analysis.BiomassCompositionProduction.getResourceRequirements(sim);
            bmComp = bmComp * mass.dryWeightFractionLipid;
            bmProd = bmProd * mass.dryWeightFractionLipid;
            byProd = byProd * mass.dryWeightFractionLipid;
            idxs = find(any(bmComp, 2) | any(bmProd, 2) | any(byProd, 2));
            content = [content;
                {0 'Lipid' '' sum(sum(bmComp(idxs, :))) sum(sum(bmProd(idxs, :))) sum(sum(byProd(idxs, :)))}
                num2cell(ones(size(idxs))) m.wholeCellModelIDs(idxs) m.names(idxs) num2cell(sum(bmComp(idxs, :), 2)) num2cell(sum(bmProd(idxs, :), 2)) num2cell(sum(byProd(idxs, :), 2))
                ];
            
            %Polyamine
            [bmComp, bmProd, byProd] = edu.stanford.covert.cell.sim.analysis.BiomassCompositionProduction.getResourceRequirements(sim);
            bmComp = bmComp * mass.dryWeightFractionPolyamine;
            bmProd = bmProd * mass.dryWeightFractionPolyamine;
            byProd = byProd * mass.dryWeightFractionPolyamine;
            idxs = find(any(bmComp, 2) | any(bmProd, 2) | any(byProd, 2));
            content = [content;
                {0 'Polyamine' '' sum(sum(bmComp(idxs, :))) sum(sum(bmProd(idxs, :))) sum(sum(byProd(idxs, :)))}
                num2cell(ones(size(idxs))) m.wholeCellModelIDs(idxs) m.names(idxs) num2cell(sum(bmComp(idxs, :), 2)) num2cell(sum(bmProd(idxs, :), 2)) num2cell(sum(byProd(idxs, :), 2))
                ];
            
            %Carbohydrate
            [bmComp, bmProd, byProd] = edu.stanford.covert.cell.sim.analysis.BiomassCompositionProduction.getResourceRequirements(sim);
            bmComp = bmComp * mass.dryWeightFractionCarbohydrate;
            bmProd = bmProd * mass.dryWeightFractionCarbohydrate;
            byProd = byProd * mass.dryWeightFractionCarbohydrate;
            idxs = find(any(bmComp, 2) | any(bmProd, 2) | any(byProd, 2));
            content = [content;
                {0 'Carbohydrate' '' sum(sum(bmComp(idxs, :))) sum(sum(bmProd(idxs, :))) sum(sum(byProd(idxs, :)))}
                num2cell(ones(size(idxs))) m.wholeCellModelIDs(idxs) m.names(idxs) num2cell(sum(bmComp(idxs, :), 2)) num2cell(sum(bmProd(idxs, :), 2)) num2cell(sum(byProd(idxs, :), 2))
                ];
            
            %Vitamin
            [bmComp, bmProd, byProd] = edu.stanford.covert.cell.sim.analysis.BiomassCompositionProduction.getResourceRequirements(sim);
            bmComp = bmComp * mass.dryWeightFractionVitamin;
            bmProd = bmProd * mass.dryWeightFractionVitamin;
            byProd = byProd * mass.dryWeightFractionVitamin;
            idxs = find(any(bmComp, 2) | any(bmProd, 2) | any(byProd, 2));
            content = [content;
                {0 'Vitamin' '' sum(sum(bmComp(idxs, :))) sum(sum(bmProd(idxs, :))) sum(sum(byProd(idxs, :)))}
                num2cell(ones(size(idxs))) m.wholeCellModelIDs(idxs) m.names(idxs) num2cell(sum(bmComp(idxs, :), 2)) num2cell(sum(bmProd(idxs, :), 2)) num2cell(sum(byProd(idxs, :), 2))
                ];
            
            %Ion
            [bmComp, bmProd, byProd] = edu.stanford.covert.cell.sim.analysis.BiomassCompositionProduction.getResourceRequirements(sim);
            bmComp = bmComp * mass.dryWeightFractionIon;
            bmProd = bmProd * mass.dryWeightFractionIon;
            byProd = byProd * mass.dryWeightFractionIon;
            idxs = find(any(bmComp, 2) | any(bmProd, 2) | any(byProd, 2));
            content = [content;
                {0 'Ion' '' sum(sum(bmComp(idxs, :))) sum(sum(bmProd(idxs, :))) sum(sum(byProd(idxs, :)))}
                num2cell(ones(size(idxs))) m.wholeCellModelIDs(idxs) m.names(idxs) num2cell(sum(bmComp(idxs, :), 2)) num2cell(sum(bmProd(idxs, :), 2)) num2cell(sum(byProd(idxs, :), 2))
                ];
            
            %format output
            indentation = cell2mat(content(:, 1));
            content = content(:, 2:end);
        end
    end
    
    methods (Static)
        function plotBiomassComposition(sim, axesHandle, bmComp)
            axesHandle1 = subplot(3, 3, 1:6);
            axesHandle2 = subplot(3, 3, 7);
            axesHandle3 = subplot(3, 3, 8);
            axesHandle4 = subplot(3, 3, 9);
            
            pos = get(axesHandle1, 'Position');
            set(axesHandle1, 'Position', [(pos(1:2) + pos(3:4)) - 1.1 * pos(3:4)  1.1 * pos(3:4)]);
            
            pos = get(axesHandle2, 'Position');
            set(axesHandle2, 'Position', [(pos(1:2) + pos(3:4)) - 1.5 * pos(3:4)  1.5 * pos(3:4)]);
            
            pos = get(axesHandle3, 'Position');
            set(axesHandle3, 'Position', [(pos(1:2) + pos(3:4)) - 1.5 * pos(3:4)  1.5 * pos(3:4)]);
            
            pos = get(axesHandle4, 'Position');
            set(axesHandle4, 'Position', [(pos(1:2) + pos(3:4)) - 1.5 * pos(3:4)  1.5 * pos(3:4)]);            
            
            %% weight fractions
            m = sim.state('Mass');
            values = [
                m.dryWeightFractionCarbohydrate
                m.dryWeightFractionDNA
                m.dryWeightFractionIon
                m.dryWeightFractionLipid
                m.dryWeightFractionPolyamine
                m.dryWeightFractionProtein
                m.dryWeightFractionRNA
                m.dryWeightFractionVitamin
                m.dryWeightFractionNucleotide];
            values = values / sum(values);
            labels = {
                'Carbohydrate'
                'DNA'
                'Ion'
                'Lipid'
                'Polyamine'
                'Protein'
                'RNA'
                'Vitamin'
                'Nucleotide'
                };
            
            %pie chart
            h = pie(axesHandle1, values, ...
                cellfun(@(value, label) {label sprintf('(%.1f%%)', value * 100)}, num2cell(values), labels, 'UniformOutput', false));
            set(h(2:2:end), 'FontSize', 12);
            
            %turn off labels of smallest weight fractions
            idx = find(values < 0.02) * 2;
            for i = 1:numel(idx)
                set(h(idx(i)), 'Visible', 'off')
            end
            
            %% dNMP
            m = sim.state('Metabolite');
            idxs = [m.dnmpIndexs; m.getIndexs({'m6dAMP'})];
            values = sum(bmComp(idxs, :), 2);
            values = values / sum(values);
            labels = m.wholeCellModelIDs(idxs);
            labels = {
                'dAMP'
                'dCMP'
                'dGMP'
                'dTMP'
                'm^6dAMP'
                };
            
            %pie chart
            h = pie(axesHandle2, values, labels);
            set(h(2:2:end), 'FontSize', 12);
            
            %% NMP
            m = sim.state('Metabolite');
            idxs = m.nmpIndexs;
            values = sum(bmComp(idxs, :), 2);
            values = values / sum(values);
            labels = m.names(idxs);
            
            %pie chart
            h = pie(axesHandle3, values, labels);
            set(h(2:2:end), 'FontSize', 12);
            
            %% AA
            m = sim.state('Metabolite');
            idxs = m.aminoAcidIndexs(1:20);
            values = sum(bmComp(idxs, :), 2);
            values = values / sum(values);
            labels = num2cell(edu.stanford.covert.cell.kb.ProteinMonomer.bases);
            
            %pie chart
            h = pie(axesHandle4, values, labels);
            set(h(2:2:end), 'FontSize', 12);
        end
        
        function plotWeightFractions(sim, axesHandle)
            m = sim.state('Mass');
            
            values = [
                m.dryWeightFractionCarbohydrate
                m.dryWeightFractionDNA
                m.dryWeightFractionIon
                m.dryWeightFractionLipid
                m.dryWeightFractionPolyamine
                m.dryWeightFractionProtein
                m.dryWeightFractionRNA
                m.dryWeightFractionVitamin
                m.dryWeightFractionNucleotide
                ];
            labels = {
                'Carbohydrate'
                'DNA'
                'Ion'
                'Lipid'
                'Polyamine'
                'Protein'
                'RNA'
                'Vitamin'
                'Nucleotide'
                };
            
            %pie chart
            h = pie(axesHandle, values, ...
                cellfun(@(value, label) {label sprintf('(%.1f%%)', value * 100)}, num2cell(values), labels, 'UniformOutput', false));
            
            %turn off labels of smallest weight fractions
            idx = find(values < 0.009) * 2;
            for i = 1:numel(idx)
                set(h(idx(i)), 'Visible', 'off')
            end
        end
        
        function plotDNMPComposition(sim, axesHandle, bmComp)
            m = sim.state('Metabolite');
            idxs = [m.dnmpIndexs; m.getIndexs({'m6dAMP'})];
            values = sum(bmComp(idxs, :), 2);
            values = values / sum(values);
            labels = m.names(idxs);
            
            %pie chart
            pie(axesHandle, values, ...
                cellfun(@(value, label) sprintf('%s\n(%.1f%%)', label, value * 100), num2cell(values), labels, 'UniformOutput', false));
        end
        
        function plotNMPComposition(sim, axesHandle, bmComp)
            m = sim.state('Metabolite');
            idxs = m.nmpIndexs;
            values = sum(bmComp(idxs, :), 2);
            values = values / sum(values);
            labels = m.names(idxs);
            
            %pie chart
            pie(axesHandle, values, ...
                cellfun(@(value, label) sprintf('%s\n(%.1f%%)', label, value * 100), num2cell(values), labels, 'UniformOutput', false));
        end
        
        function plotAAComposition(sim, axesHandle, bmComp)
            m = sim.state('Metabolite');
            idxs = m.aminoAcidIndexs(1:20);
            values = sum(bmComp(idxs, :), 2);
            values = values / sum(values);
            labels = m.names(idxs);
            
            %pie chart
            pie(axesHandle, values, ...
                cellfun(@(value, label) sprintf('%s\n(%.1f%%)', label, value * 100), num2cell(values), labels, 'UniformOutput', false));
        end
    end
    
    methods (Static)
        function [bmComp, bmProd, byProd] = getResourceRequirements(sim)
            %references
            m = sim.state('Metabolite');
            
            %calculate biomass composition, maintenance energy, and byproducts
            bmComp = m.biomassComposition;
            bmProd = m.biomassProduction;
            byProd = m.byproducts;
        end
    end
end