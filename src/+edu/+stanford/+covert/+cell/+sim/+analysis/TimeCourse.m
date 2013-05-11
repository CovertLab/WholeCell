%TimeCourse
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/12/2010
classdef TimeCourse
    methods (Static)
        function run(sim, fileName)
            import edu.stanford.covert.cell.sim.analysis.TimeCourse;
            
            % plots
            if nargin == 1
                TimeCourse.plotRNAs(sim, PlotUtil.newAxesHandle(), 1, cIdx);
                TimeCourse.plotrRNAs(sim, PlotUtil.newAxesHandle(), 1, cIdx);
                TimeCourse.plotsRNAs(sim, PlotUtil.newAxesHandle(), 1, cIdx);
                TimeCourse.plottRNAs(sim, PlotUtil.newAxesHandle(), 1, cIdx);
                TimeCourse.plotProteins(sim, PlotUtil.newAxesHandle(), 1, compIdxs);
                TimeCourse.plotWeights(sim, PlotUtil.newAxesHandle(), 1, compIdxs);
                TimeCourse.plotRNA_Weight_Fractions(sim, PlotUtil.newAxesHandle(), 1, compIdxs);
                TimeCourse.plotNTPs(sim, PlotUtil.newAxesHandle(), 1, cIdx);
                TimeCourse.plotAAs(sim, PlotUtil.newAxesHandle(), 1, cIdx);
            else
                [axesHandle, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
                
                cIdx = sim.compartment.cytosolIndexs;
                compIdxs = 1:sim.compartment.count;
                
                cla(axesHandle);
                TimeCourse.plotRNAs(sim, axesHandle, 1, cIdx);
                saveas(figHandle, [fileName '-RNAs.pdf']);
                
                cla(axesHandle);
                TimeCourse.plotrRNAs(sim, axesHandle, 1, cIdx);
                saveas(figHandle, [fileName '-rRNAs.pdf']);
                
                cla(axesHandle);
                TimeCourse.plotsRNAs(sim, axesHandle, 1, cIdx);
                saveas(figHandle, [fileName '-sRNAs.pdf']);
                
                cla(axesHandle);
                TimeCourse.plottRNAs(sim, axesHandle, 1, cIdx);
                saveas(figHandle, [fileName '-tRNAs.pdf']);
                
                cla(axesHandle);
                TimeCourse.plotProteins(sim, axesHandle, 1, compIdxs);
                saveas(figHandle, [fileName '-Proteins.pdf']);
                
                cla(axesHandle);
                TimeCourse.plotWeights(sim, axesHandle, 1, compIdxs);
                saveas(figHandle, [fileName '-Weights.pdf']);
                
                cla(axesHandle);
                TimeCourse.plotRNA_Weight_Fractions(sim, axesHandle, 1, compIdxs);
                saveas(figHandle, [fileName '-Weight Fractions.pdf']);
                
                cla(axesHandle);
                TimeCourse.plotNTPs(sim, axesHandle, 1, cIdx);
                saveas(figHandle, [fileName '-NTPs.pdf']);
                
                cla(axesHandle);
                TimeCourse.plotAAs(sim, axesHandle, 1, cIdx);
                saveas(figHandle, [fileName '-AAs.pdf']);
                
                close(figHandle);
            end
        end
    end
    
    %plots
    methods (Static)
        % plot RNA
        function plotRNAs(sim, axesHandle, time, compartments)
            r = sim.state('Rna');
            
            values = edu.stanford.covert.cell.sim.util.PlotUtil.selectCompartmentsForPlot(sim, compartments, [
                sum(r.counts, 1);
                sum(r.counts(r.matureIndexs(r.matureMRNAIndexs), :, :), 1);
                sum(r.counts(r.matureIndexs(r.matureRRNAIndexs), :, :), 1);
                sum(r.counts(r.matureIndexs(r.matureSRNAIndexs), :, :), 1);
                sum(r.counts(r.matureIndexs(r.matureTRNAIndexs), :, :), 1)]);
            
            plotHandles = plot(axesHandle, time/3600, values);
            sim.labelCompartmentPlot(plotHandles, compartments, {'Total','mRNA','rRNA','sRNA','tRNA'});
            xlabel(axesHandle, 'Time (h)', 'fontSize', 8);
            ylabel(axesHandle, 'RNA', 'fontSize', 8);
            xlim(axesHandle, [0 max(time)] / 3600);
        end
        
        function plotrRNAs(sim, axesHandle, time, compartments)
            r = sim.state('Rna');
            
            values = edu.stanford.covert.cell.sim.util.PlotUtil.selectCompartmentsForPlot(...
                sim, compartments, r.counts(r.matureIndexs(r.matureRRNAIndexs), :, :));
            plotHandles = plot(axesHandle, time/3600, values);
            sim.labelCompartmentPlot(plotHandles, compartments, r.wholeCellModelIDs(r.matureIndexs(r.matureRRNAIndexs)));
            xlabel(axesHandle, 'Time (h)', 'fontSize', 8);
            ylabel(axesHandle, 'rRNA', 'fontSize', 8);
            xlim(axesHandle, [0 max(time)] / 3600);
        end
        
        function plotsRNAs(sim, axesHandle, time, compartments)
            r = sim.state('Rna');
            
            values = edu.stanford.covert.cell.sim.util.PlotUtil.selectCompartmentsForPlot(...
                sim, compartments, r.counts(r.matureIndexs(r.matureSRNAIndexs), :, :));
            plotHandles = plot(axesHandle, time/3600, values);
            sim.labelCompartmentPlot(plotHandles, compartments, r.wholeCellModelIDs(r.matureIndexs(r.matureSRNAIndexs)));
            xlabel(axesHandle, 'Time (h)', 'fontSize', 8);
            ylabel(axesHandle, 'sRNA', 'fontSize', 8);
            xlim(axesHandle, [0 max(time)] / 3600);
        end
        
        function plottRNAs(sim, axesHandle, time, compartments)
            r = sim.state('Rna');
            
            values = edu.stanford.covert.cell.sim.util.PlotUtil.selectCompartmentsForPlot(...
                sim, compartments, r.counts(r.matureIndexs(r.matureTRNAIndexs), :, :));
            plotHandles = plot(axesHandle, time/3600, values);
            sim.labelCompartmentPlot(plotHandles, compartments);
            xlabel(axesHandle, 'Time (h)', 'fontSize', 8);
            ylabel(axesHandle, 'tRNA', 'fontSize', 8);
            xlim(axesHandle, [0 max(time)] / 3600);
        end
        
        % plot proteins
        function plotProteins(sim, axesHandle, time, compartments)
            m = sim.state('ProteinMonomer');
            c = sim.state('ProteinComplex');
            
            values = edu.stanford.covert.cell.sim.util.PlotUtil.selectCompartmentsForPlot(...
                sim, compartments, [sum(m.counts); sum(c.counts)]);
            plotHandles = plot(axesHandle, time/3600, values);
            sim.labelCompartmentPlot(plotHandles, compartments, {'Monomers','Complexes'});
            xlabel(axesHandle, 'Time (h)', 'fontSize', 8);
            ylabel(axesHandle, 'Protein', 'fontSize', 8);
            xlim(axesHandle, [0 max(time)] / 3600);
        end
        
        % plot weights
        function plotWeights(sim, axesHandle, time, compartments)
            mass = sim.state('Mass');
            values = edu.stanford.covert.cell.sim.util.PlotUtil.selectCompartmentsForPlot(sim, compartments, [
                mass.cell;
                mass.rnaWt;
                mass.proteinWt]) / (mass.cellInitialDryWeight * (1 + 1/mass.fractionWetWeight));
            plotHandles = plot(axesHandle, time/3600, values);
            sim.labelCompartmentPlot(plotHandles, compartments, {'Total', 'RNA', 'Protein'});
            xlabel(axesHandle, 'Time (h)', 'fontSize', 8);
            ylabel(axesHandle, 'Weight (Cells)', 'fontSize', 8);
            xlim(axesHandle, [0 max(time)] / 3600);
        end
        
        function plotRNA_Weight_Fractions(sim, axesHandle, time, compartments)
            r = sim.state('Rna');
            
            weightFractions = [...
                multiprod(r.molecularWeights(r.matureIndexs(r.matureMRNAIndexs))', ...
                r.counts(r.matureIndexs(r.matureMRNAIndexs), :, :), [1 2], [1 2]);...
                repmat(r.molecularWeights(r.matureIndexs(r.matureRibosomalRRNAIndexs)),...
                [1 numel(sim.compartment.wholeCellModelIDs) size(r.counts,3)]).*...
                r.counts(r.matureIndexs(r.matureRibosomalRRNAIndexs),:,:);...
                multiprod(r.molecularWeights(r.matureIndexs([r.matureSRNAIndexs;r.matureTRNAIndexs]))', ...
                r.counts(r.matureIndexs([r.matureSRNAIndexs;r.matureTRNAIndexs]), :, :),[1 2],[1 2])];
            weightFractions = weightFractions./repmat(sum(weightFractions(:)),[size(weightFractions,1),size(weightFractions,2), 1]);
            values = edu.stanford.covert.cell.sim.util.PlotUtil.selectCompartmentsForPlot(sim, compartments, weightFractions);
            
            plotHandles = plot(axesHandle, time/3600, values);
            sim.labelCompartmentPlot(plotHandles, compartments, {'mRNA','rRNA 5S','rRNA 16S','rRNA 23S','sRNA, tRNA'});
            xlabel(axesHandle, 'Time (h)', 'fontSize', 8);
            ylabel(axesHandle, 'Weight Fraction', 'fontSize', 8);
            xlim(axesHandle, [0 max(time)] / 3600);
        end
        
        function plotNTPs(sim, axesHandle, time, compartments)
            r = sim.state('Rna');
            
            ntpIncorporation = multiprod(r.baseCounts(:,sim.state('Metabolite').nmpIndexs)',r.counts,[1 2],[1 2]);
            ntpIncorporation = ntpIncorporation./repmat(sum(sum(ntpIncorporation,1),2),[size(ntpIncorporation,1), size(ntpIncorporation,2), 1]);
            values = edu.stanford.covert.cell.sim.util.PlotUtil.selectCompartmentsForPlot(sim, compartments, ntpIncorporation);
            plotHandles = plot(axesHandle, time/3600, values);
            sim.labelCompartmentPlot(plotHandles, compartments, sim.state('Metabolite').wholeCellModelIDs(sim.state('Metabolite').ntpIndexs));
            xlabel(axesHandle, 'Time (h)', 'fontSize',8);
            ylabel(axesHandle, 'NMPs', 'fontSize',8);
            xlim(axesHandle, [0 max(time)] / 3600);
        end
        
        function plotAAs(sim, axesHandle, time, compartments)
            m = sim.state('ProteinMonomer');
            c = sim.state('ProteinComplex');
            
            aaIncorporation=...
                multiprod(m.baseCounts(:,sim.state('Metabolite').aminoAcidIndexs)',m.counts,[1 2],[1 2]) + ...
                multiprod(c.naseCounts(:,sim.state('Metabolite').aminoAcidIndexs)',c.counts,[1 2],[1 2]);
            
            aaIncorporation=aaIncorporation./repmat(sum(sum(aaIncorporation,1),2), [size(aaIncorporation,1), size(aaIncorporation,2), 1]);
            values=edu.stanford.covert.cell.sim.util.PlotUtil.selectCompartmentsForPlot(sim, compartments, aaIncorporation);
            plotHandles=plot(axesHandle, time/3600, values);
            
            sim.labelCompartmentPlot(plotHandles, compartments, sim.state('Metabolite').wholeCellModelIDs(sim.state('Metabolite').aminoAcidIndexs));
            xlabel(axesHandle,'Time (h)','fontSize',8);
            ylabel(axesHandle,'AAs','fontSize',8);
            xlim(axesHandle,[0 max(time)]/3600);
        end
    end
end