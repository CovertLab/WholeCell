%Translation
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/6/2011
classdef Translation
    %plotting
    methods (Static = true)
        function plotRibosome_State_Occupancies(sim, axesHandle, time, ~)
            s = sim.state('Ribosome');
            
            plotHandles = plot(axesHandle, time / 3600, permute(s.stateOccupancies, [3 1 2]));
            xlabel(axesHandle, 'Time (h)', 'fontSize',8);
            ylabel(axesHandle, 'Occupancy', 'fontSize',8);
            xlim(axesHandle, [0 max(1, time(end))] / 3600);
            
            labels = cell(3, 1);
            labels{s.activeIndex} = 'Active';
            labels{s.notExistIndex} = 'Not Exist';
            labels{s.stalledIndex} = 'Stalled';
            legend(plotHandles, labels);
        end
        
        function plotRibosomes(sim, axesHandle, time, ~)
            s = sim.state('Ribosome');
            
            ribosomeGenomeStates = zeros(size(s.states));
            for i = 1:size(s.states, 1)
                for k = 1:size(s.states, 3)
                    if s.states(i, 1, k) == s.activeValue
                        if s.gene.strands(s.boundMRNAs(i, 1, k)) == 1
                            ribosomeGenomeStates(i, 1, k) = ...
                                s.gene.startCoordinates(s.boundMRNAs(i, 1, k)) + ...
                                s.polypeptide.nascentMonomerLengths(i, 1, k);
                        else
                            ribosomeGenomeStates(i, 1, k) = ...
                                s.gene.startCoordinates(s.boundMRNAs(i, 1, k)) + ...
                                s.gene.lengths(s.boundMRNAs(i, 1, k)) - 1 - ...
                                s.polypeptide.nascentMonomerLengths(i, 1, k);
                        end
                    end
                end
            end
            
            plot(axesHandle, time/3600, permute(ribosomeGenomeStates,[3 1 2]), '.');
            xlabel(axesHandle,'Time (h)','fontSize',8);
            ylabel(axesHandle,'Genome Coordinate','fontSize',8);
            xlim(axesHandle, [0 max(1, time(end))]/3600);
            ylim(axesHandle, [s.stalledValue-1 s.chromosome.sequenceLen + 1]);
        end
        
        function plotRibosome_States(sim, axesHandle, time, ~)
            import edu.stanford.covert.util.ConstantUtil;
            
            s = sim.state('Ribosome');
            
            colors = zeros([size(s.states) 3]);
            for i = 1:size(s.states, 1)
                for k = 1:size(s.states, 3)
                    switch s.states(i, 1, k)
                        case s.activeValue
                            colors(i, 1, k, 2) = s.polypeptide.nascentMonomerLengths(i, 1, k) / s.polypeptide.monomerLengths(s.boundMRNAs(i, 1, k));
                        case s.stalledValue
                            colors(i, 1, k, 1) = s.polypeptide.proteolysisTagLengths(i, 1, k) / s.proteolysisTagLength;
                        case s.notExistValue
                            colors(i, 1, k, 3) = 1;
                    end
                end
            end
            
            colors = reshape(colors, [], size(s.states, 3), 3);
            image(colors, 'Parent', axesHandle);
            
            xlabel(axesHandle, 'Time (h)', 'fontSize', 8);
            ylabel(axesHandle, 'Ribosome', 'fontSize', 8);
            
            xTick = get(axesHandle, 'XTick');
            xTickLabel = round((xTick - 0.5) * max(time) / ConstantUtil.secondsPerHour / max(xTick));
            set(axesHandle, 'XTickLabel', xTickLabel);
        end
        
        function plotRibosome_Bound_Genes(sim, axesHandle, time, ~)
            import edu.stanford.covert.util.ConstantUtil;
            
            s = sim.state('Ribosome');
            
            boundGenes = zeros(length(s.gene.wholeCellModelIDs), 1, length(time));
            for i = 1:size(s.boundMRNAs, 1)
                for k = 1:size(s.boundMRNAs, 3)
                    if s.boundMRNAs(i, 1, k) <= 0
                        continue;
                    end
                    boundGenes(s.boundMRNAs(i,1,k), 1, k) = ...
                        boundGenes(s.boundMRNAs(i, 1, k), 1, k) + 1;
                end
            end
            
            boundGenes = permute(boundGenes, [1 3 2]);
            colors = repmat(boundGenes / (max(max(boundGenes))), [1 1 3]);
            
            image(colors, 'Parent', axesHandle);
            xlabel(axesHandle, 'Time (h)', 'fontSize',8);
            ylabel(axesHandle, 'Gene', 'fontSize',8);
            
            xTick = get(axesHandle, 'XTick');
            xTickLabel = round((xTick - 0.5) * max(time) / ConstantUtil.secondsPerHour / max(xTick));
            set(axesHandle, 'XTickLabel', xTickLabel);
        end
        
        function plotBound_Translation_Factors(sim, axesHandle, time, ~)
            p = sim.process('Translation');
            
            translationFactors      = p.enzymes(p.enzymeIndexs_translationFactors, :, :);
            boundTranslationFactors = p.boundEnzymes(p.enzymeIndexs_translationFactors, :, :);
            
            plotHandles = plot(axesHandle, time / 3600, permute(boundTranslationFactors ./ (translationFactors + boundTranslationFactors), [3 1 2]));
            legend(plotHandles, p.enzymeNames(p.enzymeIndexs_translationFactors));
            xlabel(axesHandle, 'Time (h)', 'fontSize', 8);
            ylabel(axesHandle, 'Fraction Bound', 'fontSize', 8);
            xlim(axesHandle, [0 max(1, time(end))] / 3600);
            ylim(axesHandle, [0 1]);
        end
        
        function plotAAs(sim, axesHandle, time, ~)
            p = sim.process('tRNAAminoacylation');
            
            plotHandles = plot(axesHandle,time/3600, permute(p.substrates(p.substrateIndexs_aminoAcids,:,:),[3 1 2]));
            legend(plotHandles, p.substrateWholeCellModelIDs(p.substrateIndexs_aminoAcids));
            xlabel(axesHandle,'Time (h)','fontSize',8);
            ylabel(axesHandle,'AAs','fontSize',8);
            xlim(axesHandle,[0 max(1, time(end))]/3600);
        end
        
        function plottRNA_Synthetases(sim, axesHandle, time, ~)
            p = sim.process('tRNAAminoacylation');
            
            plotHandles = plot(axesHandle, time/3600, permute(p.tRNASynthetases,[3 1 2]));
            legend(plotHandles, p.enzymeNames(p.enzymeIndexs_tRNASynthetases));
            xlabel(axesHandle,'Time (h)','fontSize',8);
            ylabel(axesHandle,'Complexs','fontSize',8);
            xlim(axesHandle,[0 max(1, time(end))]/3600);
        end
        
        function plottRNA_Transferases(sim, axesHandle, time, ~)
            p = sim.process('tRNAAminoacylation');
            
            plotHandles = plot(axesHandle,time/3600, permute(p.tRNATransferases,[3 1 2]));
            legend(plotHandles, p.enzymeNames(p.enzymeIndexs_tRNATransferases));
            xlabel(axesHandle,'Time (h)','fontSize',8);
            ylabel(axesHandle,'Complexs','fontSize',8);
            xlim(axesHandle,[0 max(1, time(end))]/3600);
        end
    end
end