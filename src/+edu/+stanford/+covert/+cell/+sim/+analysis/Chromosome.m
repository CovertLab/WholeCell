%Chromosome
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/6/2011
classdef Chromosome
    %Chromosome
    methods (Static = true)
        %Plots locations of DNA damages on circular chromosome.
        %(+) strand damages extend outward
        %(-) strand damages extend inward
        %
        %Feature                     Color
        %==========                  =========
        %chromosome                  grey
        %gap sites                   red
        %abasic sites                green
        %damaged sugar-phosphates    blue
        %damaged bases               cyan
        %intrastrand cross links     magenta
        %strand breaks               orange
        %holliday junctions          yellow
        function plotCircular_Chromosome(sim, axesHandle, ~, ~)
            c = sim.state('Chromosome');
            
            %parameters
            damageLength = 0.1;
            strandSeparation = 0.02;
            
            %genome
            genomeLength = size(c.gapSites,1);
            [x,y]=pol2cart((0:0.01:1)*2*pi, repmat(1+strandSeparation,1,101));
            h = line(x, y, 'Parent', axesHandle); %(+) strand
            set(h,'Color',[0.6 0.6 0.6],'LineWidth',2);
            
            [x,y]=pol2cart((0:0.01:1)*2*pi, repmat(1-strandSeparation,1,101));
            h = line(x, y, 'Parent', axesHandle); %(-) strand
            set(h,'Color',[0.6 0.6 0.6],'LineWidth',2);
            
            %damage
            damageTypes={'gapSites','abasicSites','damagedSugarPhosphates','damagedBases','intrastrandCrossLinks','strandBreaks','hollidayJunctions'};
            damageNames={'Gap Sites','Abasic Sites','Damaged Sugar Phosphates','Damaged Bases','Intrastrand Cross Links','Strand Breaks','Holliday Junctions'};
            colors = [
                255 0 0;
                0 255 0;
                0 0 255;
                0, 183, 235;
                255, 0, 255;
                255, 127, 0;
                255, 205, 0]/255;
            lineHandles = zeros(0, 1);
            lineHandleIndexs = zeros(size(damageTypes));
            damageLabels = cell(size(damageTypes));
            for i=1:length(damageTypes)
                damagesTensor = sum(c.(damageTypes{i}), 3);
                if size(damagesTensor,3)>1
                    damagesTensor = sum(damagesTensor,3);
                end
                damages = find(damagesTensor);
                if size(damagesTensor,2)==1
                    damages = [damages; damages + repmat([0 1], size(damages,1), 1)]; %#ok<AGROW>
                end
                
                [x,y]=pol2cart((repmat(damages(:,1),1,2)*2*pi/genomeLength)', ...
                    (1 + repmat(2*(damages(:,2)-1)-1,1,2).*repmat([strandSeparation damageLength],size(damages,1),1))');
                
                tmpHandles = line(x,y,'Parent',axesHandle, 'Color', colors(i,:), 'LineWidth',0.1);
                if ~isempty(tmpHandles)
                    lineHandleIndexs(i) = size(lineHandles,1)+1;
                end
                lineHandles = [lineHandles; tmpHandles]; %#ok<AGROW>
                
                damageLabels{i}=sprintf('%s (%d)',damageNames{i}, nnz(c.(damageTypes{i})));
            end
            
            xlim((1+2*damageLength)*[-1 1])
            ylim((1+2*damageLength)*[-1 1])
            set(axesHandle,'Box','off','Color','none','Visible','off','DataAspectRatio',[1 1 1],'XTick',[],'YTick',[],'XMinorTick','off','YMinorTick','off','GridLineStyle','none')
            legend(lineHandles(lineHandleIndexs(lineHandleIndexs ~= 0)), damageLabels{lineHandleIndexs ~= 0}, 'Location', 'NorthEastOutside');
            legend(axesHandle,'show');
        end
        
        %Plots locations of DNA damages on chromosome displayed as array of
        %lines
        %(+) strand damages extend upward
        %(-) strand damages extend downward
        %
        %Feature                     Color
        %==========                  =========
        %chromosome                  grey
        %gap sites                   red
        %abasic sites                green
        %damaged sugar-phosphates    blue
        %damaged bases               cyan
        %intrastrand cross links     magenta
        %strand breaks               orange
        %holliday junctions          yellow
        function plotPlanar_Chromosome(sim, axesHandle, ~, ~)
            c = sim.state('Chromosome');
            
            %parameters
            damageLength = 0.3;
            strandSeparation = 0.05;
            rows = 10;
            
            %genome
            genomeLength = size(c.gapSites,1);
            line(repmat([0 1],rows,1)', repmat((1:rows)'-0.5+strandSeparation,1,2)','Parent',axesHandle,'Color',[0.6 0.6 0.6],'LineWidth',2); %(+) strand
            line(repmat([0 1],rows,1)', repmat((1:rows)'-0.5-strandSeparation,1,2)','Parent',axesHandle,'Color',[0.6 0.6 0.6],'LineWidth',2); %(-) strand
            
            %damage
            damageTypes={'gapSites','abasicSites','damagedSugarPhosphates','damagedBases','intrastrandCrossLinks','strandBreaks','hollidayJunctions'};
            damageNames={'Gap Sites','Abasic Sites','Damaged Sugar Phosphates','Damaged Bases','Intrastrand Cross Links','Strand Breaks','Holliday Junctions'};
            colors = [
                255 0 0;
                0 255 0;
                0 0 255;
                0, 183, 235;
                255, 0, 255;
                255, 127, 0;
                255, 205, 0]/255;
            lineHandles = zeros(0, 1);
            lineHandleIndexs = zeros(size(damageTypes));
            damageLabels = cell(size(damageTypes));
            for i=1:length(damageTypes)
                damagesTensor = sum(c.(damageTypes{i}), 3);
                if size(damagesTensor,3)>1
                    damagesTensor = sum(damagesTensor,3);
                end
                damages = find(damagesTensor);
                if size(damagesTensor,2)==1
                    damages = [damages; damages + repmat([0 1], size(damages,1), 1)]; %#ok<AGROW>
                end
                
                tmpHandles = line(...
                    repmat(mod(damages(:,1)/genomeLength*rows,1),1,2)', ...
                    (repmat(rows-0.5-floor(damages(:,1)/genomeLength*rows) + ...
                    strandSeparation * 2*(damages(:,2)-1)-1,1,2) + ...
                    [zeros(size(damages,1),1) damageLength*2*(damages(:,2)-1)-1])',...
                    'Parent',axesHandle,'Color', colors(i,:), 'LineWidth',1);
                if ~isempty(tmpHandles)
                    lineHandleIndexs(i) = size(lineHandles,1)+1;
                end
                lineHandles = [lineHandles; tmpHandles]; %#ok<AGROW>
                
                damageLabels{i}=sprintf('%s (%d)',damageNames{i}, nnz(c.(damageTypes{i})));
            end
            xlim(axesHandle,[0 1]);
            ylim(axesHandle,[0 rows]);
            set(axesHandle,'Box','on','XTick',[],'YTick',[],'XMinorTick','off','YMinorTick','off')
            legend(lineHandles(lineHandleIndexs(lineHandleIndexs ~= 0)), damageLabels{lineHandleIndexs ~= 0}, 'Location', 'NorthEastOutside');
            legend(axesHandle,'show');
        end
        
        %Plots state of DisA DNA integrity scanning system.
        function plotDisA_DNA_Integrity_Scanning_Status(sim, axesHandle, time, ~)
            p = sim.process('DNARepair');
            
            free  = sum(sum(p.enzymes(p.enzymeIndexs_DisA, :, :), 3), 2);
            bound = sum(sum(p.boundEnzymes(p.enzymeIndexs_DisA,  :, :), 3), 2);
            
            fracBound = bound ./ (bound + free);
            
            plot(axesHandle, time/3600, permute(fracBound, [3 1 2]));
            xlabel(axesHandle, 'Time (h)', 'fontSize',8);
            ylabel(axesHandle, 'Fraction Bound', 'fontSize',8);
            xlim(axesHandle, [0 max(1, time(end))] / 3600);
            minVal = min(fracBound);
            maxVal = max(fracBound);
            if minVal == maxVal
                maxVal = min(1, minVal + 0.05);
                minVal = max(0, minVal - 0.05);
            end
            if isfinite(minVal) && isfinite(maxVal)
                ylim([minVal - (maxVal - minVal) * 0.1   maxVal + (maxVal - minVal) * 0.1])
            end
        end
        
        %Plots state of restriction/modification system on circular chromosome.
        %(+) strand damages extend outward
        %(-) strand damages extend inward
        %
        %Feature                     Color
        %==========                  =========
        %chromosome                  grey
        %MunI region                 black
        %damaged region              red
        %unmethylated                magenta
        %hemimethylated              blue
        %methylayed                  green
        %cleaved                     cyan
        function plotRestriction_Modification_Status(sim, axesHandle, ~, ~)
            p = sim.process('DNARepair');
            
            %parameters
            damageLength = 0.1;
            strandSeparation = 0.02;
            
            %genome
            genomeLength = size(p.chromosome.gapSites,1);
            [x, y] = pol2cart((0:0.01:1) * 2 * pi, repmat(1+strandSeparation, 1, 101));
            h = line(x, y, 'Parent', axesHandle); %(+) strand
            set(h, 'Color', [0.6 0.6 0.6], 'LineWidth', 2);
            
            [x, y] = pol2cart((0:0.01:1) * 2 * pi, repmat(1-strandSeparation, 1, 101));
            h = line(x, y, 'Parent', axesHandle); %(-) strand
            set(h, 'Color', [0.6 0.6 0.6], 'LineWidth', 2);
            
            %get status of MunI sites
            [unmethylatedSites, hemimethylatedSites, methylatedSites, cleavedSites, inaccessibleRegions] = ...
                p.chromosome.rmStatus(p.RM_MunI_RecognitionSites, p.RM_MunI_MethylatedPositions, p.RM_MunI_RestrictionPositions, ...
                [], p.enzymeGlobalIndexs(p.enzymeIndexs_RM_typeII));
            sites = p.RM_MunI_RecognitionSites;
            unmethylatedSites   = find(collapse(unmethylatedSites, -1));
            hemimethylatedSites = find(collapse(hemimethylatedSites, -1));
            methylatedSites     = find(collapse(methylatedSites, -1));
            cleavedSites        = find(collapse(cleavedSites, -1));
            inaccessibleRegions = find(collapse(inaccessibleRegions, -1));
            
            %MunI sites
            [x, y] = pol2cart(sites(:, [1 end])' / genomeLength * 2 * pi, repmat(1 + strandSeparation, 2, size(sites, 1)));
            line(x, y, 'Parent',axesHandle, 'Color',[0.6 0.6 0.6], 'LineWidth',0.1); %(+) strand
            
            [x, y] = pol2cart(sites(:, [1 end])' / genomeLength * 2 * pi, repmat(1 - strandSeparation, 2, size(sites, 1)));
            tmpHandle = line(x, y, 'Parent', axesHandle, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.1); %(-) strand
            
            h = tmpHandle(1);
            legendLabels = {sprintf('MunI Sites (%d)', size(sites, 1))};
            
            %unmethylated sites
            if ~isempty(unmethylatedSites)
                [x, y] = pol2cart(sites(unmethylatedSites, [1 end])' / genomeLength * 2 * pi, repmat(1 + strandSeparation, 2, size(unmethylatedSites, 1)));
                line(x, y, 'Parent', axesHandle, 'Color', [1 0 1], 'LineWidth', 0.1); %(+) strand
                
                [x, y] = pol2cart(sites(unmethylatedSites, [1 end])' / genomeLength * 2 * pi, repmat(1 - strandSeparation, 2, size(unmethylatedSites, 1)));
                tmpHandle = line(x, y, 'Parent', axesHandle, 'Color', [1 0 1], 'LineWidth', 0.1); %(-) strand
                
                h(end + 1) = tmpHandle(1);
                legendLabels{end + 1} = sprintf('Unmethylated (%d)', numel(unmethylatedSites));
            end
            
            %inaccessible sites
            if ~isempty(inaccessibleRegions)
                [x, y] = pol2cart(sites(inaccessibleRegions, [1 end])' / genomeLength * 2 * pi, repmat(1 + strandSeparation, 2, size(inaccessibleRegions, 1)));
                line(x, y, 'Parent',axesHandle, 'Color',[1 0 0], 'LineWidth',0.1); %(+) strand
                
                [x, y] = pol2cart(sites(inaccessibleRegions, [1 end])' / genomeLength * 2 * pi, repmat(1 - strandSeparation, 2, size(inaccessibleRegions, 1)));
                tmpHandle = line(x, y, 'Parent', axesHandle, 'Color', [1 0 0], 'LineWidth', 0.1); %(-) strand
                
                h(end + 1) = tmpHandle(1);
                legendLabels{end + 1} = sprintf('Inaccessible (%d)', numel(inaccessibleRegions));
            end
            
            %hemimethylated sites
            if ~isempty(hemimethylatedSites)
                strand1 = find(collapse(p.chromosome.damagedBases(sites(hemimethylatedSites, methylatedPositions(1)), 1, :) == p.substrateGlobalIndexs(p.substrateIndexs_m6AD), -1));
                strand2 = find(collapse(p.chromosome.damagedBases(sites(hemimethylatedSites, methylatedPositions(2)), 2, :) == p.substrateGlobalIndexs(p.substrateIndexs_m6AD), -1));
                
                if ~isempty(strand1)
                    [x, y] = pol2cart(repmat(sites(hemimethylatedSites(strand1), methylatedPositions(1)) / genomeLength * 2 * pi, 1, 2),...
                        repmat(1 + strandSeparation + [0 damageLength], length(strand1), 1));
                    tmpHandle = line(x', y', 'Parent', axesHandle, 'Color', [0 0 1], 'LineWidth', 0.1);
                    
                    h(end + 1) = tmpHandle(1);
                    legendLabels{end + 1} = sprintf('Hemimethylated (%d)', numel(strand1) + numel(strand2));
                end
                
                if ~isempty(strand2)
                    [x, y] = pol2cart(repmat(sites(hemimethylatedSites(strand2), methylatedPositions(2)) / genomeLength * 2 * pi, 1, 2),...
                        repmat(1 - (strandSeparation + [0 damageLength]), length(strand2), 1));
                    tmpHandle = line(x', y', 'Parent', axesHandle, 'Color', [0 0 1], 'LineWidth', 0.1);
                    
                    if isempty(strand1)
                        h(end + 1) = tmpHandle(1);
                        legendLabels{end + 1} = sprintf('Hemimethylated (%d)', numel(strand2));
                    end
                end
            end
            
            %methylated sites
            if ~isempty(methylatedSites)
                [x, y] = pol2cart(repmat(sites(methylatedSites, methylatedSites(1))/genomeLength*2*pi, 1, 2),...
                    repmat(1+strandSeparation + [0 damageLength], length(methylatedSites), 1));
                line(x', y', 'Parent', axesHandle, 'Color', [0 1 0], 'LineWidth', 0.1);
                
                [x, y] = pol2cart(repmat(sites(methylatedSites, methylatedSites(2))/genomeLength*2*pi, 1, 2),...
                    repmat(1-(strandSeparation + [0 damageLength]), length(methylatedSites), 1));
                tmpHandle = line(x', y', 'Parent', axesHandle, 'Color', [0 1 0], 'LineWidth', 0.1);
                
                h(end + 1) = tmpHandle(1);
                legendLabels{end + 1} = sprintf('Methylated (%d)', numel(methylatedSites));
            end
            
            %cleaved sites
            if ~isempty(cleavedSites)
                [x, y] = pol2cart(repmat(sites(cleavedSites, methylatedPositions(1)) / genomeLength * 2 * pi, 1, 2),...
                    repmat(1 + strandSeparation + [0 damageLength], length(cleavedSites), 1));
                line(x', y', 'Parent', axesHandle, 'Color', [0, 183, 235] / 255, 'LineWidth', 0.1);
                
                [x, y] = pol2cart(repmat(sites(cleavedSites,methylatedPositions(2)) / genomeLength * 2 * pi, 1, 2),...
                    repmat(1 - (strandSeparation + [0 damageLength]), length(cleavedSites), 1));
                tmpHandle = line(x', y', 'Parent', axesHandle, 'Color', [0, 183, 235] / 255, 'LineWidth', 0.1);
                
                h(end + 1) = tmpHandle(1);
                legendLabels{end+1} = sprintf('Cleaved (%d)', numel(cleavedSites));
            end
            
            %axes properties
            xlim((1 + 2 * damageLength) * [-1 1]);
            ylim((1 + 2 * damageLength) * [-1 1]);
            set(axesHandle, ...
                'Box', 'off', ...
                'Color', 'none', ...
                'Visible', 'off', ...
                'DataAspectRatio', [1 1 1], ...
                'XTick', [], ...
                'YTick', [], ...
                'XMinorTick', 'off', ...
                'YMinorTick', 'off', ...
                'GridLineStyle', 'none')
            
            %legend
            lgd = legend(h, legendLabels{:});
            position = get(lgd, 'Position');
            set(lgd, 'Position', [(1 - position(3:4)) / 2 position(3:4)]);
            legend(axesHandle, 'show');
        end
    end
end