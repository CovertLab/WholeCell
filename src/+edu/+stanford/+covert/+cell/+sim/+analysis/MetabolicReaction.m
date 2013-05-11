%MetabolicReaction
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/6/2011
classdef MetabolicReaction
    %plotting metabolic reaction
    methods (Static = true)
        % plot growth
        function plotGrowth_Rate(sim, axesHandle, time, ~)
            %import classes
            import edu.stanford.covert.util.ConstantUtil;
            
            s = sim.state('MetabolicReaction');
            
            plot(axesHandle, time/3600, permute(s.growth,[3 1 2]) * ConstantUtil.secondsPerHour);
            xlabel(axesHandle, 'Time (h)','fontSize',8);
            ylabel(axesHandle, 'Growth (cell/h)','fontSize',8);
            xlim(axesHandle, [0 max(1, time(end))]/3600);
        end
        
        % plot doubling time
        function plotDoubling_Time(sim, axesHandle, time, ~)
            %import classes
            import edu.stanford.covert.util.ConstantUtil;
            
            s = sim.state('MetabolicReaction');
            
            plot(axesHandle, time/3600, permute(s.doublingTime,[3 1 2]) / ConstantUtil.secondsPerHour);
            xlabel(axesHandle, 'Time (h)', 'fontSize', 8);
            ylabel(axesHandle, 'Doubling Time (h)', 'fontSize', 8);
            xlim(axesHandle, [0 max(1, time(end))] / 3600);
        end
        
        %Plots fluxes of process's reactions as line graph with 1 line per
        %reaction. Each line shows the flux of the reaction versus time
        function plotReactionFluxes(sim, axesHandle, time, ~)
            s = sim.state('MetabolicReaction');
            
            plotHandles = plot(axesHandle, time/3600, permute(s.fluxs, [3 1 2]));
            if length(s.reactionWholeCellModelIDs) <= 10
                legend(plotHandles, s.reactionNames);
            end
            xlabel(axesHandle, 'Time (h)','fontSize',8);
            ylabel(axesHandle, 'Reaction Fluxes','fontSize',8);
            xlim(axesHandle, [0 max(1, time(end))] / 3600);
        end
    end
end