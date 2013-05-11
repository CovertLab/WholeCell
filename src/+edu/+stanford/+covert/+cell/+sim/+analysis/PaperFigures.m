% Figures for first whole cell paper
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 12/19/2011
classdef PaperFigures
    methods (Static)
        function run()
            %import classes
            import edu.stanford.covert.cell.sim.analysis.Figures34;
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            
            %base directory for output
            baseDir = 'documentation/paper/figures/';
            
            %% Main paper
            %figure 1: designed with Adobe Illustrator
            %- documentation/paper/figures/model overview.ai
            %- documentation/paper/figures/algorithm overview with process-state network.ai
            
            %figure 2, 3, 5: plotted entirely with MATLAB
            Figures34.run('documentation/paper/figures');
            
            %figure 4: combination of MATLAB, illustrator
            Figures34.run('documentation/paper/figures');
            
            %figure 6: plotted entirely with MATLAB
            SingleGeneDeletions.run();                %run analysis and cache results
            SingleGeneDeletions.plotOverview(...      %plot overview figure from cached results
                [baseDir filesep 'figure6.pdf'], ...
                [baseDir filesep 'figure6.xls']);
            
            %figure 7: combination of MATLAB, illustrator
            
            
            %% Supplement
            %state, process schematics: designed with Adobe Illustrator
        end
    end
end