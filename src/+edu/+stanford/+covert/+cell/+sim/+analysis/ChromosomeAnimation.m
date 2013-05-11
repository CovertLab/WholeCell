%ChromosomeAnimation
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 7/20/2011
classdef ChromosomeAnimation < edu.stanford.covert.cell.sim.analysis.Animation
    properties (SetAccess = protected)
        title = 'Whole Cell DNA Animation'
        description = 'Animation of dynamics of DNA in a single Mycoplasma genitalium whole-cell simulation.'
        author = 'Jonathan Karr, Covert Lab, Department of Bioengineering, Stanford University'        
    end
    
    methods
        function this = ChromosomeAnimation(varargin)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            this = this@edu.stanford.covert.cell.sim.analysis.Animation(varargin{:});
            
            this.movieWidth = round((600/90 * this.movieDPI) / 2) * 2;
            this.movieHeight = round((300/90 * this.movieDPI) / 2) * 2;
        end
        
        function [times, states] = loadSimulationData(this)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            stateNames = {
                'Time'       'values'
                'Chromosome' 'polymerizedRegions'
                'Chromosome' 'monomerBoundSites'
                'Chromosome' 'complexBoundSites'
                };
            simDir = SimulationDiskUtil.getSimulation([this.simGroupId filesep num2str(this.simId)]);
            states = DiskLogger.load(simDir, stateNames, [], [], [], 'extract');
            times = states.Time.values;
        end
        
        function svg = drawFrame(this, timeIdx)
            chr = this.simulation.state('Chromosome');
            
            %open svg
            svg = [];
            svg = [svg sprintf('<?xml version="1.0" standalone="no"?>\n')];
            svg = [svg sprintf('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n')];
            svg = [svg sprintf('<svg width="%d" height="%d" version="1.1" xmlns="http://www.w3.org/2000/svg" viewBox="%d %d %d %d">\n', 600, 300, 0, 0, 600, 300)];
            
            %white background
            svg = [svg sprintf('<rect x="%d" y="%d" width="%d" height="%d" style="fill:white; stroke:none"/>\n', ...
                0, 0, 600, 300)];
            
            %polymerized Regions
            [posStrnds, lens] = find(this.states.Chromosome.polymerizedRegions);
            lens = lens(posStrnds(:, 3) == timeIdx, :);
            posStrnds = posStrnds(posStrnds(:, 3) == timeIdx, 1:2);
            svg = [svg sprintf(this.drawArc(posStrnds, lens, [95 105], 'black', 2))];
            
            %bound proteins
            [monPosStrnds, monIdxs] = find(this.states.Chromosome.monomerBoundSites);
            [cpxPosStrnds, cpxIdxs] = find(this.states.Chromosome.complexBoundSites);
            posStrnds = [
                monPosStrnds(monPosStrnds(:, 3) == timeIdx, 1:2)
                cpxPosStrnds(cpxPosStrnds(:, 3) == timeIdx, 1:2)
                ];
            lens = [
                chr.monomerDNAFootprints(monIdxs(monPosStrnds(:, 3) == timeIdx, :))
                chr.complexDNAFootprints(cpxIdxs(cpxPosStrnds(:, 3) == timeIdx, :))
                ];
            svg = [svg sprintf(this.drawArc(posStrnds, lens, [90 110], 'red', 1))];
            
            %close svg
            svg = [svg sprintf('</svg>\n')];
        end
        
        function svg = drawArc(this, posStrnds, lens, r0, strokeColor, strokeWidth)
            chr = this.simulation.state('Chromosome');
            
            svg = [];
            for j = 1:size(posStrnds, 1)
                theta = posStrnds(j, 1) / chr.sequenceLen * 2 * pi - pi/2;
                dTheta = lens(j) / chr.sequenceLen * 2 * pi;
                largeArc = dTheta >= pi;
                r = r0(iseven(posStrnds(j, 2)) + 1);
                if posStrnds(j, 2) <= 2
                    cx = 50 + 100;
                else
                    cx = 3*50 + 3*100;
                end
                cy = 50 + 100;
                if lens(j) == chr.sequenceLen
                    svg = [svg sprintf('<circle cx="%d" cy="%d" r="%d" stroke="%s" stroke-width="%f" fill="none"/>\n', ...
                        cx, cy, r, ...
                        strokeColor, strokeWidth)]; %#ok<*AGROW>
                else
                    svg = [svg sprintf('<path d="M%f,%f A%f,%f %d %d,%d %f %f" stroke="%s" stroke-width="%d" fill="none"/>\n', ...
                        cx + r*cos(theta), cy + r*sin(theta), ...
                        r, r, 0, largeArc, 1, cx + r*cos(theta+dTheta), cy + r*sin(theta+dTheta), ...
                        strokeColor, strokeWidth)];
                end
            end
        end
    end
end