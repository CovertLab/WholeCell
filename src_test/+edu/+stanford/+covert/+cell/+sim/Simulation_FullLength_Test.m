%Full Length simulation test
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/12/2011
classdef Simulation_FullLength_Test < TestCase
    methods
        function this = Simulation_FullLength_Test(name)
            this = this@TestCase(name);
        end
    end
    
    methods
        function testSimulation(~)
            runSimulation('outDir', 'output/runSimulationTests', 'logToDisk', true, 'logToDb', false);
        end
    end
end