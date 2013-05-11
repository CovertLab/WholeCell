% Gene test cases
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 9/5/2010
classdef Gene_Test < TestCase
    %properties
    properties
        gene
    end
    
    %constructor
    methods
        function this = Gene_Test(name)
            this = this@TestCase(name);
        end
    end
    
    %setup
    methods
        function setUp(this)
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], true);
            this.gene = sim.gene;
        end
    end
    
    %tests
    methods
        function testGetIndexsByPosition(this)
            g = this.gene;
            assertEqual(10, g.getIndexsByPosition(g.startCoordinates(10)));
        end
    end
end
