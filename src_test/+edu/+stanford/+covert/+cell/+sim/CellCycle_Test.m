%CellCycle medium test
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 11/23/2010
classdef CellCycle_Test < TestCase
    properties
        simulation
    end

    methods
        function this = CellCycle_Test(name)
            this = this@TestCase(name);
        end

        function setUp(this)
            this.simulation = edu.stanford.covert.cell.sim.SimulationFixture.load();
        end
        
        %Test that the expected lengths of the 3 phases (replication initation,
        %replication, cytokinesis) of the cell cycle sum to the expected cell
        %cycle length.
        %
        %The expected lengths of the 3 phases of the cell cycle have been
        %constructed as follows:
        %1. The expected length of the replication phase was calculated
        %   analytically using experimentally observed values for the
        %   replication constants.
        %2. The expected length of the cytokinesis phase was estimated by
        %   simulation. This is implemented in the cytokinesis medium test.
        %3. The expected length of the replication initiation phase was set to
        %   be the difference of the total cell cycle length and that of the
        %   replication and cytokinesis phases, and the replication initation
        %   constants (in particular the site and state cooperativity constants)
        %   were fit to achieve this expected replication initation length. This
        %   is implemented in the replication initation process test class.
        function testCellCycleLength(this)
            sim = this.simulation;
            
            %% lengths of individual phases of cell cycle
            %1. Replication initiation: replication initiation test asserts that
            %   replication initiation takes the expected time 
            %2. Replication: replication unit test asserts that replication
            %   takes the expected time
            %3. Cytokinesis: Cytokinesis medium test asserts that cytokinesis
            %   takes the expected time
            
            %% total length of cell cycle
            time = sim.state('Time');
            assertEqual(time.cellCycleLength, ...
                + time.replicationInitiationDuration + ...
                + time.replicationDuration + ...
                + time.cytokinesisDuration);
        end
    end
end
