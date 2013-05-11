%Time
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 1/5/2011
classdef Time < edu.stanford.covert.cell.sim.CellState
    %Annotations
    properties (Constant)
        optionNames             = {   %names of properties that are options
            'verbosity';
            'seed';
            };
        fixedConstantNames      = {   %names of process properties that are considered fixed constants
            'cellCycleLength';            
            'replicationDuration';
            'cytokinesisDuration';
            };
        fittedConstantNames     = {   %names of process properties that are considered fitted constants, and should be stored with the simulation as such
            'replicationInitiationDuration';
            };
        stateNames              = {   %names of properties which are part of the simulation's state
            'values'};
        dependentStateNames     = {}; %names of properties which can be calculated from the simulation's state
    end
    
    %constants
    properties
        cellCycleLength                  %length of cell cycle (s) [PUB_0094, PUB_0556]
        replicationInitiationDuration    %Duration of the replication initiation phase of the cell cycle (s)
        replicationDuration              %Duration of the replication phase of the cell cycle (s)
        cytokinesisDuration              %Duration of the cytokinesis phase of the cell cycle (s)
    end
    
    %state
    properties
        values
    end
    
    %alternative interface to state
    properties (Constant)
        dryWeight = 0
    end
    
    %constructor
    methods
        function this = Time(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.CellState(wholeCellModelID, name);
        end
    end
    
    methods
        %allocate memory for state
        function allocateMemory(this, numTimePoints)
            this.values = zeros(1, 1, numTimePoints);
        end
        
        %initialize state
        function initialize(this)
            this.allocateMemory(1);
        end
    end
end
