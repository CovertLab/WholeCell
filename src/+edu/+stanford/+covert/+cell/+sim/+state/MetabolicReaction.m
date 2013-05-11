%MetabolicReaction
%
% @wholeCellModelID State_MetabolicReaction
% @name             Metabolic reaction
% @description
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/5/2011
classdef MetabolicReaction < edu.stanford.covert.cell.sim.CellState
    %property annotations
    properties (Constant)
        optionNames             = {   %names of properties that are options
            'verbosity';
            'seed';
            };
        fixedConstantNames      = {}; %names of process properties that are considered fixed constants
        fittedConstantNames     = {}; %names of process properties that are considered fitted constants, and should be stored with the simulation as such
        stateNames              = {   %names of properties which are part of the simulation's state
            'growth';
            'fluxs';
            };
        dependentStateNames     = {   %names of properties which can be calculated from the simulation's state
            'doublingTime'
            };
    end
    
    %constants
    properties
        reactionWholeCellModelIDs  %reaction whole cell model ids
        growth0                    %target growth rate
        initialGrowthFilterWidth   %tolerance of initial growth rate
        meanInitialGrowthRate      %mean initial growth rate
    end
    
    %state
    properties
        growth  %growth rate (cell/s)
        fluxs   %flux of each reaction, size: [length(reactionWholeCellModelIDs) X 1 X time]
    end
    
    %dependent state
    properties (Dependent = true, SetAccess = protected)
        doublingTime  %doubling time (s)
    end
    
    %dependent state
    properties (Constant)
        dryWeight = 0; %dry weight of this class' state properties
    end
    
    %constructor
    methods
        function this = MetabolicReaction(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.CellState(wholeCellModelID, name);
        end
    end
    
    methods
        function initializeConstants(this, knowledgeBase, simulation)
            this.initializeConstants@edu.stanford.covert.cell.sim.CellState(knowledgeBase, simulation);
            
            m = findobj(knowledgeBase.processes, 'wholeCellModelID', 'Process_Metabolism');
            this.reactionWholeCellModelIDs = {m.reactions.wholeCellModelID}';
        end
    end
    
    %allocate memory for state
    methods
        function allocateMemory(this, numTimePoints)
            this.growth = zeros(1, 1, numTimePoints);
            this.fluxs  = zeros(numel(this.reactionWholeCellModelIDs), 1, numTimePoints);
        end
    end
    
    %initialize state
    methods
        function initialize(~)
        end
    end
    
    %getters
    methods
        %doubling time (seconds/cell)
        function value = get.doublingTime(this)
            value = 1 ./ this.growth;
        end
    end
end
