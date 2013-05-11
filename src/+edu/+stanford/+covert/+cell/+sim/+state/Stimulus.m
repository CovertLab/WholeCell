%Stimulus
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 1/5/2011
classdef Stimulus < edu.stanford.covert.cell.sim.CellState
    %Constants
    properties (Constant)
        optionNames             = {   %names of properties that are options
            'verbosity';
            'seed';
            };
        fixedConstantNames      = {}; %names of process properties that are considered fixed constants
        fittedConstantNames     = {}; %names of process properties that are considered fitted constants, and should be stored with the simulation as such
        stateNames              = {   %names of properties which are part of the simulation's state
            'values'};
        dependentStateNames     = {}; %names of properties which can be calculated from the simulation's state
    end
    
    %constants
    properties
        wholeCellModelIDs
        names
        setValues
        
        compartment
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
        function this = Stimulus(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.CellState(wholeCellModelID, name);
        end
    end
    
    methods
        function storeObjectReferences(this, simulation)
            this.compartment = simulation.compartment;
        end
        
        function initializeConstants(this, knowledgeBase, simulation)
            import edu.stanford.covert.cell.sim.constant.Condition;
            
            this.initializeConstants@edu.stanford.covert.cell.sim.CellState(knowledgeBase, simulation);            
            
            this.wholeCellModelIDs  = {knowledgeBase.stimulis.wholeCellModelID}';
            this.names  = {knowledgeBase.stimulis.name}';
            this.setValues = zeros(0, 6);
            
            for i = 1:knowledgeBase.numStimulis
                value = zeros(length(knowledgeBase.stimulis(i).compartments), 6);
                value(:, Condition.objectIndexs)      = knowledgeBase.stimulis(i).idx;
                value(:, Condition.compartmentIndexs) = [knowledgeBase.stimulis(i).compartments.idx];
                value(:, Condition.valueIndexs)       = knowledgeBase.stimulis(i).values;
                value(:, Condition.initialTimeIndexs) = knowledgeBase.stimulis(i).initialTimes;
                value(:, Condition.finalTimeIndexs)   = knowledgeBase.stimulis(i).finalTimes;
                value(:, Condition.objectCompartmentIndexs) = sub2ind(...
                    [knowledgeBase.numStimulis knowledgeBase.numCompartments],...
                    value(:, Condition.objectIndexs),...
                    value(:, Condition.compartmentIndexs));
                
                this.setValues = [this.setValues; value];
            end
        end
        
        %allocate memory for state
        function allocateMemory(this, numTimePoints)
            this.values = zeros([numel(this.wholeCellModelIDs), this.compartment.count, numTimePoints]);
        end
                
        %initialize state
        function initialize(this)
            this.allocateMemory(1);
        end
    end
    
    %helper methods
    methods
        function value = getIndexs(this, wholeCellModelIDs)
            [~, value] = ismember(wholeCellModelIDs, this.wholeCellModelIDs);
        end
    end
end
