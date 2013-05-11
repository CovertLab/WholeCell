%MoleculeCountState
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 1/5/2011
classdef MoleculeCountState < edu.stanford.covert.cell.sim.CellState
    %Constants
    properties (Constant, Abstract)
        fixedConstantNames            %names of process properties that are considered fixed constants
        fittedConstantNames           %names of process properties that are considered fitted constants, and should be stored with the simulation as such
        stateNames                    %names of properties which are part of the simulation's state
    end
    
    %constants
    properties
        wholeCellModelIDs
        names
        molecularWeights
        baseCounts
        lengths
        halfLives
        compartments
        
        compartment
    end
    
    properties (Dependent)
        decayRates
    end
    
    %state
    properties
        counts
    end
    
    %alternative interface to state
    properties (Dependent = true, SetAccess = protected)
        dryWeight
    end   
    
    %constructor
    methods
        function this = MoleculeCountState(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.CellState(wholeCellModelID, name);
        end
    end
    
    methods
        function storeObjectReferences(this, simulation)
            this.compartment = simulation.compartment;
        end       
        
        %allocate memory for state
        function allocateMemory(this, numTimePoints)
            this.counts = zeros([numel(this.wholeCellModelIDs), this.compartment.count, numTimePoints]);
        end
                
        %initialize state
        function initialize(this)
            this.allocateMemory(1);
        end
    end
    
    %helper functions
    methods
        function value = getIndexs(this, wholeCellModelIDs)
            [~, value] = ismember(wholeCellModelIDs, this.wholeCellModelIDs);
        end
    end
    
    %getters
    methods            
        %decay rates of monomers (molecules/s)
        function value = get.decayRates(this)
            value = log(2) ./ this.halfLives;
        end
        
        function set.decayRates(this, value)
            this.halfLives = log(2) ./ value;
        end
        
        function value = get.dryWeight(this)
            value = this.calcDryWeight();
        end
        
        function value = calcDryWeight(this)
            if size(this.counts, 3) == 1
                value = this.molecularWeights' * this.counts;
            else
                value = multiprod(this.molecularWeights', this.counts, [1 2], [1 2]);
            end
            value = value / edu.stanford.covert.util.ConstantUtil.nAvogadro;
        end
    end
end
