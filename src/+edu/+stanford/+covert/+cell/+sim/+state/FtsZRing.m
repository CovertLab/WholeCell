%FtsZRing
%
% @wholeCellModelID State_FtsZRing
% @name             FtsZ ring
% @description
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 11/30/2010
classdef FtsZRing < edu.stanford.covert.cell.sim.CellState
    %property annotations
    properties (Constant)
        optionNames             = {   %names of properties that are options
            'verbosity';
            'seed';
            };
        fixedConstantNames      = {   %names of process properties that are considered fixed constants
            'numFtsZSubunitsPerNm';
            'numFtsZSubunitsPerFilament';
            'filamentLengthInNm';
            };
        fittedConstantNames     = {}; %names of process properties that are considered fitted constants, and should be stored with the simulation as such
        stateNames              = {   %names of properties which are part of the simulation's state
            'numEdgesOneStraight';
            'numEdgesTwoStraight';
            'numEdgesTwoBent';
            'numResidualBent'
            };
        dependentStateNames     = {   %names of properties which can be calculated from the simulation's state
            'numEdges'
            };
    end
    
    properties (Constant)
        dryWeight = 0; %dry weight of this class' state properties
    end
    
    %fixed biological constants
    properties
        numFtsZSubunitsPerNm            %Linear density of FtsZ monomers in filaments (1/nm) [0.23 1/nm; Anderson 2004]
        numFtsZSubunitsPerFilament      %FtsZ filament length (no. monomer subunits) [9; Anderson 2004]
        filamentLengthInNm              %FtsZ filament length (nm) [40 nm; Anderson 2004]
    end
    
    %state
    properties
        numEdgesOneStraight      %number of edges currently bound by one straight ftsZ polymer
        numEdgesTwoStraight      %number of edges currently bound by two straight ftsZ polymers
        numEdgesTwoBent          %number of edges currently bound by two bent ftsZ polymers
        numResidualBent          %number of residual (singly) bound ftsZ polymers from previous cycle
    end
    
    properties (Dependent = true, SetAccess = protected)
        numEdges
    end
    
    %references to other parts of cell state
    properties
        geometry
    end
    
    %constructor
    methods
        function this = FtsZRing(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.CellState(wholeCellModelID, name);
        end
    end
    
    methods
        function storeObjectReferences(this, simulation)
            this.geometry = simulation.state('Geometry');
        end
    end
    
    %communication between process/simulation
    methods
        function initializeConstants(this, knowledgeBase, simulation)
            this.initializeConstants@edu.stanford.covert.cell.sim.CellState(knowledgeBase, simulation);
            this.filamentLengthInNm = this.numFtsZSubunitsPerFilament / this.numFtsZSubunitsPerNm;
        end
    end
    
    %allocate memory for state
    methods
        function allocateMemory(this, numTimePoints)
            this.numEdgesOneStraight = zeros(1, 1, numTimePoints);
            this.numEdgesTwoStraight = zeros(1, 1, numTimePoints);
            this.numEdgesTwoBent     = zeros(1, 1, numTimePoints);
            this.numResidualBent     = zeros(1, 1, numTimePoints);
        end
    end
    
    %initialization
    methods
        %initialize to no bound FtsZ rings
        function initialize(this)
            this.allocateMemory(1);
        end
    end
    
    methods
        %nFtsZ = [nFtsZGTP; nFtsZGDP]
        function notUpdatingFtsZ = releaseFtsZ(~, nFtsZ)
            notUpdatingFtsZ = nFtsZ;
            if any(nFtsZ)
                warning('WholeCell:warning', 'FtsZ not decayed');
            end
        end
    end
    
    %getters
    methods
        function result = get.numEdges(this)
            g = this.geometry;
            if g.pinchedDiameter == 0
                result = 0;
                return;
            end
            
            result = this.calcNumEdges(g.pinchedDiameter, this.filamentLengthInNm);
        end
    end
    
    methods (Static)
        function result = calcNumEdges(pinchedDiameter, filamentLengthInNm)
            result = floor(pi / asin(filamentLengthInNm * 1e-9 / pinchedDiameter));
        end
    end
end
