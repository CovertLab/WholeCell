%SimulationStateSideEffectItem
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/12/2010
classdef SimulationStateSideEffectItem
    properties (SetAccess = protected)
        stateName          %name of simulation state to be updated
        propertyName       %name of state property to be updated
        indexName          %name of index over state property
        componentIndex     %component (row) index within state property/property index to be updated
        compartmentIndex   %compartment (column) index within state property to be updated
        deltaValue         %change in value of state property element
    end
    
    %constructor
    methods
        function this = SimulationStateSideEffectItem(stateName, propertyName, indexName, componentIndex, compartmentIndex, deltaValue)
            this.stateName        = stateName;
            this.propertyName     = propertyName;
            this.indexName        = indexName;
            this.componentIndex   = componentIndex;
            this.compartmentIndex = compartmentIndex;
            this.deltaValue       = deltaValue;
        end
    end
    
    %method to update simulation state
    methods
        function simulation = updateSimulationState(this, simulation)
            deltaMass = 0;
            
            for i = 1:numel(this)
                state = simulation.state(this(i).stateName);
                
                if this(i).indexName
                    iComponent = state.(this(i).indexName)(this(i).componentIndex);
                else
                    iComponent = this(i).componentIndex;
                end
                
                iCompartment = this(i).compartmentIndex;
                state.(this(i).propertyName)(iComponent, iCompartment) = ...
                    state.(this(i).propertyName)(iComponent, iCompartment) + ...
                    this(i).deltaValue;
                deltaMass = deltaMass + ...
                    state.molecularWeights(iComponent) * this(i).deltaValue;
            end
            
            if deltaMass ~= 0
                warning('WholeCell:warning',[...
                    'Simulation state side effects should preserve mass (eg. '...
                    'mass should not be passed from processes state properties to '...
                    'SimulationStateSideEffect objects, and each SimulationStateSideEffectItem '...
                    'array should represent a mass-balanced change to the state of '...
                    'the simulation)'])
            end
        end
    end
    
    %print to standard output
    methods 
        function disp(this, indent)
            if nargin < 2
                indent = 0;
            end
            
            sizCellArr = cellfun(@num2str, num2cell(size(this)), 'UniformOutput', false);
            fprintf([repmat('  ',1, indent) '%s SimulationStateSideEffectItem\n'], ...
                strjoin('x', sizCellArr{:}));
            stateNameMaxLength = max(cellfun(@length, {this.stateName}));
            propertyNameMaxLength = max(cellfun(@length, {this.propertyName}));
            indexNameMaxLength  = max(cellfun(@length, {this.indexName}));
            for i = 1:numel(this)
                fprintf([repmat('  ',1, indent+1) '%' num2str(stateNameMaxLength) 's\t%' num2str(propertyNameMaxLength) 's\t%' num2str(indexNameMaxLength) 's\t%8d\t%8d\t%4d\n'], ...
                    this(i).stateName, this(i).propertyName, this(i).indexName, this(i).componentIndex, this(i).compartmentIndex, this(i).deltaValue);
            end
        end
    end
end