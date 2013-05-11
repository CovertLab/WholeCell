%SimulationStateSideEffect
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/12/2010
classdef SimulationStateSideEffect
    properties
        items  %SimulationStateSideEffectItem array
    end
    
    %constructor
    methods
        function this = SimulationStateSideEffect(items)
            import edu.stanford.covert.cell.sim.SimulationStateSideEffectItem;
            
            if nargin == 0; return; end; %to permit allocating arrays
            
            if ~isa(items, 'SimulationStateSideEffectItem')
                throw(MException('SimulationStateSideEffect:invalidInput', 'item must be a SimulationStateSideEffectItem array'));
            end
            
            this.items = items;
        end
    end
    
    %method to update simulation state
    methods
        function simulation = updateSimulationState(this, simulation)
            for i = 1:numel(this)
                simulation = this(i).items.updateSimulationState(simulation);
            end
        end
    end
    
    %print to standard output
    methods
        function disp(this)
            sizCellArr = cellfun(@num2str, num2cell(size(this)), 'UniformOutput', false);
            fprintf('%s SimulationStateSideEffect\n', strjoin('x', sizCellArr{:}));
            for i = 1:numel(this)
                disp(this(i).items, 1);
            end
        end
    end
end