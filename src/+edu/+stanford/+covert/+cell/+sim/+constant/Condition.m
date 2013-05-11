%Condition
%
%  Author: Jonathan Karr
%  Affilitation: Covert Lab, Department of Bioengineering, Stanford University
%  Last updated: 1/10/2011
classdef Condition
    %indices of specific states within state vectors
    properties (Constant)
        objectIndexs            = 1;
        compartmentIndexs       = 2;
        valueIndexs             = 3;
        initialTimeIndexs       = 4;
        finalTimeIndexs         = 5;
        objectCompartmentIndexs = 6;
    end
    
    methods (Static)
        %Applies conditions at current time point
        function values = applyConditions(values, conditions, time)
            import edu.stanford.covert.cell.sim.constant.Condition;
            
            conditions = conditions(...
                time >= conditions(:, Condition.initialTimeIndexs) & ...
                time <= conditions(:, Condition.finalTimeIndexs), :);
            
            values(conditions(:, Condition.objectCompartmentIndexs)) = ...
                conditions(:, Condition.valueIndexs);
        end
    end
end