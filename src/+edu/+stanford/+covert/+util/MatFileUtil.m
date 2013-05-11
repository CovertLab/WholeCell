% Mat file utility functions
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 11/23/2010
classdef MatFileUtil
    methods (Static)
        function diff(file1, file2)
            struct1 = load(file1);
            struct2 = load(file2);
            
            edu.stanford.covert.util.StructUtil.diff(struct1, struct2);
        end
    end
end