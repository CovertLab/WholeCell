classdef SetUtil
    %unique
    %- inputs should be finite
    %- these methods will behave differently from the builtin unique method if
    %  the inputs are infinite
    methods (Static)
        function x = unique(x)
            x = sort(x(:));
            x = x([diff(x) ~= 0; true], 1);
        end
        
        function tf = isunique(x)
            tf = all(diff(sort(x(:))));
        end
        
        function x = unique_rows(x)
            x = sortrows(x);
            x = x([any(diff(x, 1), 2); true], :);
        end
        
        function tf = isunique_rows(x)
            tf = all(any(diff(sortrows(x), 1), 2));
        end
    end
    
    %ismember
    methods (Static)
        function ismember
        end
    end
end