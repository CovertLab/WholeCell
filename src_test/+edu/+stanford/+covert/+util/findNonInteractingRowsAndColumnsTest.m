classdef findNonInteractingRowsAndColumnsTest < TestCase
    methods
        function this = findNonInteractingRowsAndColumnsTest(methodName)
            this = this@TestCase(methodName);
        end
        
        function test1(~)
            [rowAssignments, colAssignments, blocks] = ...
                edu.stanford.covert.util.findNonInteractingRowsAndColumns([
                    1 0 1;
                    0 0 0;
                    1 0 0]);
            assertEqual([2 1 2], rowAssignments');
            assertEqual([2 1 2], colAssignments');
            assertEqual({0; [1 1; 1 0]}, blocks);
        end
        
        function test2(~)
            [rowAssignments, colAssignments, blocks] = ...
                edu.stanford.covert.util.findNonInteractingRowsAndColumns([
                    0 0 1;
                    0 1 0;
                    1 0 0]);
            assertEqual([1 2 3], rowAssignments');
            assertEqual([3 2 1], colAssignments');
            assertEqual({1; 1; 1}, blocks);
        end
        
        function test3(~)
            [rowAssignments, colAssignments, blocks] = ...
                edu.stanford.covert.util.findNonInteractingRowsAndColumns([
                    1 0 1;
                    0 1 0;
                    1 0 0]);
            assertEqual([1 2 1], rowAssignments');
            assertEqual([1 2 1], colAssignments');
            assertEqual({[1 1; 1 0]; 1}, blocks);
        end
    end
end
