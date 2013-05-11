classdef parseDoublesTest < TestCase    
    methods
        function this = parseDoublesTest(methodName)
            this = this@TestCase(methodName);
        end
        
        function testNoValues(~)
            assertEqual([], parse(',', ''));
        end
        
        function testOneValue(~)
            assertEqual(3.14, parse(',', '3.14'));
        end
        
        function testOneValueTrailingDelimiter(~)
            assertEqual(3.14, parse(',', '3.14,'));
        end
        
        function testThreeValues(~)
            assertEqual([1 2 -1], parse(',', '1,2,-1'));
        end
        
        function testThreeValuesTrailingDelimiter(~)
            assertEqual([1 2 -1], parse(',', '1,2,-1,'));
        end
        
        function testLongDelimiter(~)
            assertEqual([2e9 -10.2], parse(',_.', '2e9,_.-10.2'));
        end
        
        function testLongTrailingDelimiter(~)
            assertEqual([2e9 -10.2], parse(',_.', '2e9,_.-10.2,_.'));
        end
    end
end

function result = parse(delimiter, s)
    result = edu.stanford.covert.util.parseDoubles(delimiter, s);
end
