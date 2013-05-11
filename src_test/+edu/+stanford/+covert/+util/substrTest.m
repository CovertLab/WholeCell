classdef substrTest < TestCase    
    methods
        function this = substrTest(methodName)
            this = this@TestCase(methodName);
        end
        
        function testEmptyString(~)
            s = '';
            assertEqual(s(1:0), substr('', 1, 0));
            assertExceptionThrown(@() substr('', 1, 1), 'MATLAB:badsubscript');
            assertExceptionThrown(@() substr('', 0, 0), 'MATLAB:badsubscript');
        end
        
        function testWholeString(~)
            assertEqual('abc', substr('abc', 1, 0));
        end
        
        function testPartOfString(~)
            assertEqual('bc', substr('abcd', 2, 3));
        end
        
        function testPartOfStringFromEnd(~)
            assertEqual('bc', substr('abcd', 2, -1));
        end
    end
end

function result = substr(s, i, j)
    result = edu.stanford.covert.util.substr(s, i, j);
end
