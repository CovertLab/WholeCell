classdef SetUtil_Test < TestCase
    methods
        function this = SetUtil_Test(name)
            this = this@TestCase(name);
        end
    end
    
    methods
        function testIsUnique(~)
            import edu.stanford.covert.util.SetUtil;
            
            assertFalse(SetUtil.isunique([0 1 0 0 2 3]));
            assertFalse(SetUtil.isunique([1 0 0 2 3]));
            assertFalse(SetUtil.isunique([0 1 0 0 2 3]'));
            assertFalse(SetUtil.isunique([0 1 0; 0 2 3]));
            assertTrue(SetUtil.isunique([3 1 4 -1 NaN]));
            assertTrue(SetUtil.isunique([3 1 4 -1 NaN NaN]));
        end
    end
end