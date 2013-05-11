classdef revisionTest < TestCase
    methods
        function this = revisionTest(name)
            this = this@TestCase(name);
        end
        
        function testRevision(~)
            [rev, differences] = edu.stanford.covert.util.revision();            
            
            assertTrue(isa(rev, 'double'));
            assertIn(rev, [1 inf]);
            assertEqual(fix(rev), rev);
            
            assertTrue(isa(differences, 'char'));
        end
    end
end
