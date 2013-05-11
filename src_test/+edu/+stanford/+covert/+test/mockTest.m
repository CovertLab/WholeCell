classdef mockTest < TestCase
    methods
        function this = mockTest(name)
            this = this@TestCase(name);
        end
        
        function testConstructionAndPropertyAccess(~)
            m = edu.stanford.covert.test.mock('a', [1 2 3], 'b', 'foo');
            assertEqual([1 2 3], m.a);
            assertEqual('foo', m.b);
        end

        function testNestedStructures(~)
            m = edu.stanford.covert.test.mock('a', struct('d', [1 2 3]));
            assertEqual([1 2 3], m.a.d);
            assertEqual(struct('d', [1 2 3]), struct(m.a));
        end

        function testNewPropertyByAssignment(~)
            m = edu.stanford.covert.test.mock(struct('a', struct));
            m.b = [1 2 3];
            m.a.d.e = 'USA';
            assertEqual([1 2 3], m.b);
            assertEqual('USA', m.a.d.e);
        end

        function testPassByReference(~)
            m = edu.stanford.covert.test.mock();
            mutate(m, 'foo', 'bar');
            assertEqual('bar', m.foo);
        end
        
        function testCastToStructure(~)
            assertEqual(...
                struct('a', 'A'),...
                struct(edu.stanford.covert.test.mock('a', 'A')));
            assertEqual(...
                struct('a', struct('b', 0)),...
                struct(edu.stanford.covert.test.mock('a', struct('b', 0))));
        end
        
        function testParenthesesInSubsAsgn(~)
            m = edu.stanford.covert.test.mock('a', []);
            m.a(1) = 3;
            assertEqual(3, m.a(1));
        end
        
        function testCurlyBracesInSubsAsgn(~)
            m = edu.stanford.covert.test.mock();
            m.a = {struct};
            m.a{1}.d = {struct};
            m.a{1}.d{1}.e = 'USA';
            assertEqual('USA', m.a{1}.d{1}.e);
        end
        
        function testBindingOfFunctionHandles(~)
            function mock = f(mock)
                mock.calls = [mock.calls 'f'];
            end
            function result = g(mock, three, four)
                mock.calls = [mock.calls 'g'];
                assertEqual(3, three);
                assertEqual(4, four);
                result = 'G';
            end
            m = edu.stanford.covert.test.mock('calls', '', 'foo', @f);
            m.bar.baz = @g;
            m.foo();
            assertEqual('edu.stanford.covert.test.mock', class(m.foo()));
            assertEqual('G', m.bar.baz(3, 4));
            assertEqual('ffg', m.calls);
        end
    end
end

function mutate(o, name, val)
    o.(name) = val;
end
