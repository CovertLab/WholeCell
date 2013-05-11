classdef jsonParseTest < TestCase    
    methods
        function this = jsonParseTest(name)
            this = this@TestCase(name);
        end

        function testChar(~)
            assertEqual('', parse('""'));
            assertEqual('a', parse('"a"'));
            assertEqual(sprintf('-\t-\n'), parse('"-\t-\n"'));
            assertEqual('-\t-\n', parse('"-\\t-\\n"'));
            assertEqual(reshape('',[0 1]), parse('[[0,1],""]'));
            assertEqual(reshape('',[3 1 0]), parse('[[3,1,0],""]'));
            assertEqual(['a1';'b2'], parse('[[2,2],"ab12"]'));
        end

        function testLogical(~)
            assertEqual(true, parse('true'));
            assertEqual(false, parse('false'));
            assertEqual(logical.empty, parse('[true,[0,0]]'));
            assertEqual(true(3,1,0), parse('[true,[3,1,0]]'));
            assertEqual(logical([1 0; 1 0]), parse('[true,[2,2],1,1,0,0]'));
        end
        
        function testDouble(~)
            assertEqual(0, parse('0.0'));
            assertEqual(1, parse('1.0'));
            assertEqual(pi, parse('3.141592653589793'));
            assertEqual(double.empty, parse('[]'));
            assertEqual([3 1 0], parse('[3.0,1.0,0.0]'));
            assertEqual(ones(3,1,0), parse('[[3,1,0]]'));
            assertEqual([-1 -realmax; 0 realmax],...
                parse('[[2,2],-1.0,0.0,-1.7976931348623157E308,1.7976931348623157E308]'));
            assertExceptionThrown(@() parse('NaN'), 'jsonParse:InvalidJSON');
            assertExceptionThrown(@() parse('[[2,2],0,0,NaN,0]'),...
                'jsonParse:InvalidJSON');
        end

        function testSingle(~)
            assertEqual(single(0), parse('[1,[1,1],0.0]'));
            assertEqual(single(pi), parse('[1,[1,1],3.1415927]'));
            assertEqual(single.empty, parse('[1,[0,0]]'));
            assertEqual(ones(3,1,0,'single'), parse('[1,[3,1,0]]'));
            assertEqual(single([-1 -realmax('single'); 0 realmax('single')]),...
                parse('[1,[2,2],-1.0,0.0,-3.4028235E38,3.4028235E38]'));
            assertExceptionThrown(@() parse('[1,[2,2],0,0,NaN,0]'),...
                'jsonParse:InvalidJSON');
        end

        function testInt8(~)
            assertEqual(intmin('int8'), parse('[7,[1,1],-128]'));
            assertEqual(intmax('int8'), parse('[7,[1,1],127]'));
            assertEqual(int8.empty, parse('[7,[0,0]]'));
            assertEqual(ones(3,1,0,'int8'), parse('[7,[3,1,0]]'));
            assertEqual(int8([2 23; -123 0]), parse('[7,[2,2],2,-123,23,0]'));
        end

        function testInt16(~)
            assertEqual(intmin('int16'), parse('[15,[1,1],-32768]'));
            assertEqual(intmax('int16'), parse('[15,[1,1],32767]'));
            assertEqual(int16.empty, parse('[15,[0,0]]'));
            assertEqual(zeros(3,1,0,'int16'), parse('[15,[3,1,0]]'));
            assertEqual(int16([2 23; -234 0]), parse('[15,[2,2],2,-234,23,0]'));
        end

        function testInt32(~)
            assertEqual(intmin('int32'), parse('[31,[1,1],-2147483648]'));
            assertEqual(intmax('int32'), parse('[31,[1,1],2147483647]'));
            assertEqual(int32.empty, parse('[31,[0,0]]'));
            assertEqual(zeros(3,1,0,'int32'), parse('[31,[3,1,0]]'));
            assertEqual(int32([2 23; -234 0]), parse('[31,[2,2],2,-234,23,0]'));
        end

        function testInt64(~)
            assertEqual(intmin('int64'), parse('[63,[1,1],-9223372036854775808]'));
            assertEqual(intmax('int64'), parse('[63,[1,1],9223372036854775807]'));
            assertEqual(int64.empty, parse('[63,[0,0]]'));
            assertEqual(zeros(3,1,0,'int64'), parse('[63,[3,1,0]]'));
            assertEqual(int64([2 23; -234 0]), parse('[63,[2,2],2,-234,23,0]'));
        end

        function testUint8(~)
            assertEqual(intmin('uint8'), parse('[8,[1,1],0]'));
            assertEqual(intmax('uint8'), parse('[8,[1,1],255]'));
            assertEqual(uint8.empty, parse('[8,[0,0]]'));
            assertEqual(zeros(3,1,0,'uint8'), parse('[8,[3,1,0]]'));
            assertEqual(uint8([2 23; 234 0]), parse('[8,[2,2],2,234,23,0]'));
        end

        function testUint16(~)
            assertEqual(intmin('uint16'), parse('[16,[1,1],0]'));
            assertEqual(intmax('uint16'), parse('[16,[1,1],65535]'));
            assertEqual(uint16.empty, parse('[16,[0,0]]'));
            assertEqual(zeros(3,1,0,'uint16'), parse('[16,[3,1,0]]'));
            assertEqual(uint16([2 23; 234 0]), parse('[16,[2,2],2,234,23,0]'));
        end

        function testUint32(~)
            assertEqual(intmin('uint32'), parse('[32,[1,1],0]'));
            assertEqual(intmax('uint32'), parse('[32,[1,1],4294967295]'));
            assertEqual(uint32.empty, parse('[32,[0,0]]'));
            assertEqual(zeros(3,1,0,'uint32'), parse('[32,[3,1,0]]'));
            assertEqual(uint32([2 23; 234 0]), parse('[32,[2,2],2,234,23,0]'));
        end

        function testUint64(~)
            assertEqual(intmin('uint64'), parse('[64,[1,1],0]'));
            assertEqual(intmax('uint64'), parse('[64,[1,1],18446744073709551615]'));
            assertEqual(uint64.empty, parse('[64,[0,0]]'));
            assertEqual(zeros(3,1,0,'uint64'), parse('[64,[3,1,0]]'));
            assertEqual(uint64([2 23; 234 0]), parse('[64,[2,2],2,234,23,0]'));
        end
        
        function testCellArray(~)
            assertEqual(cell.empty, parse('[[]]'));
            assertEqual(reshape({'jared',0;'jacobs',1}, [2 1 2]),...
                parse('[[2,1,2],["jared","jacobs",0.0,1.0]]'));
            assertEqual({'jane' [] struct('a',0); 432 false -1},...
                parse('[[2,3],["jane",432.0,[],false,{"a":0.0},-1.0]]'));
        end

        function testStruct(~)
            assertEqual(struct.empty, parse('[{},[0,0],[]]'));
            assertEqual(struct('jane', [], 'b', 432),...
                parse('{"b":432.0,"jane":[]}'));
            s = struct;
            s.p = cell(0);
            assertEqual([struct('p', 0) struct('p',true); s struct('p','')],...
                parse('[{},[2,2],["p"],0.0,[[]],true,""]'));
        end

        function testClass(~)
            import edu.stanford.covert.util.SparseMat;
            assertEqual(SparseMat(), parse([
                '["edu.stanford.covert.util.SparseMat",'...
                 '{"siz":[0.0,2.0],"subs":[[0,2]],"vals":[[0,1]]}]']));
            assertEqual(SparseMat([0 0 0; 0 0 0; 0 1 0; 0 0 -1]), parse([
                '["edu.stanford.covert.util.SparseMat",'...
                 '{"siz":[4.0,3.0],"subs":[[2,2],3.0,4.0,2.0,3.0],"vals":[[2,1],1.0,-1.0]}]']));
        end
    end
end

function result = parse(jsonString)
    result = edu.stanford.covert.io.jsonParse(jsonString);
end
