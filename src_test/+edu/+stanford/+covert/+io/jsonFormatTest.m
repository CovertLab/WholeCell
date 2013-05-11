classdef jsonFormatTest < TestCase
    
    methods
        function this = jsonFormatTest(name)
            this = this@TestCase(name);
        end
        
        function testChar(~)
            assertEqual('""', json(''));
            assertEqual('"a"', json('a'));
            assertEqual('"-\t-\n"', json(sprintf('-\t-\n')));
            assertEqual('"-\\t-\\n"', json('-\t-\n'));
            assertEqual('[[0,1],""]', json(reshape('',[0 1])));
            assertEqual('[[3,1,0],""]', json(reshape('',[3 1 0])));
            assertEqual('[[2,2],"ab12"]', json(['a1';'b2']));
        end

        function testLogical(~)
            assertEqual('true', json(true));
            assertEqual('false', json(false));
            assertEqual('[true,[0,0]]', json(logical.empty));
            assertEqual('[true,[3,1,0]]', json(true(3,1,0)));
            assertEqual('[true,[2,2],1,1,0,0]', json(logical([1 0; 1 0])));
        end
        
        function testDouble(~)
            assertEqual('0.0', json(0));
            assertEqual('1.0', json(1));
            assertEqual('3.141592653589793', json(pi));
            assertEqual('[]', json(double.empty));
            assertEqual('[3.0,1.0,0.0]', json([3 1 0]));
            assertEqual('[[3,1,0]]', json(zeros(3,1,0)));
            assertEqual('[[2,2],-1.0,0.0,-1.7976931348623157E308,1.7976931348623157E308]',...
                json([-1 -inf; 0 inf]));
            assertExceptionThrown(@() json([0 NaN; 0 0]), 'jsonFormat:NaN');
        end

        function testSingle(~)
            assertEqual('[1,[1,1],0.0]', json(single(0)));
            assertEqual('[1,[1,1],3.1415927]', json(single(pi)));
            assertEqual('[1,[0,0]]', json(single.empty));
            assertEqual('[1,[3,1,0]]', json(zeros(3,1,0,'single')));
            assertEqual('[1,[2,2],-1.0,0.0,-3.4028235E38,3.4028235E38]',...
                json(single([-1 -inf; 0 inf])));
            assertExceptionThrown(@() json(single([0 NaN; 0 0])), 'jsonFormat:NaN');
        end

        function testInt8(~)
            assertEqual('[7,[1,1],-128]', json(intmin('int8')));
            assertEqual('[7,[1,1],127]', json(intmax('int8')));
            assertEqual('[7,[0,0]]', json(int8.empty));
            assertEqual('[7,[3,1,0]]', json(zeros(3,1,0,'int8')));
            assertEqual('[7,[2,2],2,-123,23,0]', json(int8([2 23; -123 0])));
        end

        function testInt16(~)
            assertEqual('[15,[1,1],-32768]', json(intmin('int16')));
            assertEqual('[15,[1,1],32767]', json(intmax('int16')));
            assertEqual('[15,[0,0]]', json(int16.empty));
            assertEqual('[15,[3,1,0]]', json(zeros(3,1,0,'int16')));
            assertEqual('[15,[2,2],2,-234,23,0]', json(int16([2 23; -234 0])));
        end

        function testInt32(~)
            assertEqual('[31,[1,1],-2147483648]', json(intmin('int32')));
            assertEqual('[31,[1,1],2147483647]', json(intmax('int32')));
            assertEqual('[31,[0,0]]', json(int32.empty));
            assertEqual('[31,[3,1,0]]', json(zeros(3,1,0,'int32')));
            assertEqual('[31,[2,2],2,-234,23,0]', json(int32([2 23; -234 0])));
        end

        function testInt64(~)
            assertEqual('[63,[1,1],-9223372036854775808]', json(intmin('int64')));
            assertEqual('[63,[1,1],9223372036854775807]', json(intmax('int64')));
            assertEqual('[63,[0,0]]', json(int64.empty));
            assertEqual('[63,[3,1,0]]', json(zeros(3,1,0,'int64')));
            assertEqual('[63,[2,2],2,-234,23,0]', json(int64([2 23; -234 0])));
        end

        function testUint8(~)
            assertEqual('[8,[1,1],0]', json(intmin('uint8')));
            assertEqual('[8,[1,1],255]', json(intmax('uint8')));
            assertEqual('[8,[0,0]]', json(uint8.empty));
            assertEqual('[8,[3,1,0]]', json(zeros(3,1,0,'uint8')));
            assertEqual('[8,[2,2],2,234,23,0]', json(uint8([2 23; 234 0])));
        end

        function testUint16(~)
            assertEqual('[16,[1,1],0]', json(intmin('uint16')));
            assertEqual('[16,[1,1],65535]', json(intmax('uint16')));
            assertEqual('[16,[0,0]]', json(uint16.empty));
            assertEqual('[16,[3,1,0]]', json(zeros(3,1,0,'uint16')));
            assertEqual('[16,[2,2],2,234,23,0]', json(uint16([2 23; 234 0])));
        end

        function testUint32(~)
            assertEqual('[32,[1,1],0]', json(intmin('uint32')));
            assertEqual('[32,[1,1],4294967295]', json(intmax('uint32')));
            assertEqual('[32,[0,0]]', json(uint32.empty));
            assertEqual('[32,[3,1,0]]', json(zeros(3,1,0,'uint32')));
            assertEqual('[32,[2,2],2,234,23,0]', json(uint32([2 23; 234 0])));
        end

        function testUint64(~)
            assertEqual('[64,[1,1],0]', json(intmin('uint64')));
            if verLessThan('matlab', '7.11')  %num2str loses precision before 2010b
                assertEqual(['[64,[1,1],' num2str(intmax('uint64')) ']'], json(intmax('uint64')));
            else
                assertEqual('[64,[1,1],18446744073709551615]', json(intmax('uint64')));
            end
            assertEqual('[64,[0,0]]', json(uint64.empty));
            assertEqual('[64,[3,1,0]]', json(zeros(3,1,0,'uint64')));
            assertEqual('[64,[2,2],2,234,23,0]', json(uint64([2 23; 234 0])));
        end
        
        function testCellArray(~)
            assertEqual('[[]]', json(cell.empty));
            assertEqual('[[2,1,2],["jared","jacobs",0.0,1.0]]',...
                json(reshape({'jared',0;'jacobs',1}, [2 1 2])));
            assertEqual('[[2,3],["jane",432.0,[],false,{"a":0.0},-1.0]]',...
                json({'jane' [] struct('a',0); 432 false -1}));
        end

        function testStruct(~)
            assertEqual('[{},[0,0],[]]', json(struct.empty));
            assertEqual('{"b":432.0,"jane":[]}',...
                json(struct('jane', [], 'b', 432)));
            s = struct;
            s.p = cell(0);
            assertEqual('[{},[2,2],["p"],0.0,[[]],true,""]',...
                json([struct('p', 0) struct('p',true); s struct('p','')]));
        end

        function testClass(~)
            import edu.stanford.covert.util.SparseMat;
            assertEqual(...
                ['["edu.stanford.covert.util.SparseMat",'...
                 '{"siz":[0.0,2.0],"subs":[[0,2]],"vals":[[0,1]]}]'],...
                json(SparseMat()));
            assertEqual(...
                ['["edu.stanford.covert.util.SparseMat",'...
                 '{"siz":[4.0,3.0],"subs":[[2,2],3.0,4.0,2.0,3.0],"vals":[[2,1],1.0,-1.0]}]'],...
                json(SparseMat([0 0 0; 0 0 0; 0 1 0; 0 0 -1])));
        end
    end
end

function result = json(value)
    result = edu.stanford.covert.io.jsonFormat(value);
end
