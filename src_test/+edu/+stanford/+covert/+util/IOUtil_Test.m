% IOUtil test cases
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 9/5/2010
classdef IOUtil_Test < TestCase
    methods
        function this = IOUtil_Test(name)
            this = this@TestCase(name);
        end
    end
    
    %test saving/loading data as binary file
    methods
        function testStringArray(~)
            import edu.stanford.covert.util.IOUtil;
            tmp=char((2^16-1)*rand(100,100));
            IOUtil.writeBinary('tmp/tmp.bin',tmp);
            assertEqual(tmp,IOUtil.readBinary('tmp/tmp.bin'));
        end
        
        function testCellStringArray(~)
            import edu.stanford.covert.util.IOUtil;
            tmp=mat2cell(char((2^16-1)*rand(100,100)), ones(100,1), ones(100,1));
            IOUtil.writeBinary('tmp/tmp.bin',tmp);
            assertEqual(tmp,IOUtil.readBinary('tmp/tmp.bin'));
        end
        
        function testUint8Array(~)
            import edu.stanford.covert.util.IOUtil;
            tmp=uint8(ceil((2^8-1)*rand(4,3)));
            IOUtil.writeBinary('tmp/tmp.bin',tmp);
            assertEqual(tmp,IOUtil.readBinary('tmp/tmp.bin'));
        end
        
        function testUint16Array(~)
            import edu.stanford.covert.util.IOUtil;
            tmp=uint16(ceil((2^16-1)*rand(4,3)));
            IOUtil.writeBinary('tmp/tmp.bin',tmp);
            assertEqual(tmp,IOUtil.readBinary('tmp/tmp.bin'));
        end
        
        function testUint32Array(~)
            import edu.stanford.covert.util.IOUtil;
            tmp=uint32(ceil((2^32-1)*rand(4,3)));
            IOUtil.writeBinary('tmp/tmp.bin',tmp);
            assertEqual(tmp,IOUtil.readBinary('tmp/tmp.bin'));
        end
        
        function testInt8Array(~)
            import edu.stanford.covert.util.IOUtil;
            tmp=int8(ceil((2^7-1)*rand(4,3)));
            IOUtil.writeBinary('tmp/tmp.bin',tmp);
            assertEqual(tmp,IOUtil.readBinary('tmp/tmp.bin'));
        end
        
        function testInt16Array(~)
            import edu.stanford.covert.util.IOUtil;
            tmp=int16(ceil((2^15-1)*rand(4,3)));
            IOUtil.writeBinary('tmp/tmp.bin',tmp);
            assertEqual(tmp,IOUtil.readBinary('tmp/tmp.bin'));
        end
        
        function testInt32Array(~)
            import edu.stanford.covert.util.IOUtil;
            tmp=int32(ceil((2^31-1)*rand(4,3)));
            IOUtil.writeBinary('tmp/tmp.bin',tmp);
            assertEqual(tmp,IOUtil.readBinary('tmp/tmp.bin'));
        end
        
        function testSingleArray(~)
            import edu.stanford.covert.util.IOUtil;
            tmp=single(rand(4,3));
            IOUtil.writeBinary('tmp/tmp.bin',tmp);
            assertEqual(tmp,IOUtil.readBinary('tmp/tmp.bin'));
        end
        
        function testDoubleArray(~)
            import edu.stanford.covert.util.IOUtil;
            tmp=double(rand(4,3));
            IOUtil.writeBinary('tmp/tmp.bin',tmp);
            assertEqual(tmp,IOUtil.readBinary('tmp/tmp.bin'));
        end
        
        function testBinaryArray(~)
            import edu.stanford.covert.util.IOUtil;
            tmp=logical(rand(4,3)>0.5);
            IOUtil.writeBinary('tmp/tmp.bin',tmp);
            assertEqual(tmp,IOUtil.readBinary('tmp/tmp.bin'));
        end
    end
    
    methods
        function testMemmapFileReadWrite(~)
            import edu.stanford.covert.util.IOUtil;
            
            %write memmap file
            data = rand(20, 30, 40);
            fileName = 'output/runSmallTests/test.dat';
            IOUtil.writeMemmapFile(data, fileName);
            
            %check header written correctly
            assertEqual({class(data) size(data) IOUtil.memmapfileDataField}, IOUtil.readBinary(fileName));
            
            %check correct file size -- header + data
            fileInfo = dir(fileName);
            assertEqual(IOUtil.memmapfileHeaderSize + numel(data) * 8, fileInfo.bytes);
            
            %check data
            m = IOUtil.readMemmapFile(fileName);
            assertEqual(data, m.Data(1).(IOUtil.memmapfileDataField))
            
            %cleanup
            clear m;
            delete(fileName);
        end
    end
end
