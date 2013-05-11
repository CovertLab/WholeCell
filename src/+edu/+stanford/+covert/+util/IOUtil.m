% IO utility functions to
% - save data to files.
% - read, edit mat file timestamps
%
% Binary Format for matrices and cell arrays (with matrix and cell array
% elements) inspired by Arvid Bottiger's write_matrix_bin,
% 1. 1*uint32               Indicator of cell array
% 2. 1*uint32               Dimensions of array
% 3. dimensions*uint32      Lengths of dimensions
% 4. a. Cell Array => recurse
%    b. Matrix
%       1*uint32            Indicator of character array
%       uint32 or float32   Dat%
%
%
%References:
%===============
%- write_matrix_bin by Arvid Bottiger
%   http://www.mathworks.com/matlabcentral/fileexchange/24483
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/14/2009
classdef IOUtil
    properties (Constant = true)
        memmapfileHeaderSize = 1024
        memmapfileDataField = 'data'
    end
    
    methods (Static)
        %Writes data to disk in binary format
        function writeBinary(filename, data)
            import edu.stanford.covert.util.IOUtil;
            
            %open file
            fid = fopen(filename, 'w');
            
            %write data
            IOUtil.writeDataBinary(fid, data);
            
            %close file
            fclose(fid);
        end
        
        %Reads data from disk in binary format
        function data = readBinary(filename)
            import edu.stanford.covert.util.IOUtil;
            
            %open file
            fid = fopen(filename, 'r');
            
            %read data
            data = IOUtil.readDataBinary(fid);
            
            %close file
            fclose(fid);
        end
    end
    
    methods (Static, Access = protected)
        %Writes data to disk in binary format
        function writeDataBinary(fid, data)
            import edu.stanford.covert.util.IOUtil;
            
            %cell array, number of dimensions, size
            dataIsCell = iscell(data);
            dataSize = size(data);
            dimensions = length(dataSize);
            fwrite(fid, [dataIsCell dimensions dataSize], 'uint32');
            
            %reshape data
            data = reshape(data,1,[]);
            
            %write data
            if dataIsCell
                for i = 1:length(data)
                    IOUtil.writeDataBinary(fid, data{i});
                end
            else
                switch class(data)
                    case 'char',    dataType = 1;
                    case 'uint8',   dataType = 2;
                    case 'uint16',  dataType = 3;
                    case 'uint32',  dataType = 4;
                    case 'uint64',  dataType = 5;
                    case 'int8',    dataType = 6;
                    case 'int16',   dataType = 7;
                    case 'int32',   dataType = 8;
                    case 'int64',   dataType = 9;
                    case 'single',  dataType = 10;
                    case 'double',  dataType = 11;
                    case 'logical', dataType = 12;
                    otherwise
                        throw(MException('IOUtil:unsupportedDatatype', 'Data type ''%s'' not supported', class(data)));
                end
                
                %character, write no. dimensions, size
                fwrite(fid,  dataType, 'uint32');
                
                %write data
                switch dataType
                    case 1,  fwrite(fid, data+0, 'uint16');
                    case 2,  fwrite(fid, data, 'uint16');  %Note not using uint8
                    case 3,  fwrite(fid, data, 'uint16');
                    case 4,  fwrite(fid, data, 'uint32');
                    case 5,  fwrite(fid, data, 'float64'); %Note not using uint64
                    case 6,  fwrite(fid, data, 'int16');   %Note not using int8
                    case 7,  fwrite(fid, data, 'int16');
                    case 8,  fwrite(fid, data, 'int32');
                    case 9,  fwrite(fid, data, 'float64'); %Note not using int64
                    case 10, fwrite(fid, data, 'float32');
                    case 11, fwrite(fid, data, 'float64');
                    case 12, fwrite(fid, data, 'uint16');  %Note not using ubit1
                    otherwise
                        throw(MException('IOUtil:unsupportedDatatype', 'Data type ''%s'' not supported', dataType));
                end
            end
        end
        
        %Reads data from disk in binary format
        function data = readDataBinary(fid)
            import edu.stanford.covert.util.IOUtil;
            
            %cell array, number of dimensions, size
            dataIsCell = fread(fid, 1, 'uint32');
            dimensions = fread(fid, 1, 'uint32');
            dataSize = fread(fid, dimensions, 'uint32')';
            count = prod(dataSize);
            
            %read data
            if dataIsCell
                data = cell(count, 1);
                for i = 1:count
                    data{i} = IOUtil.readDataBinary(fid);
                end
            else
                dataType = fread(fid, 1, 'uint32');
                switch dataType
                    case 1,  data = char(fread(fid, count, 'uint16')); %#ok<FREAD>
                    case 2,  data = uint8(fread(fid, count, 'uint16'));
                    case 3,  data = uint16(fread(fid, count, 'uint16'));
                    case 4,  data = uint32(fread(fid, count, 'uint32'));
                    case 5,  data = uint64(fread(fid, count, 'uint64'));
                    case 6,  data = int8(fread(fid, count, 'int16'));
                    case 7,  data = int16(fread(fid, count, 'int16'));
                    case 8,  data = int32(fread(fid, count, 'int32'));
                    case 9,  data = int64(fread(fid, count, 'int64'));
                    case 10, data = single(fread(fid, count, 'float32'));
                    case 11, data = double(fread(fid, count, 'float64'));
                    case 12, data = logical(fread(fid, count, 'uint16'));
                    otherwise
                        throw(MException('IOUtil:unsupportedDatatype', 'Data type ''%d'' not supported', dataType));
                end
            end
            
            if dimensions > 1
                data = reshape(data, dataSize);
            end
        end
    end
    
    methods (Static)
        function writeMemmapFile(data, fileName)
            import edu.stanford.covert.util.IOUtil;
            
            %write header
            format = {class(data) size(data) IOUtil.memmapfileDataField};
            IOUtil.writeBinary(fileName, format);
            
            %write spacing between header and data, and data
            fileInfo = dir(fileName);
            if fileInfo.bytes > IOUtil.memmapfileHeaderSize
                throw(MException('IOUtil:invalidHeader', 'Header size cannot be greater than %d bytes', IOUtil.memmapfileHeaderSize));
            end
            
            fid = fopen(fileName, 'a');
            fwrite(fid, zeros(IOUtil.memmapfileHeaderSize - fileInfo.bytes, 1, 'uint8'), 'uint8'); %spacing
            fwrite(fid, data, class(data)); %data
            fclose(fid);
        end
        
        function result = readMemmapFile(fileName)
            import edu.stanford.covert.util.IOUtil;
            
            format = IOUtil.readBinary(fileName);
            if all(format{2})
                result = memmapfile(fileName, ...
                    'Offset', IOUtil.memmapfileHeaderSize, ...
                    'Format', format, ...
                    'Repeat', 1, ...
                    'Writable', false);
            else
                result = struct('Data', struct(format{3}, []));
                result.Data.(format{3}) = zeros(format{2}, format{1});
            end
        end
    end
    
    methods (Static)
        function directories = getDirectoryNamesRecursively(directory)
            import edu.stanford.covert.util.IOUtil;
            if directory(end) == '/' || directory(end) == '\'
                directory = directory(1:end-1);
            end
            directories = {directory};
            files = dir(directory);
			files = files([files.isdir]);
            for i = 1:numel(files)
                if files(i).name(1) ~= '.'
                    directories = [directories;
                        IOUtil.getDirectoryNamesRecursively([directory filesep files(i).name])]; %#ok<AGROW>
                end
            end
        end
    end
end