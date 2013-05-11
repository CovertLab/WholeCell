function [timestamp, platform, version] = readMatFileHeader(matFile)
%READMATFILEHEADER Reads timestamp of mat file as serial date number.
%
%References:
%===============
%MAT-File Format Version 7
%   http://www.serc.iisc.ernet.in/ComputingFacilities/software/matfile_format.pdf

fid = fopen(matFile, 'r');
tokens = regexp(fread(fid, [1 116], '*char'), ...
    'MATLAB (.*?) MAT-file, Platform: (.*?), Created on: (\w{3,3} \w{3,3}  ?\d{1,2} \d{2,2}:\d{2,2}:\d{2,2} \d{4,4})', ...
    'tokens');
fclose(fid);

if isempty(tokens) || numel(tokens) > 1
    throw(MException('readMatFileHeader:invalidFile', 'file is not valid mat file'));
end

version = tokens{1}{1};
platform = tokens{1}{2};
timestamp = datenum(tokens{1}{3}, 'ddd mmm dd HH:MM:SS yyyy');
end
