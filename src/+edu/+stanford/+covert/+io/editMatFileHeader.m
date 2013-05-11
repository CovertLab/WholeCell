function editMatFileHeader(matFile, timestamp, platform, version)
%EDITMATFILEHEADER Edits timestamp of .mat file.
%
% timestamp = the serial date number (optional)
%
%References:
%===============
%MAT-File Format Version 7
%   http://www.serc.iisc.ernet.in/ComputingFacilities/software/matfile_format.pdf

if ~exist('timestamp', 'var')
    timestamp = 730486;  %(serial date number)
end
if ~exist('platform', 'var')
    platform = 'PCWIN';
end
if ~exist('version', 'var')
    version = 5.0;
end
if ~ischar(version)
    version = num2str(version);
end

%check that file is valid mat file
edu.stanford.covert.io.readMatFileHeader(matFile);

%overwrite header
fid = fopen(matFile, 'r+');
header = sprintf('MATLAB %s MAT-file, Platform: %s, Created on: %s',...
    version, platform, datestr(timestamp, 'ddd mmm dd HH:MM:SS yyyy'));
fwrite(fid, [header repmat(' ', 1, 116 - length(header))], 'char');
fclose(fid);
end
