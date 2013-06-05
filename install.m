function install
%Installs whole-cell model
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 2/2/2012

% check MATLAB version
if verLessThan('matlab', '7.9')
    throw(MException('install:newerMatlabVersionRequired', 'The whole-cell model software requires MATLAB R2009b (7.9) or newer'));
end

%% set warnings and path
setWarnings()
setPath();

%% configure stored simulation results path
outputPath = input(strrep(sprintf('Enter an absolute file path where simulated data will be stored [%s%soutput%srunSimulation]: ', pwd, filesep, filesep), '\','\\'), 's');
if ~isdir(outputPath)
    mkdir(outputPath);
end
if isempty(outputPath)
    outputPath = sprintf('%s%soutput%srunSimulation', pwd, filesep, filesep);
else
    outputPath = absolutepath(outputPath);
end

%% set server, s3cmd, bitmill-bash configurations
reply = ' ';
while ~(isequal(upper(reply), 'Y') || isequal(upper(reply), 'N') || isempty(reply))
    reply = input('Would you like to setup your server configuration? Y/N [N]: ', 's');
end
if isequal(upper(reply), 'N') || isempty(reply)
    return;
end

%prompt user
hostName = input('Enter knowledge base server hostname: ', 's');
schema = input('Enter knowledge base schema: ', 's');
userName = input('Enter knowledge base username: ', 's');
password = input('Enter knowledge base password: ', 's');
s3cmdPath = input('Enter path to parent of s3cmd (e.g. /usr/bin): ', 's');
bitmillBashPath = input('Enter path to bitmill-bash (e.g. /home/<user_name>/bitmill-bash): ', 's');

%edit knowledge base configuration
fid = fopen('getConfig.m', 'w');
fprintf(fid, 'function config = getConfig()\n');
fprintf(fid, '%%GETCONFIG Returns configuration values.\n');
fprintf(fid, '%%\n');
fprintf(fid, '%% Author: Jonathan Karr, jkarr@stanford.edu\n');
fprintf(fid, '%% Affilitation: Covert Lab, Department of Bioengineering, Stanford University\n');
fprintf(fid, '%% Last updated: 9/21/2010\n');
fprintf(fid, '\n');
fprintf(fid, '%% knowledgebase\n');
fprintf(fid, 'config.hostName = ''%s'';\n', strrep(hostName, '''', ''''''));
fprintf(fid, 'config.schema   = ''%s'';\n', strrep(schema, '''', ''''''));
fprintf(fid, 'config.userName = ''%s'';\n', strrep(userName, '''', ''''''));
fprintf(fid, 'config.password = ''%s'';\n', strrep(password, '''', ''''''));
fprintf(fid, 'config.outputPath = ''%s'';\n', strrep(outputPath, '''', ''''''));
fprintf(fid, 'config.s3cmdPath = ''%s'';\n', strrep(s3cmdPath, '''', ''''''));
fprintf(fid, 'config.bitmillBashPath = ''%s'';\n', strrep(bitmillBashPath, '''', ''''''));
fclose(fid);

%% create output folders
folders = {'bin', 'doc', 'output', 'tmp'};
for i = 1:numel(folders)
    if ~exists(folders{i}, 'dir')
        mkdir(folders{i});
    end
end

%% rebuild fixtures for current MATLAB version
generateTestFixtures(false);