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
dataFilePath = input(strrep(sprintf('Enter an absolute file path where simulated data will be stored [%s%soutput%srunSimulation]: ', pwd, filesep, filesep), '\','\\'), 's');
if ~isdir(dataFilePath)
    mkdir(dataFilePath);
end
if isempty(dataFilePath)
    dataFilePath = sprintf('%s%soutput%srunSimulation', pwd, filesep, filesep);
else
    dataFilePath = absolutepath(dataFilePath);
end

%edit stored simulation results configuration
fid = fopen('src/+edu/+stanford/+covert/+cell/+sim/+util/SimulationDiskUtil.m', 'r');
str = [];
while ~feof(fid)
    str = [str fgetl(fid) sprintf('\n')]; %#ok<AGROW>
end
fclose(fid);

str = strrep(str, '% value = ''/absolute_path/to/output/directory'';', ...
    sprintf('value = ''%s'';', strrep(dataFilePath, '''', '''''')));

fid = fopen('src/+edu/+stanford/+covert/+cell/+sim/+util/SimulationDiskUtil.m', 'w');
fwrite(fid, str);
fclose(fid);

%% set server configuration
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

%edit knowledge base configuration
fid = fopen('config.m', 'w');
fprintf(fid, 'function dbConnectionParameters = config()\n');
fprintf(fid, '%%CONFIG Returns configuration values.\n');
fprintf(fid, '%%\n');
fprintf(fid, '%% Author: Jonathan Karr, jkarr@stanford.edu\n');
fprintf(fid, '%% Affilitation: Covert Lab, Department of Bioengineering, Stanford University\n');
fprintf(fid, '%% Last updated: 9/21/2010\n');
fprintf(fid, '\n');
fprintf(fid, '%% knowledgebase\n');
fprintf(fid, 'dbConnectionParameters.hostName = ''%s'';\n', strrep(hostName, '''', ''''''));
fprintf(fid, 'dbConnectionParameters.schema   = ''%s'';\n', strrep(schema, '''', ''''''));
fprintf(fid, 'dbConnectionParameters.userName = ''%s'';\n', strrep(userName, '''', ''''''));
fprintf(fid, 'dbConnectionParameters.password = ''%s'';\n', strrep(password, '''', ''''''));
fclose(fid);