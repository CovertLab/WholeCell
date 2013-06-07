%postCloudSimulation
% Posts job to BitMill queue. Parameters and results are saved to several
% S3 files:
% - <bucketURL>/<simName>.parameters.{mat|xml}: file with specified parameter values
% - <bucketURL>/<simName>.mat: average in silico experiments over the simulated cell population
% - <bucketURL>/<simName>.out: concatentation of stdout from the individual simulations and averaging calculation
% - <bucketURL>/<simName>.err: concatentation of stderr from the individual simulations and averaging calculation
%
% Input:
% - parameterVals [struct]: struct containing desired values of simulation
%   parameters. Initialize using the getAllParameters method of the
%   simulation class to get default parameter values. Overwrite struct
%   elements to set parameter values.
% - parameterValsPath [.mat/.xml file path]: .mat file path for stored
%   struct of parameter values or .xml file describing parameter values.
%   Use http://wholecell.stanford.edu/simulation/runSimulations.php to
%   generate XML file.
% - bucketURL [char]: S3 bucket URL e.g. s3://dream-sims
% - simName [char]: textual name of simulation
%
% Output:
% - jobId [char]: BitMill ID of new job
% - status [double]: 0==>success, otherwise==>failure
% - errMsg [char]: error message from BitMill
%
% Example:
%   >> [jobId, status, errMsg] = postCloudSimulation('simName', 'test20', 'bucketUrl', 's3://<your_bucket>')
%      jobId = 82ddfe87-6263-4b40-a981-99bd17b8a68c
%      status = 0
%      errMsg = []
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/12/2013
function [jobId, status, errMsg] = postCloudSimulation(varargin)
import com.numerate.bitmill.BitMill;
import com.numerate.bitmill.s3cmd;

%% parse inputs
ip = inputParser;

ip.addParamValue('parameterVals', [], @(x) isstruct(x));
ip.addParamValue('parameterValsPath', [], @(x) exist(x, 'file'));
ip.addParamValue('bucketUrl', '', @(x) ischar(x));
ip.addParamValue('simName', '', @(x) ischar(x));

ip.parse(varargin{:});

parameterVals = ip.Results.parameterVals;
parameterValsPath = ip.Results.parameterValsPath;
bucketUrl = ip.Results.bucketUrl;
simName = ip.Results.simName;

if ~isempty(parameterVals) && ~isempty(parameterValsPath)
    throw(MException('postCloudSimulation:error', 'Only 1 of parameterValsPath and parameterValsPath can be specified'))
end

validateattributes(bucketUrl, {'char'}, {'nonempty'});
validateattributes(simName, {'char'}, {'nonempty'});

%% Chose S3 URLs
remoteParameterValsPath = sprintf('%s/%s.parameters.mat', bucketUrl, simName);
remotePredictionsPath = sprintf('%s/%s.predictions.mat', bucketUrl, simName);
remoteDistancesPath = sprintf('%s/%s.distances.mat', bucketUrl, simName);
remoteStdoutPath = sprintf('%s/%s.out', bucketUrl, simName);
remoteStderrPath = sprintf('%s/%s.err', bucketUrl, simName);

%% save parameter values to S3
%check that user doesn't already have file
[isFile, status] = s3cmd.exist(remoteParameterValsPath); if status ~= 0 ||  isFile, throw(MException('postCloudSimulation:error', 'File already exists: %s. Choose a unique simName.', remoteParameterValsPath)), end
[isFile, status] = s3cmd.exist(remoteDistancesPath);     if status ~= 0 ||  isFile, throw(MException('postCloudSimulation:error', 'File already exists: %s. Choose a unique simName.', remoteDistancesPath)),     end
[isFile, status] = s3cmd.exist(remoteStdoutPath);        if status ~= 0 ||  isFile, throw(MException('postCloudSimulation:error', 'File already exists: %s. Choose a unique simName.', remoteStdoutPath)),        end
[isFile, status] = s3cmd.exist(remoteStderrPath);        if status ~= 0 ||  isFile, throw(MException('postCloudSimulation:error', 'File already exists: %s. Choose a unique simName.', remoteStderrPath)),        end

%save to temporary file if necessary
saveParameterVals = isempty(parameterValsPath);
if saveParameterVals
    if isempty(parameterVals)
        parameterVals = struct(); %#ok<NASGU>
    end
    parameterValsPath = sprintf('%s.mat', tempname);
    save(parameterValsPath, '-struct', 'parameterVals');
end

%copy to S3 and grant permissions
[status, errMsg] = s3cmd.put(parameterValsPath, remoteParameterValsPath);
if status ~= 0
    throw(MException('postCloudSimulation:error', 'Unable to upload file: %s', errMsg));
end
s3cmd.grantAclRead(remoteParameterValsPath, BitMill.s3Account);

%cleanup temporary file
if saveParameterVals
    delete(parameterValsPath);
end

%% Post job
type = 'dream_sim';
parameters = [
    struct('name', 'notification_email', 'value', BitMill.getNotificationEmail())
    ]; %#ok<NBRAK>
inputs = [
    struct('name', 'parameters', 'url', remoteParameterValsPath)
    ]; %#ok<NBRAK>
outputs = [
    struct('name', 'predictions', 'url', remotePredictionsPath)
    struct('name', 'distances', 'url', remoteDistancesPath)    
    struct('name', 'stdout', 'url', remoteStdoutPath)
    struct('name', 'stderr', 'url', remoteStderrPath)
    ];

[jobId, status, errMsg] = BitMill.post(...
    'type', type, ...
    'parameters', parameters, ...
    'inputs', inputs, ...
    'outputs', outputs ...
    );
if isempty(jobId) || status ~= 0
    throw(MException('postCloudSimulationStatus:error', 'Unable to post simulation job to BitMill: %s', errMsg))
end
