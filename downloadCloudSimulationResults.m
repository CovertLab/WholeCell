%downloadCloudSimulationResults
% Downloads all simulation results from S3
% - <localFolder>/<simName>.parameters.mat
% - <localFolder>/<simName>.mat
% - <localFolder>/<simName>.out
% - <localFolder>/<simName>.err
%
% Input
% - simName [char]: prefix of simulation file names
% - bucketUrl [char]: S3 bucket where simulation stored
% - localFolder [char]: location where simulation files should be
%   downloaded to
%
% Ouput
% - status
% - errMsg
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/12/2013
function [status, errMsg] = downloadCloudSimulationResults(varargin)
%% parse inputs
ip = inputParser;

ip.addParamValue('simName', '', @(x) ischar(x));
ip.addParamValue('bucketUrl', '', @(x) ischar(x));
ip.addParamValue('localFolder', [], @(x) exist(x, 'dir'));

ip.parse(varargin{:});

simName = ip.Results.simName;
bucketUrl = ip.Results.bucketUrl;
localFolder = ip.Results.localFolder;

%% download
[status, errMsg] = com.numerate.bitmill.s3cmd.get(...
    sprintf('%s/%s.*', bucketUrl, simName), localFolder);