%getCloudSimulationStatus
% Returns BitMill job statuses. If jobId is specified returns status
% of just the 1 job. Otherwise returns statuses of all jobs.
%
% Optional input:
% - jobId [string]: bitmill job id
%
% Output:
% - result [struct]: struct containing job id, pool, submission time, status, and details
% - status [double]: 0==>success, otherwise==>failure
% - errMsg [char]: error message from BitMill
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/12/2013
function [result, status, errMsg] = getCloudSimulationStatus(varargin)
[result, status, errMsg] = com.numerate.bitmill.BitMill.getStatus(varargin{:});
