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
% Examples
%   >> getCloudSimulationStatus()
%      ans = 34x1 struct array with fields:
%   >> getCloudSimulationStatus('82ddfe87-6263-4b40-a981-99bd17b8a68c')
%            id: '82ddfe87-6263-4b40-a981-99bd17b8a68c'
%          pool: 'ZscfWK4Bw8EAPJRC'
%     submitted: '2013-05-25T15:52:52.743-0700'
%        status: 'cancelled'
%       details: [1x1 struct]
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/12/2013
function [result, status, errMsg] = getCloudSimulationStatus(varargin)
[result, status, errMsg] = com.numerate.bitmill.BitMill.getStatus(varargin{:});
